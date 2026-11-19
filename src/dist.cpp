#define USE_FC_LEN_T
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <errno.h>
#include <limits>
#include <sys/mman.h>
#include <unistd.h>
#include <vector>

#ifndef R_NO_REMAP
#  define R_NO_REMAP
#endif
#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>

#ifdef length
#undef length
#endif
#ifdef error
#undef error
#endif
#ifndef FCONE
#define FCONE
#endif

#include "ProgressReporter.h"
#include "tgstat.h"

// Returns the "tidy" data frame of the distances that are less or equal than the threshold.
// res holds the distances in the order of the 'dist' object.
static SEXP tidy_dist(const double *res, uint64_t num_points, double threshold, SEXP _rrownames)
{
    enum { ROW1, ROW2, DIST, NUM_COLS };
    const char *COL_NAMES[NUM_COLS] = { "row1", "row2", "dist" };
    uint64_t res_size = num_points * (num_points - 1) / 2;
    uint64_t answer_size = 0;
    SEXP answer, rcol1, rcol2, rdist, rrownames, rcolnames;

    for (uint64_t i = 0; i < res_size; ++i) {
        if (res[i] <= threshold)
            ++answer_size;
    }

    rprotect(answer = RSaneAllocVector(VECSXP, NUM_COLS));
    rprotect(rcol1 = RSaneAllocVector(INTSXP, answer_size));
    rprotect(rcol2 = RSaneAllocVector(INTSXP, answer_size));
    rprotect(rdist = RSaneAllocVector(REALSXP, answer_size));
    rprotect(rcolnames = RSaneAllocVector(STRSXP, NUM_COLS));
    rprotect(rrownames = RSaneAllocVector(INTSXP, answer_size));

    for (uint64_t i = 0; i < NUM_COLS; i++)
        SET_STRING_ELT(rcolnames, i, Rf_mkChar(COL_NAMES[i]));

    uint64_t idx1 = 0;
    uint64_t idx2 = 0;

    for (uint64_t ipoint1 = 0; ipoint1 < num_points; ++ipoint1) {
        for (uint64_t ipoint2 = ipoint1 + 1; ipoint2 < num_points; ++ipoint2) {
            if (res[idx2] <= threshold) {
                INTEGER(rcol1)[idx1] = ipoint1 + 1;
                INTEGER(rcol2)[idx1] = ipoint2 + 1;
                REAL(rdist)[idx1] = res[idx2];
                INTEGER(rrownames)[idx1] = idx1 + 1;
                ++idx1;
            }
            ++idx2;
        }
    }

    if (_rrownames != R_NilValue) {
        Rf_setAttrib(rcol1, R_LevelsSymbol, _rrownames);
        Rf_setAttrib(rcol1, R_ClassSymbol, Rf_mkString("factor"));
        Rf_setAttrib(rcol2, R_LevelsSymbol, _rrownames);
        Rf_setAttrib(rcol2, R_ClassSymbol, Rf_mkString("factor"));
    }

    SET_VECTOR_ELT(answer, ROW1, rcol1);
    SET_VECTOR_ELT(answer, ROW2, rcol2);
    SET_VECTOR_ELT(answer, DIST, rdist);

    Rf_setAttrib(answer, R_NamesSymbol, rcolnames);
    Rf_setAttrib(answer, R_ClassSymbol, Rf_mkString("data.frame"));
    Rf_setAttrib(answer, R_RowNamesSymbol, rrownames);
    return answer;
}

static void set_attrs(SEXP answer, SEXP _attrs)
{
    SEXP names = Rf_getAttrib(_attrs, R_NamesSymbol);
    for (int i = 0; i < LENGTH(_attrs); i++)
        Rf_setAttrib(answer, Rf_install(Rf_translateChar(STRING_ELT(names, i))), VECTOR_ELT(_attrs, i));
}

static SEXP get_rownames(SEXP rx)
{
    SEXP rdimnames = Rf_getAttrib(rx, R_DimNamesSymbol);
    return Rf_isNull(rdimnames) ? R_NilValue : VECTOR_ELT(rdimnames, 0);
}

static void set_factor(SEXP rcol, SEXP rlevels)
{
    if (rlevels != R_NilValue) {
        Rf_setAttrib(rcol, R_LevelsSymbol, rlevels);
        Rf_setAttrib(rcol, R_ClassSymbol, Rf_mkString("factor"));
    }
}

// Checks that x and y (unless it is NULL) are numeric matrices with the same number of columns
static void check_points(SEXP _x, SEXP _y)
{
    SEXP rmat[2] = { _x, _y };
    const char *mat_name[2] = { "x", "y" };

    for (int k = 0; k < (Rf_isNull(_y) ? 1 : 2); ++k) {
        SEXP rdim;

        if ((!Rf_isReal(rmat[k]) && !Rf_isInteger(rmat[k])) || Rf_xlength(rmat[k]) < 1 ||
            !Rf_isInteger(rdim = Rf_getAttrib(rmat[k], R_DimSymbol)) || Rf_xlength(rdim) != 2)
            verror("\"%s\" argument must be a matrix of numeric values", mat_name[k]);
    }

    if (!Rf_isNull(_y) && Rf_ncols(_x) != Rf_ncols(_y))
        verror("\"x\" and \"y\" matrices have different number of columns");
}

static uint64_t check_knn(SEXP _knn, SEXP _threshold)
{
    if ((!Rf_isReal(_knn) && !Rf_isInteger(_knn)) || Rf_xlength(_knn) != 1 || Rf_asReal(_knn) < 1 || Rf_asReal(_knn) != (double)Rf_asInteger(_knn))
        verror("\"knn\" argument must be a positive integer");

    if ((!Rf_isReal(_threshold) && !Rf_isInteger(_threshold)) || Rf_xlength(_threshold) != 1)
        verror("\"threshold\" argument must be a numeric value");

    return Rf_asInteger(_knn);
}

// Copies the rows of the matrix one after the other (v[ipoint * num_dims + idim]), with NaN for missing and
// non-finite values, and flags the rows that have them
static void read_points(SEXP rx, vector<double> &v, vector<bool> &nan_in_point)
{
    uint64_t num_points = Rf_nrows(rx);
    uint64_t num_dims = Rf_ncols(rx);
    const double *xreal = Rf_isReal(rx) ? REAL(rx) : NULL;
    const int *xint = Rf_isReal(rx) ? NULL : INTEGER(rx);

    v.resize(num_points * num_dims);
    nan_in_point.assign(num_points, false);

    for (uint64_t idim = 0, i = 0; idim < num_dims; ++idim) {
        for (uint64_t ipoint = 0; ipoint < num_points; ++ipoint, ++i) {
            if (xreal ? !R_FINITE(xreal[i]) : xint[i] == NA_INTEGER) {
                v[ipoint * num_dims + idim] = numeric_limits<double>::quiet_NaN();
                nan_in_point[ipoint] = true;
            } else
                v[ipoint * num_dims + idim] = xreal ? xreal[i] : xint[i];
        }
    }
}

// Copies the values of the matrix (column-major) into memory aligned for BLAS, with NaN for missing and non-finite
// values. Returns true if there are any.
static bool read_blas_points(SEXP rx, double *&m)
{
    uint64_t num_vals = Rf_xlength(rx);
    const double *xreal = Rf_isReal(rx) ? REAL(rx) : NULL;
    const int *xint = Rf_isReal(rx) ? NULL : INTEGER(rx);
    bool nan_in_vals = false;

    // some BLAS implementations ask to align double arrays to 64 for improved efficiency
    if (posix_memalign((void **)&m, 64, sizeof(double) * num_vals))
        verror("%s", strerror(errno));

    for (uint64_t i = 0; i < num_vals; ++i) {
        if (xreal ? !R_FINITE(xreal[i]) : xint[i] == NA_INTEGER) {
            m[i] = numeric_limits<double>::quiet_NaN();
            nan_in_vals = true;
        } else
            m[i] = xreal ? xreal[i] : xint[i];
    }
    return nan_in_vals;
}

// The distance between two points as in tgs_dist: the squared distance over the dimensions where both values are not
// NaN is scaled by num_dims / count, and the distance is NA if there are no such dimensions. nan: whether one of the
// points has a NaN.
static double point_dist(const double *p1, const double *p2, uint64_t num_dims, bool nan)
{
    double dist = 0;

    if (!nan) {
        // four partial sums, which are added in parallel
        double sums[4] = { 0, 0, 0, 0 };
        uint64_t idim = 0;

        for (; idim + 4 <= num_dims; idim += 4) {
            for (int k = 0; k < 4; ++k) {
                double dev = p1[idim + k] - p2[idim + k];
                sums[k] += dev * dev;
            }
        }
        for (; idim < num_dims; ++idim) {
            double dev = p1[idim] - p2[idim];
            sums[0] += dev * dev;
        }
        return sqrt((sums[0] + sums[1]) + (sums[2] + sums[3]));
    }

    uint64_t count = 0;

    for (uint64_t idim = 0; idim < num_dims; ++idim) {
        if (!std::isnan(p1[idim]) && !std::isnan(p2[idim])) {
            double dev = p1[idim] - p2[idim];
            dist += dev * dev;
            ++count;
        }
    }

    if (!count)
        return NA_REAL;
    return count == num_dims ? sqrt(dist) : sqrt(dist * num_dims / (double)count);
}

// The squared distance between two rows of column-major matrices with num_points1 and num_points2 rows, over the
// dimensions where both masks are 1 (all the dimensions if mask1 is NULL)
static double direct_sq_dist(const double *p1, const double *mask1, uint64_t num_points1, const double *p2, const double *mask2,
                             uint64_t num_points2, uint64_t num_dims)
{
    double dist = 0;

    for (uint64_t idim = 0; idim < num_dims; ++idim) {
        if (!mask1 || (mask1[idim * num_points1] && mask2[idim * num_points2])) {
            double dev = p1[idim * num_points1] - p2[idim * num_points2];
            dist += dev * dev;
        }
    }
    return dist;
}

// The sum of squares of each row of a column-major matrix
static vector<double> row_sums_square(const double *m, uint64_t num_points, uint64_t num_dims)
{
    vector<double> res(num_points, 0.);

    for (uint64_t idim = 0, i = 0; idim < num_dims; ++idim) {
        for (uint64_t ipoint = 0; ipoint < num_points; ++ipoint, ++i)
            res[ipoint] += m[i] * m[i];
    }
    return res;
}

// With BLAS the squared distance between x and y is computed as |x|^2 + |y|^2 - 2 x.y, whose rounding error is up to
// about num_dims * DBL_EPSILON * (|x|^2 + |y|^2). When the result is below this error times SQ_DIST_TOL the squared
// distance is computed directly instead, so that identical points are at distance 0.
static const double SQ_DIST_TOL = 4;

// Returns the "tidy" data frame of the distances between the rows of x and the rows of y that are less or equal
// than the threshold. res is the num_points1 x num_points2 matrix of the distances.
static SEXP tidy_cross_dist(const double *res, uint64_t num_points1, uint64_t num_points2, double threshold, SEXP rrownames1, SEXP rrownames2)
{
    enum { ROW1, ROW2, DIST, NUM_COLS };
    const char *COL_NAMES[NUM_COLS] = { "row1", "row2", "dist" };
    uint64_t res_size = num_points1 * num_points2;
    uint64_t answer_size = 0;
    SEXP answer, rcol1, rcol2, rdist, rrownames, rcolnames;

    for (uint64_t i = 0; i < res_size; ++i) {
        if (res[i] <= threshold)
            ++answer_size;
    }

    rprotect(answer = RSaneAllocVector(VECSXP, NUM_COLS));
    rprotect(rcol1 = RSaneAllocVector(INTSXP, answer_size));
    rprotect(rcol2 = RSaneAllocVector(INTSXP, answer_size));
    rprotect(rdist = RSaneAllocVector(REALSXP, answer_size));
    rprotect(rcolnames = RSaneAllocVector(STRSXP, NUM_COLS));
    rprotect(rrownames = RSaneAllocVector(INTSXP, answer_size));

    for (uint64_t i = 0; i < NUM_COLS; i++)
        SET_STRING_ELT(rcolnames, i, Rf_mkChar(COL_NAMES[i]));

    for (uint64_t ipoint1 = 0, idx = 0; ipoint1 < num_points1; ++ipoint1) {
        for (uint64_t ipoint2 = 0; ipoint2 < num_points2; ++ipoint2) {
            double dist = res[ipoint1 + ipoint2 * num_points1];

            if (dist <= threshold) {
                INTEGER(rcol1)[idx] = ipoint1 + 1;
                INTEGER(rcol2)[idx] = ipoint2 + 1;
                REAL(rdist)[idx] = dist;
                INTEGER(rrownames)[idx] = idx + 1;
                ++idx;
            }
        }
    }

    set_factor(rcol1, rrownames1);
    set_factor(rcol2, rrownames2);

    SET_VECTOR_ELT(answer, ROW1, rcol1);
    SET_VECTOR_ELT(answer, ROW2, rcol2);
    SET_VECTOR_ELT(answer, DIST, rdist);

    Rf_setAttrib(answer, R_NamesSymbol, rcolnames);
    Rf_setAttrib(answer, R_ClassSymbol, Rf_mkString("data.frame"));
    Rf_setAttrib(answer, R_RowNamesSymbol, rrownames);
    return answer;
}

static void set_cross_dist_attrs(SEXP answer, uint64_t num_points1, uint64_t num_points2, SEXP rrownames1, SEXP rrownames2)
{
    SEXP rdim;

    rprotect(rdim = RSaneAllocVector(INTSXP, 2));
    INTEGER(rdim)[0] = num_points1;
    INTEGER(rdim)[1] = num_points2;
    Rf_setAttrib(answer, R_DimSymbol, rdim);

    if (rrownames1 != R_NilValue || rrownames2 != R_NilValue) {
        SEXP rdimnames;

        rprotect(rdimnames = RSaneAllocVector(VECSXP, 2));
        SET_VECTOR_ELT(rdimnames, 0, rrownames1);
        SET_VECTOR_ELT(rdimnames, 1, rrownames2);
        Rf_setAttrib(answer, R_DimNamesSymbol, rdimnames);
    }
}

struct Neighbour {
    double dist;
    int    row;
};

// nearer first, ties broken by the lower row
static bool nearer(const Neighbour &n1, const Neighbour &n2)
{
    return n1.dist < n2.dist || (n1.dist == n2.dist && n1.row < n2.row);
}

// Adds a neighbour to the heap of the knn nearest neighbours found so far, which has the farthest at its front
static void add_neighbour(vector<Neighbour> &heap, uint64_t knn, const Neighbour &neighbour)
{
    if (heap.size() == knn)
        pop_heap(heap.begin(), heap.end(), nearer);
    else
        heap.resize(heap.size() + 1);
    heap.back() = neighbour;
    push_heap(heap.begin(), heap.end(), nearer);
}

// Returns the data frame of the neighbours: neighbours[ipoint * knn] holds the num_neighbours[ipoint] rows of y nearest
// to row ipoint of x, nearest first
static SEXP knn_answer(const Neighbour *neighbours, const uint64_t *num_neighbours, uint64_t num_points, uint64_t knn, SEXP rrownames1, SEXP rrownames2)
{
    enum { ROW1, ROW2, DIST, RANK, NUM_COLS };
    const char *COL_NAMES[NUM_COLS] = { "row1", "row2", "dist", "rank" };
    uint64_t answer_size = 0;
    SEXP answer, rcol1, rcol2, rdist, rrank, rrownames, rcolnames;

    for (uint64_t ipoint = 0; ipoint < num_points; ++ipoint)
        answer_size += num_neighbours[ipoint];

    rprotect(answer = RSaneAllocVector(VECSXP, NUM_COLS));
    rprotect(rcol1 = RSaneAllocVector(INTSXP, answer_size));
    rprotect(rcol2 = RSaneAllocVector(INTSXP, answer_size));
    rprotect(rdist = RSaneAllocVector(REALSXP, answer_size));
    rprotect(rrank = RSaneAllocVector(INTSXP, answer_size));
    rprotect(rcolnames = RSaneAllocVector(STRSXP, NUM_COLS));
    rprotect(rrownames = RSaneAllocVector(INTSXP, answer_size));

    for (int i = 0; i < NUM_COLS; i++)
        SET_STRING_ELT(rcolnames, i, Rf_mkChar(COL_NAMES[i]));

    for (uint64_t ipoint = 0, idx = 0; ipoint < num_points; ++ipoint) {
        for (uint64_t i = 0; i < num_neighbours[ipoint]; ++i, ++idx) {
            INTEGER(rcol1)[idx] = ipoint + 1;
            INTEGER(rcol2)[idx] = neighbours[ipoint * knn + i].row + 1;
            REAL(rdist)[idx] = neighbours[ipoint * knn + i].dist;
            INTEGER(rrank)[idx] = i + 1;
            INTEGER(rrownames)[idx] = idx + 1;
        }
    }

    set_factor(rcol1, rrownames1);
    set_factor(rcol2, rrownames2);

    SET_VECTOR_ELT(answer, ROW1, rcol1);
    SET_VECTOR_ELT(answer, ROW2, rcol2);
    SET_VECTOR_ELT(answer, DIST, rdist);
    SET_VECTOR_ELT(answer, RANK, rrank);

    Rf_setAttrib(answer, R_NamesSymbol, rcolnames);
    Rf_setAttrib(answer, R_ClassSymbol, Rf_mkString("data.frame"));
    Rf_setAttrib(answer, R_RowNamesSymbol, rrownames);
    return answer;
}

extern "C" {

SEXP tgs_dist(SEXP _x, SEXP _attrs, SEXP _tidy, SEXP _threshold, SEXP _rrownames, SEXP _envir)
{
    SEXP answer = R_NilValue;
    void *shm = (double *)MAP_FAILED;
    uint64_t shm_sizeof = 0;

	try {
        TGStat tgstat(_envir);

		if ((!Rf_isReal(_x) && !Rf_isInteger(_x)) || Rf_xlength(_x) < 1)
			verror("\"x\" argument must be a matrix of numeric values");

        if (!Rf_isLogical(_tidy) || Rf_xlength(_tidy) != 1)
            verror("\"tidy\" argument must be a logical value");

        if ((!Rf_isReal(_threshold) && !Rf_isInteger(_threshold)) || Rf_xlength(_threshold) != 1)
            verror("\"threshold\" argument must be a numeric value");

        SEXP rdim = Rf_getAttrib(_x, R_DimSymbol);

        if (!Rf_isInteger(rdim) || Rf_xlength(rdim) != 2)
            verror("\"x\" argument must be a matrix of numeric values");

        bool tidy = Rf_asLogical(_tidy);
        double threshold = fabs(Rf_asReal(_threshold));

        uint64_t num_points = Rf_nrows(_x);
        uint64_t num_dims = Rf_ncols(_x);

        if (num_points < 1 || num_dims < 1)
            verror("\"x\" argument must be a matrix of numeric values");

        uint64_t num_vals = num_points * num_dims;
        vector<double> vals;

        vals.reserve(num_vals);

        for (uint64_t i = 0; i < num_vals; ++i) {
            if ((Rf_isReal(_x) && !R_FINITE(REAL(_x)[i])) || (Rf_isInteger(_x) && INTEGER(_x)[i] == NA_INTEGER))
                vals.push_back(numeric_limits<double>::quiet_NaN());
            else
                vals.push_back(Rf_isReal(_x) ? REAL(_x)[i] : INTEGER(_x)[i]);
        }

        // the distances are stored in the order of the 'dist' object: (2,1), (3,1), ..., (N,1), (3,2), ...
        vdebug("Allocating shared memory for results\n");
        uint64_t res_size = num_points * (num_points - 1) / 2;
        shm_sizeof = sizeof(double) * max(res_size, (uint64_t)1);
        shm = mmap(NULL, shm_sizeof, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);

        if (shm == (double *)MAP_FAILED)
            verror("Failed to allocate shared memory: %s", strerror(errno));

        double *res = (double *)shm;

        int num_processes = (int)min(num_points / 2, (uint64_t)g_tgstat->num_processes());

        ProgressReporter progress;
        progress.init(num_points * num_points / 2 - num_points, 1);

        vdebug("num_processes: %d\n", num_processes);
        TGStat::prepare4multitasking();

        for (int iprocess = 0; iprocess < num_processes; ++iprocess) {
            if (!TGStat::launch_process()) {     // child process
                // two parts per process, one from each end, to balance the triangle; integer bounds so that
                // the parts cover all the rows: floor(i * N / (2P)) and N - ceil(i * N / (2P))
                uint64_t num_parts = 2 * num_processes;
                uint64_t srow[2] = { iprocess * num_points / num_parts, num_points - ((iprocess + 1) * num_points + num_parts - 1) / num_parts };
                uint64_t erow[2] = { (iprocess + 1) * num_points / num_parts, num_points - (iprocess * num_points + num_parts - 1) / num_parts };
                uint64_t itr_idx = 0;

                for (int ipart = 0; ipart < 2; ipart++) {
                    for (uint64_t irow1 = 0; irow1 < num_points; ++irow1) {
                        uint64_t start = max(srow[ipart], irow1 + 1);
                        if (start >= erow[ipart])
                            break;

                        uint64_t idx = irow1 * num_points - irow1 * (irow1 + 1) / 2 + start - irow1 - 1;

                        for (uint64_t irow2 = start; irow2 < erow[ipart]; ++irow2) {
                            double dist = 0;
                            double dev;
                            uint64_t count = 0;
                            uint64_t idx1 = irow1;
                            uint64_t idx2 = irow2;

                            for (uint64_t icol = 0; icol < num_dims; ++icol) {
                                if (!std::isnan(vals[idx1]) && !std::isnan(vals[idx2])) {
                                    dev = vals[idx1] - vals[idx2];
                                    dist += dev * dev;
                                    ++count;
                                }
                                idx1 += num_points;
                                idx2 += num_points;
                            }

                            if (count) {
                                if (count == num_dims)
                                    res[idx] = sqrt(dist);
                                else
                                    res[idx] = sqrt(dist * num_dims / (double)count);
                            } else
                                res[idx] = NA_REAL;

                            ++idx;
                            TGStat::itr_idx(++itr_idx);
                        }
                    }
                }
                rexit();
            }
        }

        while (TGStat::wait_for_kids(3000))
            progress.report(TGStat::itr_idx_sum() - progress.get_elapsed_steps());

        progress.report_last();

        // assemble the answer
        if (tidy)
            answer = tidy_dist(res, num_points, threshold, _rrownames);
        else {
            rprotect(answer = RSaneAllocVector(REALSXP, res_size));

            // copy the distances in chunks and free the pages of each copied chunk, so that the distances are not
            // held twice in memory. (munmap of a part of a shared mapping would not free its pages.)
            const uint64_t chunk_size = ((uint64_t)64 << 20) / sizeof(double);

            for (uint64_t i = 0; i < res_size; i += chunk_size) {
                uint64_t n = min(chunk_size, res_size - i);
                memcpy(REAL(answer) + i, res + i, n * sizeof(double));
#ifdef MADV_REMOVE
                madvise((char *)(res + i), n * sizeof(double), MADV_REMOVE);
#endif
            }

            set_attrs(answer, _attrs);
        }
    } catch (TGLException &e) {
        if (!TGStat::is_kid() && shm != (double *)MAP_FAILED) {
            munmap((char *)shm, shm_sizeof);     // needs to be char * for some versions of Solaris
            shm = (double *)MAP_FAILED;
        }
		rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        if (!TGStat::is_kid() && shm != (double *)MAP_FAILED) {
            munmap((char *)shm, shm_sizeof);     // needs to be char * for some versions of Solaris
            shm = (double *)MAP_FAILED;
        }
        rerror("Out of memory");
    }

    if (!TGStat::is_kid() && shm != (double *)MAP_FAILED) {
        munmap((char *)shm, shm_sizeof);     // needs to be char * for some versions of Solaris
        shm = (double *)MAP_FAILED;
    }
	rreturn(answer);
}

SEXP tgs_dist_blas(SEXP _x, SEXP _attrs, SEXP _tidy, SEXP _threshold, SEXP _rrownames, SEXP _envir)
{
    SEXP answer = R_NilValue;

	try {
        struct Mem {
            double *m{NULL};
            double *m2{NULL};
            double *mask{NULL};
            double *n{NULL};
            double *res{NULL};
            double *d{NULL};
            ~Mem() { free(m); free(m2); free(mask); free(n); free(res); free(d); }
        } mem;

        TGStat tgstat(_envir);

		if ((!Rf_isReal(_x) && !Rf_isInteger(_x)) || Rf_xlength(_x) < 1)
			verror("\"x\" argument must be a matrix of numeric values");

        if (!Rf_isLogical(_tidy) || Rf_xlength(_tidy) != 1)
            verror("\"tidy\" argument must be a logical value");

        if ((!Rf_isReal(_threshold) && !Rf_isInteger(_threshold)) || Rf_xlength(_threshold) != 1)
            verror("\"threshold\" argument must be a numeric value");

        SEXP rdim = Rf_getAttrib(_x, R_DimSymbol);

        if (!Rf_isInteger(rdim) || Rf_xlength(rdim) != 2)
            verror("\"x\" argument must be a matrix of numeric values");

        bool tidy = Rf_asLogical(_tidy);
        double threshold = fabs(Rf_asReal(_threshold));

        uint64_t num_points = Rf_nrows(_x);
        uint64_t num_dims = Rf_ncols(_x);
        int num_points32 = (int)num_points;
        int num_dims32 = (int)num_dims;

        if (num_points < 1 || num_dims < 1)
            verror("\"x\" argument must be a matrix of numeric values");

        uint64_t num_vals = num_points * num_dims;
        uint64_t res_size = num_points * (num_points - 1) / 2;
        bool nan_in_vals = false;

        // some BLAS implementations ask to align double arrays to 64 for improved efficiency
        if (posix_memalign((void **)&mem.m, 64, sizeof(double) * num_vals))
            verror("%s", strerror(errno));

        if (posix_memalign((void **)&mem.mask, 64, sizeof(double) * num_vals))
            verror("%s", strerror(errno));

        for (uint64_t i = 0; i < num_vals; ++i) {
            if ((Rf_isReal(_x) && !R_FINITE(REAL(_x)[i])) || (Rf_isInteger(_x) && INTEGER(_x)[i] == NA_INTEGER)) {
                mem.m[i] = mem.mask[i] = 0.;
                nan_in_vals = true;
            } else {
                double val = Rf_isReal(_x) ? REAL(_x)[i] : INTEGER(_x)[i];
                mem.m[i] = val;
                mem.mask[i] = 1.;
            }
        }

        // distance without nans is calculated as:
        //     sqrt(s_x^2 + s_y^2 - 2 * (m %*% t(m)))
        // distance with nans is calculated as:
        //     sqrt((m^2 %*% t(mask) + mask %*% t(m^2) - 2 * (m %*% t(m))) * num_dims / (mask %*% t(mask)))
        //
        // The lower triangle of these N x N matrices is computed in blocks of columns, and the distances are
        // written in the order of the 'dist' object directly into the answer (or into a buffer if tidy).

        vector<double> s_v2 = row_sums_square(mem.m, num_points, num_dims);
        double tol = SQ_DIST_TOL * num_dims * DBL_EPSILON;

        if (nan_in_vals) {
            if (posix_memalign((void **)&mem.m2, 64, sizeof(double) * num_vals))
                verror("%s", strerror(errno));

            for (uint64_t i = 0; i < num_vals; ++i)
                mem.m2[i] = mem.m[i] * mem.m[i];
        }

        double *d;

        if (tidy) {
            if (posix_memalign((void **)&mem.d, 64, sizeof(double) * max(res_size, (uint64_t)1)))
                verror("%s", strerror(errno));
            d = mem.d;
        } else {
            rprotect(answer = RSaneAllocVector(REALSXP, res_size));
            d = REAL(answer);
        }

        // block of at most 2^23 values (64 MB)
        uint64_t block_cols = min(num_points, max((uint64_t)1, ((uint64_t)1 << 23) / num_points));

        if (posix_memalign((void **)&mem.res, 64, sizeof(double) * num_points * block_cols))
            verror("%s", strerror(errno));

        if (nan_in_vals && posix_memalign((void **)&mem.n, 64, sizeof(double) * num_points * block_cols))
            verror("%s", strerror(errno));

        ProgressReporter progress;
        progress.init(num_points, 1);

        for (uint64_t col1 = 0, idx = 0; col1 + 1 < num_points; col1 += block_cols) {
            // rows [col1, num_points) x columns [col1, col2)
            uint64_t col2 = min(col1 + block_cols, num_points);
            int num_rows32 = (int)(num_points - col1);
            int num_cols32 = (int)(col2 - col1);
            char transa = 'N';
            char transb = 'T';
            double alpha = -2;
            double beta = 0;

            // res <- -2 * (m * t(m))
            F77_NAME(dgemm)(&transa, &transb, &num_rows32, &num_cols32, &num_dims32, &alpha, mem.m + col1, &num_points32, mem.m + col1, &num_points32, &beta, mem.res, &num_rows32 FCONE FCONE);

            if (nan_in_vals) {
                // res <- num_dims * m^2 %*% t(mask) + num_dims * mask %*% t(m^2) + num_dims * res
                alpha = num_dims;
                beta = num_dims;
                F77_NAME(dgemm)(&transa, &transb, &num_rows32, &num_cols32, &num_dims32, &alpha, mem.m2 + col1, &num_points32, mem.mask + col1, &num_points32, &beta, mem.res, &num_rows32 FCONE FCONE);
                beta = 1;
                F77_NAME(dgemm)(&transa, &transb, &num_rows32, &num_cols32, &num_dims32, &alpha, mem.mask + col1, &num_points32, mem.m2 + col1, &num_points32, &beta, mem.res, &num_rows32 FCONE FCONE);

                // n <- mask %*% t(mask)
                alpha = 1;
                beta = 0;
                F77_NAME(dgemm)(&transa, &transb, &num_rows32, &num_cols32, &num_dims32, &alpha, mem.mask + col1, &num_points32, mem.mask + col1, &num_points32, &beta, mem.n, &num_rows32 FCONE FCONE);

                // sqrt(res / n), computed directly if it is within the rounding error (see SQ_DIST_TOL); NA if the points
                // share no dimension (as without BLAS)
                for (uint64_t ipoint1 = col1; ipoint1 < col2; ++ipoint1) {
                    uint64_t offset = (ipoint1 - col1) * (num_points - col1) - col1;
                    for (uint64_t ipoint2 = ipoint1 + 1; ipoint2 < num_points; ++ipoint2) {
                        double n = mem.n[offset + ipoint2];
                        double sq_dist = mem.res[offset + ipoint2];     // times num_dims

                        if (n && sq_dist < num_dims * tol * (s_v2[ipoint1] + s_v2[ipoint2]))
                            sq_dist = num_dims * direct_sq_dist(mem.m + ipoint1, mem.mask + ipoint1, num_points, mem.m + ipoint2, mem.mask + ipoint2, num_points, num_dims);
                        d[idx++] = n ? sqrt(max(sq_dist / n, 0.)) : NA_REAL;
                    }
                }
            } else {
                // sqrt(res + s_x^2 + s_y^2)
                for (uint64_t ipoint1 = col1; ipoint1 < col2; ++ipoint1) {
                    uint64_t offset = (ipoint1 - col1) * (num_points - col1) - col1;
                    for (uint64_t ipoint2 = ipoint1 + 1; ipoint2 < num_points; ++ipoint2) {
                        double s_v2_sum = s_v2[ipoint1] + s_v2[ipoint2];
                        double sq_dist = mem.res[offset + ipoint2] + s_v2_sum;

                        if (sq_dist < tol * s_v2_sum)
                            sq_dist = direct_sq_dist(mem.m + ipoint1, NULL, num_points, mem.m + ipoint2, NULL, num_points, num_dims);
                        d[idx++] = sqrt(max(sq_dist, 0.));
                    }
                }
            }

            check_interrupt();
            progress.report(col2 - col1);
        }
        progress.report_last();

        // assemble the answer
        if (tidy)
            answer = tidy_dist(d, num_points, threshold, _rrownames);
        else
            set_attrs(answer, _attrs);
    } catch (TGLException &e) {
		rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

    rreturn(answer);
}

// The distances between the rows of x and the rows of y
SEXP tgs_cross_dist(SEXP _x, SEXP _y, SEXP _tidy, SEXP _threshold, SEXP _envir)
{
    SEXP answer = R_NilValue;
    double *res = (double *)MAP_FAILED;
    uint64_t res_sizeof = 0;

    try {
        TGStat tgstat(_envir);

        if (Rf_isNull(_y))
            verror("\"y\" argument must be a matrix of numeric values");

        check_points(_x, _y);

        if (!Rf_isLogical(_tidy) || Rf_xlength(_tidy) != 1)
            verror("\"tidy\" argument must be a logical value");

        if ((!Rf_isReal(_threshold) && !Rf_isInteger(_threshold)) || Rf_xlength(_threshold) != 1)
            verror("\"threshold\" argument must be a numeric value");

        bool tidy = Rf_asLogical(_tidy);
        double threshold = fabs(Rf_asReal(_threshold));
        uint64_t num_points[2] = { (uint64_t)Rf_nrows(_x), (uint64_t)Rf_nrows(_y) };
        uint64_t num_dims = Rf_ncols(_x);
        uint64_t res_size = num_points[0] * num_points[1];
        vector<double> vals[2];
        vector<bool> nan_in_point[2];

        read_points(_x, vals[0], nan_in_point[0]);
        read_points(_y, vals[1], nan_in_point[1]);

        // the distances are stored as the num_points[0] x num_points[1] answer matrix
        vdebug("Allocating shared memory for results\n");
        res_sizeof = sizeof(double) * res_size;
        res = (double *)mmap(NULL, res_sizeof, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);

        if (res == (double *)MAP_FAILED)
            verror("Failed to allocate shared memory: %s", strerror(errno));

        int num_processes = (int)min(num_points[0], (uint64_t)g_tgstat->num_processes());

        ProgressReporter progress;
        progress.init(res_size, 1);

        vdebug("num_processes: %d\n", num_processes);
        TGStat::prepare4multitasking();

        for (int iprocess = 0; iprocess < num_processes; ++iprocess) {
            if (!TGStat::launch_process()) {     // child process
                // rows [spoint, epoint) of x
                uint64_t spoint = iprocess * num_points[0] / num_processes;
                uint64_t epoint = (iprocess + 1) * num_points[0] / num_processes;

                for (uint64_t ipoint2 = 0; ipoint2 < num_points[1]; ++ipoint2) {
                    const double *p2 = &vals[1].front() + ipoint2 * num_dims;
                    double *dists = res + ipoint2 * num_points[0];

                    for (uint64_t ipoint1 = spoint; ipoint1 < epoint; ++ipoint1)
                        dists[ipoint1] = point_dist(&vals[0].front() + ipoint1 * num_dims, p2, num_dims, nan_in_point[0][ipoint1] || nan_in_point[1][ipoint2]);

                    TGStat::itr_idx((ipoint2 + 1) * (epoint - spoint));
                }
                rexit();
            }
        }

        while (TGStat::wait_for_kids(3000))
            progress.report(TGStat::itr_idx_sum() - progress.get_elapsed_steps());

        progress.report_last();

        // assemble the answer
        SEXP rrownames[2] = { get_rownames(_x), get_rownames(_y) };

        if (tidy)
            answer = tidy_cross_dist(res, num_points[0], num_points[1], threshold, rrownames[0], rrownames[1]);
        else {
            rprotect(answer = RSaneAllocVector(REALSXP, res_size));

            // copy the distances in chunks and free the pages of each copied chunk, so that the distances are not
            // held twice in memory. (munmap of a part of a shared mapping would not free its pages.)
            const uint64_t chunk_size = ((uint64_t)64 << 20) / sizeof(double);

            for (uint64_t i = 0; i < res_size; i += chunk_size) {
                uint64_t n = min(chunk_size, res_size - i);
                memcpy(REAL(answer) + i, res + i, n * sizeof(double));
#ifdef MADV_REMOVE
                madvise((char *)(res + i), n * sizeof(double), MADV_REMOVE);
#endif
            }

            set_cross_dist_attrs(answer, num_points[0], num_points[1], rrownames[0], rrownames[1]);
        }
    } catch (TGLException &e) {
        if (!TGStat::is_kid() && res != (double *)MAP_FAILED) {
            munmap((char *)res, res_sizeof);     // needs to be char * for some versions of Solaris
            res = (double *)MAP_FAILED;
        }
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        if (!TGStat::is_kid() && res != (double *)MAP_FAILED) {
            munmap((char *)res, res_sizeof);     // needs to be char * for some versions of Solaris
            res = (double *)MAP_FAILED;
        }
        rerror("Out of memory");
    }

    if (!TGStat::is_kid() && res != (double *)MAP_FAILED) {
        munmap((char *)res, res_sizeof);     // needs to be char * for some versions of Solaris
        res = (double *)MAP_FAILED;
    }
    rreturn(answer);
}

SEXP tgs_cross_dist_blas(SEXP _x, SEXP _y, SEXP _tidy, SEXP _threshold, SEXP _envir)
{
    SEXP answer = R_NilValue;

    try {
        struct Mem {
            double *m[2]{NULL, NULL};
            double *m2[2]{NULL, NULL};
            double *mask[2]{NULL, NULL};
            double *n{NULL};
            double *d{NULL};
            ~Mem() {
                for (int k = 0; k < 2; ++k) {
                    free(m[k]);
                    free(m2[k]);
                    free(mask[k]);
                }
                free(n);
                free(d);
            }
        } mem;

        TGStat tgstat(_envir);

        if (Rf_isNull(_y))
            verror("\"y\" argument must be a matrix of numeric values");

        check_points(_x, _y);

        if (!Rf_isLogical(_tidy) || Rf_xlength(_tidy) != 1)
            verror("\"tidy\" argument must be a logical value");

        if ((!Rf_isReal(_threshold) && !Rf_isInteger(_threshold)) || Rf_xlength(_threshold) != 1)
            verror("\"threshold\" argument must be a numeric value");

        bool tidy = Rf_asLogical(_tidy);
        double threshold = fabs(Rf_asReal(_threshold));
        SEXP rmat[2] = { _x, _y };
        uint64_t num_points[2] = { (uint64_t)Rf_nrows(_x), (uint64_t)Rf_nrows(_y) };
        uint64_t num_dims = Rf_ncols(_x);
        int num_points32[2] = { (int)num_points[0], (int)num_points[1] };
        int num_dims32 = (int)num_dims;
        uint64_t res_size = num_points[0] * num_points[1];
        bool nan_in_vals = false;

        for (int k = 0; k < 2; ++k) {
            if (read_blas_points(rmat[k], mem.m[k]))
                nan_in_vals = true;
        }

        // with missing values: m and mask are 0 where a value is missing, mask is 1 elsewhere
        if (nan_in_vals) {
            for (int k = 0; k < 2; ++k) {
                uint64_t num_vals = num_points[k] * num_dims;

                if (posix_memalign((void **)&mem.mask[k], 64, sizeof(double) * num_vals) ||
                    posix_memalign((void **)&mem.m2[k], 64, sizeof(double) * num_vals))
                    verror("%s", strerror(errno));

                for (uint64_t i = 0; i < num_vals; ++i) {
                    if (std::isnan(mem.m[k][i]))
                        mem.m[k][i] = mem.mask[k][i] = 0.;
                    else
                        mem.mask[k][i] = 1.;
                    mem.m2[k][i] = mem.m[k][i] * mem.m[k][i];
                }
            }
        }

        vector<double> s_v2[2] = { row_sums_square(mem.m[0], num_points[0], num_dims), row_sums_square(mem.m[1], num_points[1], num_dims) };

        // distance without missing values is calculated as:
        //     sqrt(s_x^2 + s_y^2 - 2 * (mx %*% t(my)))
        // distance with missing values is calculated as:
        //     sqrt((mx^2 %*% t(masky) + maskx %*% t(my^2) - 2 * (mx %*% t(my))) * num_dims / (maskx %*% t(masky)))
        //
        // The columns of the answer (the rows of y) are computed in blocks, directly in the answer unless it is tidy.
        double *d;

        if (tidy) {
            if (posix_memalign((void **)&mem.d, 64, sizeof(double) * res_size))
                verror("%s", strerror(errno));
            d = mem.d;
        } else {
            rprotect(answer = RSaneAllocVector(REALSXP, res_size));
            d = REAL(answer);
        }

        // block of at most 2^23 values (64 MB)
        uint64_t block_cols = min(num_points[1], max((uint64_t)1, ((uint64_t)1 << 23) / num_points[0]));

        if (nan_in_vals && posix_memalign((void **)&mem.n, 64, sizeof(double) * num_points[0] * block_cols))
            verror("%s", strerror(errno));

        double tol = SQ_DIST_TOL * num_dims * DBL_EPSILON;
        ProgressReporter progress;
        progress.init(num_points[1], 1);

        for (uint64_t col1 = 0; col1 < num_points[1]; col1 += block_cols) {
            uint64_t col2 = min(col1 + block_cols, num_points[1]);
            int num_cols32 = (int)(col2 - col1);
            double *res = d + col1 * num_points[0];
            char transa = 'N';
            char transb = 'T';
            double alpha = -2;
            double beta = 0;

            // res <- -2 * mx %*% t(my[col1:(col2 - 1), ])
            F77_NAME(dgemm)(&transa, &transb, &num_points32[0], &num_cols32, &num_dims32, &alpha, mem.m[0], &num_points32[0], mem.m[1] + col1, &num_points32[1], &beta, res, &num_points32[0] FCONE FCONE);

            if (nan_in_vals) {
                // res <- res + mx^2 %*% t(masky) + maskx %*% t(my^2)
                alpha = 1;
                beta = 1;
                F77_NAME(dgemm)(&transa, &transb, &num_points32[0], &num_cols32, &num_dims32, &alpha, mem.m2[0], &num_points32[0], mem.mask[1] + col1, &num_points32[1], &beta, res, &num_points32[0] FCONE FCONE);
                F77_NAME(dgemm)(&transa, &transb, &num_points32[0], &num_cols32, &num_dims32, &alpha, mem.mask[0], &num_points32[0], mem.m2[1] + col1, &num_points32[1], &beta, res, &num_points32[0] FCONE FCONE);

                // n <- maskx %*% t(masky)
                beta = 0;
                F77_NAME(dgemm)(&transa, &transb, &num_points32[0], &num_cols32, &num_dims32, &alpha, mem.mask[0], &num_points32[0], mem.mask[1] + col1, &num_points32[1], &beta, mem.n, &num_points32[0] FCONE FCONE);
            }

            for (uint64_t ipoint2 = col1, idx = 0; ipoint2 < col2; ++ipoint2) {
                for (uint64_t ipoint1 = 0; ipoint1 < num_points[0]; ++ipoint1, ++idx) {
                    double n = nan_in_vals ? mem.n[idx] : num_dims;
                    double s_v2_sum = s_v2[0][ipoint1] + s_v2[1][ipoint2];
                    double sq_dist = nan_in_vals ? res[idx] : res[idx] + s_v2_sum;

                    // NA if the points share no dimension, as without BLAS
                    if (!n) {
                        res[idx] = NA_REAL;
                        continue;
                    }

                    if (sq_dist < tol * s_v2_sum)
                        sq_dist = direct_sq_dist(mem.m[0] + ipoint1, nan_in_vals ? mem.mask[0] + ipoint1 : NULL, num_points[0],
                                                 mem.m[1] + ipoint2, nan_in_vals ? mem.mask[1] + ipoint2 : NULL, num_points[1], num_dims);

                    sq_dist = max(sq_dist, 0.);
                    res[idx] = n == num_dims ? sqrt(sq_dist) : sqrt(sq_dist * num_dims / n);
                }
            }

            check_interrupt();
            progress.report(col2 - col1);
        }
        progress.report_last();

        // assemble the answer
        SEXP rrownames[2] = { get_rownames(_x), get_rownames(_y) };

        if (tidy)
            answer = tidy_cross_dist(d, num_points[0], num_points[1], threshold, rrownames[0], rrownames[1]);
        else
            set_cross_dist_attrs(answer, num_points[0], num_points[1], rrownames[0], rrownames[1]);
    } catch (TGLException &e) {
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

    rreturn(answer);
}

// The knn nearest rows of y to each row of x (of x itself, excluding the row, if y is NULL)
SEXP tgs_dist_knn(SEXP _x, SEXP _y, SEXP _knn, SEXP _threshold, SEXP _envir)
{
    SEXP answer = R_NilValue;
    char *res = (char *)MAP_FAILED;
    uint64_t res_sizeof = 0;

    try {
        TGStat tgstat(_envir);

        check_points(_x, _y);

        bool cross = !Rf_isNull(_y);
        uint64_t knn = check_knn(_knn, _threshold);
        double threshold = fabs(Rf_asReal(_threshold));
        uint64_t num_points[2] = { (uint64_t)Rf_nrows(_x), (uint64_t)Rf_nrows(cross ? _y : _x) };
        uint64_t num_dims = Rf_ncols(_x);

        if (!cross && num_points[0] < 2)
            verror("\"x\" argument must have at least two rows");

        knn = min(knn, cross ? num_points[1] : num_points[0] - 1);

        vector<double> vals[2];
        vector<bool> nan_in_point[2];

        read_points(_x, vals[0], nan_in_point[0]);
        if (cross)
            read_points(_y, vals[1], nan_in_point[1]);

        const double *vals2 = &vals[cross ? 1 : 0].front();
        const vector<bool> &nan_in_point2 = nan_in_point[cross ? 1 : 0];

        // shared memory: the neighbours of each row of x, nearest first, followed by the number of neighbours of each row
        vdebug("Allocating shared memory for results\n");
        res_sizeof = (sizeof(Neighbour) * knn + sizeof(uint64_t)) * num_points[0];
        res = (char *)mmap(NULL, res_sizeof, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);

        if (res == (char *)MAP_FAILED)
            verror("Failed to allocate shared memory: %s", strerror(errno));

        Neighbour *neighbours = (Neighbour *)res;
        uint64_t *num_neighbours = (uint64_t *)(res + sizeof(Neighbour) * knn * num_points[0]);
        int num_processes = (int)min(num_points[0], (uint64_t)g_tgstat->num_processes());

        ProgressReporter progress;
        progress.init(num_points[0] * num_points[1], 1);

        vdebug("num_processes: %d\n", num_processes);
        TGStat::prepare4multitasking();

        for (int iprocess = 0; iprocess < num_processes; ++iprocess) {
            if (!TGStat::launch_process()) {     // child process
                // rows [spoint, epoint) of x
                uint64_t spoint = iprocess * num_points[0] / num_processes;
                uint64_t epoint = (iprocess + 1) * num_points[0] / num_processes;
                vector<Neighbour> heap;

                heap.reserve(knn);

                for (uint64_t ipoint1 = spoint; ipoint1 < epoint; ++ipoint1) {
                    const double *p1 = &vals[0].front() + ipoint1 * num_dims;

                    heap.clear();
                    for (uint64_t ipoint2 = 0; ipoint2 < num_points[1]; ++ipoint2) {
                        if (!cross && ipoint2 == ipoint1)
                            continue;

                        Neighbour neighbour = { point_dist(p1, vals2 + ipoint2 * num_dims, num_dims, nan_in_point[0][ipoint1] || nan_in_point2[ipoint2]), (int)ipoint2 };

                        // NA distances fail the comparison with the threshold
                        if (neighbour.dist <= threshold && (heap.size() < knn || nearer(neighbour, heap.front())))
                            add_neighbour(heap, knn, neighbour);
                    }

                    sort_heap(heap.begin(), heap.end(), nearer);
                    copy(heap.begin(), heap.end(), neighbours + ipoint1 * knn);
                    num_neighbours[ipoint1] = heap.size();
                    TGStat::itr_idx((ipoint1 - spoint + 1) * num_points[1]);
                }
                rexit();
            }
        }

        while (TGStat::wait_for_kids(3000))
            progress.report(TGStat::itr_idx_sum() - progress.get_elapsed_steps());

        progress.report_last();

        answer = knn_answer(neighbours, num_neighbours, num_points[0], knn, get_rownames(_x), get_rownames(cross ? _y : _x));
    } catch (TGLException &e) {
        if (!TGStat::is_kid() && res != (char *)MAP_FAILED) {
            munmap(res, res_sizeof);
            res = (char *)MAP_FAILED;
        }
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        if (!TGStat::is_kid() && res != (char *)MAP_FAILED) {
            munmap(res, res_sizeof);
            res = (char *)MAP_FAILED;
        }
        rerror("Out of memory");
    }

    if (!TGStat::is_kid() && res != (char *)MAP_FAILED) {
        munmap(res, res_sizeof);
        res = (char *)MAP_FAILED;
    }
    rreturn(answer);
}

// The same as tgs_dist_knn for data without missing values. The squared distances of a block of rows of x to all the
// rows of y are computed as s_x^2 + s_y^2 - 2 * (x %*% t(y)), with one matrix product.
SEXP tgs_dist_knn_blas(SEXP _x, SEXP _y, SEXP _knn, SEXP _threshold, SEXP _envir)
{
    SEXP answer = R_NilValue;

    try {
        struct Mem {
            double *m[2]{NULL, NULL};
            double *res{NULL};
            ~Mem() { free(m[0]); free(m[1]); free(res); }
        } mem;

        TGStat tgstat(_envir);

        check_points(_x, _y);

        bool cross = !Rf_isNull(_y);
        int num_mats = cross ? 2 : 1;
        SEXP rmat[2] = { _x, _y };
        uint64_t knn = check_knn(_knn, _threshold);
        double threshold = fabs(Rf_asReal(_threshold));
        uint64_t num_points[2] = { (uint64_t)Rf_nrows(_x), (uint64_t)Rf_nrows(cross ? _y : _x) };
        uint64_t num_dims = Rf_ncols(_x);
        int num_points32[2] = { (int)num_points[0], (int)num_points[1] };
        int num_dims32 = (int)num_dims;

        if (!cross && num_points[0] < 2)
            verror("\"x\" argument must have at least two rows");

        knn = min(knn, cross ? num_points[1] : num_points[0] - 1);

        for (int k = 0; k < num_mats; ++k) {
            if (read_blas_points(rmat[k], mem.m[k]))
                verror("BLAS implementation of tgs_dist_knn does not support missing values");
        }

        const double *mx = mem.m[0];
        const double *my = mem.m[num_mats - 1];
        vector<double> s_x = row_sums_square(mx, num_points[0], num_dims);
        vector<double> s_y = cross ? row_sums_square(my, num_points[1], num_dims) : s_x;
        double tol = SQ_DIST_TOL * num_dims * DBL_EPSILON;

        // block of at most 2^23 values (64 MB)
        uint64_t block_rows = min(num_points[0], max((uint64_t)1, ((uint64_t)1 << 23) / num_points[1]));

        if (posix_memalign((void **)&mem.res, 64, sizeof(double) * num_points[1] * block_rows))
            verror("%s", strerror(errno));

        vector<Neighbour> heap;
        vector<Neighbour> neighbours(num_points[0] * knn);
        vector<uint64_t> num_neighbours(num_points[0], 0);

        heap.reserve(knn);

        ProgressReporter progress;
        progress.init(num_points[0], 1);

        for (uint64_t row1 = 0; row1 < num_points[0]; row1 += block_rows) {
            // res <- -2 * my %*% t(mx[row1:(row2 - 1), ])
            uint64_t row2 = min(row1 + block_rows, num_points[0]);
            int num_rows32 = (int)(row2 - row1);
            char transa = 'N';
            char transb = 'T';
            double alpha = -2;
            double beta = 0;

            F77_NAME(dgemm)(&transa, &transb, &num_points32[1], &num_rows32, &num_dims32, &alpha, my, &num_points32[1], mx + row1, &num_points32[0], &beta, mem.res, &num_points32[1] FCONE FCONE);

            for (uint64_t ipoint1 = row1; ipoint1 < row2; ++ipoint1) {
                const double *res = mem.res + (ipoint1 - row1) * num_points[1];
                const double *s_y2 = &s_y.front();
                double s_x1 = s_x[ipoint1];

                // the squared distance to row ipoint2 of y
                auto sq_dist = [&](uint64_t ipoint2) {
                    double s_v2_sum = s_x1 + s_y2[ipoint2];
                    double sq_dist = res[ipoint2] + s_v2_sum;

                    if (sq_dist < tol * s_v2_sum)
                        sq_dist = direct_sq_dist(mx + ipoint1, NULL, num_points[0], my + ipoint2, NULL, num_points[1], num_dims);
                    return max(sq_dist, 0.);
                };

                // the squared distance of the farthest neighbour once the heap is full. The rows of y are visited in
                // increasing order, so a row that is not nearer than the farthest neighbour is not a neighbour.
                double bound = numeric_limits<double>::infinity();

                heap.clear();
                for (uint64_t ipoint2 = 0; ipoint2 < num_points[1]; ++ipoint2) {
                    // most rows are not nearer than the farthest neighbour: when the formula is above its rounding
                    // error, it is the squared distance and is compared to the bound directly
                    double s_v2_sum = s_x1 + s_y2[ipoint2];
                    double sq = res[ipoint2] + s_v2_sum;

                    if (sq >= bound && sq >= tol * s_v2_sum)
                        continue;

                    if (!cross && ipoint2 == ipoint1)
                        continue;

                    sq = sq_dist(ipoint2);

                    if (!(sq < bound))
                        continue;

                    Neighbour neighbour = { sqrt(sq), (int)ipoint2 };

                    if (neighbour.dist <= threshold && (heap.size() < knn || nearer(neighbour, heap.front()))) {
                        add_neighbour(heap, knn, neighbour);
                        if (heap.size() == knn)
                            bound = sq_dist(heap.front().row);
                    }
                }

                sort_heap(heap.begin(), heap.end(), nearer);
                copy(heap.begin(), heap.end(), neighbours.begin() + ipoint1 * knn);
                num_neighbours[ipoint1] = heap.size();
            }

            check_interrupt();
            progress.report(row2 - row1);
        }
        progress.report_last();

        answer = knn_answer(&neighbours.front(), &num_neighbours.front(), num_points[0], knn, get_rownames(_x), get_rownames(cross ? _y : _x));
    } catch (TGLException &e) {
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

    rreturn(answer);
}

}

