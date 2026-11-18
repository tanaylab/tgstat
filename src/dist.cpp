#define USE_FC_LEN_T
#include <algorithm>
#include <cmath>
#include <errno.h>
#include <limits>
#include <sys/mman.h>
#include <unistd.h>

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

        vector<double> s_v2;

        if (nan_in_vals) {
            if (posix_memalign((void **)&mem.m2, 64, sizeof(double) * num_vals))
                verror("%s", strerror(errno));

            for (uint64_t i = 0; i < num_vals; ++i)
                mem.m2[i] = mem.m[i] * mem.m[i];
        } else {
            s_v2.resize(num_points, 0.);

            for (uint64_t idim = 0, idx = 0; idim < num_dims; ++idim) {
                for (uint64_t ipoint = 0; ipoint < num_points; ++ipoint) {
                    s_v2[ipoint] += mem.m[idx] * mem.m[idx];
                    ++idx;
                }
            }
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

                // sqrt(res / n)
                for (uint64_t ipoint1 = col1; ipoint1 < col2; ++ipoint1) {
                    uint64_t offset = (ipoint1 - col1) * (num_points - col1) - col1;
                    for (uint64_t ipoint2 = ipoint1 + 1; ipoint2 < num_points; ++ipoint2)
                        d[idx++] = sqrt(mem.res[offset + ipoint2] / mem.n[offset + ipoint2]);
                }
            } else {
                // sqrt(res + s_x^2 + s_y^2)
                for (uint64_t ipoint1 = col1; ipoint1 < col2; ++ipoint1) {
                    uint64_t offset = (ipoint1 - col1) * (num_points - col1) - col1;
                    for (uint64_t ipoint2 = ipoint1 + 1; ipoint2 < num_points; ++ipoint2)
                        d[idx++] = sqrt(max(mem.res[offset + ipoint2] + s_v2[ipoint1] + s_v2[ipoint2], 0.));
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

}

