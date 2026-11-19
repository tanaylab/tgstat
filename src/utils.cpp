#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <errno.h>
#include <sys/mman.h>
#include <vector>
#ifndef R_NO_REMAP
#  define R_NO_REMAP
#endif
#include <R.h>
#include <Rinternals.h>

#ifdef length
#undef length
#endif
#ifdef error
#undef error
#endif

#include "tgstat.h"

void flag_constant_cols(const double *vals, uint64_t num_rows, uint64_t num_cols, vector<bool> &flags)
{
    for (uint64_t icol = 0; icol < num_cols; ++icol) {
        const double *col = vals + icol * num_rows;

        if (!flags[icol] && all_of(col, col + num_rows, [col](double val) { return val == col[0]; }))
            flags[icol] = true;
    }
}

vector<bool> constant_cols(SEXP x)
{
    uint64_t num_rows = Rf_nrows(x);
    uint64_t num_cols = Rf_ncols(x);
    const double *xreal = Rf_isReal(x) ? REAL(x) : NULL;
    const int *xint = Rf_isReal(x) ? NULL : INTEGER(x);
    vector<bool> res(num_cols, true);

    for (uint64_t icol = 0, i = 0; icol < num_cols; ++icol) {
        bool found = false;
        double first = 0;

        for (uint64_t irow = 0; irow < num_rows; ++irow, ++i) {
            if (xreal ? !R_FINITE(xreal[i]) : xint[i] == NA_INTEGER)
                continue;

            double val = xreal ? xreal[i] : xint[i];

            if (!found) {
                found = true;
                first = val;
            } else if (val != first)
                res[icol] = false;
        }
    }
    return res;
}

bool constant_in_pair(SEXP x, uint64_t xcol, SEXP y, uint64_t ycol)
{
    uint64_t num_rows = Rf_nrows(x);
    SEXP mats[2] = { x, y };
    uint64_t cols[2] = { xcol, ycol };
    const double *reals[2];
    const int *ints[2];
    bool found = false;
    bool constant[2] = { true, true };
    double first[2] = { 0, 0 };

    for (int k = 0; k < 2; ++k) {
        reals[k] = Rf_isReal(mats[k]) ? REAL(mats[k]) + cols[k] * num_rows : NULL;
        ints[k] = Rf_isReal(mats[k]) ? NULL : INTEGER(mats[k]) + cols[k] * num_rows;
    }

    for (uint64_t irow = 0; irow < num_rows; ++irow) {
        double vals[2];

        for (int k = 0; k < 2; ++k)
            vals[k] = reals[k] ? reals[k][irow] : (ints[k][irow] == NA_INTEGER ? NA_REAL : ints[k][irow]);

        if (!R_FINITE(vals[0]) || !R_FINITE(vals[1]))
            continue;

        if (!found) {
            found = true;
            first[0] = vals[0];
            first[1] = vals[1];
        } else {
            constant[0] = constant[0] && vals[0] == first[0];
            constant[1] = constant[1] && vals[1] == first[1];
            if (!constant[0] && !constant[1])
                return false;
        }
    }
    return found && (constant[0] || constant[1]);
}

bool standardize(double *v, uint64_t n, bool spearman, vector<pair<double, uint64_t>> &sorted)
{
    // constant values are checked first: after centering they are not exactly zero if their mean is rounded
    if (all_of(v, v + n, [v](double val) { return val == v[0]; }))
        return false;

    if (spearman) {
        // average ranks for ties, as in tgs_cor
        sorted.resize(n);
        for (uint64_t i = 0; i < n; ++i)
            sorted[i] = { v[i], i };
        sort(sorted.begin(), sorted.end());
        for (uint64_t start = 0, end; start < n; start = end) {
            for (end = start + 1; end < n && sorted[end].first == sorted[start].first; ++end)
                ;
            double rank = (start + end + 1) / 2.;
            for (uint64_t k = start; k < end; ++k)
                v[sorted[k].second] = rank;
        }
    }

    double mean = 0;
    for (uint64_t i = 0; i < n; ++i)
        mean += v[i];
    mean /= n;

    double sum_sq = 0;
    for (uint64_t i = 0; i < n; ++i) {
        v[i] -= mean;
        sum_sq += v[i] * v[i];
    }

    if (!(sum_sq > 0) || !R_FINITE(sum_sq))
        return false;

    double scale = 1 / sqrt(sum_sq);
    for (uint64_t i = 0; i < n; ++i)
        v[i] *= scale;
    return true;
}

bool all_finite(SEXP x)
{
    uint64_t len = Rf_xlength(x);

    if (Rf_isReal(x)) {
        const double *vals = REAL(x);
        return all_of(vals, vals + len, [](double val) { return R_FINITE(val); });
    }

    const int *vals = INTEGER(x);
    return all_of(vals, vals + len, [](int val) { return val != NA_INTEGER; });
}

void standardize_cols(SEXP rmat[], int num_mats, bool spearman, double *z[], uint64_t z_sizeof[], vector<bool> flagged[])
{
    uint64_t num_rows = Rf_nrows(rmat[0]);
    uint64_t num_cols[2] = { 0, 0 };
    uint64_t num_vals = 0;

    for (int k = 0; k < num_mats; ++k) {
        num_cols[k] = Rf_ncols(rmat[k]);
        z_sizeof[k] = sizeof(double) * max(num_rows * num_cols[k], (uint64_t)1);
        z[k] = (double *)mmap(NULL, z_sizeof[k], PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);

        if (z[k] == (double *)MAP_FAILED)
            verror("Failed to allocate shared memory: %s", strerror(errno));

        num_vals += num_rows * num_cols[k];
    }

    // one byte per column, set by the process that standardizes the column
    struct Flags {
        unsigned char *p{(unsigned char *)MAP_FAILED};
        uint64_t size{0};
        ~Flags() { if (!TGStat::is_kid() && p != (unsigned char *)MAP_FAILED) munmap((char *)p, size); }
    } flags;

    flags.size = max(num_cols[0] + num_cols[1], (uint64_t)1);
    flags.p = (unsigned char *)mmap(NULL, flags.size, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);

    if (flags.p == (unsigned char *)MAP_FAILED)
        verror("Failed to allocate shared memory: %s", strerror(errno));

    auto standardize_part = [&](uint64_t ipart, uint64_t num_parts) {
        vector<pair<double, uint64_t>> sorted;

        for (int k = 0; k < num_mats; ++k) {
            const double *xreal = Rf_isReal(rmat[k]) ? REAL(rmat[k]) : NULL;
            const int *xint = Rf_isReal(rmat[k]) ? NULL : INTEGER(rmat[k]);
            unsigned char *col_flags = flags.p + (k ? num_cols[0] : 0);
            uint64_t scol = ipart * num_cols[k] / num_parts;
            uint64_t ecol = (ipart + 1) * num_cols[k] / num_parts;

            for (uint64_t icol = scol; icol < ecol; ++icol) {
                double *v = z[k] + icol * num_rows;
                uint64_t offset = icol * num_rows;
                bool finite = true;

                for (uint64_t irow = 0; irow < num_rows && finite; ++irow) {
                    if (xreal ? !R_FINITE(xreal[offset + irow]) : xint[offset + irow] == NA_INTEGER)
                        finite = false;
                    else
                        v[irow] = xreal ? xreal[offset + irow] : xint[offset + irow];
                }

                if (!finite || !standardize(v, num_rows, spearman, sorted)) {
                    memset(v, 0, sizeof(double) * num_rows);
                    col_flags[icol] = 1;
                }
            }
        }
    };

    // a process per 500K values to rank or 4M values to center and scale: forking takes longer than a smaller part
    int num_processes = (int)min(num_vals / (spearman ? 500000 : 4000000), (uint64_t)g_tgstat->num_processes());

    if (num_processes > 1) {
        TGStat::prepare4multitasking();

        for (int iprocess = 0; iprocess < num_processes; ++iprocess) {
            if (!TGStat::launch_process()) {     // child process
                standardize_part(iprocess, num_processes);
                rexit();
            }
        }

        while (TGStat::wait_for_kids(3000))
            ;
    } else
        standardize_part(0, 1);

    for (int k = 0; k < num_mats; ++k) {
        const unsigned char *col_flags = flags.p + (k ? num_cols[0] : 0);
        flagged[k].assign(col_flags, col_flags + num_cols[k]);
    }
}

extern "C" {

SEXP tgs_finite(SEXP _x, SEXP _envir)
{
	try {
        TGStat tgstat(_envir);
        uint64_t len = Rf_xlength(_x);

        if (!Rf_isReal(_x) && !Rf_isInteger(_x))
            verror("\"x\" argument must be numeric or integer");

        for (uint64_t i = 0; i < len; ++i) {
            if ((Rf_isReal(_x) && !R_FINITE(REAL(_x)[i])) || (Rf_isInteger(_x) && INTEGER(_x)[i] == NA_INTEGER))
                rreturn(Rf_ScalarLogical(false));
        }
    } catch (TGLException &e) {
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

    rreturn(Rf_ScalarLogical(true));
}

// Replaces the values of each column of x by their ranks (if spearman is true), then centers each column and
// scales it to unit norm, so that crossprod() of the result is the correlation matrix of x.
// Returns NULL if x has fewer than 2 rows, non-finite values or constant columns.
SEXP tgs_standardize(SEXP _x, SEXP _spearman, SEXP _envir)
{
    SEXP answer = R_NilValue;
    double *shm = (double *)MAP_FAILED;
    uint64_t shm_sizeof = 0;

	try {
        TGStat tgstat(_envir);

		if ((!Rf_isReal(_x) && !Rf_isInteger(_x)) || Rf_xlength(_x) < 1)
			verror("\"x\" argument must be a matrix of numeric values");

        SEXP rdim = Rf_getAttrib(_x, R_DimSymbol);

        if (!Rf_isInteger(rdim) || Rf_xlength(rdim) != 2)
            verror("\"x\" argument must be a matrix of numeric values");

        if (!Rf_isLogical(_spearman) || Rf_xlength(_spearman) != 1)
            verror("\"spearman\" argument must be a logical value");

        bool spearman = Rf_asLogical(_spearman);
        uint64_t num_rows = Rf_nrows(_x);
        uint64_t num_cols = Rf_ncols(_x);
        uint64_t num_vals = num_rows * num_cols;

        if (num_rows < 2)
            rreturn(R_NilValue);

        const double *xreal = Rf_isReal(_x) ? REAL(_x) : NULL;
        const int *xint = Rf_isReal(_x) ? NULL : INTEGER(_x);

        rprotect(answer = RSaneAllocVector(REALSXP, num_vals));
        double *res = REAL(answer);

        for (uint64_t icol = 0, i = 0; icol < num_cols; ++icol) {
            bool constant = true;

            for (uint64_t irow = 0; irow < num_rows; ++irow, ++i) {
                if (xreal ? !R_FINITE(xreal[i]) : xint[i] == NA_INTEGER)
                    rreturn(R_NilValue);
                res[i] = xreal ? xreal[i] : xint[i];
                constant = constant && res[i] == res[icol * num_rows];
            }

            if (constant)
                rreturn(R_NilValue);
        }

        // the ranking is done in parallel by child processes that write into shared memory, with at least
        // 500K values per process, since forking takes longer than ranking a smaller part
        int num_processes = spearman ? (int)min(num_vals / 500000, (uint64_t)g_tgstat->num_processes()) : 1;

        if (num_processes > 1) {
            shm_sizeof = sizeof(double) * num_vals;
            shm = (double *)mmap(NULL, shm_sizeof, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);

            if (shm == (double *)MAP_FAILED)
                verror("Failed to allocate shared memory: %s", strerror(errno));

            TGStat::prepare4multitasking();

            for (int iprocess = 0; iprocess < num_processes; ++iprocess) {
                if (!TGStat::launch_process()) {     // child process
                    vector<pair<double, uint64_t>> sorted;
                    uint64_t scol = iprocess * num_cols / num_processes;
                    uint64_t ecol = (iprocess + 1) * num_cols / num_processes;

                    memcpy(shm + scol * num_rows, res + scol * num_rows, sizeof(double) * (ecol - scol) * num_rows);
                    for (uint64_t icol = scol; icol < ecol; ++icol)
                        standardize(shm + icol * num_rows, num_rows, true, sorted);
                    rexit();
                }
            }

            while (TGStat::wait_for_kids(3000))
                ;

            memcpy(res, shm, shm_sizeof);
        } else {
            vector<pair<double, uint64_t>> sorted;

            for (uint64_t icol = 0; icol < num_cols; ++icol) {
                if (!standardize(res + icol * num_rows, num_rows, spearman, sorted))
                    rreturn(R_NilValue);
            }
        }

        Rf_setAttrib(answer, R_DimSymbol, rdim);
        Rf_setAttrib(answer, R_DimNamesSymbol, Rf_getAttrib(_x, R_DimNamesSymbol));
    } catch (TGLException &e) {
        if (!TGStat::is_kid() && shm != (double *)MAP_FAILED) {
            munmap((char *)shm, shm_sizeof);
            shm = (double *)MAP_FAILED;
        }
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        if (!TGStat::is_kid() && shm != (double *)MAP_FAILED) {
            munmap((char *)shm, shm_sizeof);
            shm = (double *)MAP_FAILED;
        }
        rerror("Out of memory");
    }

    if (!TGStat::is_kid() && shm != (double *)MAP_FAILED)
        munmap((char *)shm, shm_sizeof);

    rreturn(answer);
}

}
