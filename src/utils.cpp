#include <algorithm>
#include <cmath>
#include <cstdint>
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

        if (num_rows < 2)
            rreturn(R_NilValue);

        rprotect(answer = RSaneAllocVector(REALSXP, num_rows * num_cols));
        double *res = REAL(answer);
        const double *xreal = Rf_isReal(_x) ? REAL(_x) : NULL;
        const int *xint = Rf_isReal(_x) ? NULL : INTEGER(_x);
        vector<pair<double, uint64_t>> sorted(spearman ? num_rows : 0);

        for (uint64_t icol = 0; icol < num_cols; ++icol) {
            double *v = res + icol * num_rows;

            for (uint64_t irow = 0, i = icol * num_rows; irow < num_rows; ++irow, ++i) {
                if (xreal ? !R_FINITE(xreal[i]) : xint[i] == NA_INTEGER)
                    rreturn(R_NilValue);
                v[irow] = xreal ? xreal[i] : xint[i];
            }

            if (spearman) {
                // average ranks for ties, as in tgs_cor
                for (uint64_t irow = 0; irow < num_rows; ++irow)
                    sorted[irow] = { v[irow], irow };
                sort(sorted.begin(), sorted.end());
                for (uint64_t start = 0, end; start < num_rows; start = end) {
                    for (end = start + 1; end < num_rows && sorted[end].first == sorted[start].first; ++end)
                        ;
                    double rank = (start + end + 1) / 2.;
                    for (uint64_t k = start; k < end; ++k)
                        v[sorted[k].second] = rank;
                }
            }

            double mean = 0;
            for (uint64_t irow = 0; irow < num_rows; ++irow)
                mean += v[irow];
            mean /= num_rows;

            double sum_sq = 0;
            for (uint64_t irow = 0; irow < num_rows; ++irow) {
                v[irow] -= mean;
                sum_sq += v[irow] * v[irow];
            }

            // constant column (or overflow)
            if (!(sum_sq > 0) || !R_FINITE(sum_sq))
                rreturn(R_NilValue);

            double scale = 1 / sqrt(sum_sq);
            for (uint64_t irow = 0; irow < num_rows; ++irow)
                v[irow] *= scale;
        }

        Rf_setAttrib(answer, R_DimSymbol, rdim);
        Rf_setAttrib(answer, R_DimNamesSymbol, Rf_getAttrib(_x, R_DimNamesSymbol));
    } catch (TGLException &e) {
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

    rreturn(answer);
}

}
