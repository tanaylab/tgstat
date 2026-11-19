#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <errno.h>
#include <new>
#include <sys/mman.h>
#include <unistd.h>
#include <vector>

#ifndef R_NO_REMAP
#  define R_NO_REMAP
#endif
#include <R.h>
#include <Rinternals.h>

#include "ProgressReporter.h"
#include "tgstat.h"

// included last: its macros (beta, sign, ...) must not reach the standard headers
#include <Rmath.h>

enum Test { KS, WILCOX };
enum Alternative { TWO_SIDED, LESS, GREATER };

// Upper tail of the limit distribution of the two-sided Kolmogorov statistic: P(K > x). Port of K2l in R's
// src/library/stats/src/ks.c (lower = FALSE), with the tolerance that ks.test uses for the asymptotic p-value.
static double pkolmogorov_two_limit_upper(double x)
{
    const double tol = 1e-6;

    if (x <= 0)
        return 1;

    if (x < 1) {
        int k_max = (int)sqrt(2 - log(tol));
        double z = -(M_PI_2 * M_PI_4) / (x * x);
        double w = log(x);
        double s = 0;

        for (int k = 1; k < k_max; k += 2)
            s += exp(k * k * z - w);
        return 1 - s / M_1_SQRT_2PI;
    }

    double z = -2 * x * x;
    double s = -1;
    int k = 2;
    double old = 0;
    double p = 2 * exp(z);

    while (fabs(old - p) > tol) {
        old = p;
        p += 2 * s * exp(z * k * k);
        s *= -1;
        k++;
    }
    return p;
}

// Computes the statistics of column icol of x: res[icol + i * num_cols] is the i-th statistic. Both samples are
// sorted and walked together, one group of tied values at a time.
static void two_sample_test(Test test, Alternative alternative, const double *xreal, const int *xint, const int *labels,
                            uint64_t num_rows, uint64_t num_cols, uint64_t icol, vector<double> &s1, vector<double> &s2,
                            double *res)
{
    uint64_t offset = icol * num_rows;

    s1.clear();
    s2.clear();
    for (uint64_t irow = 0; irow < num_rows; ++irow) {
        if (labels[irow] == NA_LOGICAL)
            continue;

        double v;
        if (xreal)
            v = xreal[offset + irow];
        else if (xint[offset + irow] == NA_INTEGER)
            continue;
        else
            v = xint[offset + irow];

        if (std::isnan(v))
            continue;

        if (labels[irow])
            s1.push_back(v);
        else
            s2.push_back(v);
    }

    int num_stats = test == KS ? 2 : 3;
    uint64_t n1 = s1.size();
    uint64_t n2 = s2.size();

    if (!n1 || !n2) {
        for (int i = 0; i < num_stats; ++i)
            res[icol + i * num_cols] = NA_REAL;
        return;
    }

    sort(s1.begin(), s1.end());
    sort(s2.begin(), s2.end());

    uint64_t i1 = 0;
    uint64_t i2 = 0;
    int64_t max_diff = 0;       // max and min over the groups of n1 * n2 * (F1 - F2), F being the empirical CDF
    int64_t min_diff = 0;
    double w = 0;               // number of pairs (v1, v2) with v1 > v2, ties counted as 1/2
    double ties = 0;            // sum of t^3 - t over the groups of t tied values

    while (i1 < n1 || i2 < n2) {
        double v = i1 < n1 && (i2 == n2 || s1[i1] <= s2[i2]) ? s1[i1] : s2[i2];
        uint64_t c1 = 0;
        uint64_t c2 = 0;

        while (i1 < n1 && s1[i1] == v) {
            ++i1;
            ++c1;
        }
        while (i2 < n2 && s2[i2] == v) {
            ++i2;
            ++c2;
        }

        if (test == KS) {
            int64_t diff = (int64_t)(i1 * n2) - (int64_t)(i2 * n1);
            max_diff = max(max_diff, diff);
            min_diff = min(min_diff, diff);
        } else {
            double t = c1 + c2;
            w += c1 * ((i2 - c2) + 0.5 * c2);
            ties += t * t * t - t;
        }
    }

    double nx = n1;
    double ny = n2;
    double pval;

    if (test == KS) {
        // ks.test(s1, s2, alternative = alternative, exact = FALSE)
        int64_t d = alternative == TWO_SIDED ? max(max_diff, -min_diff) : alternative == GREATER ? max_diff : -min_diff;
        double D = d / (nx * ny);
        double n = nx * ny / (nx + ny);

        if (alternative == TWO_SIDED)
            pval = pkolmogorov_two_limit_upper(sqrt(n) * D);
        else
            pval = 1 - (-expm1(-2 * n * (D * D)));

        res[icol] = D;
        res[icol + num_cols] = min(1., max(0., pval));
    } else {
        // wilcox.test(s1, s2, alternative = alternative, exact = FALSE, correct = TRUE)
        double N = nx + ny;
        double z = w - nx * ny / 2;
        double sigma = sqrt((nx * ny / 12) * ((N + 1) - ties / (N * (N - 1))));
        double correction = alternative == TWO_SIDED ? (z > 0) - (z < 0) : alternative == GREATER ? 1 : -1;

        z = (z - 0.5 * correction) / sigma;
        if (alternative == LESS)
            pval = Rf_pnorm5(z, 0, 1, 1, 0);
        else if (alternative == GREATER)
            pval = Rf_pnorm5(z, 0, 1, 0, 0);
        else
            pval = 2 * min(Rf_pnorm5(z, 0, 1, 1, 0), Rf_pnorm5(z, 0, 1, 0, 0));

        res[icol] = w;
        res[icol + num_cols] = w / (nx * ny);
        res[icol + 2 * num_cols] = std::isnan(pval) ? NA_REAL : pval;
    }
}

static SEXP two_sample_tests(Test test, SEXP _x, SEXP _labels, SEXP _alternative, SEXP _envir)
{
    SEXP answer = R_NilValue;
    double *res = (double *)MAP_FAILED;
    uint64_t res_sizeof = 0;

    try {
        TGStat tgstat(_envir);

        SEXP _rdims;
        if ((!Rf_isReal(_x) && !Rf_isInteger(_x)) || Rf_xlength(_x) < 1 ||
            !Rf_isInteger(_rdims = Rf_getAttrib(_x, R_DimSymbol)) || Rf_xlength(_rdims) != 2)
            verror("\"x\" argument must be a non-empty numeric matrix");

        uint64_t num_rows = INTEGER(_rdims)[0];
        uint64_t num_cols = INTEGER(_rdims)[1];

        if (!Rf_isLogical(_labels) || (uint64_t)Rf_xlength(_labels) != num_rows)
            verror("\"labels\" argument must be a logical vector of length nrow(x)");

        if (!Rf_isString(_alternative) || Rf_xlength(_alternative) != 1)
            verror("\"alternative\" argument must be a character string");

        Alternative alternative = TWO_SIDED;
        const char *alt = CHAR(STRING_ELT(_alternative, 0));
        if (!strcmp(alt, "two.sided"))
            alternative = TWO_SIDED;
        else if (!strcmp(alt, "less"))
            alternative = LESS;
        else if (!strcmp(alt, "greater"))
            alternative = GREATER;
        else
            verror("\"alternative\" argument must be one of \"two.sided\", \"less\" or \"greater\"");

        const double *xreal = Rf_isReal(_x) ? REAL(_x) : NULL;
        const int *xint = Rf_isReal(_x) ? NULL : INTEGER(_x);
        const int *labels = LOGICAL(_labels);
        int num_stats = test == KS ? 2 : 3;

        res_sizeof = sizeof(double) * num_stats * num_cols;
        res = (double *)mmap(NULL, res_sizeof, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);

        if (res == (double *)MAP_FAILED)
            verror("Failed to allocate shared memory: %s", strerror(errno));

        // at least 500K values per process, since forking takes longer than sorting a smaller part
        int num_processes = (int)min(min(num_rows * num_cols / 500000, num_cols), (uint64_t)g_tgstat->num_processes());
        vector<double> s1, s2;

        s1.reserve(num_rows);
        s2.reserve(num_rows);

        ProgressReporter progress;
        progress.init(num_cols, 1);

        vdebug("num_processes: %d\n", num_processes);

        if (num_processes > 1) {
            TGStat::prepare4multitasking();

            for (int iprocess = 0; iprocess < num_processes; ++iprocess) {
                if (!TGStat::launch_process()) {     // child process
                    uint64_t scol = iprocess * num_cols / num_processes;
                    uint64_t ecol = (iprocess + 1) * num_cols / num_processes;

                    for (uint64_t icol = scol; icol < ecol; ++icol) {
                        two_sample_test(test, alternative, xreal, xint, labels, num_rows, num_cols, icol, s1, s2, res);
                        TGStat::itr_idx(icol - scol + 1);
                    }
                    rexit();
                }
            }

            while (TGStat::wait_for_kids(3000))
                progress.report(TGStat::itr_idx_sum() - progress.get_elapsed_steps());
        } else {
            for (uint64_t icol = 0; icol < num_cols; ++icol) {
                two_sample_test(test, alternative, xreal, xint, labels, num_rows, num_cols, icol, s1, s2, res);
                check_interrupt();
                progress.report(1);
            }
        }
        progress.report_last();

        // assemble the answer: num_cols x num_stats matrix
        static const char *STAT_NAMES[2][3] = { { "D", "pval" }, { "W", "auc", "pval" } };
        SEXP dim, dimnames, colnames;
        SEXP _xdimnames = Rf_getAttrib(_x, R_DimNamesSymbol);

        rprotect(answer = RSaneAllocVector(REALSXP, num_stats * num_cols));
        memcpy(REAL(answer), res, res_sizeof);

        rprotect(dim = RSaneAllocVector(INTSXP, 2));
        INTEGER(dim)[0] = num_cols;
        INTEGER(dim)[1] = num_stats;
        Rf_setAttrib(answer, R_DimSymbol, dim);

        rprotect(dimnames = RSaneAllocVector(VECSXP, 2));
        rprotect(colnames = RSaneAllocVector(STRSXP, num_stats));
        for (int i = 0; i < num_stats; ++i)
            SET_STRING_ELT(colnames, i, Rf_mkChar(STAT_NAMES[test][i]));
        SET_VECTOR_ELT(dimnames, 0, Rf_isNull(_xdimnames) ? R_NilValue : VECTOR_ELT(_xdimnames, 1));
        SET_VECTOR_ELT(dimnames, 1, colnames);
        Rf_setAttrib(answer, R_DimNamesSymbol, dimnames);
    } catch (TGLException &e) {
        if (!TGStat::is_kid() && res != (double *)MAP_FAILED) {
            munmap((char *)res, res_sizeof);
            res = (double *)MAP_FAILED;
        }
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

    if (!TGStat::is_kid() && res != (double *)MAP_FAILED) {
        munmap((char *)res, res_sizeof);
        res = (double *)MAP_FAILED;
    }

    rreturn(answer);
}

extern "C" {

SEXP tgs_ks(SEXP _x, SEXP _labels, SEXP _alternative, SEXP _envir)
{
    return two_sample_tests(KS, _x, _labels, _alternative, _envir);
}

SEXP tgs_wilcox(SEXP _x, SEXP _labels, SEXP _alternative, SEXP _envir)
{
    return two_sample_tests(WILCOX, _x, _labels, _alternative, _envir);
}

}
