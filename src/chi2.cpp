#include <algorithm>
#include <cmath>
#include <errno.h>
#include <sys/mman.h>
#include <unistd.h>
#include <vector>
#include <cstdint>

#ifndef R_NO_REMAP
#  define R_NO_REMAP
#endif
#include <R.h>
#include <Rinternals.h>
// Note: For df=1, pchisq(x, 1, lower.tail=FALSE) = erfc(sqrt(x/2))
// This avoids per-element R API calls and is much faster.

#include "ProgressReporter.h"
#include "tgstat.h"

extern "C" {

SEXP tgs_chi2(SEXP _x, SEXP _yates, SEXP _envir)
{
    SEXP answer = R_NilValue;
    double *res = (double *)MAP_FAILED;
    uint64_t res_sizeof = 0;

    try {
        TGStat tgstat(_envir);

        if (!Rf_isLogical(_yates) || Rf_xlength(_yates) != 1)
            verror("\"yates\" argument must be a logical value");

        int yates_val = Rf_asLogical(_yates);
        if (yates_val == NA_LOGICAL)
            verror("\"yates\" argument must not be NA");

        bool yates = yates_val;

        SEXP _rdims = R_NilValue;
        SEXP _xclass = Rf_getAttrib(_x, R_ClassSymbol);
        SEXP _xx = Rf_getAttrib(_x, Rf_install("x"));   // non zero values of sparse matrix
        SEXP _xi = Rf_getAttrib(_x, Rf_install("i"));   // row number within sparse matrix
        SEXP _xp = Rf_getAttrib(_x, Rf_install("p"));   // index offset of the column within sparse matrix

        if ((!Rf_isReal(_x) && !Rf_isInteger(_x) && !Rf_isObject(_x)) ||
            ((Rf_isReal(_x) || Rf_isInteger(_x)) && (Rf_xlength(_x) < 1 || !Rf_isInteger(_rdims = Rf_getAttrib(_x, R_DimSymbol)) || Rf_xlength(_rdims) != 2)) ||
            (Rf_isObject(_x) && (!Rf_isString(_xclass) || Rf_xlength(_xclass) < 1 || strcmp(CHAR(STRING_ELT(_xclass, 0)), "dgCMatrix") ||
                              !Rf_isInteger(_rdims = Rf_getAttrib(_x, Rf_install("Dim"))) || Rf_xlength(_rdims) != 2 ||
                              (Rf_xlength(_xx) > 0 && !Rf_isInteger(_xx) && !Rf_isReal(_xx)) ||
                              !Rf_isInteger(_xi) || Rf_xlength(_xi) != Rf_xlength(_xx) ||
                              !Rf_isInteger(_xp) || Rf_xlength(_xp) != INTEGER(_rdims)[1] + 1)))
            verror("\"x\" argument must be a matrix of numeric values");

        uint64_t num_rows = INTEGER(_rdims)[0];
        uint64_t num_cols = INTEGER(_rdims)[1];

        if (num_cols != 2)
            verror("\"x\" argument must have exactly 2 columns");

        if (num_rows < 1)
            verror("\"x\" argument must have at least 1 row");

        bool is_sparse = _xclass != R_NilValue;
        bool is_real = is_sparse ? (Rf_xlength(_xx) > 0 ? Rf_isReal(_xx) : true) : Rf_isReal(_x);

        SEXP _xdimnames = is_sparse ? Rf_getAttrib(_x, Rf_install("Dimnames")) : Rf_getAttrib(_x, R_DimNamesSymbol);

        // Compute column sums
        double colSum0 = 0.0;
        double colSum1 = 0.0;

        if (!is_sparse) {
            if (is_real) {
                double *vals = REAL(_x);
                for (uint64_t i = 0; i < num_rows; ++i) {
                    double v0 = vals[i];
                    double v1 = vals[num_rows + i];
                    if (std::isnan(v0) || std::isnan(v1))
                        continue;
                    colSum0 += v0;
                    colSum1 += v1;
                }
            } else {
                int *vals = INTEGER(_x);
                for (uint64_t i = 0; i < num_rows; ++i) {
                    if (vals[i] == NA_INTEGER || vals[num_rows + i] == NA_INTEGER)
                        continue;
                    colSum0 += vals[i];
                    colSum1 += vals[num_rows + i];
                }
            }
        } else {
            int *xp = INTEGER(_xp);
            if (is_real) {
                double *xx = REAL(_xx);
                for (int j = xp[0]; j < xp[1]; ++j)
                    if (!std::isnan(xx[j])) colSum0 += xx[j];
                for (int j = xp[1]; j < xp[2]; ++j)
                    if (!std::isnan(xx[j])) colSum1 += xx[j];
            } else {
                int *xx = INTEGER(_xx);
                for (int j = xp[0]; j < xp[1]; ++j)
                    colSum0 += xx[j];
                for (int j = xp[1]; j < xp[2]; ++j)
                    colSum1 += xx[j];
            }
        }

        // Allocate shared memory for results: num_rows chi2 values + num_rows pval values
        res_sizeof = sizeof(double) * 2 * num_rows;
        res = (double *)mmap(NULL, res_sizeof, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);

        if (res == (double *)MAP_FAILED)
            verror("Failed to allocate shared memory: %s", strerror(errno));

        int num_processes;
        if (num_rows < 500) {
            num_processes = 0;  // no forking for small matrices
        } else {
            num_processes = (int)min((uint64_t)5, min(num_rows / 10, (uint64_t)(g_tgstat->num_processes() / 2)));
            if (num_rows > 0)
                num_processes = max(num_processes, 1);
        }

        vdebug("num_processes: %d\n", num_processes);

        ProgressReporter progress;
        progress.init(num_rows, 1);

        if (num_processes == 0) {
            // Compute directly in parent process without forking
            if (!is_sparse) {
                for (uint64_t irow = 0; irow < num_rows; ++irow) {
                    double a, c_val;

                    if (is_real) {
                        a = REAL(_x)[irow];
                        c_val = REAL(_x)[num_rows + irow];
                    } else {
                        int ia = INTEGER(_x)[irow];
                        int ic = INTEGER(_x)[num_rows + irow];
                        if (ia == NA_INTEGER || ic == NA_INTEGER) {
                            res[irow] = NA_REAL;
                            res[num_rows + irow] = NA_REAL;
                            continue;
                        }
                        a = (double)ia;
                        c_val = (double)ic;
                    }

                    if (is_real && (std::isnan(a) || std::isnan(c_val))) {
                        res[irow] = NA_REAL;
                        res[num_rows + irow] = NA_REAL;
                        continue;
                    }

                    double b = colSum0 - a;
                    double d = colSum1 - c_val;
                    double N = colSum0 + colSum1;
                    double r1 = colSum0;
                    double r2 = colSum1;
                    double c1 = a + c_val;
                    double c2 = b + d;
                    double denom = (double)r1 * r2 * c1 * c2;

                    if (denom == 0.0) {
                        res[irow] = 0.0;
                        res[num_rows + irow] = 1.0;
                    } else if (yates) {
                        double ad_bc = fabs(a * d - b * c_val);
                        double chi2_val = ad_bc > N / 2.0 ? (ad_bc - N / 2.0) * (ad_bc - N / 2.0) * N / denom : 0.0;
                        res[irow] = chi2_val;
                        res[num_rows + irow] = erfc(sqrt(chi2_val / 2.0));
                    } else {
                        double diff = a * d - b * c_val;
                        double chi2_val = diff * diff * N / denom;
                        res[irow] = chi2_val;
                        res[num_rows + irow] = erfc(sqrt(chi2_val / 2.0));
                    }
                }
            } else {
                // Sparse (dgCMatrix) path in parent
                int *xi = INTEGER(_xi);
                int *xp = INTEGER(_xp);

                int off0 = xp[0];
                int end0 = xp[1];
                int off1 = xp[1];
                int end1 = xp[2];

                for (uint64_t irow = 0; irow < num_rows; ++irow) {
                    double a = 0.0;
                    double c_val = 0.0;

                    if (off0 < end0 && xi[off0] == (int)irow) {
                        a = is_real ? REAL(_xx)[off0] : (double)INTEGER(_xx)[off0];
                        ++off0;
                    }
                    if (off1 < end1 && xi[off1] == (int)irow) {
                        c_val = is_real ? REAL(_xx)[off1] : (double)INTEGER(_xx)[off1];
                        ++off1;
                    }

                    if (std::isnan(a) || std::isnan(c_val)) {
                        res[irow] = NA_REAL;
                        res[num_rows + irow] = NA_REAL;
                        continue;
                    }

                    double b = colSum0 - a;
                    double d = colSum1 - c_val;
                    double N = colSum0 + colSum1;
                    double r1 = colSum0;
                    double r2 = colSum1;
                    double c1 = a + c_val;
                    double c2 = b + d;
                    double denom = (double)r1 * r2 * c1 * c2;

                    if (denom == 0.0) {
                        res[irow] = 0.0;
                        res[num_rows + irow] = 1.0;
                    } else if (yates) {
                        double ad_bc = fabs(a * d - b * c_val);
                        double chi2_val = ad_bc > N / 2.0 ? (ad_bc - N / 2.0) * (ad_bc - N / 2.0) * N / denom : 0.0;
                        res[irow] = chi2_val;
                        res[num_rows + irow] = erfc(sqrt(chi2_val / 2.0));
                    } else {
                        double diff = a * d - b * c_val;
                        double chi2_val = diff * diff * N / denom;
                        res[irow] = chi2_val;
                        res[num_rows + irow] = erfc(sqrt(chi2_val / 2.0));
                    }
                }
            }

            progress.report(num_rows);
            progress.report_last();
        } else {
            // Multi-process path
            TGStat::prepare4multitasking();

            for (int iprocess = 0; iprocess < num_processes; ++iprocess) {
                if (!TGStat::launch_process()) {     // child process
                    uint64_t srow = num_rows * (iprocess / (double)num_processes);
                    uint64_t erow = num_rows * ((iprocess + 1) / (double)num_processes);

                    if (!is_sparse) {
                        // Dense matrix path
                        for (uint64_t irow = srow; irow < erow; ++irow) {
                            double a, c_val;

                            if (is_real) {
                                a = REAL(_x)[irow];
                                c_val = REAL(_x)[num_rows + irow];
                            } else {
                                int ia = INTEGER(_x)[irow];
                                int ic = INTEGER(_x)[num_rows + irow];
                                if (ia == NA_INTEGER || ic == NA_INTEGER) {
                                    res[irow] = NA_REAL;
                                    res[num_rows + irow] = NA_REAL;
                                    TGStat::itr_idx(irow - srow + 1);
                                    continue;
                                }
                                a = (double)ia;
                                c_val = (double)ic;
                            }

                            if (is_real && (std::isnan(a) || std::isnan(c_val))) {
                                res[irow] = NA_REAL;
                                res[num_rows + irow] = NA_REAL;
                                TGStat::itr_idx(irow - srow + 1);
                                continue;
                            }

                            double b = colSum0 - a;
                            double d = colSum1 - c_val;
                            double N = colSum0 + colSum1;
                            double r1 = colSum0;
                            double r2 = colSum1;
                            double c1 = a + c_val;
                            double c2 = b + d;
                            double denom = (double)r1 * r2 * c1 * c2;

                            if (denom == 0.0) {
                                res[irow] = 0.0;
                                res[num_rows + irow] = 1.0;
                            } else if (yates) {
                                double ad_bc = fabs(a * d - b * c_val);
                                double chi2_val = ad_bc > N / 2.0 ? (ad_bc - N / 2.0) * (ad_bc - N / 2.0) * N / denom : 0.0;
                                res[irow] = chi2_val;
                                res[num_rows + irow] = erfc(sqrt(chi2_val / 2.0));
                            } else {
                                double diff = a * d - b * c_val;
                                double chi2_val = diff * diff * N / denom;
                                res[irow] = chi2_val;
                                res[num_rows + irow] = erfc(sqrt(chi2_val / 2.0));
                            }

                            TGStat::itr_idx(irow - srow + 1);
                        }
                    } else {
                        // Sparse (dgCMatrix) path
                        int *xi = INTEGER(_xi);
                        int *xp = INTEGER(_xp);

                        // Initialize walking pointers using lower_bound for each column
                        int off0 = (int)(lower_bound(xi + xp[0], xi + xp[1], (int)srow) - xi);
                        int end0 = xp[1];
                        int off1 = (int)(lower_bound(xi + xp[1], xi + xp[2], (int)srow) - xi);
                        int end1 = xp[2];

                        for (uint64_t irow = srow; irow < erow; ++irow) {
                            double a = 0.0;
                            double c_val = 0.0;

                            if (off0 < end0 && xi[off0] == (int)irow) {
                                a = is_real ? REAL(_xx)[off0] : (double)INTEGER(_xx)[off0];
                                ++off0;
                            }
                            if (off1 < end1 && xi[off1] == (int)irow) {
                                c_val = is_real ? REAL(_xx)[off1] : (double)INTEGER(_xx)[off1];
                                ++off1;
                            }

                            if (std::isnan(a) || std::isnan(c_val)) {
                                res[irow] = NA_REAL;
                                res[num_rows + irow] = NA_REAL;
                                TGStat::itr_idx(irow - srow + 1);
                                continue;
                            }

                            double b = colSum0 - a;
                            double d = colSum1 - c_val;
                            double N = colSum0 + colSum1;
                            double r1 = colSum0;
                            double r2 = colSum1;
                            double c1 = a + c_val;
                            double c2 = b + d;
                            double denom = (double)r1 * r2 * c1 * c2;

                            if (denom == 0.0) {
                                res[irow] = 0.0;
                                res[num_rows + irow] = 1.0;
                            } else if (yates) {
                                double ad_bc = fabs(a * d - b * c_val);
                                double chi2_val = ad_bc > N / 2.0 ? (ad_bc - N / 2.0) * (ad_bc - N / 2.0) * N / denom : 0.0;
                                res[irow] = chi2_val;
                                res[num_rows + irow] = erfc(sqrt(chi2_val / 2.0));
                            } else {
                                double diff = a * d - b * c_val;
                                double chi2_val = diff * diff * N / denom;
                                res[irow] = chi2_val;
                                res[num_rows + irow] = erfc(sqrt(chi2_val / 2.0));
                            }

                            TGStat::itr_idx(irow - srow + 1);
                        }
                    }

                    rexit();
                }
            }

            while (TGStat::wait_for_kids(3000))
                progress.report(TGStat::itr_idx_sum() - progress.get_elapsed_steps());

            progress.report_last();
        }

        // Assemble the answer: num_rows x 2 matrix
        SEXP dim, dimnames, colnames;

        rprotect(answer = RSaneAllocVector(REALSXP, num_rows * 2));
        memcpy(REAL(answer), res, res_sizeof);

        rprotect(dim = RSaneAllocVector(INTSXP, 2));
        INTEGER(dim)[0] = num_rows;
        INTEGER(dim)[1] = 2;
        Rf_setAttrib(answer, R_DimSymbol, dim);

        rprotect(dimnames = RSaneAllocVector(VECSXP, 2));

        // Row names from input
        SEXP row_names = R_NilValue;
        if (!Rf_isNull(_xdimnames) && Rf_xlength(_xdimnames) > 0)
            row_names = VECTOR_ELT(_xdimnames, 0);
        SET_VECTOR_ELT(dimnames, 0, row_names);

        // Column names: c("chi2", "pval")
        rprotect(colnames = RSaneAllocVector(STRSXP, 2));
        SET_STRING_ELT(colnames, 0, Rf_mkChar("chi2"));
        SET_STRING_ELT(colnames, 1, Rf_mkChar("pval"));
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

}
