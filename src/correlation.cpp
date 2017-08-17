#include <algorithm>
#include <limits>
#include <sys/mman.h>
#include <unistd.h>

#include <R.h>
#include <Rinternals.h>

#ifdef length
#undef length
#endif
#ifdef error
#undef error
#endif

#include "ProgressReporter.h"
#include "tgstat.h"

struct PSort {
    bool operator()(double *p1, double *p2) const { return *p1 < *p2 || isnan(*p2); }
};

extern "C" {

SEXP tgs_cor(SEXP _x, SEXP _pairwise_complete_obs, SEXP _spearman, SEXP _tidy, SEXP _threshold, SEXP _envir)
{
    SEXP answer = R_NilValue;
    double *res = (double *)MAP_FAILED;
    size_t res_sizeof = 0;

	try {
        TGStat tgstat(_envir);

		if (!isReal(_x) && !isInteger(_x) || Rf_length(_x) < 1)
			verror("\"x\" argument must be a matrix of numeric values");

        if (!isLogical(_pairwise_complete_obs) || Rf_length(_pairwise_complete_obs) != 1)
            verror("\"pairwise.complete.obs\" argument must be a logical value");

        if (!isLogical(_spearman) || Rf_length(_spearman) != 1)
            verror("\"spearman\" argument must be a logical value");

        if (!isLogical(_tidy) || Rf_length(_tidy) != 1)
            verror("\"tidy\" argument must be a logical value");

        if (!isReal(_threshold) && !isInteger(_threshold) || Rf_length(_threshold) != 1)
            verror("\"threshold\" argument must be a numeric value");

        SEXP rdim = getAttrib(_x, R_DimSymbol);

        if (!isInteger(rdim) || Rf_length(rdim) != 2)
            verror("\"x\" argument must be a matrix of numeric values");

        bool pairwise_complete_obs = asLogical(_pairwise_complete_obs);
        bool tidy = asLogical(_tidy);
        double threshold = fabs(asReal(_threshold));
        int num_rows = nrows(_x);
        int num_cols = ncols(_x);

        if (num_rows <= 1 || num_cols <= 1)
            verror("\"x\" argument must be a matrix of numeric values");

        size_t num_vals = num_rows * num_cols;
        vector<bool> nan_in_col(num_cols, false);
        vector<double> sums(num_cols, 0);
        vector<double> sums_square(num_cols, 0);
        vector<double> means(num_cols, 0);
        vector<double> stddevs(num_cols, 0);
        vector<double> vals;
        vector<double *> pvals;
        vector<double> col_vals1, col_vals2;
        vector<double *> col_pvals1, col_pvals2;
        vector<double *> *col_pvals[2] = { &col_pvals1, &col_pvals2 };

        vals.reserve(num_vals);
        col_vals1.reserve(num_rows);
        col_vals2.reserve(num_rows);
        col_pvals1.reserve(num_rows);
        col_pvals2.reserve(num_rows);

        for (size_t i = 0; i < num_vals; ++i) {
            if (isReal(_x) && !R_FINITE(REAL(_x)[i]) || isInteger(_x) && INTEGER(_x)[i] == NA_INTEGER) {
                nan_in_col[i / num_rows] = true;
                vals.push_back(numeric_limits<double>::quiet_NaN());
            } else
                vals.push_back(isReal(_x) ? REAL(_x)[i] : INTEGER(_x)[i]);
        }

        // replace values with ranks if spearman=T
        if (asLogical(_spearman)) {
            pvals.reserve(num_vals);
            for (size_t i = 0; i < num_vals; ++i)
                pvals.push_back(&vals[i]);

            for (int icol = 0; icol < num_cols; ++icol) {
                if (nan_in_col[icol])
                    continue;

                vector<double *>::iterator sival = pvals.begin() + icol * num_rows;
                vector<double *>::iterator eival = sival + num_rows;
                vector<double *>::iterator last_ival = sival;

                sort(sival, eival, PSort());
                for (vector<double *>::iterator ival = sival; ; ++ival) {
                    if (ival == eival || **ival != **last_ival || isnan(**ival) || isnan(**last_ival)) {
                        double rank = isnan(**last_ival) ? numeric_limits<double>::quiet_NaN() : ((ival - sival) + (last_ival - sival) - 1) / 2. + 1;

                        while (last_ival != ival) {
                            **last_ival = rank;
                            ++last_ival;
                        }

                        if (ival == eival)
                            break;
                    }
                }
            }
        }

        for (int irow = 0; irow < num_rows; ++irow) {
            size_t idx = irow;
            for (int icol = 0; icol < num_cols; ++icol) {
                if (!nan_in_col[icol]) {
                    sums[icol] += vals[idx];
                    sums_square[icol] += vals[idx] * vals[idx];
                }
                idx += num_rows;
            }
        }

        for (int icol = 0; icol < num_cols; ++icol) {
            if (!nan_in_col[icol]) {
                means[icol] = sums[icol] / num_rows;

                // we are calaculating standard deviation:
                // sqrt(sum((x-mean)^2) / N)) = sqrt(sum(x^2) / N - mean^2)
                stddevs[icol] = sqrt(sums_square[icol] / num_rows - means[icol] * means[icol]);
            }
        }

        size_t res_size = num_cols * num_cols;
        res_sizeof = sizeof(double) * res_size;
        res = (double *)mmap(NULL, res_sizeof, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);

        if (res == (double *)MAP_FAILED)
            verror("Failed to allocate shared memory: %s", strerror(errno));

        for (size_t i = 0; i < res_size; ++i)
            res[i] = numeric_limits<double>::quiet_NaN();

        int num_cores = max(1, (int)sysconf(_SC_NPROCESSORS_ONLN));
        int num_processes = min(num_cols / 2, num_cores);
        double num_cols4process = num_cols / (double)num_processes;

        ProgressReporter progress;
        progress.init(num_cols * num_cols / 2 - num_cols, 1);

        TGStat::prepare4multitasking();

        for (int iprocess = 0; iprocess < num_processes; ++iprocess) {
            if (!TGStat::launch_process()) {     // child process
                int scol[2] = { (int)(iprocess * num_cols4process / 2.), (int)(num_cols - (iprocess + 1) * num_cols4process / 2.) };
                int ecol[2] = { (int)((iprocess + 1) * num_cols4process / 2.), (int)(num_cols - iprocess * num_cols4process / 2.) };
                size_t itr_idx = 0;

                for (int ipart = 0; ipart < 2; ipart++) {
                    for (int icol1 = 0; icol1 < num_cols; ++icol1) {
                        int start = max(scol[ipart], icol1 + 1);
                        if (start >= ecol[ipart])
                            break;

                        for (int icol2 = start; icol2 < ecol[ipart]; ++icol2) {
                            size_t idx = icol1 * num_cols + icol2;

                            if ((nan_in_col[icol1] || nan_in_col[icol2]) && pairwise_complete_obs) {
                                size_t idx1 = icol1 * num_rows;
                                size_t idx2 = icol2 * num_rows;
                                double sum1 = 0;
                                double sum2 = 0;
                                double sum_square1 = 0;
                                double sum_square2 = 0;
                                double mean1, mean2, stddev1, stddev2;

                                col_vals1.clear();
                                col_vals2.clear();

                                // copy non NaN pair values to col_vals
                                for (int irow = 0; irow < num_rows; ++irow) {
                                    if (!isnan(vals[irow + idx1]) && !isnan(vals[irow + idx2])) {
                                        double val1 = vals[irow + idx1];
                                        double val2 = vals[irow + idx2];

                                        col_vals1.push_back(val1);
                                        col_vals2.push_back(val2);
                                    }
                                }

                                size_t num_col_vals = col_vals1.size();

                                if (num_col_vals) {
                                    // replace values with ranking
                                    if (asLogical(_spearman)) {
                                        col_pvals1.clear();
                                        col_pvals2.clear();

                                        for (size_t i = 0; i < num_col_vals; ++i) {
                                            col_pvals1.push_back(&col_vals1[i]);
                                            col_pvals2.push_back(&col_vals2[i]);
                                        }

                                        for (int i = 0; i < 2; ++i) {
                                            vector<double *>::iterator sival = col_pvals[i]->begin();
                                            vector<double *>::iterator eival = sival + num_col_vals;
                                            vector<double *>::iterator last_ival = sival;

                                            sort(sival, eival, PSort());
                                            for (vector<double *>::iterator ival = sival; ; ++ival) {
                                                if (ival == eival || **ival != **last_ival) {
                                                    double rank = ((ival - sival) + (last_ival - sival) - 1) / 2. + 1;
                                                    while (last_ival != ival) {
                                                        **last_ival = rank;
                                                        ++last_ival;
                                                    }

                                                    if (ival == eival)
                                                        break;
                                                }
                                            }
                                        }
                                    }

                                    for (int i = 0; i < num_col_vals; ++i) {
                                        double val1 = col_vals1[i];
                                        double val2 = col_vals2[i];

                                        sum1 += val1;
                                        sum2 += val2;
                                        sum_square1 += val1 * val1;
                                        sum_square2 += val2 * val2;
                                    }

                                    mean1 = sum1 / num_col_vals;
                                    mean2 = sum2 / num_col_vals;
                                    stddev1 = sqrt(sum_square1 / num_col_vals - mean1 * mean1);
                                    stddev2 = sqrt(sum_square2 / num_col_vals - mean2 * mean2);

                                    // calculate correlation
                                    res[idx] = 0;
                                    for (int i = 0; i < num_col_vals; ++i)
                                        res[idx] += col_vals1[i] * col_vals2[i];  // => sum(X*Y)

                                    res[idx] = res[idx] / num_col_vals;           // => mean(X*Y)
                                    res[idx] -= mean1 * mean2;                    // => covariance(X,Y)
                                    res[idx] /= stddev1 * stddev2;                // => correlation(X,Y)
                                }
                            } else if (!nan_in_col[icol1] && !nan_in_col[icol2]) {
                                size_t idx1 = icol1 * num_rows;
                                size_t idx2 = icol2 * num_rows;
                                size_t end_idx1 = idx1 + num_rows;

                                res[idx] = 0;
                                while (idx1 < end_idx1)
                                    res[idx] += vals[idx1++] * vals[idx2++];  // => sum(X*Y)

                                res[idx] = res[idx] / num_rows;               // => mean(X*Y)
                                res[idx] -= means[icol1] * means[icol2];      // => covariance(X,Y)
                                res[idx] /= stddevs[icol1] * stddevs[icol2];    // => correlation(X,Y)
                            }
                            ++itr_idx;
                            TGStat::itr_idx(itr_idx);
                        }
                    }
                }
                exit(0);
            }
        }

        while (TGStat::wait_for_kids(3000))
            progress.report(TGStat::itr_idx_sum() - progress.get_elapsed_steps());

        progress.report_last();

        // assemble the answer
        SEXP rold_dimnames = getAttrib(_x, R_DimNamesSymbol);
        SEXP rold_colnames = !isNull(rold_dimnames) && Rf_length(rold_dimnames) == 2 ? VECTOR_ELT(rold_dimnames, 1) : R_NilValue;

        if (tidy) {
            enum { COL1, COL2, COR, NUM_COLS };
            const char *COL_NAMES[NUM_COLS] = { "col1", "col2", "cor" };

            rprotect(answer = allocVector(VECSXP, NUM_COLS));

            size_t answer_size = 0;

            for (int icol1 = 0; icol1 < num_cols; ++icol1) {
                size_t idx = icol1;
                for (int icol2 = 0; icol2 < icol1; ++icol2) {
                    if (!isnan(res[idx]) && fabs(res[idx]) >= threshold)
                        ++answer_size;
                    idx += num_cols;
                }
            }

            SEXP rcol1, rcol2, rcor, rrownames, rcolnames;

            SET_VECTOR_ELT(answer, COL1, (rcol1 = allocVector(INTSXP, answer_size)));
            SET_VECTOR_ELT(answer, COL2, (rcol2 = allocVector(INTSXP, answer_size)));
            SET_VECTOR_ELT(answer, COR, (rcor = allocVector(REALSXP, answer_size)));

            if (rold_colnames != R_NilValue) {
                setAttrib(rcol1, R_LevelsSymbol, rold_colnames);
                setAttrib(rcol1, R_ClassSymbol, mkString("factor"));
                setAttrib(rcol2, R_LevelsSymbol, rold_colnames);
                setAttrib(rcol2, R_ClassSymbol, mkString("factor"));
            }

            setAttrib(answer, R_NamesSymbol, (rcolnames = allocVector(STRSXP, NUM_COLS)));
            setAttrib(answer, R_ClassSymbol, mkString("data.frame"));
            setAttrib(answer, R_RowNamesSymbol, (rrownames = allocVector(INTSXP, answer_size)));

            for (int i = 0; i < NUM_COLS; i++)
                SET_STRING_ELT(rcolnames, i, mkChar(COL_NAMES[i]));

            int i = 0;
            for (int icol1 = 0; icol1 < num_cols; ++icol1) {
                for (int icol2 = icol1 + 1; icol2 < num_cols; ++icol2) {
                    size_t idx = icol1 * num_cols + icol2;
                    if (!isnan(res[idx]) && fabs(res[idx]) >= threshold) {
                        INTEGER(rcol1)[i] = icol1 + 1;
                        INTEGER(rcol2)[i] = icol2 + 1;
                        REAL(rcor)[i] = res[idx];
                        INTEGER(rrownames)[i] = i + 1;
                        ++i;
                    }
                }
            }
        } else {
            // copy the matrix below the diagonal to the upper part (the two parts are identical since cor(X,Y)=cor(Y,X)
            for (int icol1 = 0; icol1 < num_cols; ++icol1) {
                size_t idx1 = icol1 * num_cols;
                size_t idx2 = icol1;
                for (int icol2 = 0; icol2 < icol1; ++icol2) {
                    res[idx1] = res[idx2];
                    idx1++;
                    idx2 += num_cols;
                }
                res[icol1 * (num_cols + 1)] = 1.;
            }

            for (size_t i = 0; i < res_size; ++i) {
                if (isnan(res[i]))
                    res[i] = NA_REAL;
            }

            SEXP dim;

            rprotect(answer = allocVector(REALSXP, res_size));
            memcpy(REAL(answer), res, res_sizeof);

            rprotect(dim = allocVector(INTSXP, 2));
            INTEGER(dim)[0] = num_cols;
            INTEGER(dim)[1] = num_cols;
            setAttrib(answer, R_DimSymbol, dim);

            if (rold_colnames != R_NilValue) {
                SEXP dimnames;
                rprotect(dimnames = allocVector(VECSXP, 2));
                SET_VECTOR_ELT(dimnames, 0, rold_colnames);
                SET_VECTOR_ELT(dimnames, 1, rold_colnames);
                setAttrib(answer, R_DimNamesSymbol, dimnames);
            }
        }

    } catch (TGLException &e) {
        if (!TGStat::is_kid() && res != (double *)MAP_FAILED) {
            munmap(res, res_sizeof);
            res = (double *)MAP_FAILED;
        }
		rerror("%s", e.msg());
	}

    if (!TGStat::is_kid() && res != (double *)MAP_FAILED) {
        munmap(res, res_sizeof);
        res = (double *)MAP_FAILED;
    }
	rreturn(answer);
}

}


