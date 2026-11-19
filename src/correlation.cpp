#define USE_FC_LEN_T
#include <algorithm>
#include <cmath>
#include <errno.h>
#include <limits>
#include <queue>
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
# define FCONE
#endif

#include "ProgressReporter.h"
#include "tgstat.h"

// Unmaps the shared memory of standardize_cols.
static void unmap_shared(double *z[2], uint64_t z_sizeof[2])
{
    for (int k = 0; k < 2; ++k) {
        if (!TGStat::is_kid() && z[k] != (double *)MAP_FAILED) {
            munmap((char *)z[k], z_sizeof[k]);  // needs to be char * for some versions of Solaris
            z[k] = (double *)MAP_FAILED;
        }
    }
}

extern "C" {

SEXP tgs_cor(SEXP _x, SEXP _pairwise_complete_obs, SEXP _spearman, SEXP _tidy, SEXP _threshold, SEXP _envir)
{
    SEXP answer = R_NilValue;
    double *res = (double *)MAP_FAILED;
    uint64_t res_sizeof = 0;

	try {
        TGStat tgstat(_envir);

		if ((!Rf_isReal(_x) && !Rf_isInteger(_x)) || Rf_xlength(_x) < 1)
			verror("\"x\" argument must be a matrix of numeric values");

        if (!Rf_isLogical(_pairwise_complete_obs) || Rf_xlength(_pairwise_complete_obs) != 1)
            verror("\"pairwise.complete.obs\" argument must be a logical value");

        if (!Rf_isLogical(_spearman) || Rf_xlength(_spearman) != 1)
            verror("\"spearman\" argument must be a logical value");

        if (!Rf_isLogical(_tidy) || Rf_xlength(_tidy) != 1)
            verror("\"tidy\" argument must be a logical value");

        if ((!Rf_isReal(_threshold) && !Rf_isInteger(_threshold)) || Rf_xlength(_threshold) != 1)
            verror("\"threshold\" argument must be a numeric value");

        SEXP rdim = Rf_getAttrib(_x, R_DimSymbol);

        if (!Rf_isInteger(rdim) || Rf_xlength(rdim) != 2)
            verror("\"x\" argument must be a matrix of numeric values");

        bool pairwise_complete_obs = Rf_asLogical(_pairwise_complete_obs);
        bool spearman = Rf_asLogical(_spearman);
        bool tidy = Rf_asLogical(_tidy);
        double threshold = fabs(Rf_asReal(_threshold));
        uint64_t num_rows = Rf_nrows(_x);
        uint64_t num_cols = Rf_ncols(_x);

        if (num_rows < 1 || num_cols < 1)
            verror("\"x\" argument must be a matrix of numeric values");

        uint64_t num_vals = num_rows * num_cols;
        bool nan_in_vals = false;
        vector<bool> nan_in_col(num_cols, false);
        vector<double> sums(num_cols, 0);
        vector<double> sums_square(num_cols, 0);
        vector<double> means(num_cols, 0);
        vector<double> stddevs(num_cols, 0);
        vector<double> vals;
        vector<double *> pvals;
        vector<double> col_vals1(num_rows);
        vector<double> col_vals2(num_rows);
        double *pcol_vals1 = &col_vals1.front();
        double *pcol_vals2 = &col_vals2.front();
        vector<double> *col_vals[2] = { &col_vals1, &col_vals2 };
        vector<double> nan_col(num_rows, numeric_limits<double>::quiet_NaN());

        vals.reserve(num_vals);

        for (uint64_t i = 0; i < num_vals; ++i) {
            if ((Rf_isReal(_x) && !R_FINITE(REAL(_x)[i])) || (Rf_isInteger(_x) && INTEGER(_x)[i] == NA_INTEGER)) {
                nan_in_col[i / num_rows] = true;
                nan_in_vals = true;
                vals.push_back(numeric_limits<double>::quiet_NaN());
            } else
                vals.push_back(Rf_isReal(_x) ? REAL(_x)[i] : INTEGER(_x)[i]);
        }

        // the correlations of a constant column are NA, like those of a column with missing values
        flag_constant_cols(&vals.front(), num_rows, num_cols, nan_in_col);

        // with pairwise.complete.obs and missing values: the columns that are constant in every pair
        vector<bool> constant_col;

        if (nan_in_vals && pairwise_complete_obs)
            constant_col = constant_cols(_x);

        // replace values with ranks if spearman=T
        if (spearman) {
            pvals.reserve(num_vals);
            for (uint64_t i = 0; i < num_vals; ++i)
                pvals.push_back(&vals[i]);

            for (uint64_t icol = 0; icol < num_cols; ++icol) {
                if (nan_in_col[icol] && !pairwise_complete_obs)
                    continue;

                vector<double *>::iterator sival = pvals.begin() + icol * num_rows;
                vector<double *>::iterator eival = sival + num_rows;
                vector<double *>::iterator last_ival = sival;

                if (nan_in_col[icol])
                    sort(sival, eival, [](double *p1, double *p2) { return *p1 < *p2 || (!std::isnan(*p1) && std::isnan(*p2)); });
                else
                    sort(sival, eival, [](double *p1, double *p2) { return *p1 < *p2; });

                if (!nan_in_vals || !pairwise_complete_obs) {
                    for (auto ival = sival; ; ++ival) {
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
        }

        for (uint64_t irow = 0; irow < num_rows; ++irow) {
            uint64_t idx = irow;
            for (uint64_t icol = 0; icol < num_cols; ++icol) {
                if (!nan_in_col[icol]) {
                    sums[icol] += vals[idx];
                    sums_square[icol] += vals[idx] * vals[idx];
                }
                idx += num_rows;
            }
        }

        for (uint64_t icol = 0; icol < num_cols; ++icol) {
            if (!nan_in_col[icol]) {
                means[icol] = sums[icol] / num_rows;

                // we are calaculating standard deviation:
                // sqrt(sum((x-mean)^2) / N)) = sqrt(sum(x^2) / N - mean^2)
                stddevs[icol] = sqrt(sums_square[icol] / num_rows - means[icol] * means[icol]);
            }
        }

        vdebug("Allocating shared memory for results\n");
        uint64_t res_size = num_cols * num_cols;
        res_sizeof = sizeof(double) * res_size;
        res = (double *)mmap(NULL, res_sizeof, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);

        if (res == (double *)MAP_FAILED)
            verror("Failed to allocate shared memory: %s", strerror(errno));

        for (uint64_t i = 0; i < res_size; ++i)
            res[i] = numeric_limits<double>::quiet_NaN();

        int num_processes = (int)min(max((uint64_t)1, num_cols / 10), (uint64_t)g_tgstat->num_processes());

        ProgressReporter progress;
        progress.init(num_cols * num_cols / 2 - num_cols, 1);

        vdebug("num_processes: %d\n", num_processes);
        TGStat::prepare4multitasking();

        for (int iprocess = 0; iprocess < num_processes; ++iprocess) {
            if (!TGStat::launch_process()) {     // child process
                // two parts per process, one from each end, to balance the triangle; integer bounds so that
                // the parts cover all the columns: floor(i * N / (2P)) and N - ceil(i * N / (2P))
                uint64_t num_parts = 2 * num_processes;
                uint64_t scol[2] = { iprocess * num_cols / num_parts, num_cols - ((iprocess + 1) * num_cols + num_parts - 1) / num_parts };
                uint64_t ecol[2] = { (iprocess + 1) * num_cols / num_parts, num_cols - (iprocess * num_cols + num_parts - 1) / num_parts };
                uint64_t itr_idx = 0;

                for (int ipart = 0; ipart < 2; ipart++) {
                    for (uint64_t icol1 = 0; icol1 < num_cols; ++icol1) {
                        uint64_t start = max(scol[ipart], icol1 + 1);
                        if (start >= ecol[ipart])
                            break;

                        for (uint64_t icol2 = start; icol2 < ecol[ipart]; ++icol2) {
                            uint64_t idx = icol1 * num_cols + icol2;

                            if (nan_in_vals && pairwise_complete_obs && (constant_col[icol1] || constant_col[icol2])) {
                                // a constant column has no correlations
                            } else if (nan_in_vals && pairwise_complete_obs) {
                                uint64_t idx1 = icol1 * num_rows;
                                uint64_t idx2 = icol2 * num_rows;
                                double sum1 = 0;
                                double sum2 = 0;
                                double sum_square1 = 0;
                                double sum_square2 = 0;
                                double mean1, mean2;

                                if (spearman) {
                                    uint64_t indices[2] = { idx1, idx2 };
                                    vector<double *>::iterator sivals[2] = { pvals.begin() + idx1, pvals.begin() + idx2 };
                                    vector<double *>::iterator eivals[2] = { sivals[0] + num_rows, sivals[1] + num_rows };
                                    double *spvals[2] = { &vals.front() + idx1, &vals.front() + idx2 };

                                    for (int i = 0; i < 2; ++i) {
                                        // the fastest way to set all members of col_vals to NaN
                                        memcpy(&col_vals[i]->front(), &nan_col.front(), num_rows * sizeof(double));

                                        auto last_ival = sivals[i];
                                        uint64_t num_preceeding_vals = 0;
                                        uint64_t last_num_preceeding_vals = 0;

                                        for (auto ival = sivals[i]; ; ++ival) {
                                            if (ival == eivals[i] || **ival != **last_ival || std::isnan(**ival)) {
                                                double rank = (num_preceeding_vals + last_num_preceeding_vals - 1) / 2. + 1;
                                                while (last_ival != ival) {
                                                    (*col_vals[i])[*last_ival - spvals[i]] = rank;
                                                    ++last_ival;
                                                }

                                                if (ival == eivals[i] || std::isnan(**ival))
                                                    break;

                                                last_num_preceeding_vals = num_preceeding_vals;
                                            }

                                            if (!std::isnan(vals[*ival - spvals[i] + indices[1 - i]]))
                                                ++num_preceeding_vals;
                                        }
                                    }
                                } else {
                                    pcol_vals1 = &vals.front() + idx1;
                                    pcol_vals2 = &vals.front() + idx2;
                                }

                                uint64_t num_finite_pairs = 0;
                                res[idx] = 0;
                                for (uint64_t i = 0; i < num_rows; ++i) {
                                    double val1 = pcol_vals1[i];
                                    double val2 = pcol_vals2[i];

                                    if (!std::isnan(val1) && !std::isnan(val2)) {
                                        sum1 += val1;
                                        sum2 += val2;
                                        sum_square1 += val1 * val1;
                                        sum_square2 += val2 * val2;
                                        res[idx] += val1 * val2;  // => sum(X*Y)
                                        ++num_finite_pairs;
                                    }
                                }

                                if (num_finite_pairs) {
                                    mean1 = sum1 / num_finite_pairs;
                                    mean2 = sum2 / num_finite_pairs;
                                    double var1 = sum_square1 / num_finite_pairs - mean1 * mean1;
                                    double var2 = sum_square2 / num_finite_pairs - mean2 * mean2;

                                    // calculate correlation
                                    res[idx] /= num_finite_pairs;          // => mean(X*Y)
                                    res[idx] -= mean1 * mean2;         // => covariance(X,Y)
                                    res[idx] /= sqrt(var1 * var2);     // => correlation(X,Y)
                                } else
                                    res[idx] = numeric_limits<double>::quiet_NaN();
                            } else if (!nan_in_col[icol1] && !nan_in_col[icol2]) {
                                uint64_t idx1 = icol1 * num_rows;
                                uint64_t idx2 = icol2 * num_rows;
                                uint64_t end_idx1 = idx1 + num_rows;

                                res[idx] = 0;
                                while (idx1 < end_idx1)
                                    res[idx] += vals[idx1++] * vals[idx2++];  // => sum(X*Y)

                                res[idx] = res[idx] / num_rows;               // => mean(X*Y)
                                res[idx] -= means[icol1] * means[icol2];      // => covariance(X,Y)
                                res[idx] /= stddevs[icol1] * stddevs[icol2];  // => correlation(X,Y)
                            }
                            ++itr_idx;
                            TGStat::itr_idx(itr_idx);
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
        SEXP rold_dimnames = Rf_getAttrib(_x, R_DimNamesSymbol);
        SEXP rold_colnames = !Rf_isNull(rold_dimnames) && Rf_xlength(rold_dimnames) == 2 ? VECTOR_ELT(rold_dimnames, 1) : R_NilValue;

        if (tidy) {
            enum { COL1, COL2, COR, NUM_COLS };
            const char *COL_NAMES[NUM_COLS] = { "col1", "col2", "cor" };

            rprotect(answer = RSaneAllocVector(VECSXP, NUM_COLS));

            uint64_t answer_size = 0;
            auto cmp = [&res](uint64_t idx1, uint64_t idx2) { return fabs(res[idx1]) > fabs(res[idx2]); };
            priority_queue<uint64_t, vector<uint64_t>, decltype(cmp)> q(cmp);

            for (uint64_t icol1 = 0; icol1 < num_cols; ++icol1) {
                uint64_t idx = icol1;
                for (uint64_t icol2 = 0; icol2 < icol1; ++icol2) {
                    if (!std::isnan(res[idx]) && fabs(res[idx]) >= threshold)
                        ++answer_size;
                    idx += num_cols;
                }
            }

            SEXP rcol1, rcol2, rcor, rrownames, rcolnames;

            rprotect(rcol1 = RSaneAllocVector(INTSXP, answer_size));
            rprotect(rcol2 = RSaneAllocVector(INTSXP, answer_size));
            rprotect(rcor = RSaneAllocVector(REALSXP, answer_size));
            rprotect(rcolnames = RSaneAllocVector(STRSXP, NUM_COLS));
            rprotect(rrownames = RSaneAllocVector(INTSXP, answer_size));

            for (int i = 0; i < NUM_COLS; i++)
                SET_STRING_ELT(rcolnames, i, Rf_mkChar(COL_NAMES[i]));

            if (answer_size) {
                uint64_t i = 0;
                for (uint64_t icol1 = 0; icol1 < num_cols; ++icol1) {
                    for (uint64_t icol2 = icol1 + 1; icol2 < num_cols; ++icol2) {
                        uint64_t idx = (uint64_t)icol1 * num_cols + icol2;
                        if (!std::isnan(res[idx]) && fabs(res[idx]) >= threshold) {
                            INTEGER(rcol1)[i] = icol1 + 1;
                            INTEGER(rcol2)[i] = icol2 + 1;
                            REAL(rcor)[i] = res[idx];
                            INTEGER(rrownames)[i] = i + 1;
                            ++i;
                        }
                    }
                }
            }

            if (rold_colnames != R_NilValue) {
                Rf_setAttrib(rcol1, R_LevelsSymbol, rold_colnames);
                Rf_setAttrib(rcol1, R_ClassSymbol, Rf_mkString("factor"));
                Rf_setAttrib(rcol2, R_LevelsSymbol, rold_colnames);
                Rf_setAttrib(rcol2, R_ClassSymbol, Rf_mkString("factor"));
            }

            SET_VECTOR_ELT(answer, COL1, rcol1);
            SET_VECTOR_ELT(answer, COL2, rcol2);
            SET_VECTOR_ELT(answer, COR, rcor);

            Rf_setAttrib(answer, R_NamesSymbol, rcolnames);
            Rf_setAttrib(answer, R_ClassSymbol, Rf_mkString("data.frame"));
            Rf_setAttrib(answer, R_RowNamesSymbol, rrownames);
        } else {
            // copy the matrix below the diagonal to the upper part (the two parts are identical since cor(X,Y)=cor(Y,X)
            for (uint64_t icol1 = 0; icol1 < (uint64_t)num_cols; ++icol1) {
                uint64_t idx1 = icol1 * num_cols;
                uint64_t idx2 = icol1;
                for (uint64_t icol2 = 0; icol2 < icol1; ++icol2) {
                    res[idx1] = res[idx2];
                    idx1++;
                    idx2 += num_cols;
                }
                res[icol1 * (num_cols + 1)] = num_rows > 1 ? 1. : NA_REAL;
            }

            for (uint64_t i = 0; i < res_size; ++i) {
                if (std::isnan(res[i]))
                    res[i] = NA_REAL;
            }

            SEXP dim;

            rprotect(answer = RSaneAllocVector(REALSXP, res_size));
            memcpy(REAL(answer), res, res_sizeof);

            rprotect(dim = RSaneAllocVector(INTSXP, 2));
            INTEGER(dim)[0] = num_cols;
            INTEGER(dim)[1] = num_cols;
            Rf_setAttrib(answer, R_DimSymbol, dim);

            if (rold_colnames != R_NilValue) {
                SEXP dimnames;
                rprotect(dimnames = RSaneAllocVector(VECSXP, 2));
                SET_VECTOR_ELT(dimnames, 0, rold_colnames);
                SET_VECTOR_ELT(dimnames, 1, rold_colnames);
                Rf_setAttrib(answer, R_DimNamesSymbol, dimnames);
            }
        }
    } catch (TGLException &e) {
        if (!TGStat::is_kid() && res != (double *)MAP_FAILED) {
            munmap((char *)res, res_sizeof);  // needs to be char * for some versions of Solaris
            res = (double *)MAP_FAILED;
        }
		rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

    if (!TGStat::is_kid() && res != (double *)MAP_FAILED) {
        munmap((char *)res, res_sizeof);  // needs to be char * for some versions of Solaris
        res = (double *)MAP_FAILED;
    }
	rreturn(answer);
}

SEXP tgs_cross_cor(SEXP _x, SEXP _y, SEXP _pairwise_complete_obs, SEXP _spearman, SEXP _tidy, SEXP _threshold, SEXP _envir)
{
    SEXP answer = R_NilValue;
    double *res = (double *)MAP_FAILED;
    uint64_t res_sizeof = 0;

	try {
        TGStat tgstat(_envir);

		if ((!Rf_isReal(_x) && !Rf_isInteger(_x)) || Rf_xlength(_x) < 1)
			verror("\"x\" argument must be a matrix of numeric values");

        if ((!Rf_isReal(_y) && !Rf_isInteger(_y)) || Rf_xlength(_y) < 1)
            verror("\"y\" argument must be a matrix of numeric values");

        if (!Rf_isLogical(_pairwise_complete_obs) || Rf_xlength(_pairwise_complete_obs) != 1)
            verror("\"pairwise.complete.obs\" argument must be a logical value");

        if (!Rf_isLogical(_spearman) || Rf_xlength(_spearman) != 1)
            verror("\"spearman\" argument must be a logical value");

        if (!Rf_isLogical(_tidy) || Rf_xlength(_tidy) != 1)
            verror("\"tidy\" argument must be a logical value");

        if ((!Rf_isReal(_threshold) && !Rf_isInteger(_threshold)) || Rf_xlength(_threshold) != 1)
            verror("\"threshold\" argument must be a numeric value");

        SEXP rdim1 = Rf_getAttrib(_x, R_DimSymbol);
        SEXP rdim2 = Rf_getAttrib(_y, R_DimSymbol);

        if (!Rf_isInteger(rdim1) || Rf_xlength(rdim1) != 2)
            verror("\"x\" argument must be a matrix of numeric values");

        if (!Rf_isInteger(rdim2) || Rf_xlength(rdim2) != 2)
            verror("\"y\" argument must be a matrix of numeric values");

        if (Rf_nrows(_x) != Rf_nrows(_y))
            verror("\"x\" and \"y\" matrices have different number of rows");

        bool pairwise_complete_obs = Rf_asLogical(_pairwise_complete_obs);
        bool spearman = Rf_asLogical(_spearman);
        bool tidy = Rf_asLogical(_tidy);
        double threshold = fabs(Rf_asReal(_threshold));
        uint64_t num_rows = Rf_nrows(_x);
        uint64_t num_cols[2] = { (uint64_t)Rf_ncols(_x), (uint64_t)Rf_ncols(_y) };

        if (num_rows < 1 || num_cols[0] < 1)
            verror("\"x\" argument must be a matrix of numeric values");

        if (num_cols[1] < 1)
            verror("\"y\" argument must be a matrix of numeric values");

        uint64_t num_vals[2] = { num_rows * num_cols[0], num_rows * num_cols[1] };
        bool nan_in_vals = false;
        vector<bool> nan_in_col[2];
        vector<bool> constant_col[2];     // with pairwise.complete.obs: the columns that are constant in every pair
        vector<double> sums[2];
        vector<double> sums_square[2];
        vector<double> means[2];
        vector<double> stddevs[2];
        vector<double> vals[2];
        vector<double *> pvals[2];
        vector<double> col_vals1(num_rows);
        vector<double> col_vals2(num_rows);
        double *pcol_vals1 = &col_vals1.front();
        double *pcol_vals2 = &col_vals2.front();
        vector<double> *col_vals[2] = { &col_vals1, &col_vals2 };
        vector<double> nan_col(num_rows, numeric_limits<double>::quiet_NaN());

        for (int k = 0; k < 2; ++k) {
            nan_in_col[k].resize(num_cols[k], false);
            sums[k].resize(num_cols[k], 0);
            sums_square[k].resize(num_cols[k], 0);
            means[k].resize(num_cols[k], 0);
            stddevs[k].resize(num_cols[k], 0);
            vals[k].reserve(num_vals[k]);
        }

        for (int k = 0; k < 2; ++k) {
            SEXP x = k ? _y : _x;

            for (uint64_t i = 0; i < num_vals[k]; ++i) {
                if ((Rf_isReal(x) && !R_FINITE(REAL(x)[i])) || (Rf_isInteger(x) && INTEGER(x)[i] == NA_INTEGER)) {
                    nan_in_col[k][i / num_rows] = true;
                    nan_in_vals = true;
                    vals[k].push_back(numeric_limits<double>::quiet_NaN());
                } else
                    vals[k].push_back(Rf_isReal(x) ? REAL(x)[i] : INTEGER(x)[i]);
            }

            // the correlations of a constant column are NA, like those of a column with missing values
            flag_constant_cols(&vals[k].front(), num_rows, num_cols[k], nan_in_col[k]);
            if (pairwise_complete_obs)
                constant_col[k] = constant_cols(x);
        }

        // replace values with ranks if spearman=T
        if (spearman) {
            for (int k = 0; k < 2; ++k) {
                pvals[k].reserve(num_vals[k]);
                for (uint64_t i = 0; i < num_vals[k]; ++i)
                    pvals[k].push_back(&vals[k][i]);

                for (uint64_t icol = 0; icol < num_cols[k]; ++icol) {
                    if (nan_in_col[k][icol] && !pairwise_complete_obs)
                        continue;

                    vector<double *>::iterator sival = pvals[k].begin() + icol * num_rows;
                    vector<double *>::iterator eival = sival + num_rows;
                    vector<double *>::iterator last_ival = sival;

                    if (nan_in_col[k][icol])
                        sort(sival, eival, [](double *p1, double *p2) { return *p1 < *p2 || (!std::isnan(*p1) && std::isnan(*p2)); });
                    else
                        sort(sival, eival, [](double *p1, double *p2) { return *p1 < *p2; });

                    if (!nan_in_vals || !pairwise_complete_obs) {
                        for (auto ival = sival; ; ++ival) {
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
            }
        }

        for (int k = 0; k < 2; ++k) {
            for (uint64_t irow = 0; irow < num_rows; ++irow) {
                uint64_t idx = irow;
                for (uint64_t icol = 0; icol < num_cols[k]; ++icol) {
                    if (!nan_in_col[k][icol]) {
                        sums[k][icol] += vals[k][idx];
                        sums_square[k][icol] += vals[k][idx] * vals[k][idx];
                    }
                    idx += num_rows;
                }
            }

            for (uint64_t icol = 0; icol < num_cols[k]; ++icol) {
                if (!nan_in_col[k][icol]) {
                    means[k][icol] = sums[k][icol] / num_rows;

                    // we are calaculating standard deviation:
                    // sqrt(sum((x-mean)^2) / N)) = sqrt(sum(x^2) / N - mean^2)
                    stddevs[k][icol] = sqrt(sums_square[k][icol] / num_rows - means[k][icol] * means[k][icol]);
                }
            }
        }

        vdebug("Allocating shared memory for results\n");
        uint64_t res_size = num_cols[0] * num_cols[1];
        res_sizeof = sizeof(double) * res_size;
        res = (double *)mmap(NULL, res_sizeof, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);

        if (res == (double *)MAP_FAILED)
            verror("Failed to allocate shared memory: %s", strerror(errno));

        for (uint64_t i = 0; i < res_size; ++i)
            res[i] = numeric_limits<double>::quiet_NaN();

        // We split Y columns between different processes (and not X) because when
        // child process writes its result to res, the location at res will be a contiguous
        // chunk of memory. (res is organized like that: x1y1, x2y1, ..., x1y2, x2y2, ...)
        int num_processes = (int)min(max((uint64_t)1, num_cols[1] / 10), (uint64_t)g_tgstat->num_processes());

        ProgressReporter progress;
        progress.init(res_size, 1);

        vdebug("num_processes: %d\n", num_processes);
        TGStat::prepare4multitasking();

        for (int iprocess = 0; iprocess < num_processes; ++iprocess) {
            if (!TGStat::launch_process()) {     // child process
                uint64_t scol = iprocess * num_cols[1] / num_processes;
                uint64_t ecol = (iprocess + 1) * num_cols[1] / num_processes;
                uint64_t idx = scol * num_cols[0];
                uint64_t itr_idx = 0;

                for (uint64_t icol2 = scol; icol2 < ecol; ++icol2) {
                    for (uint64_t icol1 = 0; icol1 < num_cols[0]; ++icol1, ++idx) {
                        if (nan_in_vals && pairwise_complete_obs && (constant_col[0][icol1] || constant_col[1][icol2])) {
                            // a constant column has no correlations
                        } else if (nan_in_vals && pairwise_complete_obs) {
                            uint64_t idx1 = icol1 * num_rows;
                            uint64_t idx2 = icol2 * num_rows;
                            double sum1 = 0;
                            double sum2 = 0;
                            double sum_square1 = 0;
                            double sum_square2 = 0;
                            double mean1, mean2;

                            if (spearman) {
                                uint64_t indices[2] = { idx1, idx2 };
                                vector<double *>::iterator sivals[2] = { pvals[0].begin() + idx1, pvals[1].begin() + idx2 };
                                vector<double *>::iterator eivals[2] = { sivals[0] + num_rows, sivals[1] + num_rows };
                                double *spvals[2] = { &vals[0].front() + idx1, &vals[1].front() + idx2 };

                                for (int i = 0; i < 2; ++i) {
                                    // the fastest way to set all members of col_vals to NaN
                                    memcpy(&col_vals[i]->front(), &nan_col.front(), num_rows * sizeof(double));

                                    auto last_ival = sivals[i];
                                    uint64_t num_preceeding_vals = 0;
                                    uint64_t last_num_preceeding_vals = 0;

                                    for (auto ival = sivals[i]; ; ++ival) {
                                        if (ival == eivals[i] || **ival != **last_ival || std::isnan(**ival)) {
                                            double rank = (num_preceeding_vals + last_num_preceeding_vals - 1) / 2. + 1;
                                            while (last_ival != ival) {
                                                (*col_vals[i])[*last_ival - spvals[i]] = rank;
                                                ++last_ival;
                                            }

                                            if (ival == eivals[i] || std::isnan(**ival))
                                                break;

                                            last_num_preceeding_vals = num_preceeding_vals;
                                        }

                                        if (!std::isnan(vals[1 - i][*ival - spvals[i] + indices[1 - i]]))
                                            ++num_preceeding_vals;
                                    }
                                }
                            } else {
                                pcol_vals1 = &vals[0].front() + idx1;
                                pcol_vals2 = &vals[1].front() + idx2;
                            }

                            uint64_t num_finite_pairs = 0;
                            res[idx] = 0;
                            for (uint64_t i = 0; i < num_rows; ++i) {
                                double val1 = pcol_vals1[i];
                                double val2 = pcol_vals2[i];

                                if (!std::isnan(val1) && !std::isnan(val2)) {
                                    sum1 += val1;
                                    sum2 += val2;
                                    sum_square1 += val1 * val1;
                                    sum_square2 += val2 * val2;
                                    res[idx] += val1 * val2;  // => sum(X*Y)
                                    ++num_finite_pairs;
                                }
                            }

                            if (num_finite_pairs) {
                                mean1 = sum1 / num_finite_pairs;
                                mean2 = sum2 / num_finite_pairs;
                                double var1 = sum_square1 / num_finite_pairs - mean1 * mean1;
                                double var2 = sum_square2 / num_finite_pairs - mean2 * mean2;

                                // calculate correlation
                                res[idx] /= num_finite_pairs;          // => mean(X*Y)
                                res[idx] -= mean1 * mean2;         // => covariance(X,Y)
                                res[idx] /= sqrt(var1 * var2);     // => correlation(X,Y)
                            } else
                                res[idx] = numeric_limits<double>::quiet_NaN();
                        } else if (!nan_in_col[0][icol1] && !nan_in_col[1][icol2]) {
                            uint64_t idx1 = icol1 * num_rows;
                            uint64_t idx2 = icol2 * num_rows;
                            uint64_t end_idx1 = idx1 + num_rows;

                            res[idx] = 0;
                            while (idx1 < end_idx1)
                                res[idx] += vals[0][idx1++] * vals[1][idx2++];  // => sum(X*Y)

                            res[idx] = res[idx] / num_rows;               // => mean(X*Y)
                            res[idx] -= means[0][icol1] * means[1][icol2];      // => covariance(X,Y)
                            res[idx] /= stddevs[0][icol1] * stddevs[1][icol2];  // => correlation(X,Y)
                        }
                        ++itr_idx;
                        TGStat::itr_idx(itr_idx);
                    }
                }
                rexit();
            }
        }

        while (TGStat::wait_for_kids(3000))
            progress.report(TGStat::itr_idx_sum() - progress.get_elapsed_steps());

        progress.report_last();

        // assemble the answer
        SEXP rold_dimnames[2] = { Rf_getAttrib(_x, R_DimNamesSymbol), Rf_getAttrib(_y, R_DimNamesSymbol) };
        SEXP rold_colnames[2] = {
            !Rf_isNull(rold_dimnames[0]) && Rf_xlength(rold_dimnames[0]) == 2 ? VECTOR_ELT(rold_dimnames[0], 1) : R_NilValue,
            !Rf_isNull(rold_dimnames[1]) && Rf_xlength(rold_dimnames[1]) == 2 ? VECTOR_ELT(rold_dimnames[1], 1) : R_NilValue
        };

        if (tidy) {
            enum { COL1, COL2, COR, NUM_COLS };
            const char *COL_NAMES[NUM_COLS] = { "col1", "col2", "cor" };

            rprotect(answer = RSaneAllocVector(VECSXP, NUM_COLS));

            uint64_t answer_size = 0;

            for (uint64_t i = 0; i < res_size; ++i) {
                if (!std::isnan(res[i]) && fabs(res[i]) >= threshold)
                    ++answer_size;
            }

            SEXP rcol1, rcol2, rcor, rrownames, rcolnames;

            rprotect(rcol1 = RSaneAllocVector(INTSXP, answer_size));
            rprotect(rcol2 = RSaneAllocVector(INTSXP, answer_size));
            rprotect(rcor = RSaneAllocVector(REALSXP, answer_size));
            rprotect(rcolnames = RSaneAllocVector(STRSXP, NUM_COLS));
            rprotect(rrownames = RSaneAllocVector(INTSXP, answer_size));

            for (int i = 0; i < NUM_COLS; i++)
                SET_STRING_ELT(rcolnames, i, Rf_mkChar(COL_NAMES[i]));

            if (answer_size) {
                uint64_t idx_answer = 0;
                for (uint64_t icol1 = 0; icol1 < num_cols[0]; ++icol1) {
                    uint64_t idx_m = icol1;
                    for (uint64_t icol2 = 0; icol2 < num_cols[1]; ++icol2) {
                        if (!std::isnan(res[idx_m]) && fabs(res[idx_m]) >= threshold) {
                            INTEGER(rcol1)[idx_answer] = icol1 + 1;
                            INTEGER(rcol2)[idx_answer] = icol2 + 1;
                            REAL(rcor)[idx_answer] = res[idx_m];
                            INTEGER(rrownames)[idx_answer] = idx_answer + 1;
                            ++idx_answer;
                        }
                        idx_m += num_cols[0];
                    }
                }
            }

            if (rold_colnames[0] != R_NilValue) {
                Rf_setAttrib(rcol1, R_LevelsSymbol, rold_colnames[0]);
                Rf_setAttrib(rcol1, R_ClassSymbol, Rf_mkString("factor"));
            }
            if (rold_colnames[1] != R_NilValue) {
                Rf_setAttrib(rcol2, R_LevelsSymbol, rold_colnames[1]);
                Rf_setAttrib(rcol2, R_ClassSymbol, Rf_mkString("factor"));
            }

            SET_VECTOR_ELT(answer, COL1, rcol1);
            SET_VECTOR_ELT(answer, COL2, rcol2);
            SET_VECTOR_ELT(answer, COR, rcor);

            Rf_setAttrib(answer, R_NamesSymbol, rcolnames);
            Rf_setAttrib(answer, R_ClassSymbol, Rf_mkString("data.frame"));
            Rf_setAttrib(answer, R_RowNamesSymbol, rrownames);
        } else {
            SEXP dim;
            SEXP dimnames;

            for (uint64_t i = 0; i < res_size; ++i) {
                if (std::isnan(res[i]))
                    res[i] = NA_REAL;
            }

            rprotect(answer = RSaneAllocVector(REALSXP, res_size));
            memcpy(REAL(answer), res, res_size * sizeof(double));

            rprotect(dim = RSaneAllocVector(INTSXP, 2));
            INTEGER(dim)[0] = num_cols[0];
            INTEGER(dim)[1] = num_cols[1];
            Rf_setAttrib(answer, R_DimSymbol, dim);

            rprotect(dimnames = RSaneAllocVector(VECSXP, 2));
            SET_VECTOR_ELT(dimnames, 0, rold_colnames[0]);
            SET_VECTOR_ELT(dimnames, 1, rold_colnames[1]);
            Rf_setAttrib(answer, R_DimNamesSymbol, dimnames);
        }
    } catch (TGLException &e) {
        if (!TGStat::is_kid() && res != (double *)MAP_FAILED) {
            munmap((char *)res, res_sizeof);  // needs to be char * for some versions of Solaris
            res = (double *)MAP_FAILED;
        }
		rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

    if (!TGStat::is_kid() && res != (double *)MAP_FAILED) {
        munmap((char *)res, res_sizeof);  // needs to be char * for some versions of Solaris
        res = (double *)MAP_FAILED;
    }
	rreturn(answer);
}

SEXP tgs_cor_blas(SEXP _x, SEXP _pairwise_complete_obs, SEXP _spearman, SEXP _tidy, SEXP _threshold, SEXP _envir)
{
    SEXP answer = R_NilValue;
    double *z[2] = { (double *)MAP_FAILED, (double *)MAP_FAILED };
    uint64_t z_sizeof[2] = { 0, 0 };

    try {
        struct Mem {
            double *m{NULL};
            double *mask{NULL};
            double *n{NULL};
            double *s_x{NULL};
            double *cov_n{NULL};
            double *var_n{NULL};
            double *res{NULL};
            ~Mem() { free(m); free(mask); free(n); free(s_x); free(cov_n); free(var_n); free(res); }
        } mem;

        TGStat tgstat(_envir);

        if ((!Rf_isReal(_x) && !Rf_isInteger(_x)) || Rf_xlength(_x) < 1)
            verror("\"x\" argument must be a matrix of numeric values");

        if (!Rf_isLogical(_pairwise_complete_obs) || Rf_xlength(_pairwise_complete_obs) != 1)
            verror("\"pairwise.complete.obs\" argument must be a logical value");

        if (!Rf_isLogical(_spearman) || Rf_xlength(_spearman) != 1)
            verror("\"spearman\" argument must be a logical value");

        if (!Rf_isLogical(_tidy) || Rf_xlength(_tidy) != 1)
            verror("\"tidy\" argument must be a logical value");

        if ((!Rf_isReal(_threshold) && !Rf_isInteger(_threshold)) || Rf_xlength(_threshold) != 1)
            verror("\"threshold\" argument must be a numeric value");

        SEXP rdim = Rf_getAttrib(_x, R_DimSymbol);

        if (!Rf_isInteger(rdim) || Rf_xlength(rdim) != 2)
            verror("\"x\" argument must be a matrix of numeric values");

        bool pairwise_complete_obs = Rf_asLogical(_pairwise_complete_obs);
        bool spearman = Rf_asLogical(_spearman);
        bool tidy = Rf_asLogical(_tidy);
        double threshold = fabs(Rf_asReal(_threshold));
        uint64_t num_dims = Rf_nrows(_x);
        uint64_t num_points = Rf_ncols(_x);
        int num_dims32 = (int)num_dims;
        int num_points32 = (int)num_points;

        if (num_dims < 1 || num_points < 1)
            verror("\"x\" argument must be a matrix of numeric values");

        uint64_t num_vals = num_points * num_dims;
        uint64_t res_size = num_points * num_points;
        // pairwise.complete.obs changes the correlations only if there are missing values
        bool pairwise_nan = pairwise_complete_obs && !all_finite(_x);

        if (spearman && pairwise_nan)
            verror("BLAS implementation of tgs_cor does not support spearman with pairwise.complete.obs when x contains NA / NaN / Inf");

        vdebug("START BLAS COR\n");

        // the correlations are computed directly in the answer, unless it is tidy
        double *res;

        if (tidy) {
            if (posix_memalign((void **)&mem.res, 64, sizeof(double) * res_size))
                verror("%s", strerror(errno));
            res = mem.res;
        } else {
            rprotect(answer = RSaneAllocVector(REALSXP, res_size));
            res = REAL(answer);
        }

        if (pairwise_nan) {
            // some BLAS implementations ask to align double arrays to 64 for improved efficiency
            if (posix_memalign((void **)&mem.m, 64, sizeof(double) * num_vals))
                verror("%s", strerror(errno));

            if (posix_memalign((void **)&mem.mask, 64, sizeof(double) * num_vals))
                verror("%s", strerror(errno));

            for (uint64_t i = 0; i < num_vals; ++i) {
                if ((Rf_isReal(_x) && !R_FINITE(REAL(_x)[i])) || (Rf_isInteger(_x) && INTEGER(_x)[i] == NA_INTEGER))
                    mem.m[i] = mem.mask[i] = 0.;
                else {
                    mem.m[i] = Rf_isReal(_x) ? REAL(_x)[i] : INTEGER(_x)[i];
                    mem.mask[i] = 1.;
                }
            }

            // correlation with pairwise complete obs is calculated as:
            //     n <- t(mask) %*% mask
            //     s_x <- t(mask) %*% m
            //     cov_n <- (t(m) %*% m) - (s_x * t(s_x)) / n
            //     var_n <- (t(mask) %*% (m * m)) - (s_x * s_x) / n
            //     res <- cov_n / sqrt(var_n * t(var_n))

            if (posix_memalign((void **)&mem.n, 64, sizeof(double) * res_size))
                verror("%s", strerror(errno));

            if (posix_memalign((void **)&mem.s_x, 64, sizeof(double) * res_size))
                verror("%s", strerror(errno));

            if (posix_memalign((void **)&mem.cov_n, 64, sizeof(double) * res_size))
                verror("%s", strerror(errno));

            if (posix_memalign((void **)&mem.var_n, 64, sizeof(double) * res_size))
                verror("%s", strerror(errno));

            ProgressReporter progress;
            progress.init(6, 1);

            //     n <- t(mask) %*% mask
            {
                char uplo = 'L';
                char trans = 'T';
                double alpha = 1;
                double beta = 0;
                F77_NAME(dsyrk)(&uplo, &trans, &num_points32, &num_dims32, &alpha, mem.mask, &num_dims32, &beta, mem.n, &num_points32 FCONE FCONE);
                check_interrupt();
                progress.report(1);
            }

            //     s_x <- t(mask) %*% m
            {
                char transa = 'T';
                char transb = 'N';
                double alpha = 1;
                double beta = 0;
                F77_NAME(dgemm)(&transa, &transb, &num_points32, &num_points32, &num_dims32, &alpha, mem.mask, &num_dims32, mem.m, &num_dims32, &beta, mem.s_x, &num_points32 FCONE FCONE);
                check_interrupt();
                progress.report(1);
            }

            //     cov_n <- t(m) %*% m
            {
                char uplo = 'L';
                char trans = 'T';
                double alpha = 1;
                double beta = 0;
                F77_NAME(dsyrk)(&uplo, &trans, &num_points32, &num_dims32, &alpha, mem.m, &num_dims32, &beta, mem.cov_n, &num_points32 FCONE FCONE);
                check_interrupt();
                progress.report(1);
            }

            //     cov_n <- cov_n - (s_x * t(s_x)) / n
            for (uint64_t ipoint1 = 0, idx1 = 0; ipoint1 < num_points; ++ipoint1) {
                idx1 += ipoint1 + 1;
                for (uint64_t ipoint2 = ipoint1 + 1, idx2 = ipoint2 * num_points + ipoint1; ipoint2 < num_points; ++ipoint2) {
                    mem.cov_n[idx1] -= mem.s_x[idx1] * mem.s_x[idx2] / mem.n[idx1];
                    ++idx1;
                    idx2 += num_points;
                }
            }
            check_interrupt();

            //    m <- m * m
            for (uint64_t i = 0; i < num_vals; ++i)
                mem.m[i] *= mem.m[i];

            check_interrupt();
            progress.report(1);

            //     var_n <- t(mask) %*% (m * m)
            {
                char transa = 'T';
                char transb = 'N';
                double alpha = 1;
                double beta = 0;
                F77_NAME(dgemm)(&transa, &transb, &num_points32, &num_points32, &num_dims32, &alpha, mem.mask, &num_dims32, mem.m, &num_dims32, &beta, mem.var_n, &num_points32 FCONE FCONE);
                check_interrupt();
                progress.report(1);
            }

            //     var_n <- var_n - (s_x * s_x) / n
            for (uint64_t ipoint1 = 0, idx1 = 0; ipoint1 < num_points; ++ipoint1) {
                idx1 += ipoint1 + 1;
                for (uint64_t ipoint2 = ipoint1 + 1, idx2 = ipoint2 * num_points + ipoint1; ipoint2 < num_points; ++ipoint2) {
                    mem.var_n[idx1] -= mem.s_x[idx1] * mem.s_x[idx1] / mem.n[idx1];
                    mem.var_n[idx2] -= mem.s_x[idx2] * mem.s_x[idx2] / mem.n[idx1];
                    ++idx1;
                    idx2 += num_points;
                }
            }
            check_interrupt();

            //     res <- cov_n / sqrt(var_n * t(var_n))
            vector<bool> constant_col = constant_cols(_x);

            for (uint64_t ipoint1 = 0, idx1 = 0; ipoint1 < num_points; ++ipoint1) {
                idx1 += ipoint1 + 1;
                for (uint64_t ipoint2 = ipoint1 + 1, idx2 = ipoint2 * num_points + ipoint1; ipoint2 < num_points; ++ipoint2) {
                    // a column that is constant in every pair has no correlations
                    if (constant_col[ipoint1] || constant_col[ipoint2])
                        res[idx1] = res[idx2] = NA_REAL;
                    else
                        res[idx1] = res[idx2] = mem.cov_n[idx1] / sqrt(mem.var_n[idx1] * mem.var_n[idx2]);
                    ++idx1;
                    idx2 += num_points;
                }
            }
            check_interrupt();

            for (uint64_t idx = 0; idx < res_size; idx += num_points + 1)
                res[idx] = 1.;
            check_interrupt();

            progress.report_last();
        } else {
            // The columns are standardized, so that the correlation matrix is t(z) %*% z. The columns that have no
            // correlations (a missing value or constant values) are flagged and their correlations are NA.
            ProgressReporter progress;
            progress.init(3, 1);

            vector<bool> flagged[2];

            standardize_cols(&_x, 1, spearman, z, z_sizeof, flagged);
            check_interrupt();
            progress.report(1);

            //     res <- t(z) %*% z, the lower triangle
            {
                char uplo = 'L';
                char trans = 'T';
                double alpha = 1;
                double beta = 0;
                F77_NAME(dsyrk)(&uplo, &trans, &num_points32, &num_dims32, &alpha, z[0], &num_dims32, &beta, res, &num_points32 FCONE FCONE);
            }
            unmap_shared(z, z_sizeof);
            check_interrupt();
            progress.report(1);

            vector<uint64_t> flagged_points;

            for (uint64_t ipoint = 0; ipoint < num_points; ++ipoint) {
                if (flagged[0][ipoint])
                    flagged_points.push_back(ipoint);
            }

            if (!flagged_points.empty()) {
                for (uint64_t ipoint = 0; ipoint < num_points; ++ipoint) {
                    double *col = res + ipoint * num_points;

                    if (flagged[0][ipoint])
                        fill(col + ipoint, col + num_points, NA_REAL);
                    else {
                        for (auto ipoint2 = lower_bound(flagged_points.begin(), flagged_points.end(), ipoint); ipoint2 != flagged_points.end(); ++ipoint2)
                            col[*ipoint2] = NA_REAL;
                    }
                }
            }

            for (uint64_t idx = 0; idx < res_size; idx += num_points + 1)
                res[idx] = num_dims > 1 ? 1. : NA_REAL;

            if (!tidy) {
                // copy the lower triangle to the upper one, in tiles that fit in the cache
                const uint64_t tile = 64;

                for (uint64_t col1 = 0; col1 < num_points; col1 += tile) {
                    uint64_t col2 = min(col1 + tile, num_points);

                    for (uint64_t row1 = col1; row1 < num_points; row1 += tile) {
                        uint64_t row2 = min(row1 + tile, num_points);

                        for (uint64_t icol = col1; icol < col2; ++icol) {
                            for (uint64_t irow = max(row1, icol + 1); irow < row2; ++irow)
                                res[irow * num_points + icol] = res[icol * num_points + irow];
                        }
                    }
                }
            }
            check_interrupt();
            progress.report_last();
        }
        vdebug("END BLAS COR\n");

//memcpy(res, mem.var_n, sizeof(double) * res_size);
//{
//SEXP rdims;
//rprotect(answer = RSaneAllocVector(REALSXP, (uint64_t)num_points * num_points));
//rprotect(rdims = RSaneAllocVector(INTSXP, 2));
//INTEGER(rdims)[0] = INTEGER(rdims)[1] = num_points;
//Rf_setAttrib(answer, R_DimSymbol, rdims);
//memcpy(REAL(answer), res, Rf_length(answer) * sizeof(REAL(answer)[0]));
//return answer;
//}


        // assemble the answer
        SEXP rold_dimnames = Rf_getAttrib(_x, R_DimNamesSymbol);
        SEXP rold_colnames = !Rf_isNull(rold_dimnames) && Rf_xlength(rold_dimnames) == 2 ? VECTOR_ELT(rold_dimnames, 1) : R_NilValue;

        if (tidy) {
            enum { COL1, COL2, COR, NUM_COLS };
            const char *COL_NAMES[NUM_COLS] = { "col1", "col2", "cor" };

            rprotect(answer = RSaneAllocVector(VECSXP, NUM_COLS));

            uint64_t answer_size = 0;

            for (uint64_t icol1 = 0; icol1 < num_points; ++icol1) {
                uint64_t idx = icol1;
                for (uint64_t icol2 = 0; icol2 < icol1; ++icol2) {
                    if (!std::isnan(res[idx]) && fabs(res[idx]) >= threshold)
                        ++answer_size;
                    idx += num_points;
                }
            }

            SEXP rcol1, rcol2, rcor, rrownames, rcolnames;

            rprotect(rcol1 = RSaneAllocVector(INTSXP, answer_size));
            rprotect(rcol2 = RSaneAllocVector(INTSXP, answer_size));
            rprotect(rcor = RSaneAllocVector(REALSXP, answer_size));
            rprotect(rcolnames = RSaneAllocVector(STRSXP, NUM_COLS));
            rprotect(rrownames = RSaneAllocVector(INTSXP, answer_size));

            for (int i = 0; i < NUM_COLS; i++)
                SET_STRING_ELT(rcolnames, i, Rf_mkChar(COL_NAMES[i]));

            if (answer_size) {
                uint64_t i = 0;
                for (uint64_t icol1 = 0; icol1 < num_points; ++icol1) {
                    for (uint64_t icol2 = icol1 + 1; icol2 < num_points; ++icol2) {
                        uint64_t idx = (uint64_t)icol1 * num_points + icol2;
                        if (!std::isnan(res[idx]) && fabs(res[idx]) >= threshold) {
                            INTEGER(rcol1)[i] = icol1 + 1;
                            INTEGER(rcol2)[i] = icol2 + 1;
                            REAL(rcor)[i] = res[idx];
                            INTEGER(rrownames)[i] = i + 1;
                            ++i;
                        }
                    }
                }
            }

            if (rold_colnames != R_NilValue) {
                Rf_setAttrib(rcol1, R_LevelsSymbol, rold_colnames);
                Rf_setAttrib(rcol1, R_ClassSymbol, Rf_mkString("factor"));
                Rf_setAttrib(rcol2, R_LevelsSymbol, rold_colnames);
                Rf_setAttrib(rcol2, R_ClassSymbol, Rf_mkString("factor"));
            }

            SET_VECTOR_ELT(answer, COL1, rcol1);
            SET_VECTOR_ELT(answer, COL2, rcol2);
            SET_VECTOR_ELT(answer, COR, rcor);

            Rf_setAttrib(answer, R_NamesSymbol, rcolnames);
            Rf_setAttrib(answer, R_ClassSymbol, Rf_mkString("data.frame"));
            Rf_setAttrib(answer, R_RowNamesSymbol, rrownames);
        } else {
            SEXP dim;

            if (pairwise_nan) {
                for (uint64_t i = 0; i < res_size; ++i) {
                    if (std::isnan(res[i]))
                        res[i] = NA_REAL;
                }
            }

            rprotect(dim = RSaneAllocVector(INTSXP, 2));
            INTEGER(dim)[0] = num_points;
            INTEGER(dim)[1] = num_points;
            Rf_setAttrib(answer, R_DimSymbol, dim);

            if (rold_colnames != R_NilValue) {
                SEXP dimnames;
                rprotect(dimnames = RSaneAllocVector(VECSXP, 2));
                SET_VECTOR_ELT(dimnames, 0, rold_colnames);
                SET_VECTOR_ELT(dimnames, 1, rold_colnames);
                Rf_setAttrib(answer, R_DimNamesSymbol, dimnames);
            }
        }
    } catch (TGLException &e) {
        unmap_shared(z, z_sizeof);
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        unmap_shared(z, z_sizeof);
        rerror("Out of memory");
    }

    unmap_shared(z, z_sizeof);
    rreturn(answer);
}

SEXP tgs_cross_cor_blas(SEXP _x, SEXP _y, SEXP _pairwise_complete_obs, SEXP _spearman, SEXP _tidy, SEXP _threshold, SEXP _envir)
{
    SEXP answer = R_NilValue;
    double *z[2] = { (double *)MAP_FAILED, (double *)MAP_FAILED };
    uint64_t z_sizeof[2] = { 0, 0 };

    try {
        struct Mem {
            double *m[2]{NULL};
            double *mask[2]{NULL};
            double *n{NULL};
            double *s_x{NULL};
            double *s_y{NULL};
            double *cov_n{NULL};
            double *varx_n{NULL};
            double *vary_n{NULL};
            double *res{NULL};
            ~Mem() { 
                for (int i = 0; i < 2; ++i) {
                    free(m[i]);
                    free(mask[i]);
                }
                free(n);
                free(s_x);
                free(s_y);
                free(cov_n);
                free(varx_n);
                free(vary_n);
                free(res);
            }
        } mem;

        TGStat tgstat(_envir);

        if ((!Rf_isReal(_x) && !Rf_isInteger(_x)) || Rf_xlength(_x) < 1)
            verror("\"x\" argument must be a matrix of numeric values");

        if ((!Rf_isReal(_y) && !Rf_isInteger(_y)) || Rf_xlength(_y) < 1)
            verror("\"y\" argument must be a matrix of numeric values");

        if (!Rf_isLogical(_pairwise_complete_obs) || Rf_xlength(_pairwise_complete_obs) != 1)
            verror("\"pairwise.complete.obs\" argument must be a logical value");

        if (!Rf_isLogical(_spearman) || Rf_xlength(_spearman) != 1)
            verror("\"spearman\" argument must be a logical value");

        if (!Rf_isLogical(_tidy) || Rf_xlength(_tidy) != 1)
            verror("\"tidy\" argument must be a logical value");

        if ((!Rf_isReal(_threshold) && !Rf_isInteger(_threshold)) || Rf_xlength(_threshold) != 1)
            verror("\"threshold\" argument must be a numeric value");

        SEXP rdim1 = Rf_getAttrib(_x, R_DimSymbol);
        SEXP rdim2 = Rf_getAttrib(_y, R_DimSymbol);

        if (!Rf_isInteger(rdim1) || Rf_xlength(rdim1) != 2)
            verror("\"x\" argument must be a matrix of numeric values");

        if (!Rf_isInteger(rdim2) || Rf_xlength(rdim2) != 2)
            verror("\"y\" argument must be a matrix of numeric values");

        if (Rf_nrows(_x) != Rf_nrows(_y))
            verror("\"x\" and \"y\" matrices have different number of rows");

        bool pairwise_complete_obs = Rf_asLogical(_pairwise_complete_obs);
        bool spearman = Rf_asLogical(_spearman);
        bool tidy = Rf_asLogical(_tidy);
        double threshold = fabs(Rf_asReal(_threshold));
        uint64_t num_dims = Rf_nrows(_x);
        uint64_t num_points[2] = { (uint64_t)Rf_ncols(_x), (uint64_t)Rf_ncols(_y) };
        int num_dims32 = (int)num_dims;
        int num_points32[2] = { (int)num_points[0], (int)num_points[1] };

        if (num_dims < 1 || num_points[0] < 1)
            verror("\"x\" argument must be a matrix of numeric values");

        if (num_points[1] < 1)
            verror("\"y\" argument must be a matrix of numeric values");

        uint64_t num_vals[2] = { num_points[0] * num_dims, num_points[1] * num_dims };
        uint64_t res_size = num_points[0] * num_points[1];

        // pairwise.complete.obs changes the correlations only if there are missing values
        bool pairwise_nan = pairwise_complete_obs && (!all_finite(_x) || !all_finite(_y));

        if (spearman && pairwise_nan)
            verror("BLAS implementation of tgs_cor does not support spearman with pairwise.complete.obs when x or y contain NA / NaN / Inf");

        vdebug("START BLAS COR\n");

        // the correlations are computed directly in the answer, unless it is tidy
        double *res;

        if (tidy) {
            if (posix_memalign((void **)&mem.res, 64, sizeof(double) * res_size))
                verror("%s", strerror(errno));
            res = mem.res;
        } else {
            rprotect(answer = RSaneAllocVector(REALSXP, res_size));
            res = REAL(answer);
        }

        if (pairwise_nan) {
            // some BLAS implementations ask to align double arrays to 64 for improved efficiency
            for (int k = 0; k < 2; k++) {
                if (posix_memalign((void **)&mem.m[k], 64, sizeof(double) * num_vals[k]) ||
                    posix_memalign((void **)&mem.mask[k], 64, sizeof(double) * num_vals[k]))
                    verror("%s", strerror(errno));

                SEXP x = k ? _y : _x;

                for (uint64_t i = 0; i < num_vals[k]; ++i) {
                    if ((Rf_isReal(x) && !R_FINITE(REAL(x)[i])) || (Rf_isInteger(x) && INTEGER(x)[i] == NA_INTEGER))
                        mem.m[k][i] = mem.mask[k][i] = 0.;
                    else {
                        mem.m[k][i] = Rf_isReal(x) ? REAL(x)[i] : INTEGER(x)[i];
                        mem.mask[k][i] = 1.;
                    }
                }
            }

            // given that X is a matrix of d x k and Y is a matrix of d x m
            // correlation with pairwise complete jobs is calculated as:
            //     n <- t(maskx) %*% masky
            //     s_x <- t(mx) %*% masky
            //     s_y <- t(maskx) %*% my
            //     cov_n <- (t(mx) %*% my) - (s_x * s_y) / n
            //     varx_n <- (t(mx * mx) %*% masky) - (s_x * s_x) / n
            //     vary_n <- (t(maskx) %*% (my * my)) - (s_y * s_y) / n
            //     res <- cov_n / sqrt(varx_n * vary_n)
            //
            // * all the results (e.g. n, s_x, s_y, ...) are matrices of k x m

            if (posix_memalign((void **)&mem.s_x, 64, sizeof(double) * res_size) ||
                posix_memalign((void **)&mem.s_y, 64, sizeof(double) * res_size) ||
                posix_memalign((void **)&mem.varx_n, 64, sizeof(double) * res_size) ||
                posix_memalign((void **)&mem.vary_n, 64, sizeof(double) * res_size) ||
                posix_memalign((void **)&mem.n, 64, sizeof(double) * res_size) ||
                posix_memalign((void **)&mem.cov_n, 64, sizeof(double) * res_size))
                verror("%s", strerror(errno));

            ProgressReporter progress;
            progress.init(8, 1);

            //     n <- t(maskx) %*% masky
            {
                char transa = 'T';
                char transb = 'N';
                double alpha = 1;
                double beta = 0;
                F77_NAME(dgemm)(&transa, &transb, &num_points32[0], &num_points32[1], &num_dims32, &alpha, mem.mask[0], &num_dims32, mem.mask[1], &num_dims32, &beta, mem.n, &num_points32[0] FCONE FCONE);
                check_interrupt();
                progress.report(1);
            }

            //     s_x <- t(mx) %*% masky
            {
                char transa = 'T';
                char transb = 'N';
                double alpha = 1;
                double beta = 0;
                F77_NAME(dgemm)(&transa, &transb, &num_points32[0], &num_points32[1], &num_dims32, &alpha, mem.m[0], &num_dims32, mem.mask[1], &num_dims32, &beta, mem.s_x, &num_points32[0] FCONE FCONE);
                check_interrupt();
                progress.report(1);
            }

            //     s_y <- t(maskx) %*% my
            {
                char transa = 'T';
                char transb = 'N';
                double alpha = 1;
                double beta = 0;
                F77_NAME(dgemm)(&transa, &transb, &num_points32[0], &num_points32[1], &num_dims32, &alpha, mem.mask[0], &num_dims32, mem.m[1], &num_dims32, &beta, mem.s_y, &num_points32[0] FCONE FCONE);
                check_interrupt();
                progress.report(1);
            }

            //     cov_n <- t(mx) %*% my
            {
                char transa = 'T';
                char transb = 'N';
                double alpha = 1;
                double beta = 0;
                F77_NAME(dgemm)(&transa, &transb, &num_points32[0], &num_points32[1], &num_dims32, &alpha, mem.m[0], &num_dims32, mem.m[1], &num_dims32, &beta, mem.cov_n, &num_points32[0] FCONE FCONE);
                check_interrupt();
                progress.report(1);
            }

            //     cov_n <- cov_n - (s_x * s_y) / n
            for (uint64_t i = 0; i < res_size; ++i)
                mem.cov_n[i] -= mem.s_x[i] * mem.s_y[i] / mem.n[i];
            check_interrupt();

            //    m <- m * m
            for (int k = 0; k < 2; ++k) {
                for (uint64_t i = 0; i < num_vals[k]; ++i)
                    mem.m[k][i] *= mem.m[k][i];
                check_interrupt();
            }
            progress.report(1);

            //     varx_n <- t(mx * mx) %*% masky
            {
                char transa = 'T';
                char transb = 'N';
                double alpha = 1;
                double beta = 0;
                F77_NAME(dgemm)(&transa, &transb, &num_points32[0], &num_points32[1], &num_dims32, &alpha, mem.m[0], &num_dims32, mem.mask[1], &num_dims32, &beta, mem.varx_n, &num_points32[0] FCONE FCONE);
                check_interrupt();
                progress.report(1);
            }

            //     varx_n <- varx_n - (s_x * s_x) / n
            for (uint64_t i = 0; i < res_size; ++i)
                mem.varx_n[i] -= mem.s_x[i] * mem.s_x[i] / mem.n[i];
            check_interrupt();

            //     vary_n <- maskx %*% t(my * my)
            {
                char transa = 'T';
                char transb = 'N';
                double alpha = 1;
                double beta = 0;
                F77_NAME(dgemm)(&transa, &transb, &num_points32[0], &num_points32[1], &num_dims32, &alpha, mem.mask[0], &num_dims32, mem.m[1], &num_dims32, &beta, mem.vary_n, &num_points32[0] FCONE FCONE);
                check_interrupt();
                progress.report(1);
            }

            //     vary_n <- vary_n - (s_y * s_y) / n
            for (uint64_t i = 0; i < res_size; ++i)
                mem.vary_n[i] -= mem.s_y[i] * mem.s_y[i] / mem.n[i];
            check_interrupt();

            //     res <- cov_n / sqrt(varx_n * vary_n)
            // NA for a column that is constant in every pair: it has no correlations
            vector<bool> constant_col[2] = { constant_cols(_x), constant_cols(_y) };

            for (uint64_t ipoint2 = 0, i = 0; ipoint2 < num_points[1]; ++ipoint2) {
                for (uint64_t ipoint1 = 0; ipoint1 < num_points[0]; ++ipoint1, ++i) {
                    if (constant_col[0][ipoint1] || constant_col[1][ipoint2])
                        res[i] = NA_REAL;
                    else
                        res[i] = mem.cov_n[i] / sqrt(mem.varx_n[i] * mem.vary_n[i]);
                }
            }
            check_interrupt();
            progress.report_last();
        } else {
            // The columns are standardized, so that the correlation matrix is t(zx) %*% zy. The columns that have no
            // correlations (a missing value or constant values) are flagged and their correlations are NA.
            ProgressReporter progress;
            progress.init(3, 1);

            SEXP rmat[2] = { _x, _y };
            vector<bool> flagged[2];

            standardize_cols(rmat, 2, spearman, z, z_sizeof, flagged);
            check_interrupt();
            progress.report(1);

            //     res <- t(zx) %*% zy
            {
                char transa = 'T';
                char transb = 'N';
                double alpha = 1;
                double beta = 0;
                F77_NAME(dgemm)
                (&transa, &transb, &num_points32[0], &num_points32[1], &num_dims32, &alpha, z[0], &num_dims32, z[1], &num_dims32, &beta, res, &num_points32[0] FCONE FCONE);
            }
            unmap_shared(z, z_sizeof);
            check_interrupt();
            progress.report(1);

            vector<uint64_t> flagged_points;

            for (uint64_t ipoint = 0; ipoint < num_points[0]; ++ipoint) {
                if (flagged[0][ipoint])
                    flagged_points.push_back(ipoint);
            }

            for (uint64_t ipoint = 0; ipoint < num_points[1]; ++ipoint) {
                double *col = res + ipoint * num_points[0];

                if (flagged[1][ipoint])
                    fill(col, col + num_points[0], NA_REAL);
                else {
                    for (auto ipoint1 : flagged_points)
                        col[ipoint1] = NA_REAL;
                }
            }
            check_interrupt();
            progress.report_last();
        }
        vdebug("END BLAS COR\n");

        // assemble the answer
        SEXP rold_dimnames[2] = { Rf_getAttrib(_x, R_DimNamesSymbol), Rf_getAttrib(_y, R_DimNamesSymbol) };
        SEXP rold_colnames[2] = {
            !Rf_isNull(rold_dimnames[0]) && Rf_xlength(rold_dimnames[0]) == 2 ? VECTOR_ELT(rold_dimnames[0], 1) : R_NilValue,
            !Rf_isNull(rold_dimnames[1]) && Rf_xlength(rold_dimnames[1]) == 2 ? VECTOR_ELT(rold_dimnames[1], 1) : R_NilValue
        };

        if (tidy) {
            enum { COL1, COL2, COR, NUM_COLS };
            const char *COL_NAMES[NUM_COLS] = { "col1", "col2", "cor" };

            rprotect(answer = RSaneAllocVector(VECSXP, NUM_COLS));

            uint64_t answer_size = 0;

            for (uint64_t i = 0; i < res_size; ++i) {
                if (!std::isnan(res[i]) && fabs(res[i]) >= threshold)
                    ++answer_size;
            }

            SEXP rcol1, rcol2, rcor, rrownames, rcolnames;

            rprotect(rcol1 = RSaneAllocVector(INTSXP, answer_size));
            rprotect(rcol2 = RSaneAllocVector(INTSXP, answer_size));
            rprotect(rcor = RSaneAllocVector(REALSXP, answer_size));
            rprotect(rcolnames = RSaneAllocVector(STRSXP, NUM_COLS));
            rprotect(rrownames = RSaneAllocVector(INTSXP, answer_size));

            for (int i = 0; i < NUM_COLS; i++)
                SET_STRING_ELT(rcolnames, i, Rf_mkChar(COL_NAMES[i]));

            if (answer_size) {
                uint64_t idx_answer = 0;
                for (uint64_t icol1 = 0; icol1 < num_points[0]; ++icol1) {
                    uint64_t idx_m = icol1;
                    for (uint64_t icol2 = 0; icol2 < num_points[1]; ++icol2) {
                        if (!std::isnan(res[idx_m]) && fabs(res[idx_m]) >= threshold) {
                            INTEGER(rcol1)[idx_answer] = icol1 + 1;
                            INTEGER(rcol2)[idx_answer] = icol2 + 1;
                            REAL(rcor)[idx_answer] = res[idx_m];
                            INTEGER(rrownames)[idx_answer] = idx_answer + 1;
                            ++idx_answer;
                        }
                        idx_m += num_points[0];
                    }
                }
            }

            if (rold_colnames[0] != R_NilValue) {
                Rf_setAttrib(rcol1, R_LevelsSymbol, rold_colnames[0]);
                Rf_setAttrib(rcol1, R_ClassSymbol, Rf_mkString("factor"));
            }
            if (rold_colnames[1] != R_NilValue) {
                Rf_setAttrib(rcol2, R_LevelsSymbol, rold_colnames[1]);
                Rf_setAttrib(rcol2, R_ClassSymbol, Rf_mkString("factor"));
            }

            SET_VECTOR_ELT(answer, COL1, rcol1);
            SET_VECTOR_ELT(answer, COL2, rcol2);
            SET_VECTOR_ELT(answer, COR, rcor);

            Rf_setAttrib(answer, R_NamesSymbol, rcolnames);
            Rf_setAttrib(answer, R_ClassSymbol, Rf_mkString("data.frame"));
            Rf_setAttrib(answer, R_RowNamesSymbol, rrownames);
        } else {
            SEXP dim;
            SEXP dimnames;

            if (pairwise_nan) {
                for (uint64_t i = 0; i < res_size; ++i) {
                    if (std::isnan(res[i]))
                        res[i] = NA_REAL;
                }
            }

            rprotect(dim = RSaneAllocVector(INTSXP, 2));
            INTEGER(dim)[0] = num_points[0];
            INTEGER(dim)[1] = num_points[1];
            Rf_setAttrib(answer, R_DimSymbol, dim);

            rprotect(dimnames = RSaneAllocVector(VECSXP, 2));
            SET_VECTOR_ELT(dimnames, 0, rold_colnames[0]);
            SET_VECTOR_ELT(dimnames, 1, rold_colnames[1]);
            Rf_setAttrib(answer, R_DimNamesSymbol, dimnames);
        }
    } catch (TGLException &e) {
        unmap_shared(z, z_sizeof);
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        unmap_shared(z, z_sizeof);
        rerror("Out of memory");
    }

    unmap_shared(z, z_sizeof);
    rreturn(answer);
}

}

