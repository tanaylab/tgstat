#include <algorithm>
#include <limits>
#include <sys/mman.h>
#include <unistd.h>

#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>

#ifdef length
#undef length
#endif
#ifdef error
#undef error
#endif

#include "ProgressReporter.h"
#include "tgstat.h"

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
        bool spearman = asLogical(_spearman);
        bool tidy = asLogical(_tidy);
        double threshold = fabs(asReal(_threshold));
        int num_rows = nrows(_x);
        int num_cols = ncols(_x);

        if (num_rows <= 1 || num_cols <= 1)
            verror("\"x\" argument must be a matrix of numeric values");

        size_t num_vals = num_rows * num_cols;
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

        for (size_t i = 0; i < num_vals; ++i) {
            if (isReal(_x) && !R_FINITE(REAL(_x)[i]) || isInteger(_x) && INTEGER(_x)[i] == NA_INTEGER) {
                nan_in_col[i / num_rows] = true;
                nan_in_vals = true;
                vals.push_back(numeric_limits<double>::quiet_NaN());
            } else
                vals.push_back(isReal(_x) ? REAL(_x)[i] : INTEGER(_x)[i]);
        }

        // replace values with ranks if spearman=T
        if (spearman) {
            pvals.reserve(num_vals);
            for (size_t i = 0; i < num_vals; ++i)
                pvals.push_back(&vals[i]);

            for (int icol = 0; icol < num_cols; ++icol) {
                if (nan_in_col[icol] && !pairwise_complete_obs)
                    continue;

                vector<double *>::iterator sival = pvals.begin() + icol * num_rows;
                vector<double *>::iterator eival = sival + num_rows;
                vector<double *>::iterator last_ival = sival;

                if (nan_in_col[icol])
                    sort(sival, eival, [](double *p1, double *p2) { return *p1 < *p2 || !isnan(*p1) && isnan(*p2); });
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

        vdebug("Allocating shared memory for results\n");
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

        vdebug("Num cores: %d, num_processes: %d\n", num_cores, num_processes);
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

                            if (nan_in_vals && pairwise_complete_obs) {
                                size_t idx1 = icol1 * num_rows;
                                size_t idx2 = icol2 * num_rows;
                                double sum1 = 0;
                                double sum2 = 0;
                                double sum_square1 = 0;
                                double sum_square2 = 0;
                                double mean1, mean2, stddev1, stddev2;

                                if (spearman) {
                                    size_t indices[2] = { idx1, idx2 };
                                    vector<double *>::iterator sivals[2] = { pvals.begin() + idx1, pvals.begin() + idx2 };
                                    vector<double *>::iterator eivals[2] = { sivals[0] + num_rows, sivals[1] + num_rows };
                                    double *spvals[2] = { &vals.front() + idx1, &vals.front() + idx2 };

                                    for (int i = 0; i < 2; ++i) {
                                        // the fastest way to set all members of col_vals to NaN
                                        memcpy(&col_vals[i]->front(), &nan_col.front(), num_rows * sizeof(double));

                                        auto last_ival = sivals[i];
                                        int num_preceeding_vals = 0;
                                        int last_num_preceeding_vals = 0;

                                        for (auto ival = sivals[i]; ; ++ival) {
                                            if (ival == eivals[i] || **ival != **last_ival || isnan(**ival)) {
                                                double rank = (num_preceeding_vals + last_num_preceeding_vals - 1) / 2. + 1;
                                                while (last_ival != ival) {
                                                    (*col_vals[i])[*last_ival - spvals[i]] = rank;
                                                    ++last_ival;
                                                }

                                                if (ival == eivals[i] || isnan(**ival))
                                                    break;

                                                last_num_preceeding_vals = num_preceeding_vals;
                                            }

                                            if (!isnan(vals[*ival - spvals[i] + indices[1 - i]]))
                                                ++num_preceeding_vals;
                                        }
                                    }
                                } else {
                                    pcol_vals1 = &vals.front() + idx1;
                                    pcol_vals2 = &vals.front() + idx2;
                                }

                                size_t num_finite_pairs = 0;
                                res[idx] = 0;
                                for (int i = 0; i < num_rows; ++i) {
                                    double val1 = pcol_vals1[i];
                                    double val2 = pcol_vals2[i];

                                    if (!isnan(val1) && !isnan(val2)) {
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
                                    stddev1 = sqrt(sum_square1 / num_finite_pairs - mean1 * mean1);
                                    stddev2 = sqrt(sum_square2 / num_finite_pairs - mean2 * mean2);

                                    // calculate correlation
                                    res[idx] /= num_finite_pairs;          // => mean(X*Y)
                                    res[idx] -= mean1 * mean2;         // => covariance(X,Y)
                                    res[idx] /= stddev1 * stddev2;     // => correlation(X,Y)
                                } else
                                    res[idx] = numeric_limits<double>::quiet_NaN();
                            } else if (!nan_in_col[icol1] && !nan_in_col[icol2]) {
                                size_t idx1 = icol1 * num_rows;
                                size_t idx2 = icol2 * num_rows;
                                size_t end_idx1 = idx1 + num_rows;

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

SEXP tgs_cor_blas(SEXP _x, SEXP _pairwise_complete_obs, SEXP _spearman, SEXP _tidy, SEXP _threshold, SEXP _envir)
{
    SEXP answer = R_NilValue;

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
        bool spearman = asLogical(_spearman);
        bool tidy = asLogical(_tidy);
        double threshold = fabs(asReal(_threshold));
        int num_dims = nrows(_x);
        int num_points = ncols(_x);

        if (num_dims <= 1 || num_points <= 1)
            verror("\"x\" argument must be a matrix of numeric values");

        size_t num_vals = (size_t)num_points * num_dims;
        size_t res_size = (size_t)num_points * num_points;
        bool nan_in_vals = false;
        vector<bool> nan_in_point(num_points, false);

        // some BLAS implementations ask to align double arrays to 64 for improved efficiency
        mem.m = (double *)aligned_alloc(64, sizeof(double) * num_vals);
        mem.mask = (double *)aligned_alloc(64, sizeof(double) * num_vals);
        mem.res = (double *)aligned_alloc(64, sizeof(double) * res_size);

        for (size_t i = 0; i < num_vals; ++i) {
            if (isReal(_x) && !R_FINITE(REAL(_x)[i]) || isInteger(_x) && INTEGER(_x)[i] == NA_INTEGER) {
                mem.m[i] = mem.mask[i] = 0.;
                nan_in_vals = true;
                nan_in_point[i / num_dims] = true;
            } else {
                double val = isReal(_x) ? REAL(_x)[i] : INTEGER(_x)[i];
                mem.m[i] = val;
                mem.mask[i] = 1.;
            }
        }

        if (spearman && pairwise_complete_obs && nan_in_vals)
            verror("BLAS implementation of tgs_cor does not support spearman with pairwise.complete.obs with x contains NA / NaN / Inf");

        // replace values with ranks if spearman=T
        if (spearman) {
            vector<double *> pvals;
            pvals.reserve(num_vals);
            for (size_t i = 0; i < num_vals; ++i)
                pvals.push_back(&mem.m[i]);

            for (int ipoint = 0; ipoint < num_points; ++ipoint) {
                if (nan_in_point[ipoint])
                    continue;

                vector<double *>::iterator sival = pvals.begin() + ipoint * num_dims;
                vector<double *>::iterator eival = sival + num_dims;
                vector<double *>::iterator last_ival = sival;

                sort(sival, eival, [](double *p1, double *p2) { return *p1 < *p2; });
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

        if (pairwise_complete_obs && nan_in_vals) {
            // correlation with pairwise complete obs is calculated as:
            //     n <- t(mask) %*% mask
            //     s_x <- t(mask) %*% m
            //     cov_n <- (t(m) %*% m) - (s_x * t(s_x)) / n
            //     var_n <- (t(mask) %*% (m * m)) - (s_x * s_x) / n
            //     res <- cov_n / sqrt(var_n * t(var_n))

            mem.n = (double *)aligned_alloc(64, sizeof(double) * res_size);
            mem.s_x = (double *)aligned_alloc(64, sizeof(double) * res_size);
            mem.cov_n = (double *)aligned_alloc(64, sizeof(double) * res_size);
            mem.var_n = (double *)aligned_alloc(64, sizeof(double) * res_size);

            ProgressReporter progress;
            progress.init(6, 1);

            //     n <- t(mask) %*% mask
            {
                char uplo = 'L';
                char trans = 'T';
                double alpha = 1;
                double beta = 0;
                F77_NAME(dsyrk)(&uplo, &trans, &num_points, &num_dims, &alpha, mem.mask, &num_dims, &beta, mem.n, &num_points);
                check_interrupt();
                progress.report(1);
            }

            //     s_x <- t(mask) %*% m
            {
                char transa = 'T';
                char transb = 'N';
                double alpha = 1;
                double beta = 0;
                F77_NAME(dgemm)(&transa, &transb, &num_points, &num_points, &num_dims, &alpha, mem.mask, &num_dims, mem.m, &num_dims, &beta, mem.s_x, &num_points);
                check_interrupt();
                progress.report(1);
            }

            //     cov_n <- t(m) %*% m
            {
                char uplo = 'L';
                char trans = 'T';
                double alpha = 1;
                double beta = 0;
                F77_NAME(dsyrk)(&uplo, &trans, &num_points, &num_dims, &alpha, mem.m, &num_dims, &beta, mem.cov_n, &num_points);
                check_interrupt();
                progress.report(1);
            }

            //     cov_n <- cov_n - (s_x * t(s_x)) / n
            for (size_t ipoint1 = 0, idx1 = 0; ipoint1 < num_points; ++ipoint1) {
                idx1 += ipoint1 + 1;
                for (size_t ipoint2 = ipoint1 + 1, idx2 = ipoint2 * num_points + ipoint1; ipoint2 < num_points; ++ipoint2) {
                    mem.cov_n[idx1] -= mem.s_x[idx1] * mem.s_x[idx2] / mem.n[idx1];
                    ++idx1;
                    idx2 += num_points;
                }
            }
            check_interrupt();

            //    m <- m * m
            for (size_t i = 0; i < num_vals; ++i)
                mem.m[i] *= mem.m[i];

            check_interrupt();
            progress.report(1);

            //     var_n <- t(mask) %*% (m * m)
            {
                char transa = 'T';
                char transb = 'N';
                double alpha = 1;
                double beta = 0;
                F77_NAME(dgemm)(&transa, &transb, &num_points, &num_points, &num_dims, &alpha, mem.mask, &num_dims, mem.m, &num_dims, &beta, mem.var_n, &num_points);
                check_interrupt();
                progress.report(1);
            }

            //     var_n <- var_n - (s_x * s_x) / n
            for (size_t ipoint1 = 0, idx1 = 0; ipoint1 < num_points; ++ipoint1) {
                idx1 += ipoint1 + 1;
                for (size_t ipoint2 = ipoint1 + 1, idx2 = ipoint2 * num_points + ipoint1; ipoint2 < num_points; ++ipoint2) {
                    mem.var_n[idx1] -= mem.s_x[idx1] * mem.s_x[idx1] / mem.n[idx1];
                    mem.var_n[idx2] -= mem.s_x[idx2] * mem.s_x[idx2] / mem.n[idx1];
                    ++idx1;
                    idx2 += num_points;
                }
            }
            check_interrupt();

            //     res <- cov_n / sqrt(var_n * t(var_n))
            for (size_t ipoint1 = 0, idx1 = 0; ipoint1 < num_points; ++ipoint1) {
                idx1 += ipoint1 + 1;
                for (size_t ipoint2 = ipoint1 + 1, idx2 = ipoint2 * num_points + ipoint1; ipoint2 < num_points; ++ipoint2) {
                    mem.res[idx1] = mem.res[idx2] = mem.cov_n[idx1] / sqrt(mem.var_n[idx1] * mem.var_n[idx2]);
                    ++idx1;
                    idx2 += num_points;
                }
            }
            check_interrupt();

            for (size_t idx = 0; idx < res_size; idx += num_points + 1)
                mem.res[idx] = 1.;
            check_interrupt();

            progress.report_last();
        } else {
            // correlation without pairwise complete obs or with all finite values requires matrix multiplications (only one):
            //     s_xy <- t(m) %*% m
            //     stddev[i] <- sqrt(sums_square[i] / dims - e(i)^2)
            //     res[i,j] <- (s_xy[i,j] / dims - e(i) * e(j)) / (stddev(i) * stddev(j))

            ProgressReporter progress;
            progress.init(3, 1);

            vector<double> sums(num_points, 0);
            vector<double> sums_square(num_points, 0);
            vector<double> means(num_points, 0);
            vector<double> stddevs(num_points, 0);

            for (int idim = 0; idim < num_dims; ++idim) {
                size_t idx = idim;
                for (int ipoint = 0; ipoint < num_points; ++ipoint) {
                    if (!nan_in_point[ipoint]) {
                        sums[ipoint] += mem.m[idx];
                        sums_square[ipoint] += mem.m[idx] * mem.m[idx];
                    }
                    idx += num_dims;
                }
            }
            check_interrupt();

            for (int ipoint = 0; ipoint < num_points; ++ipoint) {
                if (!nan_in_point[ipoint]) {
                    means[ipoint] = sums[ipoint] / num_dims;

                    // we are calaculating standard deviation:
                    // sqrt(sum((x-mean)^2) / N)) = sqrt(sum(x^2) / N - mean^2)
                    stddevs[ipoint] = sqrt(sums_square[ipoint] / num_dims - means[ipoint] * means[ipoint]);
                }
            }
            check_interrupt();
            progress.report(1);

            //     res <- t(m) %*% m
            {
                char uplo = 'L';
                char trans = 'T';
                double alpha = 1;
                double beta = 0;
                F77_NAME(dsyrk)(&uplo, &trans, &num_points, &num_dims, &alpha, mem.m, &num_dims, &beta, mem.res, &num_points);
                check_interrupt();
                progress.report(1);
            }

            for (size_t ipoint1 = 0, idx1 = 0; ipoint1 < num_points; ++ipoint1) {
                idx1 += ipoint1 + 1;
                for (size_t ipoint2 = ipoint1 + 1, idx2 = ipoint2 * num_points + ipoint1; ipoint2 < num_points; ++ipoint2) {
                    if (nan_in_point[ipoint1] || nan_in_point[ipoint2])
                        mem.res[idx1] = mem.res[idx2] = NA_REAL;
                    else
                        mem.res[idx1] = mem.res[idx2] = (mem.res[idx1] / num_dims - means[ipoint1] * means[ipoint2]) / (stddevs[ipoint1] * stddevs[ipoint2]);
                    ++idx1;
                    idx2 += num_points;
                }
            }
            check_interrupt();

            for (size_t idx = 0; idx < res_size; idx += num_points + 1)
                mem.res[idx] = 1.;
            check_interrupt();

            progress.report_last();
        }

//memcpy(mem.res, mem.var_n, sizeof(double) * res_size);
//{
//SEXP rdims;
//rprotect(answer = allocVector(REALSXP, (size_t)num_points * num_points));
//rprotect(rdims = allocVector(INTSXP, 2));
//INTEGER(rdims)[0] = INTEGER(rdims)[1] = num_points;
//setAttrib(answer, R_DimSymbol, rdims);
//memcpy(REAL(answer), mem.res, Rf_length(answer) * sizeof(REAL(answer)[0]));
//return answer;
//}


        // assemble the answer
        SEXP rold_dimnames = getAttrib(_x, R_DimNamesSymbol);
        SEXP rold_colnames = !isNull(rold_dimnames) && Rf_length(rold_dimnames) == 2 ? VECTOR_ELT(rold_dimnames, 1) : R_NilValue;

        if (tidy) {
            enum { COL1, COL2, COR, NUM_COLS };
            const char *COL_NAMES[NUM_COLS] = { "col1", "col2", "cor" };

            rprotect(answer = allocVector(VECSXP, NUM_COLS));

            size_t answer_size = 0;

            for (int icol1 = 0; icol1 < num_points; ++icol1) {
                size_t idx = icol1;
                for (int icol2 = 0; icol2 < icol1; ++icol2) {
                    if (!isnan(mem.res[idx]) && fabs(mem.res[idx]) >= threshold)
                        ++answer_size;
                    idx += num_points;
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
            for (int icol1 = 0; icol1 < num_points; ++icol1) {
                for (int icol2 = icol1 + 1; icol2 < num_points; ++icol2) {
                    size_t idx = icol1 * num_points + icol2;
                    if (!isnan(mem.res[idx]) && fabs(mem.res[idx]) >= threshold) {
                        INTEGER(rcol1)[i] = icol1 + 1;
                        INTEGER(rcol2)[i] = icol2 + 1;
                        REAL(rcor)[i] = mem.res[idx];
                        INTEGER(rrownames)[i] = i + 1;
                        ++i;
                    }
                }
            }
        } else {
            SEXP dim;

            rprotect(answer = allocVector(REALSXP, res_size));
            memcpy(REAL(answer), mem.res, res_size * sizeof(double));

            rprotect(dim = allocVector(INTSXP, 2));
            INTEGER(dim)[0] = num_points;
            INTEGER(dim)[1] = num_points;
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
        rerror("%s", e.msg());
    }

    rreturn(answer);
}

}


