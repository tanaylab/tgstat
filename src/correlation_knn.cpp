#include <algorithm>
#include <cmath>
#include <errno.h>
#include <limits>
#include <queue>
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

SEXP tgs_cor_knn(SEXP _x, SEXP _knn, SEXP _pairwise_complete_obs, SEXP _spearman, SEXP _threshold, SEXP _envir)
{
    SEXP answer = R_NilValue;
    char *res = (char *)MAP_FAILED;
    size_t res_sizeof = 0;

    try {
        TGStat tgstat(_envir);

        if ((!isReal(_x) && !isInteger(_x)) || xlength(_x) < 1)
            verror("\"x\" argument must be a matrix of numeric values");

        if ((!isReal(_knn) && !isInteger(_knn)) || xlength(_knn) != 1 || asReal(_knn) < 1 || asReal(_knn) != (double)asInteger(_knn))
            verror("\"knn\" argument must be a positive integer");

        if (!isLogical(_pairwise_complete_obs) || xlength(_pairwise_complete_obs) != 1)
            verror("\"pairwise.complete.obs\" argument must be a logical value");

        if (!isLogical(_spearman) || xlength(_spearman) != 1)
            verror("\"spearman\" argument must be a logical value");

        if ((!isReal(_threshold) && !isInteger(_threshold)) || xlength(_threshold) != 1)
            verror("\"threshold\" argument must be a numeric value");

        SEXP rdim = getAttrib(_x, R_DimSymbol);

        if (!isInteger(rdim) || xlength(rdim) != 2)
            verror("\"x\" argument must be a matrix of numeric values");

        bool pairwise_complete_obs = asLogical(_pairwise_complete_obs);
        bool spearman = asLogical(_spearman);
        double threshold = fabs(asReal(_threshold));
        size_t num_rows = nrows(_x);
        size_t num_cols = ncols(_x);
        size_t knn = min((size_t)asInteger(_knn), num_cols - 1);

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
            if ((isReal(_x) && !R_FINITE(REAL(_x)[i])) || (isInteger(_x) && INTEGER(_x)[i] == NA_INTEGER)) {
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

            for (size_t icol = 0; icol < num_cols; ++icol) {
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

        for (size_t irow = 0; irow < num_rows; ++irow) {
            size_t idx = irow;
            for (size_t icol = 0; icol < num_cols; ++icol) {
                if (!nan_in_col[icol]) {
                    sums[icol] += vals[idx];
                    sums_square[icol] += vals[idx] * vals[idx];
                }
                idx += num_rows;
            }
        }

        for (size_t icol = 0; icol < num_cols; ++icol) {
            if (!nan_in_col[icol]) {
                means[icol] = sums[icol] / num_rows;

                // we are calaculating standard deviation:
                // sqrt(sum((x-mean)^2) / N)) = sqrt(sum(x^2) / N - mean^2)
                stddevs[icol] = sqrt(sums_square[icol] / num_rows - means[icol] * means[icol]);
            }
        }

        vdebug("Allocating shared memory for results\n");

        // Shared memory structure:
        //    [col 0]
        //       int64_t - node with the highest correlation. int64_t used to make the following double aligned to 8.
        //       double   - highest correlation
        //       int64_t - node with the second highest correlation
        //       ...
        //    [col 1]
        //       ...
        size_t res_record_sizeof = sizeof(int64_t) + sizeof(double);
        res_sizeof = num_cols * knn * res_record_sizeof;
        res = (char *)mmap(NULL, res_sizeof, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);

        if (res == (char *)MAP_FAILED)
            verror("Failed to allocate shared memory: %s", strerror(errno));

        int num_processes = (int)min(num_cols, (size_t)g_tgstat->num_processes());
        double num_cols4process = num_cols / (double)num_processes;

        ProgressReporter progress;
        progress.init(num_cols * num_cols, 1);

        vdebug("num_processes: %d\n", num_processes);
        TGStat::prepare4multitasking();

        for (int iprocess = 0; iprocess < num_processes; ++iprocess) {
            if (!TGStat::launch_process()) {     // child process
                size_t scol = iprocess * num_cols4process;
                size_t ecol = (iprocess + 1) * num_cols4process;
                size_t itr_idx = 0;
                vector<double> cors(num_cols);
                vector<double *> pcors(num_cols);

                for (size_t icol1 = scol; icol1 < ecol; ++icol1) {
                    for (size_t icol2 = 0; icol2 < num_cols; ++icol2) {
                        double cor = 0;

                        if (icol1 != icol2) {
                            if (nan_in_vals && pairwise_complete_obs) {
                                size_t idx1 = icol1 * num_rows;
                                size_t idx2 = icol2 * num_rows;
                                double sum1 = 0;
                                double sum2 = 0;
                                double sum_square1 = 0;
                                double sum_square2 = 0;
                                double mean1, mean2;

                                if (spearman) {
                                    size_t indices[2] = { idx1, idx2 };
                                    vector<double *>::iterator sivals[2] = { pvals.begin() + idx1, pvals.begin() + idx2 };
                                    vector<double *>::iterator eivals[2] = { sivals[0] + num_rows, sivals[1] + num_rows };
                                    double *spvals[2] = { &vals.front() + idx1, &vals.front() + idx2 };

                                    for (int i = 0; i < 2; ++i) {
                                        // the fastest way to set all members of col_vals to NaN
                                        memcpy(&col_vals[i]->front(), &nan_col.front(), num_rows * sizeof(double));

                                        auto last_ival = sivals[i];
                                        size_t num_preceeding_vals = 0;
                                        size_t last_num_preceeding_vals = 0;

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

                                size_t num_finite_pairs = 0;
                                for (size_t i = 0; i < num_rows; ++i) {
                                    double val1 = pcol_vals1[i];
                                    double val2 = pcol_vals2[i];

                                    if (!std::isnan(val1) && !std::isnan(val2)) {
                                        sum1 += val1;
                                        sum2 += val2;
                                        sum_square1 += val1 * val1;
                                        sum_square2 += val2 * val2;
                                        cor += val1 * val2;  // => sum(X*Y)
                                        ++num_finite_pairs;
                                    }
                                }

                                if (num_finite_pairs) {
                                    mean1 = sum1 / num_finite_pairs;
                                    mean2 = sum2 / num_finite_pairs;
                                    double var1 = sum_square1 / num_finite_pairs - mean1 * mean1;
                                    double var2 = sum_square2 / num_finite_pairs - mean2 * mean2;

                                    // calculate correlation
                                    cor /= num_finite_pairs;      // => mean(X*Y)
                                    cor -= mean1 * mean2;         // => covariance(X,Y)
                                    cor /= sqrt(var1 * var2);     // => correlation(X,Y)
                                }
                            } else if (!nan_in_col[icol1] && !nan_in_col[icol2]) {
                                size_t idx1 = icol1 * num_rows;
                                size_t idx2 = icol2 * num_rows;
                                size_t end_idx1 = idx1 + num_rows;

                                while (idx1 < end_idx1)
                                    cor += vals[idx1++] * vals[idx2++];  // => sum(X*Y)

                                cor /= num_rows;                         // => mean(X*Y)
                                cor -= means[icol1] * means[icol2];      // => covariance(X,Y)
                                cor /= stddevs[icol1] * stddevs[icol2];  // => correlation(X,Y)
                            }
                        }

                        if (!cor || fabs(cor) < threshold)
                            cor = -2;

                        cors[icol2] = cor;
                        pcors[icol2] = &cors[icol2];

                        ++itr_idx;
                        TGStat::itr_idx(itr_idx);
                    }

                    partial_sort(pcors.begin(), pcors.begin() + knn, pcors.end(),
                                 [](double *pcor1, double *pcor2) { return *pcor1 > *pcor2 || (*pcor1 == *pcor2 && pcor1 < pcor2); });

                    // pack the results into shared memory
                    size_t offset = icol1 * res_record_sizeof * knn;
                    for (size_t i = 0; i < knn; ++i) {
                        int64_t idx = pcors[i] - &cors.front();
                        *(int64_t *)(res + offset) = idx;
                        offset += sizeof(int64_t);
                        *(double *)(res + offset) = cors[idx];
                        offset += sizeof(double);

                        if (*pcors[i] == -2)
                            break;
                    }
                }
                rexit();
            }
        }

        while (TGStat::wait_for_kids(3000))
            progress.report(TGStat::itr_idx_sum() - progress.get_elapsed_steps());

        progress.report_last();

        // assemble the answer
        enum { COL1, COL2, COR, RANK, NUM_COLS };
        const char *COL_NAMES[NUM_COLS] = { "col1", "col2", "cor", "rank" };

        size_t answer_size = 0;

        for (size_t icol = 0; icol < num_cols; ++icol) {
            size_t offset = icol * res_record_sizeof * knn;

            for (size_t i = 0; i < knn; ++i) {
                if (*(double *)(res + offset + sizeof(int64_t)) == -2)
                    break;

                ++answer_size;
                offset += res_record_sizeof;
            }
        }

        rprotect(answer = RSaneAllocVector(VECSXP, NUM_COLS));

        SEXP rcol1, rcol2, rcor, rrank, rrownames, rcolnames;
        SEXP rold_dimnames = getAttrib(_x, R_DimNamesSymbol);
        SEXP rold_colnames = !isNull(rold_dimnames) && xlength(rold_dimnames) == 2 ? VECTOR_ELT(rold_dimnames, 1) : R_NilValue;

        rprotect(rcol1 = RSaneAllocVector(INTSXP, answer_size));
        rprotect(rcol2 = RSaneAllocVector(INTSXP, answer_size));
        rprotect(rcor = RSaneAllocVector(REALSXP, answer_size));
        rprotect(rrank = RSaneAllocVector(INTSXP, answer_size));
        rprotect(rcolnames = RSaneAllocVector(STRSXP, NUM_COLS));
        rprotect(rrownames = RSaneAllocVector(INTSXP, answer_size));

        for (int i = 0; i < NUM_COLS; i++)
            SET_STRING_ELT(rcolnames, i, mkChar(COL_NAMES[i]));

        if (answer_size) {
            size_t row = 0;
            for (size_t icol = 0; icol < num_cols; ++icol) {
                size_t offset = icol * res_record_sizeof * knn;

                for (size_t i = 0; i < knn; ++i) {
                    if (*(double *)(res + offset + sizeof(int64_t)) == -2)
                        break;

                    INTEGER(rcol1)[row] = icol + 1;
                    INTEGER(rcol2)[row] = *(int64_t *)(res + offset) + 1;
                    offset += sizeof(int64_t);
                    REAL(rcor)[row] = *(double *)(res + offset);
                    offset += sizeof(double);
                    INTEGER(rrank)[row] = i + 1;
                    INTEGER(rrownames)[row] = row + 1;
                    ++row;
                }
            }
        }

        if (rold_colnames != R_NilValue) {
            setAttrib(rcol1, R_LevelsSymbol, rold_colnames);
            setAttrib(rcol1, R_ClassSymbol, mkString("factor"));
            setAttrib(rcol2, R_LevelsSymbol, rold_colnames);
            setAttrib(rcol2, R_ClassSymbol, mkString("factor"));
        }

        SET_VECTOR_ELT(answer, COL1, rcol1);
        SET_VECTOR_ELT(answer, COL2, rcol2);
        SET_VECTOR_ELT(answer, COR, rcor);
        SET_VECTOR_ELT(answer, RANK, rrank);

        setAttrib(answer, R_NamesSymbol, rcolnames);
        setAttrib(answer, R_ClassSymbol, mkString("data.frame"));
        setAttrib(answer, R_RowNamesSymbol, rrownames);
    } catch (TGLException &e) {
        if (!TGStat::is_kid() && res != (char *)MAP_FAILED) {
            munmap(res, res_sizeof);
            res = (char *)MAP_FAILED;
        }
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

    if (!TGStat::is_kid() && res != (char *)MAP_FAILED) {
        munmap(res, res_sizeof);
        res = (char *)MAP_FAILED;
    }
    rreturn(answer);
}

SEXP tgs_cross_cor_knn(SEXP _x, SEXP _y, SEXP _knn, SEXP _pairwise_complete_obs, SEXP _spearman, SEXP _threshold, SEXP _envir)
{
    SEXP answer = R_NilValue;
    char *res = (char *)MAP_FAILED;
    size_t res_sizeof = 0;

    try {
        TGStat tgstat(_envir);

        if ((!isReal(_x) && !isInteger(_x)) || xlength(_x) < 1)
            verror("\"x\" argument must be a matrix of numeric values");

        if ((!isReal(_y) && !isInteger(_y)) || xlength(_y) < 1)
            verror("\"y\" argument must be a matrix of numeric values");

        if ((!isReal(_knn) && !isInteger(_knn)) || xlength(_knn) != 1 || asReal(_knn) < 1 || asReal(_knn) != (double)asInteger(_knn))
            verror("\"knn\" argument must be a positive integer");

        if (!isLogical(_pairwise_complete_obs) || xlength(_pairwise_complete_obs) != 1)
            verror("\"pairwise.complete.obs\" argument must be a logical value");

        if (!isLogical(_spearman) || xlength(_spearman) != 1)
            verror("\"spearman\" argument must be a logical value");

        if ((!isReal(_threshold) && !isInteger(_threshold)) || xlength(_threshold) != 1)
            verror("\"threshold\" argument must be a numeric value");

        SEXP rdim1 = getAttrib(_x, R_DimSymbol);
        SEXP rdim2 = getAttrib(_y, R_DimSymbol);

        if (!isInteger(rdim1) || xlength(rdim1) != 2)
            verror("\"x\" argument must be a matrix of numeric values");

        if (!isInteger(rdim2) || xlength(rdim2) != 2)
            verror("\"y\" argument must be a matrix of numeric values");

        if (nrows(_x) != nrows(_y))
            verror("\"x\" and \"y\" matrices have different number of rows");

        bool pairwise_complete_obs = asLogical(_pairwise_complete_obs);
        bool spearman = asLogical(_spearman);
        double threshold = fabs(asReal(_threshold));
        size_t num_rows = nrows(_x);
        size_t num_cols[2] = { (size_t)ncols(_x), (size_t)ncols(_y) };
        size_t knn = min((size_t)asInteger(_knn), num_cols[1]);

        if (num_rows <= 1 || num_cols[0] <= 1)
            verror("\"x\" argument must be a matrix of numeric values");

        if (num_cols[1] <= 1)
            verror("\"y\" argument must be a matrix of numeric values");

        size_t num_vals[2] = { num_rows * num_cols[0], num_rows * num_cols[1] };
        bool nan_in_vals = false;
        vector<bool> nan_in_col[2];
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

            for (size_t i = 0; i < num_vals[k]; ++i) {
                if ((isReal(x) && !R_FINITE(REAL(x)[i])) || (isInteger(x) && INTEGER(x)[i] == NA_INTEGER)) {
                    nan_in_col[k][i / num_rows] = true;
                    nan_in_vals = true;
                    vals[k].push_back(numeric_limits<double>::quiet_NaN());
                } else
                    vals[k].push_back(isReal(x) ? REAL(x)[i] : INTEGER(x)[i]);
            }
        }

        // replace values with ranks if spearman=T
        if (spearman) {
            for (int k = 0; k < 2; ++k) {
                pvals[k].reserve(num_vals[k]);
                for (size_t i = 0; i < num_vals[k]; ++i)
                    pvals[k].push_back(&vals[k][i]);

                for (size_t icol = 0; icol < num_cols[k]; ++icol) {
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
            for (size_t irow = 0; irow < num_rows; ++irow) {
                size_t idx = irow;
                for (size_t icol = 0; icol < num_cols[k]; ++icol) {
                    if (!nan_in_col[k][icol]) {
                        sums[k][icol] += vals[k][idx];
                        sums_square[k][icol] += vals[k][idx] * vals[k][idx];
                    }
                    idx += num_rows;
                }
            }

            for (size_t icol = 0; icol < num_cols[k]; ++icol) {
                if (!nan_in_col[k][icol]) {
                    means[k][icol] = sums[k][icol] / num_rows;

                    // we are calaculating standard deviation:
                    // sqrt(sum((x-mean)^2) / N)) = sqrt(sum(x^2) / N - mean^2)
                    stddevs[k][icol] = sqrt(sums_square[k][icol] / num_rows - means[k][icol] * means[k][icol]);
                }
            }
        }

        vdebug("Allocating shared memory for results\n");

        // Shared memory structure:
        //    [col 0]
        //       int64_t - node with the highest correlation: int64_t used to make the following double aligned to 8
        //       double  - highest correlation
        //       int64_t - node with the second highest correlation
        //       ...
        //    [col 1]
        //       ...
        size_t res_record_sizeof = sizeof(int64_t) + sizeof(double);
        res_sizeof = num_cols[0] * knn * res_record_sizeof;
        res = (char *)mmap(NULL, res_sizeof, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);

        if (res == (char *)MAP_FAILED)
            verror("Failed to allocate shared memory: %s", strerror(errno));

        int num_processes = (int)min(max((size_t)1, num_cols[0] / 10), (size_t)g_tgstat->num_processes());
        double num_cols4process = num_cols[0] / (double)num_processes;

        ProgressReporter progress;
        progress.init(num_cols[0] * num_cols[1], 1);

        vdebug("num_processes: %d\n", num_processes);
        TGStat::prepare4multitasking();

        for (int iprocess = 0; iprocess < num_processes; ++iprocess) {
            if (!TGStat::launch_process()) {     // child process
                size_t scol = iprocess * num_cols4process;
                size_t ecol = (iprocess + 1) * num_cols4process;
                size_t itr_idx = 0;
                vector<double> cors(num_cols[1]);
                vector<double *> pcors(num_cols[1]);

                for (size_t icol1 = scol; icol1 < ecol; ++icol1) {
                    for (size_t icol2 = 0; icol2 < num_cols[1]; ++icol2) {
                        double cor = 0;

                        if (nan_in_vals && pairwise_complete_obs) {
                            size_t idx1 = icol1 * num_rows;
                            size_t idx2 = icol2 * num_rows;
                            double sum1 = 0;
                            double sum2 = 0;
                            double sum_square1 = 0;
                            double sum_square2 = 0;
                            double mean1, mean2;

                            if (spearman) {
                                size_t indices[2] = { idx1, idx2 };
                                vector<double *>::iterator sivals[2] = { pvals[0].begin() + idx1, pvals[1].begin() + idx2 };
                                vector<double *>::iterator eivals[2] = { sivals[0] + num_rows, sivals[1] + num_rows };
                                double *spvals[2] = { &vals[0].front() + idx1, &vals[1].front() + idx2 };

                                for (int i = 0; i < 2; ++i) {
                                    // the fastest way to set all members of col_vals to NaN
                                    memcpy(&col_vals[i]->front(), &nan_col.front(), num_rows * sizeof(double));

                                    auto last_ival = sivals[i];
                                    size_t num_preceeding_vals = 0;
                                    size_t last_num_preceeding_vals = 0;

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

                            size_t num_finite_pairs = 0;
                            for (size_t i = 0; i < num_rows; ++i) {
                                double val1 = pcol_vals1[i];
                                double val2 = pcol_vals2[i];

                                if (!std::isnan(val1) && !std::isnan(val2)) {
                                    sum1 += val1;
                                    sum2 += val2;
                                    sum_square1 += val1 * val1;
                                    sum_square2 += val2 * val2;
                                    cor += val1 * val2;  // => sum(X*Y)
                                    ++num_finite_pairs;
                                }
                            }

                            if (num_finite_pairs) {
                                mean1 = sum1 / num_finite_pairs;
                                mean2 = sum2 / num_finite_pairs;
                                double var1 = sum_square1 / num_finite_pairs - mean1 * mean1;
                                double var2 = sum_square2 / num_finite_pairs - mean2 * mean2;

                                // calculate correlation
                                cor /= num_finite_pairs;      // => mean(X*Y)
                                cor -= mean1 * mean2;         // => covariance(X,Y)
                                cor /= sqrt(var1 * var2);     // => correlation(X,Y)
                            }
                        } else if (!nan_in_col[0][icol1] && !nan_in_col[1][icol2]) {
                            size_t idx1 = icol1 * num_rows;
                            size_t idx2 = icol2 * num_rows;
                            size_t end_idx1 = idx1 + num_rows;

                            while (idx1 < end_idx1)
                                cor += vals[0][idx1++] * vals[1][idx2++];  // => sum(X*Y)

                            cor /= num_rows;                               // => mean(X*Y)
                            cor -= means[0][icol1] * means[1][icol2];      // => covariance(X,Y)
                            cor /= stddevs[0][icol1] * stddevs[1][icol2];  // => correlation(X,Y)
                        }

                        if (!cor || fabs(cor) < threshold)
                            cor = -2;

                        cors[icol2] = cor;
                        pcors[icol2] = &cors[icol2];

                        ++itr_idx;
                        TGStat::itr_idx(itr_idx);
                    }

                    partial_sort(pcors.begin(), pcors.begin() + knn, pcors.end(),
                                 [](double *pcor1, double *pcor2) { return *pcor1 > *pcor2 || (*pcor1 == *pcor2 && pcor1 < pcor2); });

                    // pack the results into shared memory
                    size_t offset = icol1 * res_record_sizeof * knn;
                    for (size_t i = 0; i < knn; ++i) {
                        int64_t idx = pcors[i] - &cors.front();
                        *(int64_t *)(res + offset) = idx;
                        offset += sizeof(int64_t);
                        *(double *)(res + offset) = cors[idx];
                        offset += sizeof(double);

                        if (*pcors[i] == -2)
                            break;
                    }
                }
                rexit();
            }
        }

        while (TGStat::wait_for_kids(3000))
            progress.report(TGStat::itr_idx_sum() - progress.get_elapsed_steps());

        progress.report_last();

        // assemble the answer
        enum { COL1, COL2, COR, RANK, NUM_COLS };
        const char *COL_NAMES[NUM_COLS] = { "col1", "col2", "cor", "rank" };

        size_t answer_size = 0;

        for (size_t icol = 0; icol < num_cols[0]; ++icol) {
            size_t offset = icol * res_record_sizeof * knn;

            for (size_t i = 0; i < knn; ++i) {
                if (*(double *)(res + offset + sizeof(int64_t)) == -2)
                    break;

                ++answer_size;
                offset += res_record_sizeof;
            }
        }

        rprotect(answer = RSaneAllocVector(VECSXP, NUM_COLS));

        SEXP rcol1, rcol2, rcor, rrank, rrownames, rcolnames;
        SEXP rold_dimnames[2] = { getAttrib(_x, R_DimNamesSymbol), getAttrib(_y, R_DimNamesSymbol) };
        SEXP rold_colnames[2] = {
            !isNull(rold_dimnames[0]) && xlength(rold_dimnames[0]) == 2 ? VECTOR_ELT(rold_dimnames[0], 1) : R_NilValue,
            !isNull(rold_dimnames[1]) && xlength(rold_dimnames[1]) == 2 ? VECTOR_ELT(rold_dimnames[1], 1) : R_NilValue
        };

        rprotect(rcol1 = RSaneAllocVector(INTSXP, answer_size));
        rprotect(rcol2 = RSaneAllocVector(INTSXP, answer_size));
        rprotect(rcor = RSaneAllocVector(REALSXP, answer_size));
        rprotect(rrank = RSaneAllocVector(INTSXP, answer_size));
        rprotect(rcolnames = RSaneAllocVector(STRSXP, NUM_COLS));
        rprotect(rrownames = RSaneAllocVector(INTSXP, answer_size));

        for (int i = 0; i < NUM_COLS; i++)
            SET_STRING_ELT(rcolnames, i, mkChar(COL_NAMES[i]));

        if (answer_size) {
            size_t row = 0;
            for (size_t icol = 0; icol < num_cols[0]; ++icol) {
                size_t offset = icol * res_record_sizeof * knn;

                for (size_t i = 0; i < knn; ++i) {
                    if (*(double *)(res + offset + sizeof(int64_t)) == -2)
                        break;

                    INTEGER(rcol1)[row] = icol + 1;
                    INTEGER(rcol2)[row] = *(int64_t *)(res + offset) + 1;
                    offset += sizeof(int64_t);
                    REAL(rcor)[row] = *(double *)(res + offset);
                    offset += sizeof(double);
                    INTEGER(rrank)[row] = i + 1;
                    INTEGER(rrownames)[row] = row + 1;
                    ++row;
                }
            }
        }

        if (rold_colnames[0] != R_NilValue) {
            setAttrib(rcol1, R_LevelsSymbol, rold_colnames[0]);
            setAttrib(rcol1, R_ClassSymbol, mkString("factor"));
        }
        if (rold_colnames[1] != R_NilValue) {
            setAttrib(rcol2, R_LevelsSymbol, rold_colnames[1]);
            setAttrib(rcol2, R_ClassSymbol, mkString("factor"));
        }

        SET_VECTOR_ELT(answer, COL1, rcol1);
        SET_VECTOR_ELT(answer, COL2, rcol2);
        SET_VECTOR_ELT(answer, COR, rcor);
        SET_VECTOR_ELT(answer, RANK, rrank);

        setAttrib(answer, R_NamesSymbol, rcolnames);
        setAttrib(answer, R_ClassSymbol, mkString("data.frame"));
        setAttrib(answer, R_RowNamesSymbol, rrownames);
    } catch (TGLException &e) {
        if (!TGStat::is_kid() && res != (char *)MAP_FAILED) {
            munmap(res, res_sizeof);
            res = (char *)MAP_FAILED;
        }
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

    if (!TGStat::is_kid() && res != (char *)MAP_FAILED) {
        munmap(res, res_sizeof);
        res = (char *)MAP_FAILED;
    }
    rreturn(answer);
}

}
