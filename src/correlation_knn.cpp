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

SEXP tgs_cor_knn_blas(SEXP _x, SEXP _knn, SEXP _pairwise_complete_obs, SEXP _spearman, SEXP _threshold, SEXP _envir)
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
        //       unsigned - node with the highest correlation
        //       double   - highest correlation
        //       unsigned - node with the second highest correlation
        //       ...
        //    [col 1]
        //       ...
        size_t res_record_sizeof = sizeof(unsigned) + sizeof(double);
        res_sizeof = num_cols * knn * res_record_sizeof;
        res = (char *)mmap(NULL, res_sizeof, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);

        if (res == (char *)MAP_FAILED)
            verror("Failed to allocate shared memory: %s", strerror(errno));

        int num_cores = max(1, (int)sysconf(_SC_NPROCESSORS_ONLN));
        int num_processes = (int)min(num_cols, (size_t)num_cores);
        double num_cols4process = num_cols / (double)num_processes;

        ProgressReporter progress;
        progress.init(num_cols * num_cols, 1);

        vdebug("Num cores: %d, num_processes: %d\n", num_cores, num_processes);
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
                        unsigned idx = pcors[i] - &cors.front();
                        *(unsigned *)(res + offset) = idx;
                        offset += sizeof(unsigned);
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
                if (*(double *)(res + offset + sizeof(unsigned)) == -2)
                    break;

                ++answer_size;
                offset += res_record_sizeof;
            }
        }

        rprotect(answer = allocVector(VECSXP, NUM_COLS));

        SEXP rcol1, rcol2, rcor, rrank, rrownames, rcolnames;
        SEXP rold_dimnames = getAttrib(_x, R_DimNamesSymbol);
        SEXP rold_colnames = !isNull(rold_dimnames) && xlength(rold_dimnames) == 2 ? VECTOR_ELT(rold_dimnames, 1) : R_NilValue;

        SET_VECTOR_ELT(answer, COL1, (rcol1 = allocVector(INTSXP, answer_size)));
        SET_VECTOR_ELT(answer, COL2, (rcol2 = allocVector(INTSXP, answer_size)));
        SET_VECTOR_ELT(answer, COR, (rcor = allocVector(REALSXP, answer_size)));
        SET_VECTOR_ELT(answer, RANK, (rrank = allocVector(INTSXP, answer_size)));

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

        if (answer_size) {
            size_t row = 0;
            for (size_t icol = 0; icol < num_cols; ++icol) {
                size_t offset = icol * res_record_sizeof * knn;

                for (size_t i = 0; i < knn; ++i) {
                    if (*(double *)(res + offset + sizeof(unsigned)) == -2)
                        break;

                    INTEGER(rcol1)[row] = icol + 1;
                    INTEGER(rcol2)[row] = *(unsigned *)(res + offset) + 1;
                    offset += sizeof(unsigned);
                    REAL(rcor)[row] = *(double *)(res + offset);
                    offset += sizeof(double);
                    INTEGER(rrank)[row] = i + 1;
                    INTEGER(rrownames)[row] = row + 1;
                    ++row;
                }
            }
        }
    } catch (TGLException &e) {
        if (!TGStat::is_kid() && res != (char *)MAP_FAILED) {
            munmap(res, res_sizeof);
            res = (char *)MAP_FAILED;
        }
        rerror("%s", e.msg());
    }

    if (!TGStat::is_kid() && res != (char *)MAP_FAILED) {
        munmap(res, res_sizeof);
        res = (char *)MAP_FAILED;
    }
    rreturn(answer);
}

}

#ifdef AAAA



SEXP tgs_cor_knn_blas(SEXP _x, SEXP _knn, SEXP _pairwise_complete_obs, SEXP _spearman, SEXP _threshold, SEXP _envir)
{
    SEXP answer = R_NilValue;

    try {
        struct Mem {
            double *m{NULL};
            double *m2{NULL};
            double *mask{NULL};
            ~Mem() { free(m); free(m2); free(mask); }
        } mem;

        TGStat tgstat(_envir);

        if ((!isReal(_x) && !isInteger(_x)) || xlength(_x) < 1)
            verror("\"x\" argument must be a matrix of numeric values");

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
        int num_rows32 = (int)num_rows;
        int num_cols32 = (int)num_cols;
        size_t num_vals = num_rows * num_cols;
vector<double> vals;
bool nan_in_vals = false;
vector<bool> nan_in_col(num_cols, false);
vector<double> col_vals1(num_rows);
vector<double> col_vals2(num_rows);
double *pcol_vals1 = &col_vals1.front();
double *pcol_vals2 = &col_vals2.front();
vector<double> *col_vals[2] = { &col_vals1, &col_vals2 };

vals.reserve(num_vals);

for (size_t i = 0; i < num_vals; ++i) {
    if ((isReal(_x) && !R_FINITE(REAL(_x)[i])) || (isInteger(_x) && INTEGER(_x)[i] == NA_INTEGER)) {
        nan_in_col[i / num_rows] = true;
        nan_in_vals = true;
        vals.push_back(numeric_limits<double>::quiet_NaN());
    } else
        vals.push_back(isReal(_x) ? REAL(_x)[i] : INTEGER(_x)[i]);
}

        if (num_rows <= 1 || num_cols <= 1)
            verror("\"x\" argument must be a matrix of numeric values");

        vdebug("START BLAS COR\n");
        // some BLAS implementations ask to align double arrays to 64 for improved efficiency
        posix_memalign((void **)&mem.m, 64, sizeof(double) * num_vals);
        posix_memalign((void **)&mem.m2, 64, sizeof(double) * num_vals);
        posix_memalign((void **)&mem.mask, 64, sizeof(double) * num_vals);

        for (size_t i = 0; i < num_vals; ++i) {
            if ((isReal(_x) && !R_FINITE(REAL(_x)[i])) || (isInteger(_x) && INTEGER(_x)[i] == NA_INTEGER)) {
                mem.m[i] = mem.m2[i] = mem.mask[i] = 0.;
//                nan_in_vals = true;
//                nan_in_point[i / num_rows] = true;
            } else {
                double val = isReal(_x) ? REAL(_x)[i] : INTEGER(_x)[i];
                mem.m[i] = val;
                mem.m2[i] = val * val;
                mem.mask[i] = 1.;
            }
        }

        if (spearman && pairwise_complete_obs && nan_in_vals)
            verror("BLAS implementation of tgs_cor does not support spearman with pairwise.complete.obs when x contains NA / NaN / Inf");

        // replace values with ranks if spearman=T
//      if (spearman) {
//          vector<double *> pvals;
//          pvals.reserve(num_vals);
//          for (size_t i = 0; i < num_vals; ++i)
//              pvals.push_back(&mem.m[i]);
//
//          for (size_t ipoint = 0; ipoint < num_cols; ++ipoint) {
//              if (nan_in_point[ipoint])
//                  continue;
//
//              vector<double *>::iterator sival = pvals.begin() + ipoint * num_rows;
//              vector<double *>::iterator eival = sival + num_rows;
//              vector<double *>::iterator last_ival = sival;
//
//              sort(sival, eival, [](double *p1, double *p2) { return *p1 < *p2; });
//              for (auto ival = sival; ; ++ival) {
//                  if (ival == eival || **ival != **last_ival) {
//                      double rank = ((ival - sival) + (last_ival - sival) - 1) / 2. + 1;
//
//                      while (last_ival != ival) {
//                          **last_ival = rank;
//                          mem.m2[*last_ival - mem.m] = rank * rank;
//                          ++last_ival;
//                      }
//
//                      if (ival == eival)
//                          break;
//                  }
//              }
//          }
//      }

        int num_cores = max(1, (int)sysconf(_SC_NPROCESSORS_ONLN));
        int num_processes = (int)min(num_cols, (size_t)num_cores);
        double num_cols4process = num_cols / (double)num_processes;

        ProgressReporter progress;
        progress.init(num_cols * num_cols - num_cols, 1);

        vdebug("Num cores: %d, num_processes: %d\n", num_cores, num_processes);
        TGStat::prepare4multitasking();

        for (int iprocess = 0; iprocess < num_processes; ++iprocess) {
            if (!TGStat::launch_process()) {     // child process
                size_t scol = { (size_t)(iprocess * num_cols4process) };
                size_t ecol = { (size_t)((iprocess + 1) * num_cols4process) };
                size_t itr_idx = 0;

                for (size_t icol1 = 0; icol1 < num_cols; ++icol1) {
                    for (size_t icol2 = scol; icol2 < ecol; ++icol2) {
                        if (icol1 == icol2)
                            continue;

                        size_t idx1 = icol1 * num_rows;
                        size_t idx2 = icol2 * num_rows;
if (spearman) {
} else {
    pcol_vals1 = &vals.front() + idx1;
    pcol_vals2 = &vals.front() + idx2;
}

                        double sum1 = 0;
                        double sum2 = 0;
                        double sum_square1 = 0;
                        double sum_square2 = 0;
                        double mean1, mean2, var1, var2;
                        double sum_xy = 0;
                        size_t num_finite_pairs = 0;

                        for (size_t i = 0; i < num_rows; ++i) {
//double val1 = pcol_vals1[i];
//double val2 = pcol_vals2[i];

//if (!std::isnan(pcol_vals1[i]) && !std::isnan(pcol_vals2[i])) {
                            if (mem.mask[idx1] && mem.mask[idx2]) {
sum1 += pcol_vals1[i];
sum2 += pcol_vals2[i];
sum_square1 += pcol_vals1[i] * pcol_vals1[i];
sum_square2 += pcol_vals2[i] * pcol_vals2[i];
sum_xy += pcol_vals1[i] * pcol_vals2[i];
//                              sum1 += mem.m[idx1];
//                              sum2 += mem.m[idx2];
//                              sum_square1 += mem.m2[idx1];
//                              sum_square2 += mem.m2[idx2];
//                              sum_xy += mem.m[idx1] * mem.m[idx2];  // => sum(X*Y)
                                ++num_finite_pairs;
                                ++idx1;
                                ++idx2;
                            }
                        }

                        if (num_finite_pairs) {
                            mean1 = sum1 / num_finite_pairs;
                            mean2 = sum2 / num_finite_pairs;
                            var1 = sum_square1 / num_finite_pairs - mean1 * mean1;
                            var2 = sum_square2 / num_finite_pairs - mean2 * mean2;

                            // calculate correlation
                            double cor = (sum_xy / num_finite_pairs - mean1 * mean2) / sqrt(var1 * var2);
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
//      {
//              size_t icol1 = 1;
//              size_t icol2 = 2;
//              double sum1 = 0;
//              double sum2 = 0;
//              double sum_square1 = 0;
//              double sum_square2 = 0;
//              double mean1, mean2, stddev1, stddev2;
//              double sum_xy = 0;
//              size_t num_finite_pairs = 0;
//
//              size_t idx1 = num_rows * icol1;
//              size_t idx2 = num_rows * icol2;
//              for (size_t i = 0; i < num_rows; ++i) {
//                  if (mem.mask[idx1] && mem.mask[idx2]) {
//                      sum1 += mem.m[idx1];
//                      sum2 += mem.m[idx2];
//                      sum_square1 += mem.m2[idx1];
//                      sum_square2 += mem.m2[idx2];
//                      sum_xy += mem.m[idx1] * mem.m[idx2];  // => sum(X*Y)
//                      ++num_finite_pairs;
//                  }
//                  ++idx1;
//                  ++idx2;
//              }
//
//              if (num_finite_pairs) {
//                  mean1 = sum1 / num_finite_pairs;
//                  mean2 = sum2 / num_finite_pairs;
//                  stddev1 = sqrt(sum_square1 / num_finite_pairs - mean1 * mean1);
//                  stddev2 = sqrt(sum_square2 / num_finite_pairs - mean2 * mean2);
//
//                  // calculate correlation
//                  double cor = (sum_xy / num_finite_pairs - mean1 * mean2) / (stddev1 * stddev2);
//                  printf("COR1 = %g\n", cor);
//
//                  cor = (sum_xy - sum1 * mean2) / //(stddev1 * stddev2 * num_finite_pairs);
//                      (sqrt((sum_square1 - sum1 * mean1) * (sum_square2 - sum2 * mean2)));
//                  printf("COR2 = %g\n", cor);
//              }
//      }

//
//int inc = 1;
//int xcol=2;
//int ycol=3;
//printf("X: ");
//for (int i = 0; i < num_rows32; ++i)
//printf("%g ", mem.m[num_rows32 * xcol + i]);
//printf("\n");
//printf("Y: ");
//for (int i = 0; i < num_rows32; ++i)
//printf("%g ", mem.m[num_rows32 * ycol + i]);
//printf("\n");
//
//printf("MASK X: ");
//for (int i = 0; i < num_rows32; ++i)
//printf("%g ", mem.mask[num_rows32 * xcol + i]);
//printf("\n");
//printf("MASK Y: ");
//for (int i = 0; i < num_rows32; ++i)
//printf("%g ", mem.mask[num_rows32 * ycol + i]);
//printf("\n");
//
//double N = F77_NAME(ddot)(&num_rows32, mem.mask + num_rows32 * xcol, &inc, mem.mask + num_rows32 * ycol, &inc);
//printf("N=%g\n", N);
//
//double sum_xy = F77_NAME(ddot)(&num_rows32, mem.m + num_rows32 * xcol, &inc, mem.m + num_rows32 * ycol, &inc);
//printf("SUM XY = %g\n", sum_xy);
//
//double sum_x = F77_NAME(ddot)(&num_rows32, mem.m + num_rows32 * xcol, &inc, mem.mask + num_rows32 * ycol, &inc);
//printf("SUM X = %g\n", sum_x);
//
//double avg_x = sum_x / N;
//
//double sum_y = F77_NAME(ddot)(&num_rows32, mem.mask + num_rows32 * xcol, &inc, mem.m + num_rows32 * ycol, &inc);
//printf("SUM Y = %g\n", sum_y);
//
//double avg_y = sum_y / N;
//
//double cov_xy = sum_xy / N - avg_x * avg_y;
//printf("COV = %g (%g / %g - %g * %g)\n", cov_xy, sum_xy, N, avg_x, avg_y);
//
//double sum_x2 = F77_NAME(ddot)(&num_rows32, mem.m2 + num_rows32 * xcol, &inc, mem.mask + num_rows32 * ycol, &inc);
//printf("SUM X^2 = %g\n", sum_x2);
//
//double sum_y2 = F77_NAME(ddot)(&num_rows32, mem.mask + num_rows32 * xcol, &inc, mem.m2 + num_rows32 * ycol, &inc);
//printf("SUM Y^2 = %g\n", sum_y2);
//
//double stddev_x = sqrt(sum_x2 / N - avg_x * avg_x);
//double stddev_y = sqrt(sum_y2 / N - avg_y * avg_y);
//printf("STDDEV X = %g\n", stddev_x);
//printf("STDDEV Y = %g\n", stddev_y);
//
//double cor_xy = cov_xy / (stddev_x * stddev_y);
//printf("COR = %g\n", cor_xy);
    } catch (TGLException &e) {
        rerror("%s", e.msg());
    }

    rreturn(answer);
}

}
#endif
