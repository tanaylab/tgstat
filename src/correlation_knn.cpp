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

//struct Mem {
//    double *m{NULL};
//    double *m2{NULL};
//    double *mask{NULL};
//    double *n{NULL};
//    double *sum_x{NULL};
//    double *sum_y{NULL};
//    double *sum_xy{NULL};
//    double *sum_x2{NULL};
//    double *sum_y2{NULL};
//    ~Mem() { free(m); free(m2); free(mask); free(n); free(sum_x); free(sum_y); free(sum_xy); free(sum_x2); free(sum_y2); }
//} mem;
//
//int num_rows32 = (int)num_rows;
//int num_cols32 = (int)num_cols;
//
//posix_memalign((void **)&mem.m, 64, sizeof(double) * num_vals);
//posix_memalign((void **)&mem.m2, 64, sizeof(double) * num_vals);
//posix_memalign((void **)&mem.mask, 64, sizeof(double) * num_vals);
//posix_memalign((void **)&mem.n, 64, sizeof(double) * num_cols);
//posix_memalign((void **)&mem.sum_x, 64, sizeof(double) * num_cols);
//posix_memalign((void **)&mem.sum_y, 64, sizeof(double) * num_cols);
//posix_memalign((void **)&mem.sum_xy, 64, sizeof(double) * num_cols);
//posix_memalign((void **)&mem.sum_x2, 64, sizeof(double) * num_cols);
//posix_memalign((void **)&mem.sum_y2, 64, sizeof(double) * num_cols);
//
//for (size_t i = 0; i < num_vals; ++i) {
//    if ((isReal(_x) && !R_FINITE(REAL(_x)[i])) || (isInteger(_x) && INTEGER(_x)[i] == NA_INTEGER))
//        mem.m[i] = mem.m2[i] = mem.mask[i] = 0.;
//    else {
//        double val = isReal(_x) ? REAL(_x)[i] : INTEGER(_x)[i];
//        mem.m[i] = val;
//        mem.m2[i] = val * val;
//        mem.mask[i] = 1.;
//    }
//}

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
//num_processes=1;
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
//{
////    char transa = 'T';
//char transa = 'T';
//    double alpha = 1;
//    double beta = 0;
//    int inc = 1;
//
//    // n = t(mask) x mask[,col]
//    F77_NAME(dgemv)(&transa, &num_rows32, &num_cols32, &alpha, mem.mask, &num_rows32, mem.mask + (num_rows * icol1), &inc, &beta, mem.n, &inc);
//
//    // sum_xy = t(m) x m[,col]
//    F77_NAME(dgemv)(&transa, &num_rows32, &num_cols32, &alpha, mem.m, &num_rows32, mem.m + (num_rows * icol1), &inc, &beta, mem.sum_xy, &inc);
//
//    // sum_x = t(m) x mask[,col]
//    F77_NAME(dgemv)(&transa, &num_rows32, &num_cols32, &alpha, mem.m, &num_rows32, mem.mask + (num_rows * icol1), &inc, &beta, mem.sum_x, &inc);
//
//    // sum_y = t(mask) x m[,col]
//    F77_NAME(dgemv)(&transa, &num_rows32, &num_cols32, &alpha, mem.mask, &num_rows32, mem.m + (num_rows * icol1), &inc, &beta, mem.sum_y, &inc);
//
//    // sum_x2 = t(m^2) x mask[,col]
//    F77_NAME(dgemv)(&transa, &num_rows32, &num_cols32, &alpha, mem.m2, &num_rows32, mem.mask + (num_rows * icol1), &inc, &beta, mem.sum_x2, &inc);
//
//    // sum_y2 = t(mask) x m^2[,col]
//    F77_NAME(dgemv)(&transa, &num_rows32, &num_cols32, &alpha, mem.mask, &num_rows32, mem.m2 + (num_rows * icol1), &inc, &beta, mem.sum_y2, &inc);
//}
//
//for (size_t i = 0; i < num_cols; ++i) {
//    double cor = 0;
//
//    if (i != icol1) {
//        double avg_x = mem.sum_x[i] / mem.n[i];
//        double avg_y = mem.sum_y[i] / mem.n[i];
//        double cov_xy = mem.sum_xy[i] / mem.n[i] - avg_x * avg_y;
//        double stddev_x = sqrt(mem.sum_x2[i] / mem.n[i] - avg_x * avg_x);
//        double stddev_y = sqrt(mem.sum_y2[i] / mem.n[i] - avg_y * avg_y);
//        cor = cov_xy / (stddev_x * stddev_y);
//    }
//
//    if (!cor || fabs(cor) < threshold)
//        cor = -2;
//
//    cors[i] = cor;
//    pcors[i] = &cors[i];
//    ++itr_idx;
//    TGStat::itr_idx(itr_idx);
//}

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

