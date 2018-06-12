#include <algorithm>
#include <limits>
#include <unordered_map>

#include <R.h>
#include <Rinternals.h>

#include "HashFunc.h"

#ifdef length
#undef length
#endif
#ifdef error
#undef error
#endif

#include "ProgressReporter.h"
#include "tgstat.h"

extern "C" {

SEXP tgs_cor_graph(SEXP _ranks, SEXP _knn, SEXP _k_expand, SEXP _k_beta, SEXP _envir)
{
	try {
        TGStat tgstat(_envir);

        int *pcol1;
        int *pcol2;
        int *prank;
        size_t num_ranks;

        if ((!isReal(_k_beta) && !isInteger(_k_beta)) || xlength(_k_beta) != 1)
            verror("\"k_beta\" argument must be a numeric value");

        {
            enum { COL1, COL2, COR, RANK, NUM_COLS };
            const char *COL_NAMES[NUM_COLS] = { "col1", "col2", "cor", "rank" };

            SEXP rnames = getAttrib(_ranks, R_NamesSymbol);

    		if (!isVector(_ranks) || xlength(_ranks) != NUM_COLS || xlength(rnames) != NUM_COLS ||
                strcmp(CHAR(STRING_ELT(rnames, COL1)), COL_NAMES[COL1]) || (!isInteger(VECTOR_ELT(_ranks, COL1)) && !isFactor(VECTOR_ELT(_ranks, COL1))) ||
                strcmp(CHAR(STRING_ELT(rnames, COL2)), COL_NAMES[COL2]) || (!isInteger(VECTOR_ELT(_ranks, COL2)) && !isFactor(VECTOR_ELT(_ranks, COL2))) ||
                xlength(VECTOR_ELT(_ranks, COL2)) != xlength(VECTOR_ELT(_ranks, COL1)) ||
                !isReal(VECTOR_ELT(_ranks, COR)) || xlength(VECTOR_ELT(_ranks, COR)) != xlength(VECTOR_ELT(_ranks, COL1)) ||
                strcmp(CHAR(STRING_ELT(rnames, RANK)), COL_NAMES[RANK]) || !isInteger(VECTOR_ELT(_ranks, RANK)) || xlength(VECTOR_ELT(_ranks, RANK)) != xlength(VECTOR_ELT(_ranks, COL1)))
    			verror("\"ranks\" argument must be in the format that is returned by tgs_knn function");

            pcol1 = INTEGER(VECTOR_ELT(_ranks, COL1));
            pcol2 = INTEGER(VECTOR_ELT(_ranks, COL2));
            prank = INTEGER(VECTOR_ELT(_ranks, RANK));
            num_ranks = xlength(VECTOR_ELT(_ranks, RANK));
        }

        if (!isNull(_knn) && ((!isReal(_knn) && !isInteger(_knn)) || xlength(_knn) != 1))
            verror("\"knn\" argument must be a numeric value");

        if (!isNull(_k_expand) && ((!isReal(_k_expand) && !isInteger(_k_expand)) || xlength(_k_expand) != 1))
            verror("\"k_expand\" argument must be a numeric value");

        double knn_d = isNull(_knn) ? 0 : asReal(_knn);
        double k_expand = asReal(_k_expand);
        double k_beta = asReal(_k_beta);

        if (knn_d < 1)
            verror("\"knn\" argument must be a positive integer");

        if (k_expand <= 0)
            verror("\"k_expand\" argument must be a positive number");

        if (k_beta <= 0)
            verror("\"k_beta\" argument must be a positive number");

        vdebug("Building the graph\n");

        size_t knn = (size_t)knn_d;
        unsigned num_points = 0;
        unordered_map<pair<unsigned, unsigned>, size_t> ij2weight;
        unordered_map<pair<unsigned, unsigned>, size_t> ij2rank;
        size_t max_weight = knn * knn * k_expand;

        vdebug("Reading ranks\n");
        ij2rank.reserve(num_ranks);
        for (size_t i = 0; i < num_ranks; ++i)
            ij2rank[{pcol1[i], pcol2[i]}] = prank[i];

        vdebug("Building edges weights\n");
        ij2weight.reserve(num_ranks * 2);
        for (const auto &r : ij2rank) {
            const auto &ij = r.first;
            num_points = max(num_points, ij.first + 1);
            num_points = max(num_points, ij.second + 1);
            if (ij.first < ij.second) {
                auto itr = ij2rank.find({ij.second, ij.first});
                if (itr != ij2rank.end()) {
                    size_t weight = r.second * itr->second;    // weight = rank[i,j] * rank[j,i]
                    if (weight <= max_weight)
                        ij2weight[ij] = ij2weight[{ij.second, ij.first}] = weight;
                }
            }
        }

        {
            decltype(ij2rank) cleaner;
            ij2rank.swap(cleaner);     // force ij2rank to release memory
        }

        struct Edge {
            unsigned node;
            size_t weight;
            Edge(unsigned _node, unsigned _weight) : node(_node), weight(_weight) {}
            bool operator<(const Edge &o) const { return weight < o.weight || (weight == o.weight && node < o.node); }
        };

        vdebug("Filter out by incoming edges\n");

        // leave max k_beta * knn incoming edges
        vector<vector<Edge>> incoming(num_points);
        for (const auto &w : ij2weight) {
            unsigned i = w.first.first;
            unsigned j = w.first.second;
            incoming[j].push_back(Edge(i, w.second));
        }

        for (auto iedges = incoming.begin(); iedges < incoming.end(); ++iedges) {
            if (iedges->size() > k_beta * knn) {
                partial_sort(iedges->begin(), iedges->begin() + k_beta * knn, iedges->end());
                for (auto iedge = iedges->begin() + k_beta * knn; iedge < iedges->end(); ++iedge)
                    ij2weight.erase({iedge->node, iedges - incoming.begin()});
            }
        }

        {
            decltype(incoming) cleaner;
            incoming.swap(cleaner);     // force incoming to release memory
        }

        vdebug("Filter out by outgoing edges\n");

        // leave max knn outgoing edges and pack the answer
        vector<vector<Edge>> outgoing(num_points);

        for (const auto &w : ij2weight) {
            unsigned i = w.first.first;
            unsigned j = w.first.second;
            outgoing[i].push_back(Edge(j, w.second));
        }

        {
            decltype(ij2weight) cleaner;
            ij2weight.swap(cleaner);     // force ij2weight to release memory
        }

        vdebug("PACKING\n");

        enum { COL1, COL2, WEIGHT, NUM_COLS };
        const char *COL_NAMES[NUM_COLS] = { "col1", "col2", "weight" };

        size_t answer_size = 0;
        SEXP ranswer, rcol1, rcol2, rweight, rrownames, rcolnames, rlevels;

        for (const auto &edges : outgoing)
            answer_size += min(edges.size(), knn);

        rprotect(ranswer = allocVector(VECSXP, NUM_COLS));
        SET_VECTOR_ELT(ranswer, COL1, (rcol1 = allocVector(INTSXP, answer_size)));
        SET_VECTOR_ELT(ranswer, COL2, (rcol2 = allocVector(INTSXP, answer_size)));
        SET_VECTOR_ELT(ranswer, WEIGHT, (rweight = allocVector(REALSXP, answer_size)));

        rlevels = getAttrib(VECTOR_ELT(_ranks, 0), R_LevelsSymbol);
        if (rlevels != R_NilValue) {
            setAttrib(rcol1, R_LevelsSymbol, rlevels);
            setAttrib(rcol1, R_ClassSymbol, mkString("factor"));
        }

        rlevels = getAttrib(VECTOR_ELT(_ranks, 1), R_LevelsSymbol);
        if (rlevels != R_NilValue) {
            setAttrib(rcol2, R_LevelsSymbol, rlevels);
            setAttrib(rcol2, R_ClassSymbol, mkString("factor"));
        }

        setAttrib(ranswer, R_NamesSymbol, (rcolnames = allocVector(STRSXP, NUM_COLS)));
        setAttrib(ranswer, R_ClassSymbol, mkString("data.frame"));
        setAttrib(ranswer, R_RowNamesSymbol, (rrownames = allocVector(INTSXP, answer_size)));

        for (int i = 0; i < NUM_COLS; i++)
            SET_STRING_ELT(rcolnames, i, mkChar(COL_NAMES[i]));

        size_t idx = 0;
        for (auto iedges = outgoing.begin(); iedges < outgoing.end(); ++iedges) {
            double rank = 0;

            if (iedges->size() <= knn)
                sort(iedges->begin(), iedges->end());
            else
                partial_sort(iedges->begin(), iedges->begin() + knn, iedges->end());

            auto iedge_end = iedges->size() <= knn ? iedges->end() : iedges->begin() + knn;

            for (auto iedge = iedges->begin(); iedge < iedge_end; ++iedge) {
                int i = iedges - outgoing.begin();
                int j = iedge->node;
                INTEGER(rcol1)[idx] = i;
                INTEGER(rcol2)[idx] = j;
                REAL(rweight)[idx] = 1. - rank / knn;
                INTEGER(rrownames)[idx] = idx + 1;
                ++idx;
                ++rank;
            }
        }

        vdebug("END\n");

        rreturn(ranswer);
    } catch (TGLException &e) {
		rerror("%s", e.msg());
	}

    rreturn(R_NilValue);
}

}





#ifdef AAAAA











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
    double *res = (double *)MAP_FAILED;
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

        size_t knn = asInteger(_knn);
        bool pairwise_complete_obs = asLogical(_pairwise_complete_obs);
        bool spearman = asLogical(_spearman);
        double threshold = fabs(asReal(_threshold));
        size_t num_rows = nrows(_x);
        size_t num_cols = ncols(_x);

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
        size_t res_size = num_cols * num_cols;
        res_sizeof = sizeof(double) * res_size;
        res = (double *)mmap(NULL, res_sizeof, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);

        if (res == (double *)MAP_FAILED)
            verror("Failed to allocate shared memory: %s", strerror(errno));

        for (size_t i = 0; i < res_size; ++i)
            res[i] = numeric_limits<double>::quiet_NaN();

        int num_cores = max(1, (int)sysconf(_SC_NPROCESSORS_ONLN));
        int num_processes = (int)min(num_cols, (size_t)num_cores);
        double num_cols4process = num_cols / (double)num_processes;

        ProgressReporter progress;
        progress.init(num_cols * (num_cols - 1), 1);

        vdebug("Num cores: %d, num_processes: %d\n", num_cores, num_processes);
        TGStat::prepare4multitasking();

        for (int iprocess = 0; iprocess < num_processes; ++iprocess) {
            if (!TGStat::launch_process()) {     // child process
                size_t scol = iprocess * num_cols4process;
                size_t ecol = (iprocess + 1) * num_cols4process;
                size_t itr_idx = 0;

                for (size_t icol1 = 0; icol1 < num_cols; ++icol1) {
                    for (size_t icol2 = scol; icol2 < ecol; ++icol2) {
                        if (icol1 == icol2)
                            continue;

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
                            } else
                                cor = numeric_limits<double>::quiet_NaN();
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










        bool pairwise_complete_obs = asLogical(_pairwise_complete_obs);
        bool spearman = asLogical(_spearman);
        double threshold = fabs(asReal(_threshold));
        size_t num_rows = nrows(_x);
        size_t num_cols = ncols(_x);
        int num_rows32 = (int)num_rows;
        int num_cols32 = (int)num_cols;

        if (num_rows <= 1 || num_cols <= 1)
            verror("\"x\" argument must be a matrix of numeric values");

        size_t num_vals = num_cols * num_rows;
        bool nan_in_vals = false;
        vector<bool> nan_in_point(num_cols, false);

        vdebug("START BLAS COR\n");
        // some BLAS implementations ask to align double arrays to 64 for improved efficiency
        posix_memalign((void **)&mem.m, 64, sizeof(double) * num_vals);
        posix_memalign((void **)&mem.m2, 64, sizeof(double) * num_vals);
        posix_memalign((void **)&mem.mask, 64, sizeof(double) * num_vals);

        for (size_t i = 0; i < num_vals; ++i) {
            if ((isReal(_x) && !R_FINITE(REAL(_x)[i])) || (isInteger(_x) && INTEGER(_x)[i] == NA_INTEGER)) {
                mem.m[i] = mem.m2[i] = mem.mask[i] = 0.;
                nan_in_vals = true;
                nan_in_point[i / num_rows] = true;
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
        if (spearman) {
            vector<double *> pvals;
            pvals.reserve(num_vals);
            for (size_t i = 0; i < num_vals; ++i)
                pvals.push_back(&mem.m[i]);

            for (size_t ipoint = 0; ipoint < num_cols; ++ipoint) {
                if (nan_in_point[ipoint])
                    continue;

                vector<double *>::iterator sival = pvals.begin() + ipoint * num_rows;
                vector<double *>::iterator eival = sival + num_rows;
                vector<double *>::iterator last_ival = sival;

                sort(sival, eival, [](double *p1, double *p2) { return *p1 < *p2; });
                for (auto ival = sival; ; ++ival) {
                    if (ival == eival || **ival != **last_ival) {
                        double rank = ((ival - sival) + (last_ival - sival) - 1) / 2. + 1;

                        while (last_ival != ival) {
                            **last_ival = rank;
                            mem.m2[*last_ival - mem.m] = rank * rank;
                            ++last_ival;
                        }

                        if (ival == eival)
                            break;
                    }
                }
            }
        }

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

                        double sum1 = 0;
                        double sum2 = 0;
                        double sum_square1 = 0;
                        double sum_square2 = 0;
                        double mean1, mean2, var1, var2;
                        double sum_xy = 0;
                        size_t num_finite_pairs = 0;

                        size_t idx1 = num_rows * icol1;
                        size_t idx2 = num_rows * icol2;
                        for (size_t i = 0; i < num_rows; ++i) {
                            if (mem.mask[idx1] && mem.mask[idx2]) {
                                sum1 += mem.m[idx1];
                                sum2 += mem.m[idx2];
                                sum_square1 += mem.m2[idx1];
                                sum_square2 += mem.m2[idx2];
                                sum_xy += mem.m[idx1] * mem.m[idx2];  // => sum(X*Y)
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
                            double cor = (sum_xy / num_finite_pairs - mean1 * mean2) / sqrt(stddev1 * stddev2);
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
//
//      progress.report_last();
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

        struct Mem {
            double *m{NULL};
            double *m2{NULL};
            double *mask{NULL};
            ~Mem() { free(m); free(m2); free(mask); }
        } mem;

#endif
