#include <algorithm>
#include <cmath>
#include <errno.h>
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

#define d(x, y) res[x + y * num_points]

//struct __attribute__((__packed__)) DistTo {
//    bool operator<(const DistTo &o) const { return dist < o.dist; }
//
//    void set(int _point, double _dist) { point = _point; dist = _dist; }
//
//    int    point;
//    double dist;
//};
//
//class Knns {
//public:
//    Knns(DistTo *_buf, int _num_points, int _k) : buf(_buf), num_points(_num_points), k(_k) {
//        size_t size = (size_t)num_points * k;
//        for (size_t i = 0; i < size; ++i)
//            buf[i].dist = numeric_limits<double>::infinity();
//    }
//
//    const DistTo &kth(size_t point) const { return buf[point * k]; }
//
//    void add(int point1, int point2, double dist) {
//        update(point1, point2, dist);
//        update(point2, point1, dist);
//    }
//
//private:
//    DistTo *buf;
//    int     num_points;
//    int     k;
//
//    void update(int point1, int point2, double dist) {
//        if (dist < kth(point1).dist) {
//            DistTo *p = buf + (size_t)point1 * k;
//            pop_heap(p, p + k);
//            p[k - 1].set(point2, dist);
//            push_heap(p, p + k);
//        }
//    }
//};
//
//SEXP knn_clustering_dist(SEXP _attrs, SEXP _rrownames, vector<double> &vals, bool tidy, double threshold, int k, int num_points, int num_dims,
//                         void *&shm, size_t &shm_sizeof, bool random_seeds)
//{
//    vdebug("Allocating shared memory for results\n");
//    size_t res_size = num_points * num_points;
//    size_t res_sizeof = res_size * sizeof(double);
//    size_t knns_size = num_points * k;
//    size_t knns_sizeof = knns_size * sizeof(DistTo);
//    size_t nns_sizeof = num_points * sizeof(int);
//
//    shm_sizeof = res_sizeof + knns_sizeof + nns_sizeof;
//
//    shm = mmap(NULL, shm_sizeof, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);
//
//    if (shm == (double *)MAP_FAILED)
//        verror("Failed to allocate shared memory: %s", strerror(errno));
//
//    double *res = (double *)shm;
//    DistTo *knns_buf = (DistTo *)((char *)res + res_sizeof);
//    Knns knns(knns_buf, num_points, k);                     // k nearest neighbors for each point: it's a heap, the highest distance is at the beginning
//    int *nns = (int *)((char *)knns_buf + knns_sizeof);     // nearest neighbor for each point
//    vector<int> anchors(num_points);                        // anchor (cluster's seed) for each point
//
//    for (size_t i = 0; i < res_size; ++i)
//        res[i] = numeric_limits<double>::quiet_NaN();
//
//    for (size_t i = 0; i < num_points; ++i)
//        nns[i] = -1;
//
//    int num_anchors = (int)sqrt(num_points);
//    vector<int> seed_candidates;
//    vector<int> seeds;
//
//    seed_candidates.reserve(num_points);
//    seeds.reserve(num_anchors);
//
//    for (int i = 0; i < num_points; ++i)
//        seed_candidates.push_back(i);
//
//    if (random_seeds) {
//        // randomly draw seeds from the existing points
//        for (int i = 0; i < num_anchors; ++i) {
//            int idx = drand48() * seed_candidates.size();
//            seeds.push_back(seed_candidates[idx]);
//            swap(seed_candidates[idx], seed_candidates.back());
//            seed_candidates.pop_back();
//        }
//
//        // calculate the distance between each point and every seed;
//        // calculate knn distance for each seed;
//        // find nearest neighbor for each point
//        for (auto seed : seeds) {
//            for (int ipoint = 0; ipoint < num_points; ++ipoint) {
//                if (seed == ipoint)
//                    continue;
//
//                double dist = 0;
//                double dev;
//                size_t idx1 = seed;
//                size_t idx2 = ipoint;
//
//                for (int idim = 0; idim < num_dims; ++idim) {
//                    dev = vals[idx1] - vals[idx2];
//                    dist += dev * dev;
//                    idx1 += num_points;
//                    idx2 += num_points;
//                }
//
//                dist = sqrt(dist);
//                d(seed, ipoint) = d(ipoint, seed) = dist;
////printf("1 [%d,%d]\n", seed + 1, ipoint + 1);
//                knns.add(seed, ipoint, dist);
//
//                if (nns[ipoint] == (size_t)-1 || dist < d(ipoint, nns[ipoint]))
//                    nns[ipoint] = anchors[ipoint] = seed;
//            }
//        }
//    } else {
//        // first seed is chosen randomly
//        int idx = drand48() * seed_candidates.size();
//        seeds.push_back(seed_candidates[idx]);
//        swap(seed_candidates[idx], seed_candidates.back());
//        seed_candidates.pop_back();
//
//        while (1) {
//            // calculate the distance between each point and every seed;
//            // calculate knn distance for each seed;
//            for (int ipoint = 0; ipoint < num_points; ++ipoint) {
//                if (seeds.back() == ipoint)
//                    continue;
//
//                double dist = 0;
//                double dev;
//                size_t idx1 = seeds.back();
//                size_t idx2 = ipoint;
//
//                for (int idim = 0; idim < num_dims; ++idim) {
//                    dev = vals[idx1] - vals[idx2];
//                    dist += dev * dev;
//                    idx1 += num_points;
//                    idx2 += num_points;
//                }
//
//                dist = sqrt(dist);
//                d(seeds.back(), ipoint) = d(ipoint, seeds.back()) = dist;
////printf("2 [%d,%d]\n", seeds.back() + 1, ipoint + 1);
//                knns.add(seeds.back(), ipoint, dist);
//            }
//
//            // find nearest neighbor for each point
//            for (auto seed : seeds) {
//                for (int ipoint = 0; ipoint < num_points; ++ipoint) {
//                    if (seed != ipoint && (nns[ipoint] == (size_t)-1 || d(seed, ipoint) < d(ipoint, nns[ipoint])))
//                        nns[ipoint] = anchors[ipoint] = seed;
//                }
//            }
//
//            if (seeds.size() == num_anchors)
//                break;
//
//            // choose the seed from 80-90 quantile of the remaining candidates when they are sorted by their distance from the anchor
//            sort(seed_candidates.begin(), seed_candidates.end(), [&](int p1, int p2) { return d(p1, anchors[p1]) < d(p2, anchors[p2]); });
////printf("Seed candidates:\n");
////for (auto sss : seed_candidates)
////printf("\tseed %d, anchor: %d, dist: %g\n", sss, anchors[sss], d(anchors[sss], sss));
//            idx = (drand48() + 8.) * seed_candidates.size() / 10.;
////printf("selected: %d\n", seed_candidates[idx]);
//            seeds.push_back(seed_candidates[idx]);
//            swap(seed_candidates[idx], seed_candidates.back());
//            seed_candidates.pop_back();
//        }
////printf("2 Seed: %d\n", seed);
//    }
//
////printf("Seeds: ");
////for (auto seed : seeds)
////printf("%d ", seed + 1);
////printf("\n\n");
////
////printf("Kth dist: ");
////for (auto seed : seeds)
////printf("%g ", knns.kth(seed).dist);
////printf("\n\n");
////
////printf("Kth point: ");
////for (auto seed : seeds)
////printf("%d ", knns.kth(seed).point + 1);
////printf("\n\n");
////
////printf("Anchors: ");
////for (int i = 0; i < num_points; i++)
////printf("%d ", anchors[i] + 1);
////printf("\n\n");
////
////printf("Distances:\n");
////for (int i = 0; i < num_points; i++) {
////for (int j = 0; j < num_points; j++)
////printf("%g\t", d(i, j));
////printf("\n");
////}
//
//    int num_cores = max(1, (int)sysconf(_SC_NPROCESSORS_ONLN));
//    int num_processes = min(num_points, num_cores);
////num_processes = 1;
//    double num_row4process = num_points / (double)num_processes;
//
//    ProgressReporter progress;
//    progress.init(num_points * (num_points - 1), 1);
//
//    vdebug("Num cores: %d, num_processes: %d\n", num_cores, num_processes);
//    TGStat::prepare4multitasking();
//
//    for (int iprocess = 0; iprocess < num_processes; ++iprocess) {
//        if (!TGStat::launch_process()) {     // child process
//            int spoint = (int)num_row4process * iprocess;
//            int epoint = (int)num_row4process * (iprocess + 1);
//            size_t itr_idx = 0;
//
//            for (int ipoint1 = 0; ipoint1 < num_points; ++ipoint1) {
//                for (int ipoint2 = spoint; ipoint2 < epoint; ++ipoint2) {
//                    if (ipoint1 == ipoint2)
//                        continue;
//
//                    bool do_dist_calc = std::isnan(d(ipoint1, ipoint2));
//
//                    if (do_dist_calc) {
////                        SemLocker sl(g_tgstat->shm_sem());
//
//                        if (nns[ipoint1] != nns[ipoint2]) {
//                            // Should we calculate the distances between two points x and y?
//                            // Triangle rule:
//                            //    d(x,y) + d(x,z) >= d(y,z)
//                            // Similarly:
//                            //    d(x,y) + d(y,z) >= d(x,z)
//                            // Therefore:
//                            //    d(x,y) >= |d(x,z) - d(y,z)|
//                            //
//                            // On the other hand d(x,y) is relevant only if:
//                            //    d(x,y) <= knn_d(x)
//                            // We don't know knn_d(x) however if we know knn_d of some other point w then we can bound it:
//                            //    knn_d(x) <= knn(w) + d(x,w)
//                            // Therefore:
//                            //    |d(x,z) - d(y,z)| <= d(x,y) <= knn(w) + d(x,w)
//                            //
//                            // So if:
//                            //    |d(x,z) - d(y,z)| > knn(w) + d(x,w)
//                            // we can skip computation of d(x,y)
//                            size_t nn1 = nns[ipoint1];
//                            size_t nn2 = nns[ipoint2];
//                            double dist1 = d(ipoint1, nn1);
//                            double dist2 = 0;
//                            double knn_d_thr1 = knns.kth(anchors[ipoint1]).dist + d(ipoint1, anchors[ipoint1]);
//                            double knn_d_thr2 = std::isinf(knns.kth(nn1).dist) ? numeric_limits<double>::infinity() : knns.kth(nn1).dist + dist1;
//                            double knn_d_thr = min(knn_d_thr1, knn_d_thr2);
//
//                            if (!std::isnan(dist2 = d(ipoint2, nn1)) && fabs(dist1 - dist2) > knn_d_thr ||
//                                !std::isnan(dist1 = d(ipoint1, nn2)) && fabs(dist1 - d(ipoint2, nn2)) > knn_d_thr ||
//                                fabs(d(ipoint1, anchors[ipoint1]) - d(ipoint2, anchors[ipoint1])) > knn_d_thr ||
//                                fabs(d(ipoint1, anchors[ipoint2]) - d(ipoint2, anchors[ipoint2])) > knn_d_thr)
//{
//                                do_dist_calc = false;
////printf("No dist for [%d, %d], nn [%d, %d], knn_d: %g (knn_d(anchor): %g, knn_d(nn): %g), knn point(nn1): %d, knn_dist(nn1): %g\n",
////ipoint1 + 1, ipoint2 + 1, nn1 + 1, nn2 + 1, knn_d_thr, knn_d_thr1, knn_d_thr2, knns.kth(nn1).point + 1, knns.kth(nn1).dist);
////if (!std::isnan(dist2 = d(ipoint2, nn1)) && fabs(dist1 - dist2) > knn_d_thr)
////printf("\tReason: 1\n");
////else if (!std::isnan(dist1 = d(ipoint1, nn2)) && fabs(dist1 - d(ipoint2, nn2)) > knn_d_thr)
////printf("\tReason: 2\n");
////else if (fabs(d(ipoint1, anchors[ipoint1]) - d(ipoint2, anchors[ipoint1])) > knn_d_thr)
////printf("\tReason: 3\n");
////else
////printf("\tReason: 4\n");
////printf("\td(%d,%d)=%g\n", ipoint1 + 1, nn1 + 1, d(ipoint1, nn1));
////printf("\td(%d,%d)=%g\n", ipoint2 + 1, nn1 + 1, d(ipoint2, nn1));
////printf("\td(%d,%d)=%g\n", ipoint1 + 1, nn2 + 1, d(ipoint1, nn2));
////printf("\td(%d,%d)=%g\n", ipoint2 + 1, nn2 + 1, d(ipoint2, nn2));
////printf("\td(%d,%d a1)=%g\n", ipoint1 + 1, anchors[ipoint1] + 1, d(ipoint1, anchors[ipoint1]));
////printf("\td(%d,%d a1)=%g\n", ipoint2 + 1, anchors[ipoint1] + 1, d(ipoint2, anchors[ipoint1]));
////printf("\td(%d,%d a2)=%g\n", ipoint1 + 1, anchors[ipoint2] + 1, d(ipoint1, anchors[ipoint2]));
////printf("\td(%d,%d a2)=%g\n", ipoint2 + 1, anchors[ipoint2] + 1, d(ipoint2, anchors[ipoint2]));
//}
////else
////printf("calc dist [%d,%d]\n", ipoint1 + 1, ipoint2 + 1);
//                        }
//                    }
//
//                    if (do_dist_calc) {
//                        double dist = 0;
//                        double dev;
//                        size_t idx1 = ipoint1;
//                        size_t idx2 = ipoint2;
//
//                        for (int icol = 0; icol < num_dims; ++icol) {
//                            if (!std::isnan(vals[idx1]) && !std::isnan(vals[idx2])) {
//                                dev = vals[idx1] - vals[idx2];
//                                dist += dev * dev;
//                            }
//                            idx1 += num_points;
//                            idx2 += num_points;
//                        }
//
//                        dist = sqrt(dist);
//                        d(ipoint1, ipoint2) = d(ipoint2, ipoint1) = dist;
//
////                        SemLocker sl(g_tgstat->shm_sem());
//
//                        if (d(ipoint1, nns[ipoint1]) > dist)
////{
////printf("1nns[%d] => %d, dist: %g (before: %d, dist %g)\n", ipoint1 + 1, ipoint2 + 1, dist, nns[ipoint1] + 1, d(ipoint1, nns[ipoint1]));
//                            nns[ipoint1] = ipoint2;
////}
//                        if (d(ipoint2, nns[ipoint2]) > dist)
////{
////printf("2nns[%d] => %d, dist: %g (before: %d, dist %g)\n", ipoint2 + 1, ipoint1 + 1, dist, nns[ipoint2] + 1, d(ipoint2, nns[ipoint2]));
//                            nns[ipoint2] = ipoint1;
////}
//
////printf("3 [%d,%d]\n", ipoint1 + 1, ipoint2 + 1);
//                        knns.add(ipoint1, ipoint2, dist);
//                    }
//
//                    TGStat::itr_idx(++itr_idx);
//                }
//            }
//            exit(0);
//        }
//    }
//
//    while (TGStat::wait_for_kids(3000))
//        progress.report(TGStat::itr_idx_sum() - progress.get_elapsed_steps());
//
//    progress.report_last();
//
//    SEXP answer;
//
//    // assemble the answer
//    if (tidy) {
//        enum { ROW1, ROW2, DIST, num_dims };
//        const char *COL_NAMES[num_dims] = { "row1", "row2", "dist" };
//
//        rprotect(answer = allocVector(VECSXP, num_dims));
//
//        size_t answer_size = 0;
//
//        for (int ipoint = 0; ipoint < num_points; ++ipoint) {
//            size_t idx = k * (size_t)ipoint;
//            sort(knns_buf + idx, knns_buf + idx + k);
//            for (int i = idx; i < idx + k; ++i) {
//                if (knns_buf[i].dist > threshold)
//                    break;
//                ++answer_size;
//            }
//        }
//
//        SEXP rcol1, rcol2, rdist, rrownames, rcolnames;
//
//        SET_VECTOR_ELT(answer, ROW1, (rcol1 = allocVector(INTSXP, answer_size)));
//        SET_VECTOR_ELT(answer, ROW2, (rcol2 = allocVector(INTSXP, answer_size)));
//        SET_VECTOR_ELT(answer, DIST, (rdist = allocVector(REALSXP, answer_size)));
//
//        if (_rrownames != R_NilValue) {
//            setAttrib(rcol1, R_LevelsSymbol, _rrownames);
//            setAttrib(rcol1, R_ClassSymbol, mkString("factor"));
//            setAttrib(rcol2, R_LevelsSymbol, _rrownames);
//            setAttrib(rcol2, R_ClassSymbol, mkString("factor"));
//        }
//
//        setAttrib(answer, R_NamesSymbol, (rcolnames = allocVector(STRSXP, num_dims)));
//        setAttrib(answer, R_ClassSymbol, mkString("data.frame"));
//        setAttrib(answer, R_RowNamesSymbol, (rrownames = allocVector(INTSXP, answer_size)));
//
//        for (int i = 0; i < num_dims; i++)
//            SET_STRING_ELT(rcolnames, i, mkChar(COL_NAMES[i]));
//
//        size_t idx1 = 0;
//        for (int ipoint = 0; ipoint < num_points; ++ipoint) {
//            size_t idx2 = k * (size_t)ipoint;
//            for (int i = idx2; i < idx2 + k; ++i) {
//                if (knns_buf[i].dist > threshold)
//                    break;
//
//                INTEGER(rcol1)[idx1] = ipoint + 1;
//                INTEGER(rcol2)[idx1] = knns_buf[i].point;
//                REAL(rdist)[idx1] = knns_buf[i].dist;
//                INTEGER(rrownames)[idx1] = idx1 + 1;
//                ++idx1;
//            }
//        }
//    } else {
//        rprotect(answer = allocVector(REALSXP, (size_t)num_points * (num_points - 1) / 2));
//
//        size_t idx1 = 0;
//        double *d = REAL(answer);
//
//        for (int ipoint1 = 0; ipoint1 < num_points; ++ipoint1) {
//            size_t idx2 = (ipoint1 + 1) * num_points + ipoint1;
//            for (int ipoint2 = ipoint1 + 1; ipoint2 < num_points; ++ipoint2) {
//                d[idx1++] = res[idx2];
//                idx2 += num_points;
//            }
//        }
//
//size_t num_nans = 0;
//for (size_t i = 0; i < res_size; ++i)
//if (std::isnan(res[i]))
//++num_nans;
//printf("Saved %ld distance calculations\n", num_nans - num_points);
//
//        SEXP names = getAttrib(_attrs, R_NamesSymbol);
//        for (int i = 0; i < LENGTH(_attrs); i++)
//            setAttrib(answer, install(translateChar(STRING_ELT(names, i))), VECTOR_ELT(_attrs, i));
//    }
//
//    return answer;
//}

extern "C" {

SEXP tgs_dist(SEXP _x, SEXP _attrs, SEXP _tidy, SEXP _threshold, SEXP _rrownames, SEXP _envir)
{
    SEXP answer = R_NilValue;
    void *shm = (double *)MAP_FAILED;
    size_t shm_sizeof = 0;

	try {
        TGStat tgstat(_envir);

		if ((!isReal(_x) && !isInteger(_x)) || xlength(_x) < 1)
			verror("\"x\" argument must be a matrix of numeric values");

        if (!isLogical(_tidy) || xlength(_tidy) != 1)
            verror("\"tidy\" argument must be a logical value");

        if ((!isReal(_threshold) && !isInteger(_threshold)) || xlength(_threshold) != 1)
            verror("\"threshold\" argument must be a numeric value");

        SEXP rdim = getAttrib(_x, R_DimSymbol);

        if (!isInteger(rdim) || xlength(rdim) != 2)
            verror("\"x\" argument must be a matrix of numeric values");

        bool tidy = asLogical(_tidy);
        double threshold = fabs(asReal(_threshold));

        size_t num_points = nrows(_x);
        size_t num_dims = ncols(_x);

        if (num_points <= 1 || num_dims <= 1)
            verror("\"x\" argument must be a matrix of numeric values");

        size_t num_vals = num_points * num_dims;
        vector<double> vals;
        bool nan_in_vals = false;

        vals.reserve(num_vals);

        for (size_t i = 0; i < num_vals; ++i) {
            if ((isReal(_x) && !R_FINITE(REAL(_x)[i])) || (isInteger(_x) && INTEGER(_x)[i] == NA_INTEGER)) {
                vals.push_back(numeric_limits<double>::quiet_NaN());
                nan_in_vals = true;
            } else
                vals.push_back(isReal(_x) ? REAL(_x)[i] : INTEGER(_x)[i]);
        }

        vdebug("Allocating shared memory for results\n");
        size_t res_size = num_points * num_points;
        shm_sizeof = sizeof(double) * res_size;
        shm = mmap(NULL, shm_sizeof, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);

        if (shm == (double *)MAP_FAILED)
            verror("Failed to allocate shared memory: %s", strerror(errno));

        double *res = (double *)shm;

        int num_cores = max(1, (int)sysconf(_SC_NPROCESSORS_ONLN));
        int num_processes = (int)min(num_points / 2, (size_t)num_cores);
        double num_row4process = num_points / (double)num_processes;

        ProgressReporter progress;
        progress.init(num_points * num_points / 2 - num_points, 1);

        TGStat::prepare4multitasking();

        for (int iprocess = 0; iprocess < num_processes; ++iprocess) {
            if (!TGStat::launch_process()) {     // child process
                size_t srow[2] = { (size_t)(iprocess * num_row4process / 2.), (size_t)(num_points - (iprocess + 1) * num_row4process / 2.) };
                size_t erow[2] = { (size_t)((iprocess + 1) * num_row4process / 2.), (size_t)(num_points - iprocess * num_row4process / 2.) };
                size_t itr_idx = 0;

                for (int ipart = 0; ipart < 2; ipart++) {
                    for (size_t irow1 = 0; irow1 < num_points; ++irow1) {
                        size_t start = max(srow[ipart], irow1 + 1);
                        if (start >= erow[ipart])
                            break;

                        size_t idx = start * num_points + irow1;

                        for (size_t irow2 = start; irow2 < erow[ipart]; ++irow2) {
                            double dist = 0;
                            double dev;
                            size_t count = 0;
                            size_t idx1 = irow1;
                            size_t idx2 = irow2;

                            for (size_t icol = 0; icol < num_dims; ++icol) {
                                if (!std::isnan(vals[idx1]) && !std::isnan(vals[idx2])) {
                                    dev = vals[idx1] - vals[idx2];
                                    dist += dev * dev;
                                    ++count;
                                }
                                idx1 += num_points;
                                idx2 += num_points;
                            }

                            if (count) {
                                if (count == num_dims)
                                    res[idx] = sqrt(dist);
                                else
                                    res[idx] = sqrt(dist * num_dims / (double)count);
                            } else
                                res[idx] = NA_REAL;

                            idx += num_points;
                            TGStat::itr_idx(++itr_idx);
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
        if (tidy) {
            enum { ROW1, ROW2, DIST, num_dims };
            const char *COL_NAMES[num_dims] = { "row1", "row2", "dist" };

            rprotect(answer = allocVector(VECSXP, num_dims));

            size_t answer_size = 0;

            for (size_t ipoint1 = 0; ipoint1 < num_points; ++ipoint1) {
                size_t idx2 = (ipoint1 + 1) * num_points + ipoint1;
                for (size_t ipoint2 = ipoint1 + 1; ipoint2 < num_points; ++ipoint2) {
                    if (res[idx2] <= threshold)
                        ++answer_size;
                    idx2 += num_points;
                }
            }

            SEXP rcol1, rcol2, rdist, rrownames, rcolnames;

            SET_VECTOR_ELT(answer, ROW1, (rcol1 = allocVector(INTSXP, answer_size)));
            SET_VECTOR_ELT(answer, ROW2, (rcol2 = allocVector(INTSXP, answer_size)));
            SET_VECTOR_ELT(answer, DIST, (rdist = allocVector(REALSXP, answer_size)));

            if (_rrownames != R_NilValue) {
                setAttrib(rcol1, R_LevelsSymbol, _rrownames);
                setAttrib(rcol1, R_ClassSymbol, mkString("factor"));
                setAttrib(rcol2, R_LevelsSymbol, _rrownames);
                setAttrib(rcol2, R_ClassSymbol, mkString("factor"));
            }

            setAttrib(answer, R_NamesSymbol, (rcolnames = allocVector(STRSXP, num_dims)));
            setAttrib(answer, R_ClassSymbol, mkString("data.frame"));
            setAttrib(answer, R_RowNamesSymbol, (rrownames = allocVector(INTSXP, answer_size)));

            for (size_t i = 0; i < num_dims; i++)
                SET_STRING_ELT(rcolnames, i, mkChar(COL_NAMES[i]));

            size_t idx1 = 0;
            for (size_t ipoint1 = 0; ipoint1 < num_points; ++ipoint1) {
                size_t idx2 = (ipoint1 + 1) * num_points + ipoint1;
                for (size_t ipoint2 = ipoint1 + 1; ipoint2 < num_points; ++ipoint2) {
                    if (res[idx2] <= threshold) {
                        INTEGER(rcol1)[idx1] = ipoint1 + 1;
                        INTEGER(rcol2)[idx1] = ipoint2 + 1;
                        REAL(rdist)[idx1] = res[idx2];
                        INTEGER(rrownames)[idx1] = idx1 + 1;
                        ++idx1;
                    }
                    idx2 += num_points;
                }
            }
        } else {
            rprotect(answer = allocVector(REALSXP, (size_t)num_points * (num_points - 1) / 2));

            size_t idx1 = 0;
            double *d = REAL(answer);

            for (size_t ipoint1 = 0; ipoint1 < num_points; ++ipoint1) {
                size_t idx2 = (ipoint1 + 1) * num_points + ipoint1;
                for (size_t ipoint2 = ipoint1 + 1; ipoint2 < num_points; ++ipoint2) {
                    d[idx1++] = res[idx2];
                    idx2 += num_points;
                }
            }

            SEXP names = getAttrib(_attrs, R_NamesSymbol);
            for (int i = 0; i < LENGTH(_attrs); i++)
                setAttrib(answer, install(translateChar(STRING_ELT(names, i))), VECTOR_ELT(_attrs, i));
        }
    } catch (TGLException &e) {
        if (!TGStat::is_kid() && shm != (double *)MAP_FAILED) {
            munmap(shm, shm_sizeof);
            shm = (double *)MAP_FAILED;
        }
		rerror("%s", e.msg());
	}

    if (!TGStat::is_kid() && shm != (double *)MAP_FAILED) {
        munmap(shm, shm_sizeof);
        shm = (double *)MAP_FAILED;
    }
	rreturn(answer);
}

SEXP tgs_dist_blas(SEXP _x, SEXP _attrs, SEXP _tidy, SEXP _threshold, SEXP _rrownames, SEXP _envir)
{
    SEXP answer = R_NilValue;

	try {
        struct Mem {
            double *m{NULL};
            double *mask{NULL};
            double *n{NULL};
            double *res{NULL};
            ~Mem() { free(m); free(mask); free(n); free(res); }
        } mem;

        TGStat tgstat(_envir);

		if ((!isReal(_x) && !isInteger(_x)) || xlength(_x) < 1)
			verror("\"x\" argument must be a matrix of numeric values");

        if (!isLogical(_tidy) || xlength(_tidy) != 1)
            verror("\"tidy\" argument must be a logical value");

        if ((!isReal(_threshold) && !isInteger(_threshold)) || xlength(_threshold) != 1)
            verror("\"threshold\" argument must be a numeric value");

        SEXP rdim = getAttrib(_x, R_DimSymbol);

        if (!isInteger(rdim) || xlength(rdim) != 2)
            verror("\"x\" argument must be a matrix of numeric values");

        bool tidy = asLogical(_tidy);
        double threshold = fabs(asReal(_threshold));

        size_t num_points = nrows(_x);
        size_t num_dims = ncols(_x);
        int num_points32 = (int)num_points;
        int num_dims32 = (int)num_dims;

        if (num_points <= 1 || num_dims <= 1)
            verror("\"x\" argument must be a matrix of numeric values");

        size_t num_vals = num_points * num_dims;
        size_t res_size = num_points * num_points;
        bool nan_in_vals = false;

        // some BLAS implementations ask to align double arrays to 64 for improved efficiency
        posix_memalign((void **)&mem.m, 64, sizeof(double) * num_vals);
        posix_memalign((void **)&mem.mask, 64, sizeof(double) * num_vals);
        posix_memalign((void **)&mem.res, 64, sizeof(double) * res_size);

        for (size_t i = 0; i < num_vals; ++i) {
            if ((isReal(_x) && !R_FINITE(REAL(_x)[i])) || (isInteger(_x) && INTEGER(_x)[i] == NA_INTEGER)) {
                mem.m[i] = mem.mask[i] = 0.;
                nan_in_vals = true;
            } else {
                double val = isReal(_x) ? REAL(_x)[i] : INTEGER(_x)[i];
                mem.m[i] = val;
                mem.mask[i] = 1.;
            }
        }

        ProgressReporter progress;
        progress.init(nan_in_vals ? 4 : 2, 1);

        // distance without nans is calculated as:
        //     sqrt(s_x^2 + s_y^2 - 2 * (m %*% t(m)))
        // distance with nans is calculated as:
        //     sqrt((m^2 %*% t(mask) + mask %*% t(m^2) - 2 * (m %*% t(m))) * num_dims / (mask %*% t(mask)))

        // res <- -2 * (m * t(m))
        {
            char uplo = 'L';
            char trans = 'N';
            double alpha = -2;
            double beta = 0;
            F77_NAME(dsyrk)(&uplo, &trans, &num_points32, &num_dims32, &alpha, mem.m, &num_points32, &beta, mem.res, &num_points32);
            check_interrupt();
            progress.report(1);
        }

        if (nan_in_vals) {
            // m <- m^2
            for (size_t i = 0; i < num_vals; ++i)
                mem.m[i] = mem.m[i] * mem.m[i];

            // res <- num_dims * m^2 %*% t(mask) + num_dims * mask %*% t(m^2) + num_dims * res
            {
                char uplo = 'L';
                char trans = 'N';
                double alpha = num_dims;
                double beta = num_dims;
                F77_NAME(dsyr2k)(&uplo, &trans, &num_points32, &num_dims32, &alpha, mem.m, &num_points32, mem.mask, &num_points32, &beta, mem.res, &num_points32);
                check_interrupt();
                progress.report(1);
            }

            // n <- mask %*% t(mask)
            {
                char uplo = 'L';
                char trans = 'N';
                double alpha = 1;
                double beta = 0;
                posix_memalign((void **)&mem.n, 64, sizeof(double) * res_size);
                F77_NAME(dsyrk)(&uplo, &trans, &num_points32, &num_dims32, &alpha, mem.mask, &num_points32, &beta, mem.n, &num_points32);
                check_interrupt();
                progress.report(1);
            }

            // sqrt(res / n)
            for (size_t ipoint1 = 0, idx = 0; ipoint1 < num_points; ++ipoint1) {
                idx += ipoint1 + 1;
                for (size_t ipoint2 = ipoint1 + 1; ipoint2 < num_points; ++ipoint2) {
                    mem.res[idx] = sqrt(mem.res[idx] / mem.n[idx]);
                    ++idx;
                }
            }
            check_interrupt();
            progress.report_last();
        } else {
            vector<double> s_v2(num_points, 0.);

            for (size_t idim = 0, idx = 0; idim < num_dims; ++idim) {
                for (size_t ipoint = 0; ipoint < num_points; ++ipoint) {
                    s_v2[ipoint] += mem.m[idx] * mem.m[idx];
                    ++idx;
                }
            }

            // sqrt(res + s_x^2 + s_y^2)
            for (size_t ipoint1 = 0, idx = 0; ipoint1 < num_points; ++ipoint1) {
                idx += ipoint1 + 1;
                for (size_t ipoint2 = ipoint1 + 1; ipoint2 < num_points; ++ipoint2) {
                    mem.res[idx] = sqrt(max(mem.res[idx] + s_v2[ipoint1] + s_v2[ipoint2], 0.));
                    ++idx;
                }
            }

            check_interrupt();
            progress.report_last();
        }
        
//{
//SEXP rdims;
//rprotect(answer = allocVector(REALSXP, num_points * num_points));
//rprotect(rdims = allocVector(INTSXP, 2));
//INTEGER(rdims)[0] = INTEGER(rdims)[1] = num_points;
//setAttrib(answer, R_DimSymbol, rdims);
//memcpy(REAL(answer), mem.res, xlength(answer) * sizeof(REAL(answer)[0]));
//return answer;
//}

        // assemble the answer
        if (tidy) {
            enum { ROW1, ROW2, DIST, num_dims };
            const char *COL_NAMES[num_dims] = { "row1", "row2", "dist" };

            rprotect(answer = allocVector(VECSXP, num_dims));

            size_t answer_size = 0;

            for (size_t ipoint1 = 0, idx = 0; ipoint1 < num_points; ++ipoint1) {
                idx += ipoint1 + 1;
                for (size_t ipoint2 = ipoint1 + 1; ipoint2 < num_points; ++ipoint2) {
                    if (mem.res[idx] <= threshold)
                        ++answer_size;
                    ++idx;
                }
            }

            SEXP rcol1, rcol2, rdist, rrownames, rcolnames;

            SET_VECTOR_ELT(answer, ROW1, (rcol1 = allocVector(INTSXP, answer_size)));
            SET_VECTOR_ELT(answer, ROW2, (rcol2 = allocVector(INTSXP, answer_size)));
            SET_VECTOR_ELT(answer, DIST, (rdist = allocVector(REALSXP, answer_size)));

            if (_rrownames != R_NilValue) {
                setAttrib(rcol1, R_LevelsSymbol, _rrownames);
                setAttrib(rcol1, R_ClassSymbol, mkString("factor"));
                setAttrib(rcol2, R_LevelsSymbol, _rrownames);
                setAttrib(rcol2, R_ClassSymbol, mkString("factor"));
            }

            setAttrib(answer, R_NamesSymbol, (rcolnames = allocVector(STRSXP, num_dims)));
            setAttrib(answer, R_ClassSymbol, mkString("data.frame"));
            setAttrib(answer, R_RowNamesSymbol, (rrownames = allocVector(INTSXP, answer_size)));

            for (size_t i = 0; i < num_dims; i++)
                SET_STRING_ELT(rcolnames, i, mkChar(COL_NAMES[i]));

            size_t idx1 = 0;
            size_t idx2 = 0;

            for (size_t ipoint1 = 0; ipoint1 < num_points; ++ipoint1) {
                idx2 += ipoint1 + 1;
                for (size_t ipoint2 = ipoint1 + 1; ipoint2 < num_points; ++ipoint2) {
                    if (mem.res[idx2] <= threshold) {
                        INTEGER(rcol1)[idx1] = ipoint1 + 1;
                        INTEGER(rcol2)[idx1] = ipoint2 + 1;
                        REAL(rdist)[idx1] = mem.res[idx2];
                        INTEGER(rrownames)[idx1] = idx1 + 1;
                        ++idx1;
                    }
                    ++idx2;
                }
            }
        } else {
            rprotect(answer = allocVector(REALSXP, (size_t)num_points * (num_points - 1) / 2));

            size_t idx1 = 0;
            size_t idx2 = 0;
            double *d = REAL(answer);

            for (size_t ipoint1 = 0; ipoint1 < num_points; ++ipoint1) {
                idx2 += ipoint1 + 1;
                for (size_t ipoint2 = ipoint1 + 1; ipoint2 < num_points; ++ipoint2)
                    d[idx1++] = mem.res[idx2++];
            }

            SEXP names = getAttrib(_attrs, R_NamesSymbol);
            for (int i = 0; i < LENGTH(_attrs); i++)
                setAttrib(answer, install(translateChar(STRING_ELT(names, i))), VECTOR_ELT(_attrs, i));
        }
    } catch (TGLException &e) {
		rerror("%s", e.msg());
	}

	rreturn(answer);
}

}

