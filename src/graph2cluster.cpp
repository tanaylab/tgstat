#include <algorithm>
#include <limits>
#include <unordered_map>
#include <unordered_set>

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

SEXP tgs_graph2cluster(SEXP _graph, SEXP _min_cluster_size, SEXP _cooling, SEXP _burn_in, SEXP _envir)
{
	try {
        TGStat tgstat(_envir);

        int *pcol1;
        int *pcol2;
        double *pweight;
        size_t num_edges;

        {
            enum { COL1, COL2, WEIGHT, NUM_COLS };
            const char *COL_NAMES[NUM_COLS] = { "col1", "col2", "weight" };

            SEXP rnames = getAttrib(_graph, R_NamesSymbol);

            if (!isVector(_graph) || xlength(_graph) != NUM_COLS || xlength(rnames) != NUM_COLS ||
                strcmp(CHAR(STRING_ELT(rnames, COL1)), COL_NAMES[COL1]) || !isInteger(VECTOR_ELT(_graph, COL1)) && !isFactor(VECTOR_ELT(_graph, COL1)) ||
                strcmp(CHAR(STRING_ELT(rnames, COL2)), COL_NAMES[COL2]) || !isInteger(VECTOR_ELT(_graph, COL2)) && !isFactor(VECTOR_ELT(_graph, COL2)) ||
                xlength(VECTOR_ELT(_graph, COL2)) != xlength(VECTOR_ELT(_graph, COL1)) ||
                strcmp(CHAR(STRING_ELT(rnames, WEIGHT)), COL_NAMES[WEIGHT]) || !isReal(VECTOR_ELT(_graph, WEIGHT)) || xlength(VECTOR_ELT(_graph, WEIGHT)) != xlength(VECTOR_ELT(_graph, COL1)))
                verror("\"graph\" argument must be in the format that is returned by tgs_cor_graph function");

            pcol1 = INTEGER(VECTOR_ELT(_graph, COL1));
            pcol2 = INTEGER(VECTOR_ELT(_graph, COL2));
            pweight = REAL(VECTOR_ELT(_graph, WEIGHT));
            num_edges = xlength(VECTOR_ELT(_graph, COL1));
        }

        if (!isInteger(_min_cluster_size) && !isReal(_min_cluster_size) || xlength(_min_cluster_size) != 1 || asInteger(_min_cluster_size) < 1)
            verror("\"min_cluster_size\" argument must be a positive integer");

        if (!isInteger(_cooling) && !isReal(_cooling) || xlength(_cooling) != 1 || asReal(_cooling) < 1)
            verror("\"cooling\" argument must be a number greater or equal than 1");

        if (!isInteger(_burn_in) && !isReal(_burn_in) || xlength(_burn_in) != 1 || asInteger(_burn_in) < 0)
            verror("\"burn_in\" argument must be a positive integer");

        vdebug("Finding seeds...\n");

        unsigned min_cluster_size = asInteger(_min_cluster_size);
        float cooling_rate = asReal(_cooling);
        unsigned burn_in = asInteger(_burn_in);

        unordered_map<pair<unsigned, unsigned>, float> ij2weight;
        unsigned num_points = 0;

        for (size_t i = 0; i < num_edges; ++i) {
            ij2weight[{pcol1[i] - 1, pcol2[i] - 1}] = pweight[i];
            num_points = max(num_points, (unsigned)pcol1[i]);
            num_points = max(num_points, (unsigned)pcol2[i]);
        }

        vector<vector<unsigned>> incoming(num_points);
        vector<vector<unsigned>> outgoing(num_points);
        vector<unsigned> node2seed(num_points, -1);
        vector<unsigned> node2cluster(num_points, -1);
        vector<unsigned> cluster_sizes;
        unsigned num_clusters = 0;
        unordered_set<unsigned> cands;                  // seed candidates 
        vector<unsigned> node2weight(num_points, 0);    // node2weight[i] = | incoming(i) - COVERED |
        size_t weights_sum = 0;
        unsigned max_weight = 0;
        unsigned min_weight = min_cluster_size;

        for (size_t i = 0; i < num_edges; ++i) {
            incoming[pcol2[i] - 1].push_back(pcol1[i] - 1);
            outgoing[pcol1[i] - 1].push_back(pcol2[i] - 1);
        }

        for (size_t i = 0; i < num_points; ++i) {
            node2weight[i] = incoming[i].size();
            max_weight = max(max_weight, node2weight[i]);
            if (node2weight[i] >= min_weight)
                cands.insert(i);
            weights_sum += node2weight[i];
        }

        // stage 1: find the seeds around which the clusters would form

        while (max_weight >= min_weight) {
            // select the next seed randomly with weight = node2weight
            float v = drand48() * weights_sum;
            float sum = 0;
            unsigned seed;

            for (auto i : cands) {
                sum += node2weight[i];
                if (sum > v) {
                    seed = i;
                    break;
                }
            }

            // update all our data structures
            cluster_sizes.push_back(0);

            for (auto i : incoming[seed]) {
                if (node2seed[i] == -1) {
                    for (auto j : outgoing[i]) {
                        --node2weight[j];
                        if (node2weight[j] < min_weight)
                            cands.erase(j);
                    }
                    cands.erase(i);
                    node2seed[i] = seed;
                    node2cluster[i] = num_clusters;
                    ++cluster_sizes.back();
                }
            }
            cands.erase(seed);
            node2seed[seed] = seed;
            node2cluster[seed] = num_clusters;
            ++cluster_sizes.back();

            max_weight = 0;
            weights_sum = 0;
            for (auto i : cands) {
                weights_sum += node2weight[i];
                max_weight = max(max_weight, node2weight[i]);
            }

            check_interrupt();
            ++num_clusters;
        }

//printf("Num clusters: %d\n", num_clusters);
//for (auto i = node2seed.begin(); i < node2seed.end(); ++i)
//printf("%d -> %d\n", i - node2seed.begin(), *i);

        // look for the nodes that have not been covered yet and assign them to the seeds
        // that achieve the maximal weight along the d2 path; nodes that cannot reach the seeds within d2 will stay unassigned
        for (auto i = node2seed.begin(); i < node2seed.end(); ++i) {
            if (*i == -1) {
                float best_weight = -1;
                unsigned best_seed = -1;
                unsigned best_cluster = -1;
                unsigned origin = i - node2seed.begin();

                for (auto j : outgoing[origin]) {
                    if (node2seed[j] != -1 &&                  // j is covered
                        ij2weight.count({j, node2seed[j]}) &&  // j is connected directly to the seed (and not via someone else)
                        best_weight < ij2weight[{origin, j}] + ij2weight[{j, node2seed[j]}])    // i->j->seed achieves the highest weight
                    {
                        best_weight = ij2weight[{origin, j}] + ij2weight[{j, node2seed[j]}];
                        best_seed = node2seed[j];
                        best_cluster = node2cluster[j];
                    }
                }

                if (best_seed != -1) {
                    *i = best_seed;
                    node2cluster[i - node2seed.begin()] = node2cluster[best_seed];
                    ++cluster_sizes[node2cluster[best_seed]];
                }
            }
        }

        check_interrupt();
        vdebug("Formed %d clusters\n", num_clusters);
//printf("Cluster sizes: %ld, num_clusters: %ld\n", cluster_sizes.size(), num_clusters);
//for (auto i = node2seed.begin(); i < node2seed.end(); ++i)
//printf("%d -> %d, cluster: %d, cluster size: %d\n", i - node2seed.begin(), *i, node2cluster[i - node2seed.begin()], node2cluster[i - node2seed.begin()] == -1 ? -1 : cluster_sizes[node2cluster[i - node2seed.begin()]]);

        // stage 2: consolidate the clusters

        struct Score {
            float votes_in;          // sum of edges weights from the node to the cluster
            float votes_out;         // sum of edges weights from the cluster to the node
//            Score *pnext;             // pointer to the next cluster where both votes_in and voutes_out are non zero (= total votes is not zero);
                                      // we use this pointer to speed up our search for the highest score (the majority of clusters are expected to get zero score)

//          Score() : votes_in(0), votes_out(0), pnext(NULL) {}
//          Score(const Score &o) : votes_in(o.votes_in), votes_out(o.votes_out), pnext(o.pnext) {}
Score() : votes_in(0), votes_out(0) {}
Score(const Score &o) : votes_in(o.votes_in), votes_out(o.votes_out) {}
        };

        vector<vector<Score>> node2cluster_score(num_points);   // given node i and cluster j the score is at [i][j]
        vector<Score *> node2nzero_score(num_points, NULL);
        vector<unsigned> reassignments(num_points, 0);          // number of reassignments of node i from one cluster to another
        vector<float> cooling_rates(num_points, 1);

//printf("0. COOLING RATES SIZE: %ld\n", cooling_rates.size());
        vdebug("Starting consolidation of nodes around the clusters...\n");
        for (unsigned i = 0; i < num_points; ++i) {
            vector<Score> &cluster2score = node2cluster_score[i];
            cluster2score.resize(num_clusters);

            for (auto j : outgoing[i]) {
                if (node2cluster[j] != -1)
                    cluster2score[node2cluster[j]].votes_in += ij2weight[{i, j}];
            }

            for (auto j : incoming[i]) {
                if (node2cluster[j] != -1)
                    cluster2score[node2cluster[j]].votes_out += ij2weight[{j, i}];
            }
        }
        check_interrupt();

        // set the pnext: the list of the non-zero scored clusters starts at node2nzero_score and then follows the list by pnext
        vdebug("Setting pnext...\n");
//      for (unsigned i = 0; i < num_points; ++i) {
//          Score **pprev_score = &node2nzero_score[i];
//          vector<Score> &cluster2score = node2cluster_score[i];
//
//          for (auto &score : cluster2score) {
//              if (score.votes_in && score.votes_out) {
//                  *pprev_score = &score;
//                  pprev_score = &score.pnext;
//              }
//          }
//      }
        vdebug("Setting pnext - DONE\n");
        check_interrupt();

//for (unsigned i = 0; i < num_points; ++i) {
//Score **pprev_score = &node2nzero_score[i];
//
//printf("point %d\n", i);
//for (Score *pscore = node2nzero_score[i]; pscore; pscore = pscore->pnext) {
//unsigned cluster = pscore - &*node2cluster_score[i].begin();
//printf("\tcluster %d, votes_in: %g, votes_out: %g\n", cluster, pscore->votes_in, pscore->votes_out);
//}
//}
        vector<unsigned> candidates(num_points);

        for (unsigned i = 0; i < num_points; ++i)
            candidates[i] = i;

size_t pnext_maintenance_cost = 0;;

        while (1) {
            unsigned num_reassignments = 0;

            random_shuffle(candidates.begin(), candidates.end());

            for (auto i : candidates) {
//printf("Candidate %d\n", i);
                // find the maximal score
                int new_cluster = -1;
                int old_cluster = node2cluster[i];
                float max_score = 0;
//double old_score;

vector<Score> &cluster2score = node2cluster_score[i];
//                for (Score *pscore = node2nzero_score[i]; pscore; pscore = pscore->pnext) {
//                    int cluster = pscore - &*node2cluster_score[i].begin();
for (int cluster = 0; cluster < num_clusters; ++cluster) {
Score *pscore = &cluster2score[cluster];
if (pscore->votes_in && pscore->votes_out){
                    float score = cluster == old_cluster ?
                        cooling_rates[i] * pscore->votes_in * pscore->votes_out / (cluster_sizes[cluster] * (float)cluster_sizes[cluster]) :
                        pscore->votes_in * pscore->votes_out / ((cluster_sizes[cluster] + 1) * (float)(cluster_sizes[cluster] + 1));

                    if (score > max_score) {
                        new_cluster = cluster;
                        max_score = score;
                    }
}
//if (cluster == old_cluster)
//old_score = score;
                }

                if (new_cluster != -1 && new_cluster != old_cluster) {
//printf("Cand %d, old cluster: %d, new cluster: %d, old score: %g, new score: %g, cooling rate: %g\n", i, old_cluster, new_cluster, old_score, max_score, cooling_rates[i]);

                    // reassign the cluster
                    for (auto j : incoming[i]) {
                        float weight = ij2weight[{j, i}];
                        vector<Score> &cluster2score = node2cluster_score[j];

                        if (old_cluster != -1) {     // some nodes might not be assigned yet to any cluster
                            cluster2score[old_cluster].votes_in -= weight;

//                            if (cluster2score[old_cluster].votes_in <= 0 && cluster2score[old_cluster].votes_out > 0) {    // from now on this cluster should be bypassed
////printf("\tnode %d, cluster %d, votes in turned zero => removed\n", j, old_cluster);
//                                int cluster;
//                                for (cluster = old_cluster - 1; cluster >= 0; --cluster) {
//pnext_maintenance_cost++;
//                                    if (cluster2score[cluster].votes_in > 0 && cluster2score[cluster].votes_out > 0) {
//                                        cluster2score[cluster].pnext = cluster2score[old_cluster].pnext;
//                                        break;
//                                    }
//                                }
//                                if (cluster < 0)
//                                    node2nzero_score[j] = cluster2score[old_cluster].pnext;
//                            }
                        }

//                        if (cluster2score[new_cluster].votes_in <= 0 && cluster2score[new_cluster].votes_out > 0) {      // this cluster is currently bypassed => add it to the list
////printf("\tnode %d, cluster %d, votes in resurrected => added\n", j, new_cluster);
//                            int cluster;
//                            for (cluster = new_cluster - 1; cluster >= 0; --cluster) {
//pnext_maintenance_cost++;
//                                if (cluster2score[cluster].votes_in > 0 && cluster2score[cluster].votes_out > 0) {
//                                    cluster2score[new_cluster].pnext = cluster2score[cluster].pnext;
//                                    cluster2score[cluster].pnext = &cluster2score[new_cluster];
//                                    break;
//                                }
//                            }
//                            if (cluster < 0) {
//                                cluster2score[new_cluster].pnext = node2nzero_score[j];
//                                node2nzero_score[j] = &cluster2score[new_cluster];
//                            }
//                        }
                        cluster2score[new_cluster].votes_in += weight;
                    }

                    for (auto j : outgoing[i]) {
                        float weight = ij2weight[{i, j}];
                        vector<Score> &cluster2score = node2cluster_score[j];

                        if (old_cluster != -1) {
                            cluster2score[old_cluster].votes_out -= weight;

//                            if (cluster2score[old_cluster].votes_out <= 0 && cluster2score[old_cluster].votes_in > 0) {    // from now on this cluster should be bypassed
////printf("\tnode %d, cluster %d, votes out turned zero => removed\n", j, old_cluster);
//                                int cluster;
//                                for (cluster = old_cluster - 1; cluster >= 0; --cluster) {
//pnext_maintenance_cost++;
//                                    if (cluster2score[cluster].votes_in > 0 && cluster2score[cluster].votes_out > 0) {
//                                        cluster2score[cluster].pnext = cluster2score[old_cluster].pnext;
//                                        break;
//                                    }
//                                }
//                                if (cluster < 0)
//                                    node2nzero_score[j] = cluster2score[old_cluster].pnext;
//                            }
                        }

//                        if (cluster2score[new_cluster].votes_out <= 0 && cluster2score[new_cluster].votes_in > 0) {      // this cluster is currently bypassed => add it to the list
////printf("\tnode %d, cluster %d, votes out resurrected => added\n", j, new_cluster);
//                            int cluster;
//                            for (cluster = new_cluster - 1; cluster >= 0; --cluster) {
//pnext_maintenance_cost++;
//                                if (cluster2score[cluster].votes_in > 0 && cluster2score[cluster].votes_out > 0) {
//                                    cluster2score[new_cluster].pnext = cluster2score[cluster].pnext;
//                                    cluster2score[cluster].pnext = &cluster2score[new_cluster];
//                                    break;
//                                }
//                            }
//                            if (cluster < 0) {
//                                cluster2score[new_cluster].pnext = node2nzero_score[j];
//                                node2nzero_score[j] = &cluster2score[new_cluster];
//                            }
//                        }
                        cluster2score[new_cluster].votes_out += weight;
                    }

                    node2cluster[i] = new_cluster;
                    ++reassignments[i];
                    ++num_reassignments;

                    if (reassignments[i] > burn_in)
                        cooling_rates[i] *= cooling_rate;

                    if (old_cluster != -1)
                        --cluster_sizes[old_cluster];
                    ++cluster_sizes[new_cluster];

//printf("\n\nAFTER reassignment:\n");
//printf("global cooling rate: %g, burn in: %d\n", cooling_rate, burn_in);
//for (unsigned i = 0; i < num_clusters; ++i) {
//printf("Cluster[%d], size: %d\n", i, cluster_sizes[i]);
//}
//
//for (unsigned i = 0; i < num_points; ++i) {
//Score **pprev_score = &node2nzero_score[i];
//printf("point %d, cluster %d, cooling rate %g, reassignments %d\n", i, node2cluster[i], cooling_rates[i], reassignments[i]);
//for (Score *pscore = node2nzero_score[i]; pscore; pscore = pscore->pnext) {
//unsigned cluster = pscore - &*node2cluster_score[i].begin();
//printf("\tcluster %d, votes_in: %g, votes_out: %g\n", cluster, pscore->votes_in, pscore->votes_out);
//}
//}

//if (num_reassignments == 2)
//rreturn(R_NilValue);
                }
            }

            check_interrupt();

            if (!num_reassignments)
                break;
        }
        vdebug("Num iterations for non-zero score list maintenance: %ld\n", pnext_maintenance_cost);
        vdebug("Consolidation - DONE\n");

        if (g_tgstat->debug()) {
            double total_weight = 0;
            double in_cluster_weight = 0;
            double unassigned_weight = 0;
            vector<double> in_cluster_weights(num_clusters, 0);
            vector<double> incoming_weights(num_clusters, 0);
            vector<double> outgoing_weights(num_clusters, 0);

            for (auto &edge_weight : ij2weight) {
                unsigned i = edge_weight.first.first;
                unsigned j = edge_weight.first.second;
                double weight = edge_weight.second;

                total_weight += weight;
                if (node2cluster[i] == -1 || node2cluster[j] == -1)
                    unassigned_weight += weight;
                else {
                    if (node2cluster[i] == node2cluster[j]) {
                        in_cluster_weight += weight;
                        in_cluster_weights[node2cluster[i]] += weight;
                    } else {
                        outgoing_weights[node2cluster[i]] += weight;
                        incoming_weights[node2cluster[j]] += weight;
                    }
                }
            }

            printf("TOTAL WEIGHTS: %g\n", total_weight);
            printf("  IN CLUSTERS: %g\n", in_cluster_weight);
            printf("  UNASSIGNED:  %g\n", unassigned_weight);
            printf("  IN/TOTAL:    %g\n", in_cluster_weight / total_weight);
            printf("\n");

            for (unsigned cluster = 0; cluster < num_clusters; ++cluster) {
                printf("CLUSTER %d\n", cluster + 1);
                printf("\tIN CLUSTER: %g\n", in_cluster_weights[cluster]);
                printf("\tINCOMING:   %g\n", incoming_weights[cluster]);
                printf("\tOUTGOING:   %g\n", outgoing_weights[cluster]);
                printf("\tIN/TOTAL:   %g\n", in_cluster_weights[cluster] / (in_cluster_weights[cluster] + incoming_weights[cluster] + outgoing_weights[cluster]));
            }
        }

//for (unsigned i = 0; i < num_points; ++i) {
//printf("== point %d, cluster %d, cooling rate %g, reassignments %d\n", i, node2cluster[i], cooling_rates[i], reassignments[i]);
//}
//for (unsigned i = 0; i < num_clusters; ++i) {
//printf("Cluster[%d], size: %d\n", i, cluster_sizes[i]);
//}

        vdebug("Packing the return value\n");

        enum { NODE, CLUSTER, NUM_COLS };
        const char *COL_NAMES[NUM_COLS] = { "node", "cluster" };

        SEXP ranswer, rnode, rcluster, rrownames, rcolnames;

        rprotect(ranswer = allocVector(VECSXP, NUM_COLS));
        SET_VECTOR_ELT(ranswer, NODE, (rnode = allocVector(INTSXP, num_points)));
        SET_VECTOR_ELT(ranswer, CLUSTER, (rcluster = allocVector(INTSXP, num_points)));

        setAttrib(ranswer, R_NamesSymbol, (rcolnames = allocVector(STRSXP, NUM_COLS)));
        setAttrib(ranswer, R_ClassSymbol, mkString("data.frame"));
        setAttrib(ranswer, R_RowNamesSymbol, (rrownames = allocVector(INTSXP, num_points)));

        for (int i = 0; i < NUM_COLS; i++)
            SET_STRING_ELT(rcolnames, i, mkChar(COL_NAMES[i]));

        for (unsigned i = 0; i < num_points; ++i) {
            INTEGER(rnode)[i] = i + 1;
            INTEGER(rcluster)[i] = node2cluster[i] + 1;
            INTEGER(rrownames)[i] = i + 1;
        }

        vdebug("Packing the return value - DONE\n");

        rreturn(ranswer);
    } catch (TGLException &e) {
		rerror("%s", e.msg());
	}

    rreturn(R_NilValue);
}

}


