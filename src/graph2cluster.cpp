#include <algorithm>
#include <errno.h>
#include <limits>
#include <strings.h>
#include <unordered_map>
#include <unordered_set>
#include <unistd.h>
#include <sys/mman.h>

#include <R.h>
#include <Rinternals.h>

#include "HashFunc.h"
#include "RandomShuffle.h"

#ifdef length
#undef length
#endif
#ifdef error
#undef error
#endif

#include "ProgressReporter.h"
#include "tgstat.h"

unsigned graph2cluster(const int *pnode1, const int *pnode2, const double *pweight, size_t num_edges,
                       unsigned min_cluster_size, float cooling_rate, unsigned burn_in, unsigned *node2cluster, size_t num_nodes)
{
    vdebug("Finding seeds...\n");

    unordered_map<pair<unsigned, unsigned>, float> ij2weight;

    ij2weight.reserve(num_edges);
    for (size_t i = 0; i < num_edges; ++i)
        ij2weight[{pnode1[i] - 1, pnode2[i] - 1}] = pweight[i];

    vector<vector<unsigned>> incoming(num_nodes);
    vector<vector<unsigned>> outgoing(num_nodes);
    vector<unsigned> node2seed(num_nodes, -1);
    vector<unsigned> cluster_sizes;
    unsigned num_clusters = 0;
    unordered_set<unsigned> cands;                 // seed candidates 
    vector<unsigned> node2weight(num_nodes, 0);    // node2weight[i] = | incoming(i) - COVERED |
    size_t weights_sum = 0;
    unsigned max_weight = 0;
    unsigned min_weight = min_cluster_size;

    for (size_t i = 0; i < num_edges; ++i) {
        incoming[pnode2[i] - 1].push_back(pnode1[i] - 1);
        outgoing[pnode1[i] - 1].push_back(pnode2[i] - 1);
    }

    cands.reserve(num_nodes);
    for (size_t i = 0; i < num_nodes; ++i) {
        node2weight[i] = incoming[i].size();
        max_weight = max(max_weight, node2weight[i]);
        if (node2weight[i] >= min_weight) {
            cands.insert(i);
            weights_sum += node2weight[i];
        }
        node2cluster[i] = -1;
    }

    // stage 1: find the seeds around which the clusters would form

    while (max_weight >= min_weight) {
        // select the next seed randomly with weight = node2weight
        float v = unif_rand() * weights_sum;
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

//REprintf("Num clusters: %d\n", num_clusters);
//for (auto i = node2seed.begin(); i < node2seed.end(); ++i)
//REprintf("\tnode2seed %d -> %d\n", i - node2seed.begin(), *i);
//REprintf("Num clusters: %d\n", num_clusters);
//for (size_t i = 0; i < num_nodes; ++i)
//REprintf("\t%d -> %d\n", i, node2cluster[i]);

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
//REprintf("Cluster sizes: %ld, num_clusters: %ld\n", cluster_sizes.size(), num_clusters);
//for (auto i = node2seed.begin(); i < node2seed.end(); ++i)
//REprintf("%d -> %d, cluster: %d, cluster size: %d\n", i - node2seed.begin(), *i, node2cluster[i - node2seed.begin()], node2cluster[i - node2seed.begin()] == -1 ? -1 : cluster_sizes[node2cluster[i - node2seed.begin()]]);

    // stage 2: consolidate the clusters

    struct Score {
        float votes_in;          // sum of edges weights from the node to the cluster
        float votes_out;         // sum of edges weights from the cluster to the node

        Score() : votes_in(0), votes_out(0) {}
        Score(const Score &o) : votes_in(o.votes_in), votes_out(o.votes_out) {}
    };

    vector<vector<Score>> node2cluster_score(num_nodes);   // given node i and cluster j the score is at [i][j]
    vector<unsigned> reassignments(num_nodes, 0);          // number of reassignments of node i from one cluster to another
    vector<float> cooling_rates(num_nodes, 1);

    vdebug("Starting consolidation of nodes around the clusters...\n");
    for (size_t i = 0; i < num_nodes; ++i) {
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

    vector<unsigned> candidates(num_nodes);

    for (size_t i = 0; i < num_nodes; ++i)
        candidates[i] = i;

    while (1) {
        unsigned num_reassignments = 0;

        tgs_random_shuffle(candidates.begin(), candidates.end(), unif_rand);

        for (auto i : candidates) {
            // find the maximal score
            int new_cluster = -1;
            int old_cluster = node2cluster[i];
            float max_score = 0;

            vector<Score> &cluster2score = node2cluster_score[i];
            for (int cluster = 0; cluster < num_clusters; ++cluster) {
                const Score &cscore = cluster2score[cluster];
                if (cscore.votes_in && cscore.votes_out) {
                    float score = cluster == old_cluster ?
                        cooling_rates[i] * cscore.votes_in * cscore.votes_out / (cluster_sizes[cluster] * (float)cluster_sizes[cluster]) :
                        cscore.votes_in * cscore.votes_out / ((cluster_sizes[cluster] + 1) * (float)(cluster_sizes[cluster] + 1));

                    if (score > max_score) {
                        new_cluster = cluster;
                        max_score = score;
                    }
                }
            }

            if (new_cluster != -1 && new_cluster != old_cluster) {
                // reassign the cluster
                for (auto j : incoming[i]) {
                    float weight = ij2weight[{j, i}];
                    vector<Score> &cluster2score = node2cluster_score[j];

                    if (old_cluster != -1)      // some nodes might not be assigned yet to any cluster
                        cluster2score[old_cluster].votes_in -= weight;

                    cluster2score[new_cluster].votes_in += weight;
                }

                for (auto j : outgoing[i]) {
                    float weight = ij2weight[{i, j}];
                    vector<Score> &cluster2score = node2cluster_score[j];

                    if (old_cluster != -1)
                        cluster2score[old_cluster].votes_out -= weight;

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
            }
        }

        check_interrupt();

        if (!num_reassignments)
            break;
    }
    vdebug("Consolidation - DONE\n");

//for (unsigned i = 0; i < num_clusters; ++i) {
//REprintf("Cluster[%d], size: %d\n", i, cluster_sizes[i]);
//}
//
    if (0 && g_tgstat->debug()) {
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

        REprintf("TOTAL WEIGHTS: %g\n", total_weight);
        REprintf("  IN CLUSTERS: %g\n", in_cluster_weight);
        REprintf("  UNASSIGNED:  %g\n", unassigned_weight);
        REprintf("  IN/TOTAL:    %g\n", in_cluster_weight / total_weight);
        REprintf("\n");

        for (unsigned cluster = 0; cluster < num_clusters; ++cluster) {
            REprintf("CLUSTER %d\n", cluster + 1);
            REprintf("\tIN CLUSTER: %g\n", in_cluster_weights[cluster]);
            REprintf("\tINCOMING:   %g\n", incoming_weights[cluster]);
            REprintf("\tOUTGOING:   %g\n", outgoing_weights[cluster]);
            REprintf("\tIN/TOTAL:   %g\n", in_cluster_weights[cluster] / (in_cluster_weights[cluster] + incoming_weights[cluster] + outgoing_weights[cluster]));
        }
    }

//for (unsigned i = 0; i < num_nodes; ++i)
//REprintf("Node %d => %d\n", i + 1, node2cluster[i] + 1);

    return num_clusters;
}

void reassign_weights(const int *pnode1, const int *pnode2, unsigned knn, size_t num_edges,
                      const vector<bool> &node_selected, vector<int> &nodes1, vector<int> &nodes2, vector<double> &weights)
{
    unsigned last_idx = 0;
    unsigned last_node = pnode1[0] - 1;
    unsigned num_node_edges = node_selected[pnode2[0] - 1] ? 1 : 0;

    for (size_t i = 1; ; ++i) {
        if (i == num_edges || last_node != pnode1[i] - 1) {
            if (node_selected[last_node] && num_node_edges) {
                double rank = 0;
                num_node_edges = min(num_node_edges, knn);

                for (size_t j = last_idx; j < i; ++j) {
                    unsigned pointed_node = pnode2[j] - 1;

                    if (node_selected[pointed_node]) {
                        nodes1.push_back(last_node + 1);
                        nodes2.push_back(pointed_node + 1);
                        weights.push_back(1. - rank / num_node_edges);
                        ++rank;
                        if (rank >= knn)
                            break;
                    }
                }
            }

            if (i == num_edges)
                break;

            last_idx = i;
            last_node = pnode1[i] - 1;
            num_node_edges = 0;
        }

        if (node_selected[pnode1[i] - 1] && node_selected[pnode2[i] - 1])
            ++num_node_edges;
    }
}

void launch_kid_hash(const int *pnode1, const int *pnode2, const double *pweight, size_t num_nodes, size_t num_edges, unsigned *res, int slot,
                     unsigned knn, double p_resamp, unsigned min_cluster_size, float cooling_rate, unsigned burn_in, uint64_t seed)
{
    vdebug("Launching a working process at slot %d\n", slot);
    if (!TGStat::launch_process()) {     // child process
        size_t num_kid_nodes = max(1., p_resamp * num_nodes);
        vector<bool> node_selected(num_nodes, false);
        vector<unsigned> nodes(num_nodes);
        for (size_t i = 0; i < num_nodes; ++i)
            nodes[i] = i;

        g_tgstat->rnd_seed(seed);
        vdebug("Random seed: %ld\n", seed);
        tgs_random_shuffle(nodes.begin(), nodes.end(), unif_rand);
        sort(nodes.begin(), nodes.begin() + num_kid_nodes);   // this helps to reduce memory cache misses later
        for (unsigned i = 0; i < num_kid_nodes; ++i)
            node_selected[nodes[i]] = true;

        vector<int> nodes1;
        vector<int> nodes2;
        vector<double> weights;

        reassign_weights(pnode1, pnode2, knn, num_edges, node_selected, nodes1, nodes2, weights);
        vdebug("num child edges = %ld, num all edges: %ld\n", nodes1.size(), num_edges);

        unsigned kid_res_sizeof = sizeof(unsigned) + sizeof(unsigned) + sizeof(unsigned) * num_nodes;
        unsigned *pready = (unsigned *)((char *)res + slot * kid_res_sizeof);
        unsigned *num_clusters = pready + 1;
        unsigned *clusters = num_clusters + 1;

        g_tgstat->rnd_seed(seed);   // this will make the results reproducible by tgs_graph_cover
        *num_clusters = graph2cluster(nodes1.data(), nodes2.data(), weights.data(), nodes1.size(), min_cluster_size, cooling_rate, burn_in, clusters, num_nodes);

        for (size_t i = num_kid_nodes; i < num_nodes; ++i)
            clusters[nodes[i]] = -2;                          // mark nodes that were not selected by -2

        *pready = 1;
        rexit();
    }
}

void launch_kid_full(const int *pnode1, const int *pnode2, const double *pweight, size_t num_nodes, size_t num_edges, unsigned *pready, unsigned short *pco_clust,
                     unsigned short *psamples, unsigned knn, double p_resamp, unsigned min_cluster_size, float cooling_rate, unsigned burn_in, uint64_t seed)
{
    if (!TGStat::launch_process()) {     // child process
        unsigned num_kid_nodes = max(1., p_resamp * num_nodes);
        vector<bool> node_selected(num_nodes, false);
        vector<unsigned> nodes(num_nodes);
        for (size_t i = 0; i < num_nodes; ++i)
            nodes[i] = i;

        g_tgstat->rnd_seed(seed);
        vdebug("Random seed: %ld\n", seed);
        tgs_random_shuffle(nodes.begin(), nodes.end(), unif_rand);
        sort(nodes.begin(), nodes.begin() + num_kid_nodes);   // this helps to reduce memory cache misses later
        for (unsigned i = 0; i < num_kid_nodes; ++i) {
            node_selected[nodes[i]] = true;
            __sync_add_and_fetch(psamples + nodes[i], 1);
        }

        vector<unsigned> node2cluster(num_nodes, -1);
        vector<int> nodes1;
        vector<int> nodes2;
        vector<double> weights;

        reassign_weights(pnode1, pnode2, knn, num_edges, node_selected, nodes1, nodes2, weights);
        vdebug("num child edges = %ld, num all edges: %ld\n", nodes1.size(), num_edges);

        g_tgstat->rnd_seed(seed);   // this will make the results reproducible by tgs_graph_cover
        unsigned num_clusters = graph2cluster(nodes1.data(), nodes2.data(), weights.data(), nodes1.size(), min_cluster_size, cooling_rate, burn_in, node2cluster.data(), num_nodes);

        vector<vector<unsigned>> cluster2nodes(num_clusters);

        for (size_t i = 0; i < num_kid_nodes; ++i) {
            if (node2cluster[nodes[i]] != -1)
                cluster2nodes[node2cluster[nodes[i]]].push_back(nodes[i]);
        }

        for (auto &cluster_nodes : cluster2nodes) {
            sort(cluster_nodes.begin(), cluster_nodes.end());
            for (auto inode1 = cluster_nodes.begin(); inode1 < cluster_nodes.end(); ++inode1) {
                // this ugly offset translates node1, node2 to the index in the buffer that stores the lower triangle (incl. diagonal) of the co_cluster matrix
                unsigned short *p = pco_clust + ((size_t)*inode1 * (2 * num_nodes - (size_t)*inode1 - 1)) / 2;
                for (auto inode2 = inode1; inode2 < cluster_nodes.end(); ++inode2)
                    __sync_add_and_fetch(p + *inode2, 1);
            }
        }

        *pready = 1;
        rexit();
    }
}

void launch_kid_edges(const int *pnode1, const int *pnode2, const double *pweight, size_t num_nodes, size_t num_edges, unsigned *pready, unsigned short *pco_clust,
                      unsigned short *psamples, unsigned knn, double p_resamp, unsigned min_cluster_size, float cooling_rate, unsigned burn_in, uint64_t seed)
{
    if (!TGStat::launch_process()) {     // child process
        unsigned num_kid_nodes = max(1., p_resamp * num_nodes);
        vector<bool> node_selected(num_nodes, false);
        vector<unsigned> nodes(num_nodes);
        for (size_t i = 0; i < num_nodes; ++i)
            nodes[i] = i;

        g_tgstat->rnd_seed(seed);
        vdebug("Random seed: %ld\n", seed);
        tgs_random_shuffle(nodes.begin(), nodes.end(), unif_rand);
        sort(nodes.begin(), nodes.begin() + num_kid_nodes);   // this helps to reduce memory cache misses later
        for (unsigned i = 0; i < num_kid_nodes; ++i) {
            node_selected[nodes[i]] = true;
            __sync_add_and_fetch(psamples + nodes[i], 1);
        }

        vector<unsigned> node2cluster(num_nodes, -1);
        vector<int> nodes1;
        vector<int> nodes2;
        vector<double> weights;

        reassign_weights(pnode1, pnode2, knn, num_edges, node_selected, nodes1, nodes2, weights);
        vdebug("num child edges = %ld, num all edges: %ld\n", nodes1.size(), num_edges);

        g_tgstat->rnd_seed(seed);   // this will make the results reproducible by tgs_graph_cover
        graph2cluster(nodes1.data(), nodes2.data(), weights.data(), nodes1.size(), min_cluster_size, cooling_rate, burn_in, node2cluster.data(), num_nodes);

        for (size_t i = 0; i < num_edges; ++i) {
            int node1 = pnode1[i] - 1;
            int node2 = pnode2[i] - 1;

            if (node_selected[node1] && node_selected[node2] && node2cluster[node1] != -1 && node2cluster[node1] == node2cluster[node2])
                __sync_add_and_fetch(pco_clust + i, 1);
        }

        *pready = 1;
        rexit();
    }
}

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
                strcmp(CHAR(STRING_ELT(rnames, COL1)), COL_NAMES[COL1]) || (!isInteger(VECTOR_ELT(_graph, COL1)) && !isFactor(VECTOR_ELT(_graph, COL1))) ||
                strcmp(CHAR(STRING_ELT(rnames, COL2)), COL_NAMES[COL2]) || (!isInteger(VECTOR_ELT(_graph, COL2)) && !isFactor(VECTOR_ELT(_graph, COL2))) ||
                xlength(VECTOR_ELT(_graph, COL2)) != xlength(VECTOR_ELT(_graph, COL1)) ||
                strcmp(CHAR(STRING_ELT(rnames, WEIGHT)), COL_NAMES[WEIGHT]) || !isReal(VECTOR_ELT(_graph, WEIGHT)) || xlength(VECTOR_ELT(_graph, WEIGHT)) != xlength(VECTOR_ELT(_graph, COL1)))
                verror("\"graph\" argument must be in the format that is returned by tgs_cor_graph function");

            pcol1 = INTEGER(VECTOR_ELT(_graph, COL1));
            pcol2 = INTEGER(VECTOR_ELT(_graph, COL2));
            pweight = REAL(VECTOR_ELT(_graph, WEIGHT));
            num_edges = xlength(VECTOR_ELT(_graph, COL1));
        }

        SEXP rlevels1 = getAttrib(VECTOR_ELT(_graph, 0), R_LevelsSymbol);
        SEXP rlevels2 = getAttrib(VECTOR_ELT(_graph, 1), R_LevelsSymbol);

        if ((rlevels1 != R_NilValue && rlevels2 == R_NilValue) || (rlevels1 == R_NilValue && rlevels2 != R_NilValue) ||
            (rlevels1 != R_NilValue && rlevels2 != R_NilValue && xlength(rlevels1) != xlength(rlevels2)))
            verror("\"graph\" argument must be in the format that is returned by tgs_graph function");

        if ((!isInteger(_min_cluster_size) && !isReal(_min_cluster_size)) || xlength(_min_cluster_size) != 1 || asInteger(_min_cluster_size) < 1)
            verror("\"min_cluster_size\" argument must be a positive integer");

        if ((!isInteger(_cooling) && !isReal(_cooling)) || xlength(_cooling) != 1 || asReal(_cooling) < 1)
            verror("\"cooling\" argument must be a number greater or equal than 1");

        if ((!isInteger(_burn_in) && !isReal(_burn_in)) || xlength(_burn_in) != 1 || asInteger(_burn_in) < 0)
            verror("\"burn_in\" argument must be a positive integer");

        unsigned min_cluster_size = asInteger(_min_cluster_size);
        float cooling_rate = asReal(_cooling);
        unsigned burn_in = asInteger(_burn_in);

        size_t num_nodes = 0;
        if (rlevels1 == R_NilValue) {
            for (size_t i = 0; i < num_edges; ++i) {
                num_nodes = max(num_nodes, (size_t)pcol1[i]);
                num_nodes = max(num_nodes, (size_t)pcol2[i]);
            }
        } else
            num_nodes = xlength(rlevels1);

        vector<unsigned> node2cluster(num_nodes, -1);

        graph2cluster(pcol1, pcol2, pweight, num_edges, min_cluster_size, cooling_rate, burn_in, node2cluster.data(), num_nodes);

        vdebug("Packing the return value\n");

        enum { NODE, CLUSTER, NUM_COLS };
        const char *COL_NAMES[NUM_COLS] = { "node", "cluster" };

        SEXP ranswer, rnode, rcluster, rrownames, rcolnames;

        rprotect(ranswer = RSaneAllocVector(VECSXP, NUM_COLS));
        rprotect(rnode = RSaneAllocVector(INTSXP, num_nodes));
        rprotect(rcluster = RSaneAllocVector(INTSXP, num_nodes));
        rprotect(rcolnames = RSaneAllocVector(STRSXP, NUM_COLS));
        rprotect(rrownames = RSaneAllocVector(INTSXP, num_nodes));

        for (int i = 0; i < NUM_COLS; i++)
            SET_STRING_ELT(rcolnames, i, mkChar(COL_NAMES[i]));

        for (size_t i = 0; i < num_nodes; ++i) {
            INTEGER(rnode)[i] = i + 1;
            INTEGER(rcluster)[i] = node2cluster[i] + 1;
            INTEGER(rrownames)[i] = i + 1;
        }

        if (rlevels1 != R_NilValue) {
            setAttrib(rnode, R_LevelsSymbol, rlevels1);
            setAttrib(rnode, R_ClassSymbol, mkString("factor"));
        }

        SET_VECTOR_ELT(ranswer, NODE, rnode);
        SET_VECTOR_ELT(ranswer, CLUSTER, rcluster);

        setAttrib(ranswer, R_NamesSymbol, rcolnames);
        setAttrib(ranswer, R_ClassSymbol, mkString("data.frame"));
        setAttrib(ranswer, R_RowNamesSymbol, rrownames);

        vdebug("Packing the return value - DONE\n");

        rreturn(ranswer);
    } catch (TGLException &e) {
		rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

    rreturn(R_NilValue);
}

SEXP tgs_graph2cluster_multi_hash(SEXP _graph, SEXP _knn, SEXP _min_cluster_size, SEXP _cooling, SEXP _burn_in, SEXP _p_resamp, SEXP _n_resamp, SEXP _envir)
{
    SEXP answer = R_NilValue;
    unsigned *res = (unsigned *)MAP_FAILED;
    size_t res_sizeof = 0;

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
                strcmp(CHAR(STRING_ELT(rnames, COL1)), COL_NAMES[COL1]) || (!isInteger(VECTOR_ELT(_graph, COL1)) && !isFactor(VECTOR_ELT(_graph, COL1))) ||
                strcmp(CHAR(STRING_ELT(rnames, COL2)), COL_NAMES[COL2]) || (!isInteger(VECTOR_ELT(_graph, COL2)) && !isFactor(VECTOR_ELT(_graph, COL2))) ||
                xlength(VECTOR_ELT(_graph, COL2)) != xlength(VECTOR_ELT(_graph, COL1)) ||
                strcmp(CHAR(STRING_ELT(rnames, WEIGHT)), COL_NAMES[WEIGHT]) || !isReal(VECTOR_ELT(_graph, WEIGHT)) || xlength(VECTOR_ELT(_graph, WEIGHT)) != xlength(VECTOR_ELT(_graph, COL1)))
                verror("\"graph\" argument must be in the format that is returned by tgs_graph function");

            pcol1 = INTEGER(VECTOR_ELT(_graph, COL1));
            pcol2 = INTEGER(VECTOR_ELT(_graph, COL2));
            pweight = REAL(VECTOR_ELT(_graph, WEIGHT));
            num_edges = xlength(VECTOR_ELT(_graph, COL1));
        }

        SEXP rlevels1 = getAttrib(VECTOR_ELT(_graph, 0), R_LevelsSymbol);
        SEXP rlevels2 = getAttrib(VECTOR_ELT(_graph, 1), R_LevelsSymbol);

        if ((rlevels1 != R_NilValue && rlevels2 == R_NilValue) || (rlevels1 == R_NilValue && rlevels2 != R_NilValue) ||
            (rlevels1 != R_NilValue && rlevels2 != R_NilValue && xlength(rlevels1) != xlength(rlevels2)))
            verror("\"graph\" argument must be in the format that is returned by tgs_graph function");

        if ((!isInteger(_min_cluster_size) && !isReal(_min_cluster_size)) || xlength(_min_cluster_size) != 1 || asInteger(_min_cluster_size) < 1)
            verror("\"min_cluster_size\" argument must be a positive integer");

        if ((!isInteger(_cooling) && !isReal(_cooling)) || xlength(_cooling) != 1 || asReal(_cooling) < 1)
            verror("\"cooling\" argument must be a number greater or equal than 1");

        if ((!isInteger(_burn_in) && !isReal(_burn_in)) || xlength(_burn_in) != 1 || asInteger(_burn_in) < 0)
            verror("\"burn_in\" argument must be a positive integer");

        if ((!isNull(_knn) && ((!isReal(_knn) && !isInteger(_knn)) || xlength(_knn) != 1)) || asInteger(_knn) < 1)
            verror("\"knn\" argument must be a positive integer");

        if ((!isInteger(_n_resamp) && !isReal(_n_resamp)) || xlength(_n_resamp) != 1 || asInteger(_n_resamp) < 1 || asInteger(_n_resamp) > 0xffff)
            verror("\"n_resamp\" argument must be a positive integer withn [1, %d] range", 0xffff);

        if ((!isInteger(_p_resamp) && !isReal(_p_resamp)) || xlength(_p_resamp) != 1 || asReal(_p_resamp) > 1 || asReal(_p_resamp) <= 0)
            verror("\"p_resamp\" argument must be a number in (0,1] range");

        unsigned min_cluster_size = asInteger(_min_cluster_size);
        float cooling_rate = asReal(_cooling);
        unsigned burn_in = asInteger(_burn_in);
        unsigned knn = asInteger(_knn);
        unsigned n_resamp = asInteger(_n_resamp);
        double p_resamp = asReal(_p_resamp);
        size_t num_nodes = 0;

        if (rlevels1 == R_NilValue) {
            for (size_t i = 0; i < num_edges; ++i) {
                num_nodes = max(num_nodes, (size_t)pcol1[i]);
                num_nodes = max(num_nodes, (size_t)pcol2[i]);
            }
        } else
            num_nodes = xlength(rlevels1);

        int max_num_kids = min((int)n_resamp, max(g_tgstat->num_processes() - 1, 1));
//max_num_kids = 1;
        int num_kids_launched = 0;
        int num_kids_finished = 0;

        unordered_map<pair<unsigned, unsigned>, unsigned short> co_cluster;
        vector<unsigned short> co_cluster_diag(num_nodes, 0);
        vector<unsigned short> node2sample_cnt(num_nodes, 0);

        vdebug("Allocating shared memory for results\n");

        // shared memory is addressed as following:
        //    [child process 0]
        //       unsigned - status: 1 when the data is ready for read, otherwise 0
        //       unsigned - number of clusters
        //       unsigned - cluster of node 0
        //       ...
        //       unsigned - cluster of node i
        //    ...
        //    [child process j]
        //       ...

        size_t kid_res_sizeof = sizeof(unsigned) + sizeof(unsigned) + sizeof(unsigned) * num_nodes;

        res_sizeof = kid_res_sizeof * max_num_kids;
        res = (unsigned *)mmap(NULL, res_sizeof, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);

        if (res == (unsigned *)MAP_FAILED)
            verror("Failed to allocate shared memory: %s", strerror(errno));

        ProgressReporter progress;
        progress.init(n_resamp, 1);

        vdebug("num_processes: %d\n", max_num_kids);
        TGStat::prepare4multitasking();

        for (int i = 0; i < max_num_kids; ++i) {
            unsigned *pready = (unsigned *)((char *)res + i * kid_res_sizeof);
            *pready = 0;
            launch_kid_hash(pcol1, pcol2, pweight, num_nodes, num_edges, res, i, knn, p_resamp, min_cluster_size, cooling_rate, burn_in, unif_rand() * 0xffff);
            ++num_kids_launched;
        }

        while (num_kids_finished < n_resamp) {
            while (TGStat::wait_for_kid(3000))
                progress.report(0);

            for (int ikid = 0; ikid < max_num_kids; ++ikid) {
                unsigned *pready = (unsigned *)((char *)res + ikid * kid_res_sizeof);

                if (*pready) {
                    vdebug("======== Result at slot %d is ready (status: %d)\n", ikid, *pready);

                    unsigned *num_clusters = pready + 1;
                    unsigned *clusters = num_clusters + 1;
                    size_t num_pairs = 0;
                    vector<vector<unsigned>> cluster2nodes(*num_clusters);

                    for (size_t i = 0; i < num_nodes; ++i) {
                        if (clusters[i] != -2) {
                            ++node2sample_cnt[i];
                            if (clusters[i] != -1) {
                                cluster2nodes[clusters[i]].push_back(i);
                                co_cluster_diag[i]++;
                            }
                        }
                    }

                    for (auto &cluster_nodes : cluster2nodes)
                        num_pairs += cluster_nodes.size() * (cluster_nodes.size() - 1) / 2;

                    vdebug("Updating co_cluster... Num pairs: %ld, num_clusters: %d\n", num_pairs, *num_clusters, co_cluster.size());
                    vdebug("co_cluster.size() = %ld, load factor: %g, max load factor: %g\n", co_cluster.size(), co_cluster.load_factor(), co_cluster.max_load_factor());

                    for (unsigned cluster = 0; cluster < *num_clusters; ++cluster) {
                        const vector<unsigned> &cluster_nodes = cluster2nodes[cluster];
                        for (auto i = cluster_nodes.begin(); i < cluster_nodes.end(); ++i) {
                            for (auto j = i + 1; j < cluster_nodes.end(); ++j) {
                                unsigned node1 = *i;
                                unsigned node2 = *j;
                                if (node1 > node2)
                                    swap(node1, node2);
                                co_cluster[{node1, node2}]++;
                            }
                        }
                    }

                    *pready = 0;
                    ++num_kids_finished;
                    progress.report(1);
                    vdebug("Num processes ended: %d\n", num_kids_finished);

                    if (num_kids_launched < n_resamp) {
                        launch_kid_hash(pcol1, pcol2, pweight, num_nodes, num_edges, res, ikid, knn, p_resamp, min_cluster_size, cooling_rate, burn_in, unif_rand() * 0xffff);
                        ++num_kids_launched;
                    }
                }
            }
        }

        while (TGStat::wait_for_kids(3000))
            progress.report(0);

        progress.report_last();

        vdebug("Packing the result...\n");

        enum { NODE1, NODE2, CNT, NUM_COLS };
        const char *COL_NAMES[NUM_COLS] = { "node1", "node2", "cnt" };

        SEXP rco_clust, rsamples, rnode1, rnode2, rcount, rrownames, rcolnames, rnames;
        unsigned co_cluster_diag_size = 0;

        for (auto cnt : co_cluster_diag) {
            if (cnt)
                ++co_cluster_diag_size;
        }

        rprotect(answer = RSaneAllocVector(VECSXP, 2));

        rprotect(rco_clust = RSaneAllocVector(VECSXP, NUM_COLS));
        rprotect(rnode1 = RSaneAllocVector(INTSXP, co_cluster_diag_size + co_cluster.size()));
        rprotect(rnode2 = RSaneAllocVector(INTSXP, co_cluster_diag_size + co_cluster.size()));
        rprotect(rcount = RSaneAllocVector(INTSXP, co_cluster_diag_size + co_cluster.size()));
        rprotect(rcolnames = RSaneAllocVector(STRSXP, NUM_COLS));
        rprotect(rrownames = RSaneAllocVector(INTSXP, co_cluster_diag_size + co_cluster.size()));

        {
            int i = 0;
            for (auto ico_cluster = co_cluster_diag.begin(); ico_cluster != co_cluster_diag.end(); ++ico_cluster) {
                if (*ico_cluster) {
                    INTEGER(rnode1)[i] = INTEGER(rnode2)[i] = ico_cluster - co_cluster_diag.begin() + 1;
                    INTEGER(rcount)[i] = *ico_cluster;
                    INTEGER(rrownames)[i] = i + 1;
                    ++i;
                }
            }
        }

        {
            int i = co_cluster_diag_size;
            for (auto ico_cluster = co_cluster.begin(); ico_cluster != co_cluster.end(); ++ico_cluster, ++i) {
                INTEGER(rnode1)[i] = ico_cluster->first.first + 1;
                INTEGER(rnode2)[i] = ico_cluster->first.second + 1;
                INTEGER(rcount)[i] = ico_cluster->second;
                INTEGER(rrownames)[i] = i + 1;
            }
        }

        for (int i = 0; i < NUM_COLS; i++)
            SET_STRING_ELT(rcolnames, i, mkChar(COL_NAMES[i]));

        rprotect(rsamples = RSaneAllocVector(INTSXP, num_nodes));

        for (size_t i = 0; i < num_nodes; ++i)
            INTEGER(rsamples)[i] = node2sample_cnt[i];

        if (rlevels1 != R_NilValue) {
            setAttrib(rnode1, R_LevelsSymbol, rlevels1);
            setAttrib(rnode1, R_ClassSymbol, mkString("factor"));
            setAttrib(rnode2, R_LevelsSymbol, rlevels1);
            setAttrib(rnode2, R_ClassSymbol, mkString("factor"));
            setAttrib(rsamples, R_NamesSymbol, rlevels1);
        }

        SET_VECTOR_ELT(rco_clust, NODE1, rnode1);
        SET_VECTOR_ELT(rco_clust, NODE2, rnode2);
        SET_VECTOR_ELT(rco_clust, CNT, rcount);

        setAttrib(rco_clust, R_NamesSymbol, rcolnames);
        setAttrib(rco_clust, R_ClassSymbol, mkString("data.frame"));
        setAttrib(rco_clust, R_RowNamesSymbol, rrownames);

        rprotect(rnames = RSaneAllocVector(STRSXP, 2));
        SET_STRING_ELT(rnames, 0, mkChar("co_cluster"));
        SET_STRING_ELT(rnames, 1, mkChar("samples"));
        setAttrib(answer, R_NamesSymbol, rnames);

        SET_VECTOR_ELT(answer, 0, rco_clust);
        SET_VECTOR_ELT(answer, 1, rsamples);
    } catch (TGLException &e) {
        if (!TGStat::is_kid() && res != (unsigned *)MAP_FAILED) {
            munmap(res, res_sizeof);
            res = (unsigned *)MAP_FAILED;
        }
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

    if (!TGStat::is_kid() && res != (unsigned *)MAP_FAILED) {
        munmap(res, res_sizeof);
        res = (unsigned *)MAP_FAILED;
    }
    rreturn(answer);
}

SEXP tgs_graph2cluster_multi_full(SEXP _graph, SEXP _knn, SEXP _min_cluster_size, SEXP _cooling, SEXP _burn_in, SEXP _p_resamp, SEXP _n_resamp, SEXP _envir)
{
    SEXP answer = R_NilValue;
    unsigned *res = (unsigned *)MAP_FAILED;
    size_t res_sizeof = 0;

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
                strcmp(CHAR(STRING_ELT(rnames, COL1)), COL_NAMES[COL1]) || (!isInteger(VECTOR_ELT(_graph, COL1)) && !isFactor(VECTOR_ELT(_graph, COL1))) ||
                strcmp(CHAR(STRING_ELT(rnames, COL2)), COL_NAMES[COL2]) || (!isInteger(VECTOR_ELT(_graph, COL2)) && !isFactor(VECTOR_ELT(_graph, COL2))) ||
                xlength(VECTOR_ELT(_graph, COL2)) != xlength(VECTOR_ELT(_graph, COL1)) ||
                strcmp(CHAR(STRING_ELT(rnames, WEIGHT)), COL_NAMES[WEIGHT]) || !isReal(VECTOR_ELT(_graph, WEIGHT)) || xlength(VECTOR_ELT(_graph, WEIGHT)) != xlength(VECTOR_ELT(_graph, COL1)))
                verror("\"graph\" argument must be in the format that is returned by tgs_cor_graph function");

            pcol1 = INTEGER(VECTOR_ELT(_graph, COL1));
            pcol2 = INTEGER(VECTOR_ELT(_graph, COL2));
            pweight = REAL(VECTOR_ELT(_graph, WEIGHT));
            num_edges = xlength(VECTOR_ELT(_graph, COL1));
        }

        SEXP rlevels1 = getAttrib(VECTOR_ELT(_graph, 0), R_LevelsSymbol);
        SEXP rlevels2 = getAttrib(VECTOR_ELT(_graph, 1), R_LevelsSymbol);

        if ((rlevels1 != R_NilValue && rlevels2 == R_NilValue) || (rlevels1 == R_NilValue && rlevels2 != R_NilValue) ||
            (rlevels1 != R_NilValue && rlevels2 != R_NilValue && xlength(rlevels1) != xlength(rlevels2)))
            verror("\"graph\" argument must be in the format that is returned by tgs_graph function");

        if ((!isInteger(_min_cluster_size) && !isReal(_min_cluster_size)) || xlength(_min_cluster_size) != 1 || asInteger(_min_cluster_size) < 1)
            verror("\"min_cluster_size\" argument must be a positive integer");

        if ((!isInteger(_cooling) && !isReal(_cooling)) || xlength(_cooling) != 1 || asReal(_cooling) < 1)
            verror("\"cooling\" argument must be a number greater or equal than 1");

        if ((!isInteger(_burn_in) && !isReal(_burn_in)) || xlength(_burn_in) != 1 || asInteger(_burn_in) < 0)
            verror("\"burn_in\" argument must be a positive integer");

        if ((!isNull(_knn) && ((!isReal(_knn) && !isInteger(_knn)) || xlength(_knn) != 1)) || asInteger(_knn) < 1)
            verror("\"knn\" argument must be a positive integer");

        if ((!isInteger(_n_resamp) && !isReal(_n_resamp)) || xlength(_n_resamp) != 1 || asInteger(_n_resamp) < 1 || asInteger(_n_resamp) > 0xffff)
            verror("\"n_resamp\" argument must be a positive integer withn [1, %d] range", 0xffff);

        if ((!isInteger(_p_resamp) && !isReal(_p_resamp)) || xlength(_p_resamp) != 1 || asReal(_p_resamp) > 1 || asReal(_p_resamp) <= 0)
            verror("\"p_resamp\" argument must be a number in (0,1] range");

        unsigned min_cluster_size = asInteger(_min_cluster_size);
        float cooling_rate = asReal(_cooling);
        unsigned burn_in = asInteger(_burn_in);
        unsigned knn = asInteger(_knn);
        unsigned n_resamp = asInteger(_n_resamp);
        double p_resamp = asReal(_p_resamp);

        size_t num_nodes = 0;
        if (rlevels1 == R_NilValue) {
            for (size_t i = 0; i < num_edges; ++i) {
                num_nodes = max(num_nodes, (size_t)pcol1[i]);
                num_nodes = max(num_nodes, (size_t)pcol2[i]);
            }
        } else
            num_nodes = xlength(rlevels1);

        int max_num_kids = min((int)n_resamp, max(g_tgstat->num_processes() - 1, 1));
//max_num_kids = 2;
        int num_kids_launched = 0;
        int num_kids_finished = 0;
        size_t co_clust_size = ((num_nodes + 1) * num_nodes) / 2;

        vdebug("Allocating shared memory for results\n");

        // shared memory is addressed as following:
        //    [child process 0]
        //       unsigned - status: 1 when the data is ready for read, otherwise 0
        //    [child process j]
        //       ...
        //    short x (num_nodes + 1) x num_nodes / 2 - lower triangle of co_clust of matrix including the diagonal
        //    short x num_nodes                       - samples vector

        res_sizeof = sizeof(unsigned) * max_num_kids + sizeof(short) * co_clust_size + sizeof(short) * num_nodes;
        res = (unsigned *)mmap(NULL, res_sizeof, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);
        memset(res, 0, res_sizeof);

        unsigned *pready = (unsigned *)res;
        unsigned short *pco_clust = (unsigned short *)((unsigned *)res + max_num_kids);
        unsigned short *psamples = pco_clust + co_clust_size;

        if (res == (unsigned *)MAP_FAILED)
            verror("Failed to allocate shared memory: %s", strerror(errno));

        ProgressReporter progress;
        progress.init(n_resamp, 1);

        vdebug("num_processes: %d\n", max_num_kids);
        TGStat::prepare4multitasking();

        for (int i = 0; i < max_num_kids; ++i) {
            pready[i] = 0;
            vdebug("Launching a working process at slot %d\n", i);
            launch_kid_full(pcol1, pcol2, pweight, num_nodes, num_edges, pready + i, pco_clust, psamples,
                            knn, p_resamp, min_cluster_size, cooling_rate, burn_in, unif_rand() * 0xffff);
            ++num_kids_launched;
        }

        while (num_kids_finished < n_resamp) {
            while (TGStat::wait_for_kid(3000))
                progress.report(0);

            for (int ikid = 0; ikid < max_num_kids; ++ikid) {
                if (pready[ikid]) {
                    vdebug("======== Result at slot %d is ready (status: %d)\n", ikid, pready[ikid]);
                    pready[ikid] = 0;
                    ++num_kids_finished;
                    progress.report(1);
                    vdebug("Num processes ended: %d\n", num_kids_finished);

                    if (num_kids_launched < n_resamp) {
                        vdebug("Launching a working process at slot %d\n", ikid);
                        launch_kid_full(pcol1, pcol2, pweight, num_nodes, num_edges, pready + ikid, pco_clust, psamples,
                                        knn, p_resamp, min_cluster_size, cooling_rate, burn_in, unif_rand() * 0xffff);
                        ++num_kids_launched;
                    }
                }
            }
        }

        while (TGStat::wait_for_kids(3000))
            progress.report(0);

        progress.report_last();

        vdebug("Packing the result...\n");

        enum { NODE1, NODE2, CNT, NUM_COLS };
        const char *COL_NAMES[NUM_COLS] = { "node1", "node2", "cnt" };

        SEXP rco_clust, rsamples, rnode1, rnode2, rcount, rrownames, rcolnames, rnames;
        size_t co_clust_nz_size = 0;

        for (size_t i = 0; i < co_clust_size; ++i) {
            if (pco_clust[i])
                ++co_clust_nz_size;
        }

        rprotect(answer = RSaneAllocVector(VECSXP, 2));

        rprotect(rco_clust = RSaneAllocVector(VECSXP, NUM_COLS));
        rprotect(rnode1 = RSaneAllocVector(INTSXP, co_clust_nz_size));
        rprotect(rnode2 = RSaneAllocVector(INTSXP, co_clust_nz_size));
        rprotect(rcount = RSaneAllocVector(INTSXP, co_clust_nz_size));
        rprotect(rcolnames = RSaneAllocVector(STRSXP, NUM_COLS));
        rprotect(rrownames = RSaneAllocVector(INTSXP, co_clust_nz_size));

        {
            size_t row = 0;

            for (size_t node1 = 0; node1 < num_nodes; ++node1) {
                unsigned short *p = pco_clust + (node1 * (2 * num_nodes - node1 - 1)) / 2;
                for (size_t node2 = node1; node2 < num_nodes; ++node2) {
                    if (p[node2]) {
                        INTEGER(rnode1)[row] = node1 + 1;
                        INTEGER(rnode2)[row] = node2 + 1;
                        INTEGER(rcount)[row] = p[node2];
                        INTEGER(rrownames)[row] = row + 1;
                        ++row;
                    }
                }
            }
        }

        for (int i = 0; i < NUM_COLS; i++)
            SET_STRING_ELT(rcolnames, i, mkChar(COL_NAMES[i]));

        SET_VECTOR_ELT(rco_clust, NODE1, rnode1);
        SET_VECTOR_ELT(rco_clust, NODE2, rnode2);
        SET_VECTOR_ELT(rco_clust, CNT, rcount);

        setAttrib(rco_clust, R_NamesSymbol, rcolnames);
        setAttrib(rco_clust, R_ClassSymbol, mkString("data.frame"));
        setAttrib(rco_clust, R_RowNamesSymbol, rrownames);

        rprotect(rsamples = RSaneAllocVector(INTSXP, num_nodes));

        for (size_t i = 0; i < num_nodes; ++i)
            INTEGER(rsamples)[i] = psamples[i];

        if (rlevels1 != R_NilValue) {
            setAttrib(rnode1, R_LevelsSymbol, rlevels1);
            setAttrib(rnode1, R_ClassSymbol, mkString("factor"));
            setAttrib(rnode2, R_LevelsSymbol, rlevels1);
            setAttrib(rnode2, R_ClassSymbol, mkString("factor"));
            setAttrib(rsamples, R_NamesSymbol, rlevels1);
        }

        rprotect(rnames = RSaneAllocVector(STRSXP, 2));
        SET_STRING_ELT(rnames, 0, mkChar("co_cluster"));
        SET_STRING_ELT(rnames, 1, mkChar("samples"));
        setAttrib(answer, R_NamesSymbol, rnames);

        SET_VECTOR_ELT(answer, 0, rco_clust);
        SET_VECTOR_ELT(answer, 1, rsamples);
    } catch (TGLException &e) {
        if (!TGStat::is_kid() && res != (unsigned *)MAP_FAILED) {
            munmap(res, res_sizeof);
            res = (unsigned *)MAP_FAILED;
        }
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

    if (!TGStat::is_kid() && res != (unsigned *)MAP_FAILED) {
        munmap(res, res_sizeof);
        res = (unsigned *)MAP_FAILED;
    }
    rreturn(answer);
}

SEXP tgs_graph2cluster_multi_edges(SEXP _graph, SEXP _knn, SEXP _min_cluster_size, SEXP _cooling, SEXP _burn_in, SEXP _p_resamp, SEXP _n_resamp, SEXP _envir)
{
    SEXP answer = R_NilValue;
    unsigned *res = (unsigned *)MAP_FAILED;
    size_t res_sizeof = 0;

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
                strcmp(CHAR(STRING_ELT(rnames, COL1)), COL_NAMES[COL1]) || (!isInteger(VECTOR_ELT(_graph, COL1)) && !isFactor(VECTOR_ELT(_graph, COL1))) ||
                strcmp(CHAR(STRING_ELT(rnames, COL2)), COL_NAMES[COL2]) || (!isInteger(VECTOR_ELT(_graph, COL2)) && !isFactor(VECTOR_ELT(_graph, COL2))) ||
                xlength(VECTOR_ELT(_graph, COL2)) != xlength(VECTOR_ELT(_graph, COL1)) ||
                strcmp(CHAR(STRING_ELT(rnames, WEIGHT)), COL_NAMES[WEIGHT]) || !isReal(VECTOR_ELT(_graph, WEIGHT)) || xlength(VECTOR_ELT(_graph, WEIGHT)) != xlength(VECTOR_ELT(_graph, COL1)))
                verror("\"graph\" argument must be in the format that is returned by tgs_cor_graph function");

            pcol1 = INTEGER(VECTOR_ELT(_graph, COL1));
            pcol2 = INTEGER(VECTOR_ELT(_graph, COL2));
            pweight = REAL(VECTOR_ELT(_graph, WEIGHT));
            num_edges = xlength(VECTOR_ELT(_graph, COL1));
        }

        SEXP rlevels1 = getAttrib(VECTOR_ELT(_graph, 0), R_LevelsSymbol);
        SEXP rlevels2 = getAttrib(VECTOR_ELT(_graph, 1), R_LevelsSymbol);

        if ((rlevels1 != R_NilValue && rlevels2 == R_NilValue) || (rlevels1 == R_NilValue && rlevels2 != R_NilValue) ||
            (rlevels1 != R_NilValue && rlevels2 != R_NilValue && xlength(rlevels1) != xlength(rlevels2)))
            verror("\"graph\" argument must be in the format that is returned by tgs_graph function");

        if ((!isInteger(_min_cluster_size) && !isReal(_min_cluster_size)) || xlength(_min_cluster_size) != 1 || asInteger(_min_cluster_size) < 1)
            verror("\"min_cluster_size\" argument must be a positive integer");

        if ((!isInteger(_cooling) && !isReal(_cooling)) || xlength(_cooling) != 1 || asReal(_cooling) < 1)
            verror("\"cooling\" argument must be a number greater or equal than 1");

        if ((!isInteger(_burn_in) && !isReal(_burn_in)) || xlength(_burn_in) != 1 || asInteger(_burn_in) < 0)
            verror("\"burn_in\" argument must be a positive integer");

        if ((!isNull(_knn) && ((!isReal(_knn) && !isInteger(_knn)) || xlength(_knn) != 1)) || asInteger(_knn) < 1)
            verror("\"knn\" argument must be a positive integer");

        if ((!isInteger(_n_resamp) && !isReal(_n_resamp)) || xlength(_n_resamp) != 1 || asInteger(_n_resamp) < 1 || asInteger(_n_resamp) > 0xffff)
            verror("\"n_resamp\" argument must be a positive integer withn [1, %d] range", 0xffff);

        if ((!isInteger(_p_resamp) && !isReal(_p_resamp)) || xlength(_p_resamp) != 1 || asReal(_p_resamp) > 1 || asReal(_p_resamp) <= 0)
            verror("\"p_resamp\" argument must be a number in (0,1] range");

        unsigned min_cluster_size = asInteger(_min_cluster_size);
        float cooling_rate = asReal(_cooling);
        unsigned burn_in = asInteger(_burn_in);
        unsigned knn = asInteger(_knn);
        unsigned n_resamp = asInteger(_n_resamp);
        double p_resamp = asReal(_p_resamp);

        size_t num_nodes = 0;
        if (rlevels1 == R_NilValue) {
            for (size_t i = 0; i < num_edges; ++i) {
                num_nodes = max(num_nodes, (size_t)pcol1[i]);
                num_nodes = max(num_nodes, (size_t)pcol2[i]);
            }
        } else
            num_nodes = xlength(rlevels1);

        int max_num_kids = min((int)n_resamp, max(g_tgstat->num_processes() - 1, 1));
//max_num_kids = 1;
        int num_kids_launched = 0;
        int num_kids_finished = 0;

        vdebug("Allocating shared memory for results\n");

        // shared memory is addressed as following:
        //    [child process 0]
        //       unsigned - status: 1 when the data is ready for read, otherwise 0
        //    [child process j]
        //       ...
        //    short x num_edges - co_clust array
        //    short x num_nodes - samples vector

        res_sizeof = sizeof(unsigned) * max_num_kids + 2 * sizeof(short) * num_edges;
        res = (unsigned *)mmap(NULL, res_sizeof, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);
        memset(res, 0, res_sizeof);

        unsigned *pready = (unsigned *)res;
        unsigned short *pco_clust = (unsigned short *)((unsigned *)res + max_num_kids);
        unsigned short *psamples = pco_clust + num_edges;

        if (res == (unsigned *)MAP_FAILED)
            verror("Failed to allocate shared memory: %s", strerror(errno));

        ProgressReporter progress;
        progress.init(n_resamp, 1);

        vdebug("num_processes: %d\n", max_num_kids);
        TGStat::prepare4multitasking();

        for (int i = 0; i < max_num_kids; ++i) {
            pready[i] = 0;
            vdebug("Launching a working process at slot %d\n", i);
            launch_kid_edges(pcol1, pcol2, pweight, num_nodes, num_edges, pready + i, pco_clust, psamples,
                            knn, p_resamp, min_cluster_size, cooling_rate, burn_in, unif_rand() * 0xffff);
            ++num_kids_launched;
        }

        while (num_kids_finished < n_resamp) {
            while (TGStat::wait_for_kid(3000))
                progress.report(0);

            for (int ikid = 0; ikid < max_num_kids; ++ikid) {
                if (pready[ikid]) {
                    vdebug("======== Result at slot %d is ready (status: %d)\n", ikid, pready[ikid]);
                    pready[ikid] = 0;
                    ++num_kids_finished;
                    progress.report(1);
                    vdebug("Num processes ended: %d\n", num_kids_finished);

                    if (num_kids_launched < n_resamp) {
                        vdebug("Launching a working process at slot %d\n", ikid);
                        launch_kid_edges(pcol1, pcol2, pweight, num_nodes, num_edges, pready + ikid, pco_clust, psamples,
                                         knn, p_resamp, min_cluster_size, cooling_rate, burn_in, unif_rand() * 0xffff);
                        ++num_kids_launched;
                    }
                }
            }
        }

        while (TGStat::wait_for_kids(3000))
            progress.report(0);

        progress.report_last();

        vdebug("Packing the result...\n");

        enum { NODE1, NODE2, CNT, NUM_COLS };
        const char *COL_NAMES[NUM_COLS] = { "node1", "node2", "cnt" };

        SEXP rco_clust, rsamples, rnode1, rnode2, rcount, rrownames, rcolnames, rnames;

        rprotect(answer = RSaneAllocVector(VECSXP, 2));

        size_t co_clust_nz_size = 0;

        for (size_t i = 0; i < num_edges; ++i) {
            if (pco_clust[i])
                ++co_clust_nz_size;
        }

        rprotect(answer = RSaneAllocVector(VECSXP, 2));

        rprotect(rco_clust = RSaneAllocVector(VECSXP, NUM_COLS));
        rprotect(rnode1 = RSaneAllocVector(INTSXP, co_clust_nz_size));
        rprotect(rnode2 = RSaneAllocVector(INTSXP, co_clust_nz_size));
        rprotect(rcount = RSaneAllocVector(INTSXP, co_clust_nz_size));
        rprotect(rcolnames = RSaneAllocVector(STRSXP, NUM_COLS));
        rprotect(rrownames = RSaneAllocVector(INTSXP, co_clust_nz_size));

        {
            size_t row = 0;

            for (size_t i = 0; i < num_edges; ++i) {
                if (pco_clust[i]) {
                    INTEGER(rnode1)[row] = pcol1[i];
                    INTEGER(rnode2)[row] = pcol2[i];
                    INTEGER(rcount)[row] = pco_clust[i];
                    INTEGER(rrownames)[row] = row + 1;
                    ++row;
                }
            }
        }

        for (int i = 0; i < NUM_COLS; i++)
            SET_STRING_ELT(rcolnames, i, mkChar(COL_NAMES[i]));

        SET_VECTOR_ELT(rco_clust, NODE1, rnode1);
        SET_VECTOR_ELT(rco_clust, NODE2, rnode2);
        SET_VECTOR_ELT(rco_clust, CNT, rcount);

        setAttrib(rco_clust, R_NamesSymbol, rcolnames);
        setAttrib(rco_clust, R_ClassSymbol, mkString("data.frame"));
        setAttrib(rco_clust, R_RowNamesSymbol, rrownames);

        rprotect(rsamples = RSaneAllocVector(INTSXP, num_nodes));

        for (size_t i = 0; i < num_nodes; ++i)
            INTEGER(rsamples)[i] = psamples[i];

        if (rlevels1 != R_NilValue) {
            setAttrib(rnode1, R_LevelsSymbol, rlevels1);
            setAttrib(rnode1, R_ClassSymbol, mkString("factor"));
            setAttrib(rnode2, R_LevelsSymbol, rlevels1);
            setAttrib(rnode2, R_ClassSymbol, mkString("factor"));
            setAttrib(rsamples, R_NamesSymbol, rlevels1);
        }

        rprotect(rnames = RSaneAllocVector(STRSXP, 2));
        SET_STRING_ELT(rnames, 0, mkChar("co_cluster"));
        SET_STRING_ELT(rnames, 1, mkChar("samples"));
        setAttrib(answer, R_NamesSymbol, rnames);

        SET_VECTOR_ELT(answer, 0, rco_clust);
        SET_VECTOR_ELT(answer, 1, rsamples);
    } catch (TGLException &e) {
        if (!TGStat::is_kid() && res != (unsigned *)MAP_FAILED) {
            munmap(res, res_sizeof);
            res = (unsigned *)MAP_FAILED;
        }
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

    if (!TGStat::is_kid() && res != (unsigned *)MAP_FAILED) {
        munmap(res, res_sizeof);
        res = (unsigned *)MAP_FAILED;
    }
    rreturn(answer);
}

}

