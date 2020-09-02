#ifndef GRAPH2CLUSTER_H_INCLUDED
#define GRAPH2CLUSTER_H_INCLUDED

#include <vector>

// node1 / node2 - 1-based
// node2cluster - if node is unassigned node2cluster[node] will be -1
// returns number of clusters created
unsigned graph2cluster(const int *pnode1, const int *pnode2, const double *pweight, uint64_t num_edges,
                       unsigned min_cluster_size, float cooling_rate, unsigned burn_in, unsigned *node2cluster, uint64_t num_nodes);

#endif

