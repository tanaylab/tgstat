\name{tgs_cor_graph}
\alias{tgs_cor_graph}
\title{
Builds directed graph of correlations
}
\description{
Builds directed graph of correlations where the nodes are the matrix columns.
}
\usage{
tgs_cor_graph(x, knn, k_expand, k_alpha = 10, k_beta = 3,
              pairwise.complete.obs = F, spearman = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{numeric matrix}
  \item{knn}{maximal node degree}
  \item{knn_expand}{see below}
  \item{knn_alpha}{see below}
  \item{knn_beta}{see below}
  \item{pairwise.complete.obs}{see 'tgs_cor'}
  \item{spearman}{see 'tgs_cor'}
}
\details{
This function builds a directed graph where the nodes are the matrix columns which are connected by the edges that represent the highest correlation achieved between different columns. Edges of the high degree nodes are pruned.

Here is how 'tgs_cor_graph' works in detail:

1. Find 'knn * k_alpha' highest correlations for 'x' by calling tgs_cor_knn(x, knn * k_alpha, pairwise.complete.objs = pairwise.complete.objs, spearman = spearman, threshold = 0). A pair of columns i, j that appears in the result implies the edge in the graph from i to j: edge(i,j). Let the rank of correlation of i and j be rank(i,j).
2. Calculate symmetrised rank of i and j: sym_rank(i,j) = rank(i,j) * rank(j,i). If one of the ranks is missing from the previous result sym_rank is set to NA.
3. Prune the edges: remove edge(i,j) if sym_rank(i,j) == NA OR  sym_rank(i,j) < knn * knn * k_expand
4. Prune excessive incoming edges: remove edge(i,j) if more than knn * k_beta edges of type edge(node,j) exist and sym_rank(i,j) is higher than sym_rank(node,j) for node != j.
5. Prune excessive outgoing edges: remove edge(i,j) if more than knn edges of type edge(i,node) exist and sym_rank(i,j) is higher than sym_rank(i,node) for node != i.
}
\value{
The graph edges are returned in a data frame, with the weight of each edge varying between 1 and 0. edge(i,j) receives weight 1 if its sym_rank is the lowest among all edges of type edge(i,node). Accordingly the edge with the highest sym_rank receives weight 0. Formally defined: weight(i,j) = 1 - (place(i,j) - 1) / (knn - 1), where place(i,j) is the location of edge(i,j) within the sorted set of edges outgoing from i, i.e. edge(i,node). The sort is done by sym_rank of the edges. If knn == 1, the weight is 1.
}

\examples{
set.seed(seed = 0)
rows <- 100
cols <- 1000
vals <- sample(1 : (rows * cols / 2), rows * cols, replace = TRUE)
m <- matrix(vals, nrow = rows, ncol = cols)
m[sample(1 : (rows * cols), rows * cols / 1000)] <- NA

r<-tgs_cor_graph(m, 3, 10, pairwise.complete.obs = F, spearman = T)
}
\keyword{ ~correlation }

