\name{tgs_graph_cover_resample}
\alias{tgs_graph_cover_resample}
\title{
Clusters directed graph multiple times with randomized sample subset
}
\description{
Clusters directed graph multiple times with randomized sample subset.
}
\usage{
tgs_graph_cover_resample(graph, knn, min_cluster_size,
                         cooling = 1.05, burn_in = 10,
                         p_resamp = 0.75, n_resamp = 500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{graph}{directed graph in the format returned by tgs_cor_graph}
  \item{knn}{maximal number of edges used per node for each sample subset}
  \item{min_cluster_size}{used to determine the candidates for seeding (= min weight)}
  \item{cooling}{factor that is used to gradually increase the chance of a node to stay in the cluster}
  \item{burn_in}{number of node reassignments after which cooling is applied}
  \item{p_resamp}{fraction of total number of nodes used in each sample subset}
  \item{n_resamp}{number iterations the clustering is run on different sample subsets}
}
\details{
The algorithm is explained in a separate document.
}
\value{
A list with two members.
The first member is a data frame with 3 columns: "node1", "node2" and "cnt". "cnt" indicates the number of times "node1" and "node2" appeared in the same cluster.
The second member of the list is a vector of [number of nodes] length reflecting how many times each node was used in the subset.
}

\examples{
set.seed(seed = 0)
rows <- 1000
cols <- 2000
vals <- sample(1 : (rows * cols / 2), rows * cols, replace = T)
m <- matrix(vals, nrow = rows, ncol = cols)

r1 <- tgs_cor_graph(m, 50, 10, spearman = F)
r2 <- tgs_graph_cover_resample(r1, 10, 1)
}
\keyword{ ~cluster }

