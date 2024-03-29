% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/knn.R
\name{tgs_graph_cover_resample}
\alias{tgs_graph_cover_resample}
\title{Clusters directed graph multiple times with randomized sample subset}
\usage{
tgs_graph_cover_resample(
  graph,
  knn,
  min_cluster_size,
  cooling = 1.05,
  burn_in = 10,
  p_resamp = 0.75,
  n_resamp = 500,
  method = "hash"
)
}
\arguments{
\item{graph}{directed graph in the format returned by tgs_graph}

\item{knn}{maximal number of edges used per node for each sample subset}

\item{min_cluster_size}{used to determine the candidates for seeding (= min
weight)}

\item{cooling}{factor that is used to gradually increase the chance of a
node to stay in the cluster}

\item{burn_in}{number of node reassignments after which cooling is applied}

\item{p_resamp}{fraction of total number of nodes used in each sample subset}

\item{n_resamp}{number iterations the clustering is run on different sample
subsets}

\item{method}{method for calculating co_cluster and co_sample; valid values:
"hash", "full", "edges"}
}
\value{
If method == "hash", a list with two members. The first member is a
data frame with 3 columns: "node1", "node2" and "cnt". "cnt" indicates the
number of times "node1" and "node2" appeared in the same cluster. The second
member of the list is a vector of \strong{number of nodes} length reflecting how
many times each node was used in the subset.

If method == "full", a list containing two matrices: co_cluster and
co_sample.

If method == "edges", a list containing two data frames: co_cluster and
co_sample.
}
\description{
Clusters directed graph multiple times with randomized sample subset.
}
\details{
The algorithm is explained in a "MetaCell: analysis of single-cell RNA-seq
data using K-nn graph partitions" paper, published in "Genome Biology" #20:
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1812-2
}
\examples{
\donttest{
# Note: all the available CPU cores might be used

set.seed(seed = 0)
rows <- 100
cols <- 200
vals <- sample(1:(rows * cols / 2), rows * cols, replace = TRUE)
m <- matrix(vals, nrow = rows, ncol = cols)

r1 <- tgs_cor(m, pairwise.complete.obs = FALSE, spearman = TRUE)
r2 <- tgs_knn(r1, knn = 20, diag = FALSE)
r3 <- tgs_graph(r2, knn = 3, k_expand = 10)
r4 <- tgs_graph_cover_resample(r3, 10, 1)
}

\dontshow{
set.seed(seed = 0)
rows <- 100
cols <- 200
vals <- sample(1:(rows * cols / 2), rows * cols, replace = TRUE)
m <- matrix(vals, nrow = rows, ncol = cols)

r1 <- tgs_cor(m, pairwise.complete.obs = FALSE, spearman = TRUE)
r2 <- tgs_knn(r1, knn = 20, diag = FALSE)
r3 <- tgs_graph(r2, knn = 3, k_expand = 10)
r4 <- tgs_graph_cover_resample(r3, 10, 1, n_resamp = 5)
}

}
\seealso{
\code{\link[=tgs_graph]{tgs_graph()}}
}
\keyword{~cluster}
