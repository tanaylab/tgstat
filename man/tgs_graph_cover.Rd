\name{tgs_graph_cover}
\alias{tgs_graph_cover}
\title{
Clusters directed graph
}
\description{
Clusters directed graph.
}
\usage{
tgs_graph_cover(graph, min_cluster_size, cooling = 1.05,
                burn_in = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{graph}{directed graph in the format returned by tgs_graph}
  \item{min_cluster_size}{used to determine the candidates for seeding (= min weight)}
  \item{cooling}{factor that is used to gradually increase the chance of a node to stay in the cluster}
  \item{burn_in}{number of node reassignments after which cooling is applied}
}
\details{
The algorithm is explained in a separate document.
}
\value{
Data frame that maps each node to its cluster.
}

\seealso{
\code{\link{tgs_graph}}
}

\examples{
set.seed(seed = 0)
rows <- 100
cols <- 1000
vals <- sample(1 : (rows * cols / 2), rows * cols, replace = TRUE)
m <- matrix(vals, nrow = rows, ncol = cols)
m[sample(1 : (rows * cols), rows * cols / 1000)] <- NA

r1 <- tgs_cor(m, pairwise.complete.obs = FALSE, spearman = TRUE)
r2 <- tgs_knn(r1, knn = 30, diag = FALSE)
r3 <- tgs_graph(r2, knn = 3, k_expand = 10)
r4 <- tgs_graph_cover(r3, 5)
}
\keyword{ ~cluster }

