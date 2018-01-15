\name{tgs_cluster}
\alias{tgs_cluster}
\title{
Clusters directed graph
}
\description{
Clusters directed graph.
}
\usage{
tgs_cluster(graph, min_cluster_size, cooling = 1.05, burn_in = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{graph}{directed graph in the format returned by tgs_cor_graph}
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

\examples{
set.seed(seed = 0)
rows <- 100
cols <- 1000
vals <- sample(1 : (rows * cols / 2), rows * cols, replace = TRUE)
m <- matrix(vals, nrow = rows, ncol = cols)
m[sample(1 : (rows * cols), rows * cols / 1000)] <- NA

r1<-tgs_cor_graph(m, 50, 10, pairwise.complete.obs = F, spearman = T)
r2<-tgs_cluster(r1, 5)
}
\keyword{ ~cluster }

