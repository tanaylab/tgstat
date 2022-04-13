#' Builds directed graph of correlations
#'
#' Builds directed graph of correlations where the nodes are the matrix
#' columns.
#'
#' This function builds a directed graph based on the edges in 'x' and their
#' ranks.
#'
#' 'x' is a data frame containing 4 columns named: 'col1', 'col2', 'val',
#' 'rank'. The third column ('val' can have a different name). The result in
#' the compatible format is returned, for example, by 'tgs_knn' function.
#'
#' 'tgs_graph' prunes some of the edges in 'x' based on the following steps:
#'
#' 1. A pair of columns i, j that appears in 'x' in 'col1', 'col2' implies the
#' edge in the graph from i to j: edge(i,j). Let the rank of i and j be
#' rank(i,j).
#'
#' 2. Calculate symmetrised rank of i and j: sym_rank(i,j) = rank(i,j) *
#' rank(j,i). If one of the ranks is missing from the previous result sym_rank
#' is set to NA.
#'
#' 3. Prune the edges: remove edge(i,j) if sym_rank(i,j) == NA OR sym_rank(i,j)
#' < knn * knn * k_expand
#'
#' 4. Prune excessive incoming edges: remove edge(i,j) if more than knn *
#' k_beta edges of type edge(node,j) exist and sym_rank(i,j) is higher than
#' sym_rank(node,j) for node != j.
#'
#' 5. Prune excessive outgoing edges: remove edge(i,j) if more than knn edges
#' of type edge(i,node) exist and sym_rank(i,j) is higher than sym_rank(i,node)
#' for node != i.
#'
#' @param x see below
#' @param knn maximal node degree
#' @param k_expand see below
#' @param k_beta see below
#' @return The graph edges are returned in a data frame, with the weight of
#' each edge. edge(i,j) receives weight 1 if its sym_rank is the lowest among
#' all edges of type edge(i,node). Formally defined: weight(i,j) = 1 -
#' (place(i,j) - 1) / knn, where place(i,j) is the location of edge(i,j) within
#' the sorted set of edges outgoing from i, i.e. edge(i,node). The sort is done
#' by sym_rank of the edges.
#' @keywords ~graph
#' @examples
#' \donttest{
#' # Note: all the available CPU cores might be used
#'
#' set.seed(seed = 1)
#' rows <- 100
#' cols <- 1000
#' vals <- sample(1:(rows * cols / 2), rows * cols, replace = TRUE)
#' m <- matrix(vals, nrow = rows, ncol = cols)
#' m[sample(1:(rows * cols), rows * cols / 1000)] <- NA
#'
#' r1 <- tgs_cor(m, pairwise.complete.obs = FALSE, spearman = TRUE)
#' r2 <- tgs_knn(r1, knn = 30, diag = FALSE)
#' r3 <- tgs_graph(r2, knn = 3, k_expand = 10)
#' }
#'
#' \dontshow{
#' options(tgs_use.blas = FALSE)
#' options(tgs_max.processes = 1)
#'
#' set.seed(seed = 1)
#' rows <- 100
#' cols <- 100
#' vals <- sample(1:(rows * cols / 2), rows * cols, replace = TRUE)
#' m <- matrix(vals, nrow = rows, ncol = cols)
#' m[sample(1:(rows * cols), rows * cols / 1000)] <- NA
#'
#' r1 <- tgs_cor(m, pairwise.complete.obs = FALSE, spearman = TRUE)
#' r2 <- tgs_knn(r1, knn = 30, diag = FALSE)
#' r3 <- tgs_graph(r2, knn = 3, k_expand = 10)
#' }
#'
#' @export tgs_graph
tgs_graph <- function(x, knn, k_expand, k_beta = 3) {
    if (missing(x) || missing(knn) || missing(k_expand)) {
        stop("Usage: tgs_graph(x, knn, k_expand, k_beta = 3)", call. = FALSE)
    }

    .Call("tgs_cor_graph", x, knn, k_expand, k_beta, new.env(parent = parent.frame()))
}



#' Clusters directed graph
#'
#' Clusters directed graph.
#'
#' The algorithm is explained in a "MetaCell: analysis of single-cell RNA-seq
#' data using K-nn graph partitions" paper, published in "Genome Biology" #20:
#' https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1812-2
#'
#' @param graph directed graph in the format returned by tgs_graph
#' @param min_cluster_size used to determine the candidates for seeding (= min
#' weight)
#' @param cooling factor that is used to gradually increase the chance of a
#' node to stay in the cluster
#' @param burn_in number of node reassignments after which cooling is applied
#' @return Data frame that maps each node to its cluster.
#' @seealso [tgs_graph()]
#' @keywords ~cluster
#' @examples
#' \donttest{
#' # Note: all the available CPU cores might be used
#'
#' set.seed(seed = 0)
#' rows <- 100
#' cols <- 1000
#' vals <- sample(1:(rows * cols / 2), rows * cols, replace = TRUE)
#' m <- matrix(vals, nrow = rows, ncol = cols)
#' m[sample(1:(rows * cols), rows * cols / 1000)] <- NA
#'
#' r1 <- tgs_cor(m, pairwise.complete.obs = FALSE, spearman = TRUE)
#' r2 <- tgs_knn(r1, knn = 30, diag = FALSE)
#' r3 <- tgs_graph(r2, knn = 3, k_expand = 10)
#' r4 <- tgs_graph_cover(r3, 5)
#' }
#'
#' \dontshow{
#' options(tgs_use.blas = FALSE)
#' options(tgs_max.processes = 1)
#'
#' set.seed(seed = 0)
#' rows <- 100
#' cols <- 100
#' vals <- sample(1:(rows * cols / 2), rows * cols, replace = TRUE)
#' m <- matrix(vals, nrow = rows, ncol = cols)
#' m[sample(1:(rows * cols), rows * cols / 1000)] <- NA
#'
#' r1 <- tgs_cor(m, pairwise.complete.obs = FALSE, spearman = TRUE)
#' r2 <- tgs_knn(r1, knn = 30, diag = FALSE)
#' r3 <- tgs_graph(r2, knn = 3, k_expand = 10)
#' r4 <- tgs_graph_cover(r3, 5)
#' }
#'
#' @export tgs_graph_cover
tgs_graph_cover <- function(graph, min_cluster_size, cooling = 1.05, burn_in = 10) {
    if (missing(graph) || missing(min_cluster_size)) {
        stop("Usage: tgs_graph_cover(graph, min_cluster_size, cooling = 1.05, burn_in = 10)", call. = FALSE)
    }

    .Call("tgs_graph2cluster", graph, min_cluster_size, cooling, burn_in, new.env(parent = parent.frame()))
}



#' Clusters directed graph multiple times with randomized sample subset
#'
#' Clusters directed graph multiple times with randomized sample subset.
#'
#' The algorithm is explained in a "MetaCell: analysis of single-cell RNA-seq
#' data using K-nn graph partitions" paper, published in "Genome Biology" #20:
#' https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1812-2
#'
#' @param graph directed graph in the format returned by tgs_graph
#' @param knn maximal number of edges used per node for each sample subset
#' @param min_cluster_size used to determine the candidates for seeding (= min
#' weight)
#' @param cooling factor that is used to gradually increase the chance of a
#' node to stay in the cluster
#' @param burn_in number of node reassignments after which cooling is applied
#' @param p_resamp fraction of total number of nodes used in each sample subset
#' @param n_resamp number iterations the clustering is run on different sample
#' subsets
#' @param method method for calculating co_cluster and co_sample; valid values:
#' "hash", "full", "edges"
#' @return If method == "hash", a list with two members. The first member is a
#' data frame with 3 columns: "node1", "node2" and "cnt". "cnt" indicates the
#' number of times "node1" and "node2" appeared in the same cluster. The second
#' member of the list is a vector of __number of nodes__ length reflecting how
#' many times each node was used in the subset.
#'
#' If method == "full", a list containing two matrices: co_cluster and
#' co_sample.
#'
#' If method == "edges", a list containing two data frames: co_cluster and
#' co_sample.
#' @seealso [tgs_graph()]
#' @keywords ~cluster
#' @examples
#' \donttest{
#' # Note: all the available CPU cores might be used
#'
#' set.seed(seed = 0)
#' rows <- 100
#' cols <- 200
#' vals <- sample(1:(rows * cols / 2), rows * cols, replace = TRUE)
#' m <- matrix(vals, nrow = rows, ncol = cols)
#'
#' r1 <- tgs_cor(m, pairwise.complete.obs = FALSE, spearman = TRUE)
#' r2 <- tgs_knn(r1, knn = 20, diag = FALSE)
#' r3 <- tgs_graph(r2, knn = 3, k_expand = 10)
#' r4 <- tgs_graph_cover_resample(r3, 10, 1)
#' }
#'
#' \dontshow{
#' set.seed(seed = 0)
#' rows <- 100
#' cols <- 200
#' vals <- sample(1:(rows * cols / 2), rows * cols, replace = TRUE)
#' m <- matrix(vals, nrow = rows, ncol = cols)
#'
#' r1 <- tgs_cor(m, pairwise.complete.obs = FALSE, spearman = TRUE)
#' r2 <- tgs_knn(r1, knn = 20, diag = FALSE)
#' r3 <- tgs_graph(r2, knn = 3, k_expand = 10)
#' r4 <- tgs_graph_cover_resample(r3, 10, 1, n_resamp = 5)
#' }
#'
#' @export tgs_graph_cover_resample
tgs_graph_cover_resample <- function(graph, knn, min_cluster_size, cooling = 1.05, burn_in = 10, p_resamp = 0.75, n_resamp = 500, method = "hash") {
    if (missing(graph) || missing(knn) || missing(min_cluster_size)) {
        stop("Usage: tgs_graph_cover_resample(graph, knn, min_cluster_size, cooling = 1.05, burn_in = 10, p_resamp = 0.75, n_resamp = 500)", call. = FALSE)
    }

    if (method == "hash") {
        .Call("tgs_graph2cluster_multi_hash", graph, knn, min_cluster_size, cooling, burn_in, p_resamp, n_resamp, new.env(parent = parent.frame()))
    } else if (method == "full") {
        .Call("tgs_graph2cluster_multi_full", graph, knn, min_cluster_size, cooling, burn_in, p_resamp, n_resamp, new.env(parent = parent.frame()))
    } else if (method == "edges") {
        .Call("tgs_graph2cluster_multi_edges", graph, knn, min_cluster_size, cooling, burn_in, p_resamp, n_resamp, new.env(parent = parent.frame()))
    } else {
        stop("\"method\" argument must be equal to \"hash\", \"full\" or \"edges\"", call. = FALSE)
    }
}



#' Returns k highest values of each column
#'
#' Returns k highest values of each column.
#'
#' 'tgs_knn' returns the highest 'knn' values of each column of 'x' (if 'x' is
#' a matrix). 'x' can be also a sparse matrix given in a data frame of 'col',
#' 'row', 'value' format.
#'
#' 'NA' and 'Inf' values are skipped as well as the values below 'threshold'.
#' If 'diag' is 'F' values of the diagonal (row == col) are skipped too.
#'
#' @param x numeric matrix or data frame (see below)
#' @param knn the number of highest values returned per column
#' @param diag if 'F' values of row 'i' and col 'j' are skipped for each i == j
#' @param threshold filter out values lower than threshold
#' @return A sparse matrix in a data frame format with 'col1', 'col2', 'val'
#' and 'rank' columns. 'col1' and 'col2' represent the column and the row
#' number of 'x'.
#' @keywords ~knn
#' @examples
#' \donttest{
#' # Note: all the available CPU cores might be used
#'
#' set.seed(seed = 1)
#' rows <- 100
#' cols <- 1000
#' vals <- sample(1:(rows * cols / 2), rows * cols, replace = TRUE)
#' m <- matrix(vals, nrow = rows, ncol = cols)
#' m[sample(1:(rows * cols), rows * cols / 1000)] <- NA
#' r <- tgs_knn(m, 3)
#' }
#'
#' \dontshow{
#' options(tgs_use.blas = FALSE)
#' options(tgs_max.processes = 1)
#'
#' set.seed(seed = 1)
#' rows <- 100
#' cols <- 100
#' vals <- sample(1:(rows * cols / 2), rows * cols, replace = TRUE)
#' m <- matrix(vals, nrow = rows, ncol = cols)
#' m[sample(1:(rows * cols), rows * cols / 1000)] <- NA
#' r <- tgs_knn(m, 3)
#' }
#'
#' @export tgs_knn
tgs_knn <- function(x, knn, diag = FALSE, threshold = 0) {
    if (missing(x) || missing(knn)) {
        stop("Usage: tgs_knn(x, knn, diag = FALSE, threshold = 0)", call. = FALSE)
    }

    if (is.integer(x)) {
        tmp <- as.double(x)
        attributes(tmp) <- attributes(x)
        x <- tmp
    } else if (is.data.frame(x) && ncol(x) >= 3 && is.integer(x[, 3])) {
        tmp <- as.double(x[, 3])
        attributes(tmp) <- attributes(x[, 3])
        x[, 3] <- tmp
    }

    .Call("tgs_knn", x, knn, diag, threshold, new.env(parent = parent.frame()))
}
