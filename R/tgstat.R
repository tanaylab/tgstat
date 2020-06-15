.onLoad <- function(lib, pkg) {
}

.tgs_call <- function(...) {
	tryCatch({ res <- .Call(...) },
			 interrupt = function(interrupt){ stop("Command interrupted!", call. = FALSE); } )
	res
}

.tgs_getOption <- function(x, default = NULL) {
    if (missing(default)) 
        return(options(x)[[1L]])
    if (x %in% names(options())) 
        options(x)[[1L]]
    else default
}

.tgs_use_blas <- function() {
    .tgs_getOption("tgs_use.blas", F)
}

tgs_cor <- function(x, y = NULL, pairwise.complete.obs = F, spearman = F, tidy = F, threshold = 0) {
    if (missing(x))
        stop("Usage: tgs_cor(x, y = NULL, pairwise.complete.obs = F, spearman = F, tidy = F, threshold = 0)", call. = F)

    if (is.null(y)) {
        if (!.tgs_use_blas() || pairwise.complete.obs && spearman && !tgs_finite(x))
            .Call("tgs_cor", x, pairwise.complete.obs, spearman, tidy, threshold, new.env(parent = parent.frame()))
        else
            .Call("tgs_cor_blas", x, pairwise.complete.obs, spearman, tidy, threshold, new.env(parent = parent.frame()))
    } else {
        if (!.tgs_use_blas() || pairwise.complete.obs && spearman && (!tgs_finite(x) || !tgs_finite(y)))
            .Call("tgs_cross_cor", x, y, pairwise.complete.obs, spearman, tidy, threshold, new.env(parent = parent.frame()))
        else
            .Call("tgs_cross_cor_blas", x, y, pairwise.complete.obs, spearman, tidy, threshold, new.env(parent = parent.frame()))
    }
}

tgs_cor_knn <- function(x, y, knn, pairwise.complete.obs = F, spearman = F, threshold = 0) {
    if (missing(x) || missing(knn))
        stop("Usage: tgs_cor_knn(x, y, knn, pairwise.complete.obs = F, spearman = F, threshold = 0)", call. = F)

    if (is.null(y))
        .Call("tgs_cor_knn", x, knn, pairwise.complete.obs, spearman, threshold, new.env(parent = parent.frame()))
    else
        .Call("tgs_cross_cor_knn", x, y, knn, pairwise.complete.obs, spearman, threshold, new.env(parent = parent.frame()))
}

tgs_dist <- function(x, diag = F, upper = F, tidy = F, threshold = Inf) {
    if (missing(x))
        stop("Usage: tgs_dist(x, diag = F, upper = F, tidy = F, threshold = Inf)", call. = F)

    attrs <- list(Size = nrow(x), Labels = dimnames(x)[[1L]], Diag = diag, 
        Upper = upper, method = "euclidian", call = match.call(), class = "dist")

    if (.tgs_use_blas())
        .Call("tgs_dist_blas", x, attrs, tidy, threshold, dimnames(x)[[1L]], new.env(parent = parent.frame()))
    else
        .Call("tgs_dist", x, attrs, tidy, threshold, dimnames(x)[[1L]], new.env(parent = parent.frame()))
}

tgs_graph <- function(x, knn, k_expand, k_beta = 3) {
    if (missing(x) || missing(knn) || missing(k_expand))
        stop("Usage: tgs_graph(x, knn, k_expand, k_beta = 3)", call. = F)

    .Call("tgs_cor_graph", x, knn, k_expand, k_beta, new.env(parent = parent.frame()))
}

tgs_graph_cover <- function(graph, min_cluster_size, cooling = 1.05, burn_in = 10) {
    if (missing(graph) || missing(min_cluster_size))
        stop("Usage: tgs_graph_cover(graph, min_cluster_size, cooling = 1.05, burn_in = 10)", call. = F)

    .Call("tgs_graph2cluster", graph, min_cluster_size, cooling, burn_in, new.env(parent = parent.frame()))
}

tgs_graph_cover_resample <- function(graph, knn, min_cluster_size, cooling = 1.05, burn_in = 10, p_resamp = 0.75, n_resamp = 500, method = "hash") {
    if (missing(graph) || missing(knn) || missing(min_cluster_size))
        stop("Usage: tgs_graph_cover_resample(graph, knn, min_cluster_size, cooling = 1.05, burn_in = 10, p_resamp = 0.75, n_resamp = 500)", call. = F)

    if (method == "hash")
        .Call("tgs_graph2cluster_multi_hash", graph, knn, min_cluster_size, cooling, burn_in, p_resamp, n_resamp, new.env(parent = parent.frame()))
    else if (method == "full")
        .Call("tgs_graph2cluster_multi_full", graph, knn, min_cluster_size, cooling, burn_in, p_resamp, n_resamp, new.env(parent = parent.frame()))
    else if (method == "edges")
        .Call("tgs_graph2cluster_multi_edges", graph, knn, min_cluster_size, cooling, burn_in, p_resamp, n_resamp, new.env(parent = parent.frame()))
    else
        stop("\"method\" argument must be equal to \"hash\", \"full\" or \"edges\"", call. = F)
}

tgs_knn <- function(x, knn, diag = F, threshold = 0) {
    if (missing(x) || missing(knn))
        stop("Usage: tgs_knn(x, knn, diag = F, threshold = 0)", call. = F)

    if (is.integer(x)) {
        tmp <- as.double(x)
        attributes(tmp) <- attributes(x)
        x <- tmp
    } else if (is.data.frame(x) && ncol(x) >= 3 && is.integer(x[,3])) {
        tmp <- as.double(x[,3])
        attributes(tmp) <- attributes(x[,3])
        x[,3] <- tmp
    }

    .Call("tgs_knn", x, knn, diag, threshold, new.env(parent = parent.frame()))
}

tgs_finite <- function(x) {
    if (missing(x))
        stop("Usage: tgs_finite(x)", call. = F)

    .Call("tgs_finite", x, new.env(parent = parent.frame()))
}

tgs_matrix_tapply <- function(x, index, fun, ...) {
    if (missing(x) || missing(index) || missing(fun))
        stop("Usage: tgs_matrix_tapply(x, index, fun)", call. = F)

    args <- as.list(substitute(list(...)))[-1L]
    if (!is.factor(index))
        index <- factor(index)
    .Call("tgs_matrix_tapply", x, index, fun, as.character(substitute(fun)), args, new.env(parent = parent.frame()))
}
