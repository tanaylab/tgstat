.onLoad <- function(lib, pkg) {
}

.onAttach <- function(lib, pkg) {
	Sys.umask("0002")

	assign(".TGS_FUNCS", getNamespaceExports("tgstat"), envir = .GlobalEnv)

	if (R.Version()$major >= 3)
	    assign(".TGS_LIBDIR", path.package("tgstat"), envir = .GlobalEnv)
	else
	    assign(".TGS_LIBDIR", .path.package("tgstat"), envir = .GlobalEnv)	
}

.onDetach <- function(lib) {
	if (exists(".TGS_FUNCS", envir = .GlobalEnv))
		remove(".TGS_FUNCS", envir = .GlobalEnv)
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

tgs_cor <- function(x, pairwise.complete.obs = F, spearman = F, tidy = F, threshold = 0) {
    if (missing(x))
        stop("Usage: tgs_cor(x, pairwise.complete.obs = F, spearman = F, tidy = F, threshold = 0)", call. = F)

    if (pairwise.complete.obs && spearman && !tgs_finite(x) || !.tgs_use_blas())
        .tgs_call("tgs_cor", x, pairwise.complete.obs, spearman, tidy, threshold, NULL, new.env(parent = parent.frame()))
    else
        .tgs_call("tgs_cor_blas", x, pairwise.complete.obs, spearman, tidy, threshold, NULL, new.env(parent = parent.frame()))
}

tgs_cor_graph <- function(x, knn, k_expand, k_alpha = 10, k_beta = 3, pairwise.complete.obs = F, spearman = F) {
    if (missing(x) || missing(knn) || missing(k_expand))
        stop("Usage: tgs_cor_graph(x, knn, k_expand, k_alpha = 10, k_beta = 3, pairwise.complete.obs = F, spearman = F)", call. = F)

    r <- tgs_cor_knn(x, knn * k_alpha, pairwise.complete.obs = pairwise.complete.obs, spearman = spearman, threshold = 0)
    .tgs_call("tgs_cor_graph", r, knn, k_expand, k_beta, colnames(x), new.env(parent = parent.frame()))
}

tgs_cor_knn <- function(x, knn, pairwise.complete.obs = F, spearman = F, threshold = 0) {
    if (missing(x) || missing(knn))
        stop("Usage: tgs_cor_knn(x, knn, pairwise.complete.obs = F, spearman = F, threshold = 0)", call. = F)

    if (pairwise.complete.obs && spearman && !tgs_finite(x) || !.tgs_use_blas())
        .tgs_call("tgs_cor", x, pairwise.complete.obs, spearman, T, threshold, knn, new.env(parent = parent.frame()))
    else
        .tgs_call("tgs_cor_blas", x, pairwise.complete.obs, spearman, T, threshold, knn, new.env(parent = parent.frame()))
}

tgs_dist <- function(x, diag = FALSE, upper = FALSE, tidy = F, threshold = Inf) {
    if (missing(x))
        stop("Usage: tgs_dist(x, diag = F, upper = F, tidy = F, threshold = Inf)", call. = F)

    attrs <- list(Size = nrow(x), Labels = dimnames(x)[[1L]], Diag = diag, 
        Upper = upper, method = "euclidian", call = match.call(), class = "dist")

    if (.tgs_use_blas())
        .tgs_call("tgs_dist_blas", x, attrs, tidy, threshold, dimnames(x)[[1L]], new.env(parent = parent.frame()))
    else
        .tgs_call("tgs_dist", x, attrs, tidy, threshold, dimnames(x)[[1L]], new.env(parent = parent.frame()))
}

tgs_graph_cover <- function(graph, min_cluster_size, cooling = 1.05, burn_in = 10) {
    if (missing(graph) || missing(min_cluster_size))
        stop("Usage: tgs_graph_cover(graph, min_cluster_size, cooling = 1.05, burn_in = 10)", call. = F)

    .tgs_call("tgs_graph2cluster", graph, min_cluster_size, cooling, burn_in, new.env(parent = parent.frame()))
}

tgs_graph_cover_resample <- function(graph, knn, min_cluster_size, cooling = 1.05, burn_in = 10, p_resamp = 0.75, n_resamp = 500, method = "hash") {
    if (missing(graph) || missing(knn) || missing(min_cluster_size))
        stop("Usage: tgs_graph_cover_resample(graph, knn, min_cluster_size, cooling = 1.05, burn_in = 10, p_resamp = 0.75, n_resamp = 500)", call. = F)

    if (method == "hash")
        .tgs_call("tgs_graph2cluster_multi_hash", graph, knn, min_cluster_size, cooling, burn_in, p_resamp, n_resamp, method, new.env(parent = parent.frame()))
    else if (method == "full")
        .tgs_call("tgs_graph2cluster_multi_full", graph, knn, min_cluster_size, cooling, burn_in, p_resamp, n_resamp, method, new.env(parent = parent.frame()))
    else if (method == "edges")
        .tgs_call("tgs_graph2cluster_multi_edges", graph, knn, min_cluster_size, cooling, burn_in, p_resamp, n_resamp, method, new.env(parent = parent.frame()))
    else
        stop("\"method\" argument must be equal to \"hash\", \"full\" or \"edges\"", call. = F)
}

tgs_finite <- function(x) {
    if (missing(x))
        stop("Usage: tgs_finite(x)", call. = F)

    .tgs_call("tgs_finite", x, new.env(parent = parent.frame()))
}

