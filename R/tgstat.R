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
	tryCatch({ res <- .Call(..., PACKAGE='tgstat') },
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

#' Calculates correlation between the matrix columns
#' 
#' Calculates correlation between the matrix columns.
#' 
#' This function is very similar to 'package:stats::cor'. Unlike the latter it
#' uses all available CPU cores to compute the correlation in a much faster
#' way. The basic implementation of 'pairwise.complete.obs' is also more
#' efficient giving overall great run-time advantage.
#' 
#' Unlike 'package:stats::cor' 'tgs_cor' computes the correlation only between
#' matrix columns and implements only two modes of treating data containing NA,
#' which are equivalent to 'use="everything"' and 'use="pairwise.complete.obs".
#' Please refer the documentation of this function for more details.
#' 
#' \itemize{
#'  \item{}{\code{tgs_cor(x, spearman = F)} is equivalent to \code{cor(x, method = "pearson")}}
#'  \item{}{\code{tgs_cor(x, pairwise.complete.obs = T, spearman = T)} is equivalent to \code{cor(x, use = "pairwise.complete.obs", method = "spearman")}}
#'  \item{}{\code{tgs_cor(x, pairwise.complete.obs = T, spearman = F)} is equivalent to \code{cor(x, use = "pairwise.complete.obs", method = "pearson")}}
#' }
#' 
#' Finally 'tgs_cor' can output its result in "tidy" format: a data frame with
#' three columns named 'col1', 'col2' and 'cor'. Only the correlation values
#' which abs are equal or above the 'threshold' are reported. Correlation of
#' column number X with column Y is reported only if X < Y.
#' 
#' @param x numeric matrix
#' @param pairwise.complete.obs see below
#' @param spearman if 'TRUE' Spearman correlation is computed, otherwise
#' Pearson
#' @param tidy if 'TRUE' data is outputed in tidy format
#' @param threshold absolute threshold above which values are outputed in tidy
#' format
#' @return If 'tidy' is 'FALSE' - the matrix of correlation values, where
#' val[X,Y] represents the correlation between columns X and Y of the input
#' matrix. If 'tidy' is 'TRUE' - a data frame, where each row represents
#' correlation between two pairs of columns.
#' @keywords ~correlation
#' @examples
#' 
#' set.seed(seed = 0)
#' rows <- 100
#' cols <- 1000
#' vals <- sample(1 : (rows * cols / 2), rows * cols, replace = T)
#' m <- matrix(vals, nrow = rows, ncol = cols)
#' m[sample(1 : (rows * cols), rows * cols / 1000)] <- NA
#' m[1:5, 1:5]
#' 
#' r1 <- tgs_cor(m, spearman = F)
#' r2 <- tgs_cor(m, pairwise.complete.obs = T, spearman = T)#' 
#' r2[1:5, 1:5]
#' 
#' @export 
tgs_cor <- function(x, pairwise.complete.obs = F, spearman = F, tidy = F, threshold = 0) {
    if (missing(x))
        stop("Usage: tgs_cor(x, pairwise.complete.obs = F, spearman = F, tidy = F, threshold = 0)", call. = F)
    .tgs_call("tgs_cor", x, pairwise.complete.obs, spearman, tidy, threshold, new.env(parent = parent.frame()))
}

#' Calculates distances between the matrix rows
#' 
#' Calculates distances between the matrix rows.
#' 
#' This function is very similar to 'package:stats::dist'. Unlike the latter it
#' uses all available CPU cores to compute the distances in a much faster way.
#' 
#' Unlike 'package:stats::dist' 'tgs_dist' uses always "euclidean" metrics (see
#' 'method' parameter of 'dist' function). Thus:
#' 
#' 'tgs_dist(x)' is equivalent to 'dist(x, method = "euclidean")'
#' 
#' 'tgs_dist' can output its result in "tidy" format: a data frame with three
#' columns named 'row1', 'row2' and 'dist'. Only the distances that are less or
#' equal than the 'threshold' are reported. Distance between row number X and Y
#' is reported only if X < Y.
#' 
#' @param x numeric matrix
#' @param diag see 'dist' documentation
#' @param upper see 'dist' documentation
#' @param tidy if 'TRUE' data is outputed in tidy format
#' @param threshold threshold below which values are outputed in tidy format
#' @return If 'tidy' is 'FALSE' - the output is similar to that of 'dist'
#' function. If 'tidy' is 'TRUE' - 'tgs_dist' returns a data frame, where each
#' row represents distances between two pairs of original rows.
#' @keywords ~distance
#' @examples
#' 
#' set.seed(seed = 0)
#' rows <- 100
#' cols <- 1000
#' vals <- sample(1 : (rows * cols / 2), rows * cols, replace = TRUE)
#' m <- matrix(vals, nrow = rows, ncol = cols)
#' m[sample(1 : (rows * cols), rows * cols / 1000)] <- NA
#' r <- tgs_dist(m)
#' 
#' @export 
tgs_dist <- function(x, diag = FALSE, upper = FALSE, tidy = F, threshold = Inf) {
    if (missing(x))
        stop("Usage: tgs_dist(x, diag = F, upper = F, tidy = F, threshold = Inf)", call. = F)

    attrs <- list(Size = nrow(x), Labels = dimnames(x)[[1L]], Diag = diag, 
        Upper = upper, method = "euclidian", call = match.call(), class = "dist")

    .tgs_call("tgs_dist", x, attrs, tidy, threshold, dimnames(x)[[1L]], new.env(parent = parent.frame()))
}
