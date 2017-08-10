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
#' 'tgs_cor(x, spearman = F)' is equivalent to 'cor(x, method = "pearson")'
#' 'tgs_cor(x, spearman = T)' is equivalent to 'cor(x, method = "spearman")'
#' 'tgs_cor(x, pairwise.complete.obs = T, spearman = T)' is equivalent to
#' 'cor(x, use = "pairwise.complete.obs", method = "spearman")' 'tgs_cor(x,
#' pairwise.complete.obs = T, spearman = F)' is equivalent to 'cor(x, use =
#' "pairwise.complete.obs", method = "pearson")'
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
#' 
#' r1 <- tgs_cor(m, spearman = F)
#' r2 <- tgs_cor(m, pairwise.complete.obs = T, spearman = T)
#' 
#' @export tgs_cor
tgs_cor <- function(x = NULL, pairwise.complete.obs = F, spearman = F, tidy = F, threshold = 0) {
    if (is.null(x))
        stop("Usage: tgs_cor(x, pairwise.complete.obs = F, spearman = F, tidy = F, threshold = 0)", call. = F)
    .tgs_call("tgs_cor", x, pairwise.complete.obs, spearman, tidy, threshold, new.env(parent = parent.frame()))
}
