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
#' is reported only if X < Y. 'diag' and 'upper' parameters are ignored when
#' the result is returned in "tidy" format.
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
#' \donttest{
#' # Note: all the available CPU cores might be used
#'
#' set.seed(seed = 0)
#' rows <- 100
#' cols <- 1000
#' vals <- sample(1:(rows * cols / 2), rows * cols, replace = TRUE)
#' m <- matrix(vals, nrow = rows, ncol = cols)
#' m[sample(1:(rows * cols), rows * cols / 1000)] <- NA
#' r <- tgs_dist(m)
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
#' r <- tgs_dist(m)
#' }
#'
#' @export tgs_dist
tgs_dist <- function(x, diag = FALSE, upper = FALSE, tidy = FALSE, threshold = Inf) {
    if (missing(x)) {
        stop("Usage: tgs_dist(x, diag = FALSE, upper = FALSE, tidy = FALSE, threshold = Inf)", call. = FALSE)
    }

    attrs <- list(
        Size = nrow(x), Labels = dimnames(x)[[1L]], Diag = diag,
        Upper = upper, method = "euclidean", call = match.call(), class = "dist"
    )

    if (.tgs_use_blas()) {
        .Call("tgs_dist_blas", x, attrs, tidy, threshold, dimnames(x)[[1L]], new.env(parent = parent.frame()))
    } else {
        .Call("tgs_dist", x, attrs, tidy, threshold, dimnames(x)[[1L]], new.env(parent = parent.frame()))
    }
}
