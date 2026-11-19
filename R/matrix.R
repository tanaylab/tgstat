#' For each matrix row apply a function over a ragged array
#'
#' For each matrix row (or column) apply a function to each cell of a ragged
#' array, that is to each (non-empty) group of values given by the levels of a
#' factor.
#'
#' 'tgs_matrix_tapply(x, index, fun)' is essentialy an efficient implementation
#' of 'apply(mat, 1, function(x) tapply(x, index, fun))'.
#' 'tgs_matrix_tapply(x, index, fun, margin = 1)' groups the rows instead of the
#' columns: it is the same as 'tgs_matrix_tapply(t(x), index, fun)' without
#' transposing 'x', and for 'sum' the same as 'rowsum(x, index)'.
#'
#' 'sum', 'mean', 'max', 'min' and 'median' are computed in C++ when 'fun' is
#' passed by name (e.g. 'fun = max', not 'fun = function(v) max(v)') and the
#' only other argument is 'na.rm' (for 'sum' also the extra values to add, for
#' 'mean' also 'trim = 0'). Any other function is called from R once for every
#' group and row (or column), which is much slower.
#'
#' The implicit zeros of a sparse matrix are values like any other, e.g. for
#' 'mean', 'median', 'max' and 'min'.
#'
#' @param x a matrix or a sparse matrix of 'dgCMatrix' type
#' @param index a factor (or a vector, coerced to a factor by 'as.factor') of the
#' same length as the number of columns of 'x' (or rows, if 'margin' is 1)
#' @param fun the function to be applied
#' @param ... optional arguments to 'fun'
#' @param margin 2 (the default) to group the columns of 'x' by 'index', 1 to
#' group the rows. Can only be passed by name.
#' @return A matrix with a row for each level of 'index' and a column for each
#' row of 'x' (or column, if 'margin' is 1). The \code{[i, j]} element is the
#' result of applying 'fun' to \code{x[j, index == levels(index)[i]]} (or to
#' \code{x[index == levels(index)[i], j]} if 'margin' is 1), and NA for a level
#' of 'index' that has no values. \cr
#' Note that the return value is a dense matrix even when \code{x} is sparse.
#' @keywords ~apply ~tapply
#' @examples
#' \donttest{
#' # Note: all the available CPU cores might be used
#'
#' set.seed(seed = 1)
#' nr <- 6
#' nc <- 10
#' mat <- matrix(sample(c(rep(0, 6), 1:3), nr * nc, replace = TRUE), nrow = nr, ncol = nc)
#' index <- factor(rep_len(1:3, ncol(mat)), levels = 0:5)
#' r1 <- apply(mat, 1, function(x) tapply(x, index, sum))
#' r2 <- tgs_matrix_tapply(mat, index, sum)
#' r3 <- tgs_matrix_tapply(t(mat), index, max, margin = 1)
#' }
#'
#' \dontshow{
#' options(tgs_use.blas = FALSE)
#' options(tgs_max.processes = 1)
#'
#' set.seed(seed = 1)
#' nr <- 6
#' nc <- 10
#' mat <- matrix(sample(c(rep(0, 6), 1:3), nr * nc, replace = TRUE), nrow = nr, ncol = nc)
#' index <- factor(rep_len(1:3, ncol(mat)), levels = 0:5)
#' r1 <- apply(mat, 1, function(x) tapply(x, index, sum))
#' r2 <- tgs_matrix_tapply(mat, index, sum)
#' r3 <- tgs_matrix_tapply(t(mat), index, max, margin = 1)
#' }
#'
#' @export tgs_matrix_tapply
tgs_matrix_tapply <- function(x, index, fun, ..., margin = 2) {
    if (missing(x) || missing(index) || missing(fun)) {
        stop("Usage: tgs_matrix_tapply(x, index, fun, ..., margin = 2)", call. = FALSE)
    }

    if (!is.numeric(margin) || length(margin) != 1 || !(margin %in% c(1, 2))) {
        stop("\"margin\" must be 1 or 2", call. = FALSE)
    }

    args <- as.list(substitute(list(...)))[-1L]
    # the same as factor(index), but much faster for an integer index
    index <- as.factor(index)
    .Call("tgs_matrix_tapply", x, index, fun, as.character(substitute(fun)), args, as.integer(margin), new.env(parent = parent.frame()))
}
