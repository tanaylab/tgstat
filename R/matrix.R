#' For each matrix row apply a function over a ragged array
#'
#' For each matrix row apply a function to each cell of a ragged array, that is
#' to each (non-empty) group of values given by a unique combination of the
#' levels of certain factors.
#'
#' 'tgs_matrix_tapply(x, index, fun)' is essentialy an efficient implementation
#' of 'apply(mat, 1, function(x) tapply(x, index, fun))'.
#'
#' @param x a matrix or a sparse matrix of 'dgCMatrix' type
#' @param index a 'list' of one or more 'factor's, each of same length as the
#' number of columns in 'x'. The elements are coerced to factors by
#' 'as.factor'.
#' @param fun the function to be applied
#' @param ... optional arguments to 'fun'
#' @return A matrix of length(index) X nrow(x) size. Each \code{[i,j]} element
#' represents the result of applying 'fun' to
#' \code{x[i,which(index==levels(index)[j])]}. \cr
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
#' }
#'
#' @export tgs_matrix_tapply
tgs_matrix_tapply <- function(x, index, fun, ...) {
    if (missing(x) || missing(index) || missing(fun)) {
        stop("Usage: tgs_matrix_tapply(x, index, fun)", call. = FALSE)
    }

    args <- as.list(substitute(list(...)))[-1L]
    if (!is.factor(index)) {
        index <- factor(index)
    }
    .Call("tgs_matrix_tapply", x, index, fun, as.character(substitute(fun)), args, new.env(parent = parent.frame()))
}
