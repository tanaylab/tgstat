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

#' Calculates distances between the correlation profiles of matrix columns
#'
#' Calculates the euclidean distances between the rows of the correlation
#' matrix of the columns of 'x'.
#'
#' 'tgs_cor_dist(x, pairwise.complete.obs, spearman)' returns the same
#' distances as 'tgs_dist(tgs_cor(x, pairwise.complete.obs =
#' pairwise.complete.obs, spearman = spearman))' up to rounding errors, for
#' well-conditioned data. (When the columns of 'x' have offsets that are large
#' relative to their spread, 'tgs_cor' loses precision while 'tgs_cor_dist',
#' which centers the columns first, does not, so the two differ.)
#'
#' When 'x' has much fewer rows (G) than columns (N), 'tgs_cor_dist' does not
#' compute the NxN correlation matrix. Instead it maps each column of 'x' to a
#' point in G dimensions, such that the euclidean distances between the points
#' equal the distances between the rows of the correlation matrix, and passes
#' these points to 'tgs_dist'. This takes time proportional to G^3 + N^2 x G
#' (an eigendecomposition of a GxG matrix and the distances) instead of
#' N^2 x G + N^3, and is used when G < N / 3 with BLAS ('tgs_use.blas') and
#' G < N / 2 without, which is where it was measured to be faster.
#'
#' It is also used only when all the values of 'x' are finite and no column of
#' 'x' is constant. Otherwise 'tgs_cor_dist' returns 'tgs_dist(tgs_cor(...))'.
#'
#' @param x numeric matrix
#' @param pairwise.complete.obs see 'tgs_cor' documentation
#' @param spearman if 'TRUE' Spearman correlation is computed, otherwise
#' Pearson
#' @return An object of class 'dist' with the distances between the columns of
#' 'x', similar to the output of 'tgs_dist'.
#' @seealso [tgs_cor()], [tgs_dist()]
#' @keywords ~distance ~correlation
#' @examples
#' \donttest{
#' # Note: all the available CPU cores might be used
#'
#' set.seed(seed = 0)
#' m <- matrix(rnorm(100 * 2000), nrow = 100)
#' d <- tgs_cor_dist(m)
#' hc <- hclust(d)
#' }
#'
#' \dontshow{
#' options(tgs_use.blas = FALSE)
#' options(tgs_max.processes = 1)
#'
#' set.seed(seed = 0)
#' m <- matrix(rnorm(20 * 200), nrow = 20)
#' d <- tgs_cor_dist(m)
#' }
#'
#' @export tgs_cor_dist
tgs_cor_dist <- function(x, pairwise.complete.obs = FALSE, spearman = FALSE) {
    if (missing(x)) {
        stop("Usage: tgs_cor_dist(x, pairwise.complete.obs = FALSE, spearman = FALSE)", call. = FALSE)
    }

    if (!is.logical(pairwise.complete.obs) || length(pairwise.complete.obs) != 1) {
        stop("\"pairwise.complete.obs\" argument must be a logical value", call. = FALSE)
    }

    # z: the columns of x (or their ranks) centered and scaled to unit norm, so that the correlation
    # matrix is crossprod(z). NULL if x has non-finite values or constant columns.
    # Using z costs G^3 + N^2 G (G = nrow(x), N = ncol(x)) and the composition N^2 G + N^3. Measured on
    # N = 2000 to 8000, z is faster for G below about N / 3 with BLAS and 0.6 N without.
    z <- NULL
    if (is.matrix(x) && nrow(x) < ncol(x) / (if (.tgs_use_blas()) 3 else 2)) {
        z <- .Call("tgs_standardize", x, spearman, new.env(parent = parent.frame()))
    }

    if (is.null(z)) {
        d <- tgs_dist(tgs_cor(x, pairwise.complete.obs = pairwise.complete.obs, spearman = spearman))
    } else {
        # The squared distance between rows i and j of crossprod(z) is
        # t(z[, i] - z[, j]) %*% tcrossprod(z) %*% (z[, i] - z[, j]). With tcrossprod(z) = V L t(V),
        # this is the squared distance between rows i and j of crossprod(z, V L^(1/2)).
        e <- eigen(tcrossprod(z), symmetric = TRUE)
        d <- tgs_dist(crossprod(z, e$vectors * rep(sqrt(pmax(e$values, 0)), each = nrow(z))))
    }
    attr(d, "call") <- match.call()
    d
}
