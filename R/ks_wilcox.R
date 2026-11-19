#' Vectorized two-sample Kolmogorov-Smirnov test
#'
#' Performs a two-sample Kolmogorov-Smirnov test on each column of a matrix,
#' comparing the rows labeled \code{TRUE} with the rows labeled \code{FALSE}.
#'
#' For each column j, the statistic and the p-value are equal to those of
#' \code{ks.test(x[labels, j], x[!labels, j], alternative = alternative,
#' exact = FALSE)}: \code{D} is the maximal difference between the empirical
#' cumulative distribution functions of the two samples (two-sided), or the
#' maximal difference in the direction of \code{alternative}, and the p-value
#' is computed from the asymptotic (Kolmogorov) distribution of the statistic,
#' also when there are ties.
#'
#' Missing values (\code{NA}, \code{NaN}) are dropped from each column
#' separately, and rows whose label is \code{NA} are ignored. Infinite values
#' are kept, as in \code{ks.test}.
#'
#' @param x a numeric matrix or a data frame of numeric columns. Rows are
#' observations and columns are features.
#' @param labels a logical vector of length \code{nrow(x)}. Rows labeled
#' \code{TRUE} form the first sample and rows labeled \code{FALSE} the
#' second.
#' @param alternative one of \code{"two.sided"}, \code{"less"} or
#' \code{"greater"}, with the same meaning as in \code{ks.test}.
#' @return A numeric matrix with \code{ncol(x)} rows and 2 columns named
#' \code{"D"} and \code{"pval"}. Row names are the column names of \code{x}.
#' Columns where one of the samples has no non-missing values are \code{NA}.
#'
#' @seealso \code{\link[stats]{ks.test}}, \code{\link{tgs_wilcox}}
#'
#' @examples
#' \donttest{
#' # Note: all the available CPU cores might be used
#'
#' set.seed(42)
#' x <- matrix(rnorm(2000 * 100), ncol = 100)
#' labels <- rep(c(TRUE, FALSE), each = 1000)
#' x[labels, 1:10] <- x[labels, 1:10] + 0.5
#' result <- tgs_ks(x, labels)
#' head(result)
#' }
#'
#' \dontshow{
#' options(tgs_max.processes = 1)
#'
#' set.seed(42)
#' x <- matrix(rnorm(200 * 10), ncol = 10)
#' labels <- rep(c(TRUE, FALSE), each = 100)
#' result <- tgs_ks(x, labels)
#' }
#'
#' @export tgs_ks
tgs_ks <- function(x, labels, alternative = c("two.sided", "less", "greater")) {
    if (missing(x) || missing(labels)) {
        stop("Usage: tgs_ks(x, labels, alternative = \"two.sided\")", call. = FALSE)
    }
    alternative <- match.arg(alternative)
    x <- .tgs_as_matrix(x)

    .Call("tgs_ks", x, labels, alternative, new.env(parent = parent.frame()))
}

#' Vectorized Wilcoxon rank sum test
#'
#' Performs a Wilcoxon rank sum (Mann-Whitney) test on each column of a
#' matrix, comparing the rows labeled \code{TRUE} with the rows labeled
#' \code{FALSE}, and computes the area under the ROC curve (AUC).
#'
#' For each column j, the statistic and the p-value are equal to those of
#' \code{wilcox.test(x[labels, j], x[!labels, j], alternative = alternative,
#' exact = FALSE, correct = TRUE)}: \code{W} is the number of pairs of a
#' \code{TRUE} and a \code{FALSE} observation where the \code{TRUE}
#' observation is higher, with ties counted as 1/2, and the p-value is
#' computed from the normal approximation with correction for ties and
#' continuity correction. \code{auc = W / (n1 * n2)}, where \code{n1} and
#' \code{n2} are the numbers of non-missing observations of the column in each
#' sample: the probability that a \code{TRUE} observation is higher than a
#' \code{FALSE} one.
#'
#' Missing values (\code{NA}, \code{NaN}) are dropped from each column
#' separately, and rows whose label is \code{NA} are ignored. Infinite values
#' are kept, as in \code{wilcox.test}.
#'
#' @inheritParams tgs_ks
#' @param alternative one of \code{"two.sided"}, \code{"less"} or
#' \code{"greater"}, with the same meaning as in \code{wilcox.test}.
#' @return A numeric matrix with \code{ncol(x)} rows and 3 columns named
#' \code{"W"}, \code{"auc"} and \code{"pval"}. Row names are the column names
#' of \code{x}. Columns where one of the samples has no non-missing values are
#' \code{NA}. When the p-value is not defined (e.g. all values of the column
#' are tied and \code{alternative} is \code{"two.sided"}), it is \code{NA}.
#'
#' @seealso \code{\link[stats]{wilcox.test}}, \code{\link{tgs_ks}}
#'
#' @examples
#' \donttest{
#' # Note: all the available CPU cores might be used
#'
#' set.seed(42)
#' x <- matrix(rnorm(2000 * 100), ncol = 100)
#' labels <- rep(c(TRUE, FALSE), each = 1000)
#' x[labels, 1:10] <- x[labels, 1:10] + 0.5
#' result <- tgs_wilcox(x, labels)
#' head(result)
#' }
#'
#' \dontshow{
#' options(tgs_max.processes = 1)
#'
#' set.seed(42)
#' x <- matrix(rnorm(200 * 10), ncol = 10)
#' labels <- rep(c(TRUE, FALSE), each = 100)
#' result <- tgs_wilcox(x, labels)
#' }
#'
#' @export tgs_wilcox
tgs_wilcox <- function(x, labels, alternative = c("two.sided", "less", "greater")) {
    if (missing(x) || missing(labels)) {
        stop("Usage: tgs_wilcox(x, labels, alternative = \"two.sided\")", call. = FALSE)
    }
    alternative <- match.arg(alternative)
    x <- .tgs_as_matrix(x)

    .Call("tgs_wilcox", x, labels, alternative, new.env(parent = parent.frame()))
}
