#' Vectorized chi-squared test for 2x2 contingency tables
#'
#' Performs a chi-squared test with optional Yates' continuity correction on
#' each row of a two-column count matrix. For each row, a 2x2 contingency
#' table is constructed using the row counts and the column sums, and the
#' chi-squared statistic and p-value are computed.
#'
#' This function is useful for differential gene expression analysis, where
#' each row represents a gene and the two columns represent UMI counts in two
#' conditions. The test determines whether the gene's proportion differs
#' significantly between conditions.
#'
#' For each row i, the contingency table is:
#' \tabular{lcc}{
#'   \tab Condition 1 \tab Condition 2 \cr
#'   Gene i \tab x\[i,1\] \tab x\[i,2\] \cr
#'   Other genes \tab colSum1 - x\[i,1\] \tab colSum2 - x\[i,2\] \cr
#' }
#'
#' @param x a numeric matrix or sparse matrix of \code{dgCMatrix} type with
#' exactly 2 columns containing non-negative counts.
#' @param yates logical; if \code{TRUE} (default), Yates' continuity
#' correction is applied to the chi-squared statistic.
#' @return A numeric matrix with \code{nrow(x)} rows and 2 columns named
#' \code{"chi2"} and \code{"pval"}. Row names from \code{x} are preserved.
#' When the denominator of the chi-squared formula is zero, the statistic
#' is set to 0 and the p-value to 1.
#'
#' @examples
#' \donttest{
#' # Note: all the available CPU cores might be used
#'
#' set.seed(42)
#' # Simulate UMI counts for 1000 genes in 2 conditions
#' mat <- matrix(rpois(2000, lambda = 100), ncol = 2)
#' rownames(mat) <- paste0("gene", 1:1000)
#' colnames(mat) <- c("condition1", "condition2")
#' result <- tgs_chi2(mat)
#' head(result)
#'
#' # Without Yates' correction
#' result_no_yates <- tgs_chi2(mat, yates = FALSE)
#'
#' # With sparse matrix
#' sparse_mat <- Matrix::Matrix(mat, sparse = TRUE)
#' result_sparse <- tgs_chi2(sparse_mat)
#' }
#'
#' \dontshow{
#' options(tgs_max.processes = 1)
#'
#' set.seed(42)
#' mat <- matrix(rpois(200, lambda = 100), ncol = 2)
#' rownames(mat) <- paste0("gene", 1:100)
#' result <- tgs_chi2(mat)
#' }
#'
#' @export tgs_chi2
tgs_chi2 <- function(x, yates = TRUE) {
    if (missing(x)) {
        stop("Usage: tgs_chi2(x, yates = TRUE)", call. = FALSE)
    }

    .Call("tgs_chi2", x, yates, new.env())
}
