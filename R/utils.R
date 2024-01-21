.tgs_call <- function(...) {
    tryCatch(
        {
            res <- .Call(...)
        },
        interrupt = function(interrupt) {
            stop("Command interrupted!", call. = FALSE)
        }
    )
    return(res)
}

.tgs_use_blas <- function() {
    opt <- getOption("tgs_use.blas")
    if (is.null(opt)) {
        opt <- .tgs_guess_use_blas()
    } else {
        opt <- as.logical(opt)
    }
    return(opt)
}

.tgs_guess_use_blas <- function() {
    blas <- basename(extSoftVersion()["BLAS"])
    if (grepl("openblas", blas, ignore.case = TRUE) || grepl("mkl", blas, ignore.case = TRUE) || grepl("atlas", blas, ignore.case = TRUE) || grepl("gsl", blas, ignore.case = TRUE)) {
        return(TRUE)
    }
    return(FALSE)
}

#' Checks whether all the elements of the vector are finite
#'
#' Checks whether all the elements of the vector are finite.
#'
#' 'tgs_finite' returns 'TRUE' if all the elements of 'x' are finite numbers.
#' (See: 'is.finite'.)
#'
#' @param x numeric or integer vector or matrix
#' @return 'TRUE' if all the elements of 'x' are finite, otherwise 'FALSE'.
#' @keywords ~finite
#' @examples
#'
#' tgs_finite(1:10)
#' tgs_finite(c(1:10, NaN))
#' tgs_finite(c(1:10, Inf))
#'
#' @export tgs_finite
tgs_finite <- function(x) {
    if (missing(x)) {
        stop("Usage: tgs_finite(x)", call. = FALSE)
    }

    .Call("tgs_finite", x, new.env(parent = parent.frame()))
}
