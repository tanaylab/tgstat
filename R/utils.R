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

.tgs_getOption <- function(x, default = NULL) {
    if (missing(default)) {
        return(options(x)[[1L]])
    }
    if (x %in% names(options())) {
        return(options(x)[[1L]])
    } else {
        return(default)
    }
}

.tgs_use_blas <- function() {
    .tgs_getOption("tgs_use.blas", F)
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
