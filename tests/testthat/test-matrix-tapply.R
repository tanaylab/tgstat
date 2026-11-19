# reference: f applied to the values of each row (margin 2) or column (margin 1) of x grouped by index, NA for empty groups
tapply_ref <- function(x, index, f, margin = 2) {
    index <- as.factor(index)
    x <- as.matrix(x)
    res <- suppressWarnings(apply(x, 3 - margin, function(v) as.numeric(tapply(v, index, f))))
    matrix(res, nrow = nlevels(index), dimnames = list(levels(index), dimnames(x)[[3 - margin]]))
}

# reference for sum and mean: the values are added one by one in increasing order of column (margin 2) or row (margin 1),
# starting from "init", and a NaN result is returned as NA
sum_loop <- function(x, index, margin = 2, na.rm = FALSE, mean = FALSE, init = 0) {
    index <- as.factor(index)
    x <- as.matrix(x)
    if (margin == 1) {
        x <- t(x)
    }
    res <- matrix(NA_real_, nlevels(index), nrow(x), dimnames = list(levels(index), rownames(x)))
    for (g in seq_len(nlevels(index))) {
        cols <- which(as.integer(index) == g)
        if (length(cols) == 0) {
            next
        }
        for (i in seq_len(nrow(x))) {
            s <- init
            n <- 0
            for (j in cols) {
                if (!na.rm || !is.na(x[i, j])) {
                    s <- s + x[i, j]
                    n <- n + 1
                }
            }
            res[g, i] <- if (is.na(s) || (mean && n == 0)) NA_real_ else if (mean) s / n else s
        }
    }
    res
}

with_processes <- function(num_processes, expr) {
    old <- options(tgs_max.processes = num_processes)
    on.exit(options(old))
    expr
}

# 60 x 50 so that both margins use several processes when allowed
dense_double <- function(nr = 60, nc = 50) {
    x <- matrix(round(rnorm(nr * nc), 1), nr, nc)
    x[sample(length(x), length(x) / 3)] <- 0
    x[sample(length(x), 30)] <- NA
    x[sample(length(x), 15)] <- NaN
    x[sample(length(x), 10)] <- Inf
    x[sample(length(x), 10)] <- -Inf
    x
}

dense_integer <- function(nr = 60, nc = 50) {
    x <- matrix(sample(-5:5, nr * nc, replace = TRUE), nr, nc)
    x[sample(length(x), 30)] <- NA
    x
}

sparse_double <- function(nr = 60, nc = 50) {
    x <- Matrix::rsparsematrix(nr, nc, 0.3, rand.x = function(n) round(rnorm(n), 1))
    x@x[sample(length(x@x), 20)] <- NA
    x@x[sample(length(x@x), 10)] <- NaN
    x@x[sample(length(x@x), 5)] <- Inf
    x@x[sample(length(x@x), 5)] <- 0 # explicit zeros
    x
}

# an index with an empty level, a level with a single element and a level whose elements are all NA in the first row (margin 2)
# or column (margin 1) of x
test_index <- function(n) {
    index <- c(rep(c("a", "b", "c", "d"), length.out = n - 3), "e", "f", "f")
    factor(index, levels = c("a", "b", "c", "d", "e", "f", "empty"))
}

na_first <- function(x, index, margin) {
    if (margin == 2) {
        x[1, index == "f"] <- NA
    } else {
        x[index == "f", 1] <- NA
    }
    x
}

test_matrices <- function() {
    res <- list(dense_double = dense_double(), dense_integer = dense_integer())
    if (requireNamespace("Matrix", quietly = TRUE)) {
        res$sparse <- sparse_double()
    }
    res
}

test_that("native max, min and median match base R and the generic path", {
    set.seed(60427)
    for (x in test_matrices()) {
        for (margin in 1:2) {
            index <- test_index(if (margin == 1) nrow(x) else ncol(x))
            x <- na_first(x, index, margin)
            for (fname in c("max", "min", "median")) {
                f <- match.fun(fname)
                for (np in c(1, 8)) {
                    with_processes(np, {
                        res <- switch(fname,
                            max = tgs_matrix_tapply(x, index, max, margin = margin),
                            min = tgs_matrix_tapply(x, index, min, margin = margin),
                            median = tgs_matrix_tapply(x, index, median, margin = margin)
                        )
                        res_na_rm <- switch(fname,
                            max = tgs_matrix_tapply(x, index, max, na.rm = TRUE, margin = margin),
                            min = tgs_matrix_tapply(x, index, min, na.rm = TRUE, margin = margin),
                            median = tgs_matrix_tapply(x, index, median, na.rm = TRUE, margin = margin)
                        )
                        expect_identical(res, tapply_ref(x, index, f, margin))
                        expect_identical(res_na_rm, tapply_ref(x, index, function(v) f(v, na.rm = TRUE), margin))
                        expect_identical(res, tgs_matrix_tapply(x, index, function(v) f(v), margin = margin))
                        expect_identical(res_na_rm, tgs_matrix_tapply(x, index, function(v) f(v, na.rm = TRUE), margin = margin))
                    })
                }
            }
        }
    }
})

test_that("all-NA and empty groups give what base R gives", {
    x <- rbind(c(NA, NA, NaN, 5), c(NA, 2, NA, 6))
    index <- factor(c("a", "a", "b", "c"), levels = c("a", "b", "c", "d"))
    expect_identical(tgs_matrix_tapply(x, index, max, na.rm = TRUE)[, 1], c(a = -Inf, b = -Inf, c = 5, d = NA))
    expect_identical(tgs_matrix_tapply(x, index, min, na.rm = TRUE)[, 1], c(a = Inf, b = Inf, c = 5, d = NA))
    expect_identical(tgs_matrix_tapply(x, index, median, na.rm = TRUE)[, 1], c(a = NA, b = NA, c = 5, d = NA))
    expect_identical(tgs_matrix_tapply(x, index, max)[, 1], c(a = NA, b = NaN, c = 5, d = NA))
    expect_identical(tgs_matrix_tapply(x, index, max)[, 2], c(a = NA, b = NA, c = 6, d = NA))
    expect_identical(tgs_matrix_tapply(x, index, median)[, 1], c(a = NA, b = NA, c = 5, d = NA))
    expect_identical(tgs_matrix_tapply(x, index, max, na.rm = TRUE)[, 2], c(a = 2, b = -Inf, c = 6, d = NA))
    for (f in c("max", "min", "median")) {
        expect_identical(
            do.call(tgs_matrix_tapply, list(t(x), index, as.name(f), na.rm = TRUE, margin = 1)),
            do.call(tgs_matrix_tapply, list(x, index, as.name(f), na.rm = TRUE))
        )
    }
})

test_that("the median of an even count is the mean of the two middle values as R computes it", {
    set.seed(60427)
    x <- matrix(runif(400 * 2) * sample(c(-1, 1), 800, replace = TRUE) * 10^sample(-5:300, 800, replace = TRUE), 400, 2)
    expect_identical(tgs_matrix_tapply(x, c(1, 1), median), tapply_ref(x, c(1, 1), median))
    # pairs whose sum overflows, infinite pairs
    x <- matrix(c(.Machine$double.xmax, .Machine$double.xmax, 1.5e308, 1.6e308, -1.7e308, -1.1e308, -Inf, Inf, Inf, 1), 2)
    expect_identical(tgs_matrix_tapply(x, c(1, 1), median, margin = 1), tapply_ref(x, c(1, 1), median, margin = 1))
    expect_identical(tgs_matrix_tapply(t(x), c(1, 1), median), tapply_ref(t(x), c(1, 1), median))
})

test_that("sum and mean add the values in the same order as before", {
    set.seed(60427)
    for (x in test_matrices()) {
        for (margin in 1:2) {
            index <- test_index(if (margin == 1) nrow(x) else ncol(x))
            x <- na_first(x, index, margin)
            for (np in c(1, 8)) {
                with_processes(np, {
                    expect_identical(tgs_matrix_tapply(x, index, sum, margin = margin), sum_loop(x, index, margin))
                    expect_identical(tgs_matrix_tapply(x, index, sum, na.rm = TRUE, margin = margin), sum_loop(x, index, margin, na.rm = TRUE))
                    expect_identical(tgs_matrix_tapply(x, index, sum, 10, na.rm = TRUE, margin = margin), sum_loop(x, index, margin, na.rm = TRUE, init = 10))
                    expect_identical(tgs_matrix_tapply(x, index, mean, margin = margin), sum_loop(x, index, margin, mean = TRUE))
                    expect_identical(tgs_matrix_tapply(x, index, mean, na.rm = TRUE, margin = margin), sum_loop(x, index, margin, na.rm = TRUE, mean = TRUE))
                })
            }
            expect_equal(tgs_matrix_tapply(x, index, sum, na.rm = TRUE, margin = margin), tapply_ref(x, index, function(v) sum(v, na.rm = TRUE), margin))
            expect_equal(tgs_matrix_tapply(x, index, mean, na.rm = TRUE, margin = margin), tapply_ref(x, index, function(v) mean(v, na.rm = TRUE), margin))
        }
    }
})

test_that("sum and mean of an integer matrix with NA are NA unless na.rm", {
    x <- matrix(1:12, 3, 4)
    x[2, 1] <- NA
    index <- factor(c(1, 1, 2, 2))
    expect_identical(tgs_matrix_tapply(x, index, sum)[, 2], c(`1` = NA, `2` = 19))
    expect_identical(tgs_matrix_tapply(x, index, sum, na.rm = TRUE)[, 2], c(`1` = 5, `2` = 19))
    expect_identical(tgs_matrix_tapply(x, index, mean)[, 2], c(`1` = NA, `2` = 9.5))
    expect_identical(tgs_matrix_tapply(x, index, mean, na.rm = TRUE)[, 2], c(`1` = 5, `2` = 9.5))
    expect_identical(tgs_matrix_tapply(x, index, function(v) max(v))[, 2], c(`1` = NA, `2` = 11))
})

test_that("margin 1 is the same as margin 2 on the transposed matrix and as rowsum", {
    set.seed(60427)
    x <- matrix(rnorm(3000 * 12), 3000, 12, dimnames = list(NULL, paste0("c", 1:12)))
    index <- sample(c("x", "y", "z", "w"), 3000, replace = TRUE)
    expect_identical(tgs_matrix_tapply(x, index, sum, margin = 1), rowsum(x, index))
    for (f in c("sum", "mean", "max", "min", "median")) {
        expect_identical(
            do.call(tgs_matrix_tapply, list(x, index, as.name(f), margin = 1)),
            do.call(tgs_matrix_tapply, list(t(x), index, as.name(f)))
        )
    }
    skip_if_not_installed("Matrix")
    s <- Matrix::rsparsematrix(3000, 12, 0.1)
    expect_identical(tgs_matrix_tapply(s, index, max, margin = 1), tgs_matrix_tapply(Matrix::t(s), index, max))
    expect_identical(tgs_matrix_tapply(s, index, sum, margin = 1), tgs_matrix_tapply(Matrix::t(s), index, sum))
})

test_that("R functions get the values in order, in both margins", {
    set.seed(60427)
    f <- function(v) sum(v * seq_along(v))
    for (x in test_matrices()) {
        for (margin in 1:2) {
            index <- test_index(if (margin == 1) nrow(x) else ncol(x))
            for (np in c(1, 8)) {
                with_processes(np, {
                    expect_identical(tgs_matrix_tapply(x, index, f, margin = margin), tapply_ref(x, index, f, margin))
                    expect_identical(tgs_matrix_tapply(x, index, function(v) v[length(v)], margin = margin), tapply_ref(x, index, function(v) v[length(v)], margin))
                })
            }
        }
    }
})

test_that("extra arguments are passed to R functions, other than margin", {
    set.seed(60427)
    x <- matrix(round(rnorm(40 * 30), 1), 40, 30)
    index <- sample(1:4, 30, replace = TRUE)
    expect_identical(tgs_matrix_tapply(x, index, function(v, a) sum(v) + a, a = 100), tapply_ref(x, index, function(v) sum(v) + 100))
    expect_identical(tgs_matrix_tapply(x, index, function(v, a) sum(v) + a, 100), tapply_ref(x, index, function(v) sum(v) + 100))
    expect_identical(tgs_matrix_tapply(x, index, quantile, 0.3), tapply_ref(x, index, function(v) quantile(v, 0.3)))
    expect_identical(tgs_matrix_tapply(x, index, mean, trim = 0.2), tapply_ref(x, index, function(v) mean(v, trim = 0.2)))
    expect_identical(tgs_matrix_tapply(x, index, max, 1), tapply_ref(x, index, function(v) max(v, 1)))
    expect_identical(tgs_matrix_tapply(x, index, median, TRUE), tapply_ref(x, index, function(v) median(v, TRUE)))
    expect_identical(tgs_matrix_tapply(x, index, function(v, ...) length(list(...)), margin = 2), tapply_ref(x, index, function(v) 0))
    expect_identical(tgs_matrix_tapply(t(x), index, function(v, a) sum(v) + a, a = 100, margin = 1), tapply_ref(x, index, function(v) sum(v) + 100))
})

test_that("tgs_matrix_tapply keeps the names and handles a non-factor index", {
    x <- matrix(1:12 + 0, 3, 4, dimnames = list(c("r1", "r2", "r3"), c("c1", "c2", "c3", "c4")))
    index <- c("b", "a", "b", "a")
    res <- tgs_matrix_tapply(x, index, sum)
    expect_identical(dimnames(res), list(c("a", "b"), c("r1", "r2", "r3")))
    expect_identical(res, tapply_ref(x, index, sum))
    res <- tgs_matrix_tapply(x, c(2, 1, 2), max, margin = 1)
    expect_identical(dimnames(res), list(c("1", "2"), c("c1", "c2", "c3", "c4")))
    expect_identical(res, tapply_ref(x, c(2, 1, 2), max, margin = 1))
    expect_null(colnames(tgs_matrix_tapply(unname(x), index, sum)))
    skip_if_not_installed("Matrix")
    s <- Matrix::Matrix(x, sparse = TRUE)
    expect_identical(tgs_matrix_tapply(s, index, median), tapply_ref(x, index, median))
    expect_identical(tgs_matrix_tapply(s, c(2, 1, 2), mean, margin = 1), tapply_ref(x, c(2, 1, 2), mean, margin = 1))
})

test_that("tgs_matrix_tapply checks margin and the length of index", {
    x <- matrix(1:12 + 0, 3, 4)
    expect_error(tgs_matrix_tapply(x, 1:4, sum, margin = 3), "margin")
    expect_error(tgs_matrix_tapply(x, 1:4, sum, margin = "1"), "margin")
    expect_error(tgs_matrix_tapply(x, 1:4, sum, margin = NA), "margin")
    expect_error(tgs_matrix_tapply(x, 1:4, sum, margin = 1:2), "margin")
    expect_error(tgs_matrix_tapply(x, 1:4, sum, margin = 1), "rows")
    expect_error(tgs_matrix_tapply(x, 1:3, sum), "columns")
})
