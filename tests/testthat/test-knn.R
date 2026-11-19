# reference implementation: the highest knn values of each column, ties broken by the row number
knn_ref <- function(m, knn, diag = FALSE, threshold = 0) {
    res <- lapply(seq_len(ncol(m)), function(j) {
        v <- m[, j]
        ok <- is.finite(v) & v >= threshold
        if (!diag && j <= nrow(m)) {
            ok[j] <- FALSE
        }
        rows <- which(ok)
        rows <- utils::head(rows[order(-v[rows], rows)], knn)
        data.frame(col1 = rep(j, length(rows)), col2 = rows, val = v[rows], rank = seq_along(rows))
    })
    do.call(rbind, res)
}

expect_knn_equal <- function(m, ...) {
    res <- tgs_knn(m, ...)
    ref <- knn_ref(m, ...)
    expect_equal(as.integer(res$col1), ref$col1)
    expect_equal(as.integer(res$col2), ref$col2)
    expect_equal(res$val, ref$val)
    expect_equal(res$rank, ref$rank)
}

test_that("tgs_knn of a matrix matches the reference", {
    set.seed(60427)
    m <- matrix(rnorm(300 * 300), 300, 300)
    expect_knn_equal(m, 10)
    expect_knn_equal(m, 10, diag = TRUE)
    expect_knn_equal(m, 10, threshold = 0.5)
    expect_knn_equal(m, 10, threshold = -Inf)
    expect_knn_equal(m, 1000, threshold = -Inf)
})

test_that("tgs_knn of a matrix handles ties, non-finite values and non-square matrices", {
    set.seed(60427)
    m <- matrix(sample(0:5, 200 * 120, replace = TRUE) + 0, 200, 120)
    m[sample(length(m), 500)] <- NA
    m[sample(length(m), 100)] <- Inf
    m[sample(length(m), 100)] <- -Inf
    m[, 7] <- NA
    expect_knn_equal(m, 15)
    expect_knn_equal(m, 15, diag = TRUE, threshold = 2)
    expect_knn_equal(t(m), 15, threshold = -Inf)
})

test_that("tgs_knn of a matrix does not depend on the number of processes", {
    set.seed(60427)
    m <- matrix(rnorm(1500 * 1500), 1500, 1500)
    res <- lapply(c(1, 3, 7), function(np) {
        old <- options(tgs_max.processes = np)
        on.exit(options(old))
        tgs_knn(m, 20)
    })
    expect_identical(res[[1]], res[[2]])
    expect_identical(res[[1]], res[[3]])
    expect_knn_equal(m, 20)
})

test_that("tgs_knn keeps the column names and returns an empty data frame when nothing passes the threshold", {
    m <- matrix(c(1, 5, 3, 2, 9, 4, 7, 8, 6), 3, 3, dimnames = list(NULL, c("a", "b", "c")))
    res <- tgs_knn(m, 1)
    expect_equal(as.character(res$col1), c("a", "b", "c"))
    expect_equal(as.character(res$col2), c("b", "c", "b"))
    expect_equal(res$val, c(5, 4, 8))
    res <- tgs_knn(m, 2, threshold = 100)
    expect_equal(nrow(res), 0)
    expect_named(res, c("col1", "col2", "val", "rank"))
})

test_that("tgs_knn of a tidy data frame matches tgs_knn of the matrix", {
    set.seed(60427)
    m <- matrix(rnorm(80 * 80), 80, 80, dimnames = list(NULL, paste0("c", 1:80)))
    tidy <- tgs_knn(m, 80, threshold = -Inf)
    res <- tgs_knn(tidy, 5, threshold = -Inf)
    ref <- tgs_knn(m, 5, threshold = -Inf)
    expect_equal(res, ref, ignore_attr = TRUE)
})
