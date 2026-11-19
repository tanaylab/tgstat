test_that("tgs_dist works correctly", {
    # Create test matrix
    set.seed(60427)
    x <- matrix(rnorm(100), nrow = 10)

    # Compare to base R dist
    expect_equal(as.matrix(tgs_dist(x)), as.matrix(dist(x)))

    # Test with diag = TRUE
    diag_res <- tgs_dist(x, diag = TRUE)
    expect_equal(diag_res, dist(x, diag = TRUE), ignore_attr = TRUE)

    # Test with upper = TRUE
    upper_res <- as.matrix(tgs_dist(x, upper = TRUE))
    expect_equal(upper_res, as.matrix(dist(x, upper = TRUE)))

    # Test tidy output
    tidy_res <- tgs_dist(x, tidy = TRUE)
    expect_s3_class(tidy_res, "data.frame")
    expect_named(tidy_res, c("row1", "row2", "dist"))

    # Test threshold
    thresh_res <- tgs_dist(x, tidy = TRUE, threshold = 1)
    expect_true(all(thresh_res$dist <= 1))

    # Test with NA values
    x[1, 1] <- NA
    expect_equal(tgs_dist(x), dist(x, method = "euclidean", diag = FALSE, upper = FALSE), ignore_attr = TRUE)

    # Test with different number of rows and columns
    y <- matrix(rnorm(50), nrow = 5, ncol = 10)
    expect_equal(as.matrix(tgs_dist(y)), as.matrix(dist(y)), ignore_attr = TRUE)

    # Test with non-numeric input
    z <- matrix(letters[1:100], nrow = 10)
    expect_error(tgs_dist(z))
})

test_that("tgs_dist matches dist with and without BLAS on more points than one block", {
    set.seed(60427)
    x <- matrix(rnorm(3000 * 5), nrow = 3000, dimnames = list(paste0("p", 1:3000), NULL))
    x_na <- x
    x_na[sample(length(x), 100)] <- NA
    old <- getOption("tgs_use.blas")
    on.exit(options(tgs_use.blas = old))

    for (use_blas in c(TRUE, FALSE)) {
        options(tgs_use.blas = use_blas)
        for (m in list(x, x_na)) {
            d <- dist(m)
            expect_equal(tgs_dist(m), d, ignore_attr = "call")
            tidy <- tgs_dist(m, tidy = TRUE, threshold = 1)
            expected <- which(as.vector(d) <= 1)
            expect_equal(nrow(tidy), length(expected))
            expect_equal(tidy$dist, as.vector(d)[expected])
            expect_equal(levels(tidy$row1), rownames(x))
        }
    }
})

test_that("tgs_dist computes all the rows with process counts that are not powers of 2", {
    set.seed(60427)
    x <- matrix(rnorm(120 * 5), nrow = 120)
    old <- options(tgs_use.blas = FALSE, tgs_max.processes = 1)
    on.exit(options(old))
    for (p in c(3, 7, 11, 13, 50)) {
        options(tgs_max.processes = p)
        for (n in c(30, 31, 120)) {
            expect_equal(tgs_dist(x[1:n, ]), dist(x[1:n, ]), ignore_attr = TRUE)
        }
    }
})

test_that("tgs_dist without BLAS copies results larger than one 64 MB chunk", {
    set.seed(60427)
    x <- matrix(rnorm(4200 * 3), nrow = 4200) # 4200 * 4199 / 2 distances > 2^23
    x[sample(length(x), 50)] <- NA
    old <- options(tgs_use.blas = FALSE)
    on.exit(options(old))
    expect_equal(tgs_dist(x), dist(x), ignore_attr = TRUE)
})

test_that("tgs_dist with BLAS and NAs gives no NaN for identical points", {
    old <- options(tgs_use.blas = TRUE)
    on.exit(options(old))
    for (seed in 1:20) {
        set.seed(seed)
        x <- matrix(runif(40 * 30) * 10, nrow = 40)
        x[2, ] <- x[1, ]
        x[4, ] <- x[3, ]
        x[sample(length(x), 50)] <- NA
        d <- tgs_dist(x)
        expect_false(any(is.nan(d)))
        expect_equal(d, dist(x), ignore_attr = TRUE)
    }

    # points that share no dimension: NA, as without BLAS
    x <- rbind(c(1, NA), c(NA, 2), c(3, 4))
    expected <- c(NA, 2 * sqrt(2), 2 * sqrt(2))
    expect_equal(as.vector(tgs_dist(x)), expected)
    expect_false(is.nan(tgs_dist(x)[1]))
    options(tgs_use.blas = FALSE)
    expect_equal(as.vector(tgs_dist(x)), expected)
    expect_false(is.nan(tgs_dist(x)[1]))
})

# reference: the distances between the rows of x and the rows of y
cross_dist_ref <- function(x, y) {
    # tgstat treats infinite values as missing
    x[!is.finite(x)] <- NA
    y[!is.finite(y)] <- NA
    unname(as.matrix(dist(rbind(x, y)))[seq_len(nrow(x)), nrow(x) + seq_len(nrow(y)), drop = FALSE])
}

with_blas <- function(use_blas, expr) {
    old <- options(tgs_use.blas = use_blas)
    on.exit(options(old))
    expr
}

test_that("tgs_dist between two matrices matches dist with and without BLAS", {
    set.seed(60427)
    x <- matrix(rnorm(30 * 6), 30)
    y <- matrix(rnorm(12 * 6), 12)
    x_na <- x
    y_na <- y
    x_na[sample(length(x), 20)] <- NA
    x_na[7, 3] <- Inf
    y_na[sample(length(y), 10)] <- NaN
    y_na[4, ] <- NA
    # rows that share no dimension
    x_na[10, ] <- c(NA, NA, NA, 1, 2, 3)
    y_na[8, ] <- c(1, 2, 3, NA, NA, NA)
    xi <- matrix(sample(0:20, 30 * 6, replace = TRUE), 30)
    yi <- matrix(sample(0:20, 12 * 6, replace = TRUE), 12)
    yi[2, 5] <- NA

    for (use_blas in c(TRUE, FALSE)) {
        for (m in list(list(x, y), list(x_na, y_na), list(xi, yi), list(x, y[1, , drop = FALSE]), list(x[1, , drop = FALSE], y))) {
            expect_equal(with_blas(use_blas, tgs_dist(m[[1]], m[[2]])), cross_dist_ref(m[[1]], m[[2]]))
        }
        res <- with_blas(use_blas, tgs_dist(x_na, y_na))
        expect_true(is.na(res[10, 8]))
        expect_false(any(is.nan(res)))
        expect_true(all(is.na(res[, 4])))
    }
})

test_that("tgs_dist between two matrices returns the row names and the tidy format", {
    set.seed(60427)
    x <- matrix(rnorm(20 * 4), 20, dimnames = list(paste0("x", 1:20), NULL))
    y <- matrix(rnorm(15 * 4), 15, dimnames = list(paste0("y", 1:15), NULL))
    y[3, 2] <- NA
    y[4, ] <- NA
    ref <- cross_dist_ref(x, y)

    for (use_blas in c(TRUE, FALSE)) {
        res <- with_blas(use_blas, tgs_dist(x, y))
        expect_equal(dimnames(res), list(rownames(x), rownames(y)))
        expect_equal(dimnames(with_blas(use_blas, tgs_dist(x, unname(y)))), list(rownames(x), NULL))
        expect_null(dimnames(with_blas(use_blas, tgs_dist(unname(x), unname(y)))))

        for (threshold in c(Inf, 2)) {
            tidy <- with_blas(use_blas, tgs_dist(x, y, tidy = TRUE, threshold = threshold))
            expect_named(tidy, c("row1", "row2", "dist"))
            expect_equal(levels(tidy$row1), rownames(x))
            expect_equal(levels(tidy$row2), rownames(y))
            # the pairs ordered by the row of x, then by the row of y
            pairs <- which(t(ref) <= threshold, arr.ind = TRUE)
            expect_equal(as.integer(tidy$row1), unname(pairs[, 2]))
            expect_equal(as.integer(tidy$row2), unname(pairs[, 1]))
            expect_equal(tidy$dist, t(ref)[pairs])
        }

        tidy <- with_blas(use_blas, tgs_dist(unname(x), unname(y), tidy = TRUE))
        expect_type(tidy$row1, "integer")
        expect_type(tidy$row2, "integer")
        expect_equal(nrow(tidy), sum(!is.na(ref)))
    }
})

test_that("tgs_dist between two matrices with BLAS on more columns than one block matches the version without BLAS", {
    set.seed(60427)
    # 3000 x 7000 distances: blocks of 2^23 / 3000 columns
    x <- matrix(rnorm(3000 * 4), 3000)
    y <- matrix(rnorm(7000 * 4), 7000)
    y_na <- y
    y_na[sample(length(y), 100)] <- NA
    for (m in list(y, y_na)) {
        res <- with_blas(TRUE, tgs_dist(x, m))
        expect_equal(res, with_blas(FALSE, tgs_dist(x, m)))
        i <- sample(nrow(x), 20)
        j <- sample(nrow(m), 20)
        expect_equal(res[cbind(i, j)], diag(cross_dist_ref(x[i, ], m[j, ])))
    }
})

test_that("tgs_dist between two matrices does not depend on the number of processes", {
    set.seed(60427)
    x <- matrix(rnorm(120 * 5), 120)
    y <- matrix(rnorm(30 * 5), 30)
    x[sample(length(x), 30)] <- NA
    old <- options(tgs_use.blas = FALSE, tgs_max.processes = 1)
    on.exit(options(old))
    expected <- tgs_dist(x, y)
    expect_equal(expected, cross_dist_ref(x, y))
    for (p in c(3, 7, 11, 13, 50)) {
        options(tgs_max.processes = p)
        expect_identical(tgs_dist(x, y), expected)
        expect_identical(tgs_dist(x[1:31, ], y), expected[1:31, ])
    }
})

test_that("tgs_dist between two matrices gives 0 for identical points", {
    for (seed in 1:4) {
        set.seed(seed)
        x <- matrix(runif(40 * 30) * 1000, 40)
        x_na <- x
        x_na[sample(length(x), 50)] <- NA
        for (use_blas in c(TRUE, FALSE)) {
            for (m in list(x, x_na)) {
                res <- with_blas(use_blas, tgs_dist(m, m[c(5, 1, 17, 17, 40), ]))
                expect_identical(res[cbind(c(5, 1, 17, 17, 40), 1:5)], rep(0, 5))
                expect_false(any(is.nan(res)))
                expect_equal(res, cross_dist_ref(m, m[c(5, 1, 17, 17, 40), ]))
            }
        }
    }
})

test_that("tgs_dist accepts data frames and rejects invalid y", {
    set.seed(60427)
    x <- matrix(rnorm(10 * 3), 10)
    y <- matrix(rnorm(5 * 3), 5)
    for (use_blas in c(TRUE, FALSE)) {
        expect_equal(with_blas(use_blas, tgs_dist(as.data.frame(x), as.data.frame(y))), cross_dist_ref(x, y))
        expect_equal(with_blas(use_blas, as.matrix(tgs_dist(as.data.frame(x)))), as.matrix(dist(x)), ignore_attr = TRUE)
        expect_error(with_blas(use_blas, tgs_dist(x, y[, 1:2])), "different number of columns")
        expect_error(with_blas(use_blas, tgs_dist(x, 1:3)), "matrix of numeric values")
        expect_error(with_blas(use_blas, tgs_dist(x, matrix(letters[1:6], 2))), "matrix of numeric values")
    }

    # "diag" used to be the second argument
    expect_error(tgs_dist(x, TRUE), "diag")
    expect_error(tgs_dist(x, FALSE, TRUE), "diag")
    expect_equal(tgs_dist(x, diag = TRUE), dist(x, diag = TRUE), ignore_attr = "call")
})
