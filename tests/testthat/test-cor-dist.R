expect_cor_dist <- function(x, pairwise.complete.obs = FALSE, spearman = FALSE) {
    old <- getOption("tgs_use.blas")
    on.exit(options(tgs_use.blas = old))
    for (use_blas in c(TRUE, FALSE)) {
        options(tgs_use.blas = use_blas)
        expect_equal(
            tgs_cor_dist(x, pairwise.complete.obs = pairwise.complete.obs, spearman = spearman),
            tgs_dist(tgs_cor(x, pairwise.complete.obs = pairwise.complete.obs, spearman = spearman)),
            ignore_attr = "call"
        )
    }
}

test_that("tgs_cor_dist equals tgs_dist(tgs_cor())", {
    set.seed(60427)
    x <- matrix(rnorm(30 * 500), nrow = 30, dimnames = list(NULL, paste0("c", 1:500)))
    expect_cor_dist(x)
    expect_cor_dist(x, spearman = TRUE)
    expect_cor_dist(x, pairwise.complete.obs = TRUE)

    # compare to base R
    expect_equal(as.vector(tgs_cor_dist(x)), as.vector(dist(cor(x))))
    expect_equal(as.vector(tgs_cor_dist(x, spearman = TRUE)), as.vector(dist(cor(x, method = "spearman"))))

    # integers with ties
    xi <- matrix(sample(1:10, 30 * 300, replace = TRUE), nrow = 30)
    expect_cor_dist(xi)
    expect_cor_dist(xi, spearman = TRUE)

    # more rows than columns
    xt <- matrix(rnorm(300 * 40), nrow = 300)
    expect_cor_dist(xt)
    expect_cor_dist(xt, spearman = TRUE)

    # rows between N / 3 and N / 2: tgs_dist(tgs_cor()) with BLAS, the fast path without
    xm <- matrix(rnorm(40 * 100), nrow = 40)
    expect_cor_dist(xm)
    expect_cor_dist(xm, spearman = TRUE)

    expect_equal(attr(tgs_cor_dist(x), "call"), quote(tgs_cor_dist(x = x)))
})

test_that("tgs_cor_dist falls back to tgs_dist(tgs_cor()) on missing values and constant columns", {
    set.seed(60427)
    x <- matrix(rnorm(20 * 100), nrow = 20) # few enough rows for the fast path with and without BLAS

    x_na <- x
    x_na[sample(length(x), 10)] <- NA
    expect_cor_dist(x_na)
    expect_cor_dist(x_na, pairwise.complete.obs = TRUE)
    expect_cor_dist(x_na, pairwise.complete.obs = TRUE, spearman = TRUE)

    x_const <- x
    x_const[, 3] <- 1
    expect_cor_dist(x_const)
    expect_cor_dist(x_const, spearman = TRUE)

    # constants whose computed mean is not exact
    for (val in c(0.1, log2(1e-5))) {
        x_const[, 3] <- val
        expect_cor_dist(x_const)
        expect_cor_dist(x_const, spearman = TRUE)
    }

    expect_error(tgs_cor_dist(letters))
    expect_error(tgs_cor_dist(x, pairwise.complete.obs = "a"), "pairwise.complete.obs")
    expect_error(tgs_cor_dist(x, spearman = "a"), "spearman")
})

test_that("tgs_cor_dist handles few columns or rows", {
    set.seed(60427)
    x <- matrix(rnorm(10 * 2), nrow = 10, dimnames = list(NULL, c("a", "b")))
    expect_cor_dist(x)
    expect_cor_dist(x[, 1, drop = FALSE])
    expect_length(tgs_cor_dist(x[, 1, drop = FALSE]), 0)
    expect_cor_dist(x[1:2, ])
    expect_cor_dist(x[1, , drop = FALSE])

    y <- matrix(rnorm(2 * 8), nrow = 2)
    expect_cor_dist(y)
    expect_cor_dist(y, spearman = TRUE)
})

test_that("tgs_standardize returns columns whose crossprod is the correlation matrix", {
    standardize <- function(x, spearman = FALSE) .Call("tgs_standardize", x, spearman, new.env(), PACKAGE = "tgstat")

    set.seed(60427)
    x <- matrix(sample(1:10, 20 * 50, replace = TRUE), nrow = 20, dimnames = list(NULL, paste0("c", 1:50)))
    expect_equal(crossprod(standardize(x)), cor(x))
    expect_equal(crossprod(standardize(x, spearman = TRUE)), cor(x, method = "spearman"))
    expect_equal(dimnames(standardize(x)), dimnames(x))

    x_na <- x
    x_na[3, 4] <- NA
    expect_null(standardize(x_na))
    x_const <- x
    x_const[, 4] <- 7L
    expect_null(standardize(x_const))
    expect_null(standardize(x[1, , drop = FALSE]))
    expect_error(standardize(matrix(letters[1:20], nrow = 2)))

    # the ranking of large matrices is split between processes
    x <- matrix(sample(1:50, 100 * 20000, replace = TRUE), nrow = 100)
    old <- getOption("tgs_max.processes")
    on.exit(options(tgs_max.processes = old))
    options(tgs_max.processes = 1)
    z1 <- standardize(x, spearman = TRUE)
    options(tgs_max.processes = 4)
    expect_identical(standardize(x, spearman = TRUE), z1)
    expect_equal(crossprod(z1[, 1:100]), cor(x[, 1:100], method = "spearman"))
})

test_that("tgs_cor_dist stays accurate for columns with large offsets", {
    set.seed(60427)
    x <- matrix(rnorm(20 * 300), nrow = 20) + 1e5
    expect_equal(as.vector(tgs_cor_dist(x)), as.vector(dist(cor(x))), tolerance = 1e-12)
    expect_equal(as.vector(tgs_cor_dist(x, spearman = TRUE)), as.vector(dist(cor(x, method = "spearman"))), tolerance = 1e-12)
})

test_that("tgs_cor_dist works with more columns than one tgs_dist BLAS block", {
    set.seed(60427)
    x <- matrix(rnorm(20 * 3000), nrow = 20) # tgs_dist on 3000 points uses blocks of 2796 columns
    expect_cor_dist(x)
})
