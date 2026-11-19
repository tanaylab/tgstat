test_that("tgs_cor works correctly", {
    # Create test matrix
    set.seed(60427)
    x <- matrix(rnorm(100), nrow = 10)

    # Compare to base R cor
    expect_equal(tgs_cor(x), cor(x))

    # Test with NAs
    x[1, 1] <- NA
    expect_equal(
        tgs_cor(x, pairwise.complete.obs = TRUE),
        cor(x, use = "pairwise.complete.obs")
    )

    # Test Spearman correlation
    expect_equal(
        tgs_cor(x, spearman = TRUE),
        cor(x, method = "spearman")
    )

    # Test tidy output
    tidy_res <- tgs_cor(x, tidy = TRUE)
    expect_s3_class(tidy_res, "data.frame")
    expect_named(tidy_res, c("col1", "col2", "cor"))

    # Test threshold
    thresh_res <- tgs_cor(x, tidy = TRUE, threshold = 0.5)
    expect_true(all(abs(thresh_res$cor) >= 0.5))

    # Test cross-correlation
    y <- matrix(rnorm(50), nrow = 10)
    expect_equal(tgs_cor(x, y), cor(x, y), ignore_attr = TRUE)
})

test_that("tgs_cor computes all the columns with process counts that are not powers of 2", {
    set.seed(60427)
    x <- matrix(rnorm(10 * 120), nrow = 10)
    y <- matrix(rnorm(10 * 50), nrow = 10)
    old <- options(tgs_use.blas = FALSE, tgs_max.processes = 1)
    on.exit(options(old))
    for (p in c(3, 7, 11, 13, 50)) {
        options(tgs_max.processes = p)
        expect_equal(tgs_cor(x), cor(x))
        expect_equal(tgs_cor(y, x), cor(y, x), ignore_attr = TRUE)
    }
})

test_that("tgs_cor returns NA for constant columns", {
    set.seed(60427)
    x <- matrix(rnorm(30 * 5), 30, 5)
    y <- matrix(rnorm(30 * 4), 30, 4)
    old <- options(tgs_use.blas = FALSE)
    on.exit(options(old))

    # the variance of these constants, computed as mean(x^2) - mean(x)^2, is not exactly zero
    for (const in c(3, 0.1, 1 / 3, 1e6 + 0.1, 0)) {
        x[, 2] <- const
        y[, 3] <- const
        ref <- suppressWarnings(stats::cor(x))
        diag(ref) <- 1
        cross_ref <- suppressWarnings(stats::cor(x, y))
        for (use_blas in c(TRUE, FALSE)) {
            options(tgs_use.blas = use_blas)
            for (pairwise in c(TRUE, FALSE)) {
                expect_equal(tgs_cor(x, pairwise.complete.obs = pairwise), ref)
                expect_equal(tgs_cor(x, y, pairwise.complete.obs = pairwise), cross_ref, ignore_attr = TRUE)
            }
            expect_true(all(is.na(tgs_cor(x, spearman = TRUE)[2, -2])))
            expect_true(all(is.na(tgs_cor(x, y, spearman = TRUE)[2, ])))
            expect_true(all(is.na(tgs_cor(x, y, spearman = TRUE)[, 3])))
        }
    }
})
