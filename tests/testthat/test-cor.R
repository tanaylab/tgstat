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
