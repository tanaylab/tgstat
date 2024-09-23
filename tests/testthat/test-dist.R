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
