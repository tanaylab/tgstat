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
