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
        ref[2, ] <- ref[, 2] <- NA
        diag(ref) <- 1
        cross_ref <- suppressWarnings(stats::cor(x, y))
        cross_ref[2, ] <- cross_ref[, 3] <- NA
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

test_that("tgs_cor with BLAS matches stats::cor and tgs_cor without BLAS", {
    set.seed(60427)
    x <- matrix(rnorm(50 * 30), 50, 30, dimnames = list(NULL, paste0("x", 1:30)))
    x[, 5] <- x[, 4] * 2 + 1
    x[, 7] <- 2
    x[3, 11] <- NA
    x[8, 12] <- Inf
    y <- matrix(rnorm(50 * 20), 50, 20, dimnames = list(NULL, paste0("y", 1:20)))
    y[, 2] <- -1
    y[9, 6] <- NA

    # tgstat treats infinite values as missing, and the diagonal is always 1
    ref <- function(x, y = NULL, method = "pearson") {
        x[!is.finite(x)] <- NA
        res <- suppressWarnings(stats::cor(x, y, method = method))
        if (is.null(y)) {
            diag(res) <- 1
        }
        res
    }

    old <- options(tgs_use.blas = TRUE)
    on.exit(options(old))

    for (spearman in c(FALSE, TRUE)) {
        method <- if (spearman) "spearman" else "pearson"
        options(tgs_use.blas = TRUE)
        res <- tgs_cor(x, spearman = spearman)
        cross <- tgs_cor(x, y, spearman = spearman)
        tidy <- tgs_cor(x, spearman = spearman, tidy = TRUE, threshold = 0.3)
        cross_tidy <- tgs_cor(x, y, spearman = spearman, tidy = TRUE, threshold = 0.3)
        expect_equal(res, ref(x, method = method))
        expect_equal(cross, ref(x, y, method = method))
        expect_true(isSymmetric(res))

        options(tgs_use.blas = FALSE)
        expect_equal(res, tgs_cor(x, spearman = spearman))
        expect_equal(cross, tgs_cor(x, y, spearman = spearman))
        expect_equal(tidy, tgs_cor(x, spearman = spearman, tidy = TRUE, threshold = 0.3))
        expect_equal(cross_tidy, tgs_cor(x, y, spearman = spearman, tidy = TRUE, threshold = 0.3))
        expect_equal(nrow(tidy), sum(abs(res[upper.tri(res)]) >= 0.3, na.rm = TRUE))
    }

    # integer matrices
    xi <- matrix(sample(0:9, 40 * 12, replace = TRUE), 40, 12)
    options(tgs_use.blas = TRUE)
    expect_equal(tgs_cor(xi), ref(xi))
    expect_equal(tgs_cor(xi, xi[, 1:5], spearman = TRUE), ref(xi, xi[, 1:5], method = "spearman"), ignore_attr = TRUE)
})

test_that("tgs_cor with BLAS does not depend on the number of processes", {
    set.seed(60427)
    # enough values to rank in several processes, and to center and scale in several processes
    xs <- matrix(rnorm(1500 * 1200), 1500, 1200)
    xs[, 17] <- 4
    xs[20, 100] <- NA
    xp <- matrix(rnorm(3000 * 3000), 3000, 3000)
    old <- options(tgs_use.blas = TRUE, tgs_max.processes = 1)
    on.exit(options(old))

    ref_s <- tgs_cor(xs, spearman = TRUE)
    ref_cross <- tgs_cor(xs, xs[, 1:300], spearman = TRUE)
    ref_p <- tgs_cor(xp)
    for (np in c(3, 8)) {
        options(tgs_max.processes = np)
        expect_identical(tgs_cor(xs, spearman = TRUE), ref_s)
        expect_identical(tgs_cor(xs, xs[, 1:300], spearman = TRUE), ref_cross)
        expect_identical(tgs_cor(xp), ref_p)
    }
    expect_equal(ref_s[1:50, 1:50], suppressWarnings(stats::cor(xs[, 1:50], method = "spearman")), ignore_attr = TRUE)
    expect_true(all(is.na(ref_s[17, -17])) && all(is.na(ref_s[100, -100])))
})

test_that("tgs_cor accepts data frames and vectors", {
    set.seed(60427)
    m <- matrix(rnorm(40 * 6), 40, 6, dimnames = list(NULL, letters[1:6]))
    v <- rnorm(40)
    df <- as.data.frame(m)
    old <- options(tgs_use.blas = TRUE)
    on.exit(options(old))

    for (use_blas in c(TRUE, FALSE)) {
        options(tgs_use.blas = use_blas)
        expect_identical(tgs_cor(df), tgs_cor(m))
        expect_identical(tgs_cor(df, df[, 1:3]), tgs_cor(m, m[, 1:3]))
        expect_equal(tgs_cor(m, v)[, 1], stats::cor(m, v)[, 1])
        expect_equal(tgs_cor(v, m)[1, ], stats::cor(v, m)[1, ])
        expect_equal(tgs_cor(v, m[, 2], spearman = TRUE)[1, 1], stats::cor(v, m[, 2], method = "spearman"))
        expect_equal(dim(tgs_cor(df, v)), c(6L, 1L))
        expect_error(tgs_cor(v))
    }
})

test_that("tgs_cor returns NA for constant columns with pairwise.complete.obs and missing values", {
    set.seed(60427)
    x <- matrix(rnorm(30 * 6), 30, 6)
    y <- matrix(rnorm(30 * 4), 30, 4)
    old <- options(tgs_use.blas = FALSE)
    on.exit(options(old))

    for (const in c(3, 0.1, 1 / 3)) {
        x[, 2] <- const
        x[c(4, 9), 2] <- NA
        x[5, 4] <- NA
        y[, 3] <- const
        y[7, 1] <- NA
        # stats::cor does not return NA for the constant columns on platforms without extended precision
        ref <- suppressWarnings(stats::cor(x, use = "pairwise.complete.obs"))
        ref[2, ] <- ref[, 2] <- NA
        diag(ref) <- 1
        cross_ref <- suppressWarnings(stats::cor(x, y, use = "pairwise.complete.obs"))
        cross_ref[2, ] <- cross_ref[, 3] <- NA
        for (use_blas in c(TRUE, FALSE)) {
            options(tgs_use.blas = use_blas)
            expect_equal(tgs_cor(x, pairwise.complete.obs = TRUE), ref)
            expect_equal(tgs_cor(x, y, pairwise.complete.obs = TRUE), cross_ref, ignore_attr = TRUE)
        }
        options(tgs_use.blas = FALSE)
        expect_true(all(is.na(tgs_cor(x, pairwise.complete.obs = TRUE, spearman = TRUE)[2, -2])))
        expect_true(all(is.na(tgs_cor(x, y, pairwise.complete.obs = TRUE, spearman = TRUE)[, 3])))
    }
})

test_that("a NULL y is not converted to a matrix", {
    # before R 4.4.0 is.atomic(NULL) is TRUE
    expect_null(.tgs_as_matrix(NULL, vector_ok = TRUE))
    expect_null(.tgs_as_matrix(NULL))
    expect_equal(dim(.tgs_as_matrix(1:3, vector_ok = TRUE)), c(3L, 1L))
})

test_that("tgs_cor returns NA for a column that is constant over the rows it shares with another column", {
    set.seed(60427)
    old <- options(tgs_use.blas = FALSE)
    on.exit(options(old))

    for (const in c(3, 0.1, 1 / 3, 1e6 + 0.1, 0)) {
        for (num_shared in c(7, 20, 30)) {
            x <- matrix(rnorm(50 * 6), 50, 6)
            y <- matrix(rnorm(50 * 4), 50, 4)
            # column 2 of x is constant over the rows in which column 5 of x and column 3 of y have values
            x[1:num_shared, 2] <- const
            x[(num_shared + 1):50, 5] <- NA
            y[(num_shared + 1):50, 3] <- NA
            x[4, 1] <- NA

            # stats::cor does not return NA for constant values on platforms without extended precision
            ref <- suppressWarnings(stats::cor(x, use = "pairwise.complete.obs"))
            ref[2, 5] <- ref[5, 2] <- NA
            cross_ref <- suppressWarnings(stats::cor(x, y, use = "pairwise.complete.obs"))
            cross_ref[2, 3] <- NA
            spearman_ref <- suppressWarnings(stats::cor(x, use = "pairwise.complete.obs", method = "spearman"))
            spearman_ref[2, 5] <- spearman_ref[5, 2] <- NA
            spearman_cross_ref <- suppressWarnings(stats::cor(x, y, use = "pairwise.complete.obs", method = "spearman"))
            spearman_cross_ref[2, 3] <- NA

            for (use_blas in c(TRUE, FALSE)) {
                options(tgs_use.blas = use_blas)
                expect_equal(tgs_cor(x, pairwise.complete.obs = TRUE), ref)
                expect_equal(tgs_cor(x, y, pairwise.complete.obs = TRUE), cross_ref, ignore_attr = TRUE)
                expect_equal(tgs_cor(y, x, pairwise.complete.obs = TRUE), t(cross_ref), ignore_attr = TRUE)
            }
            options(tgs_use.blas = FALSE)
            expect_equal(tgs_cor(x, pairwise.complete.obs = TRUE, spearman = TRUE), spearman_ref)
            expect_equal(tgs_cor(x, y, pairwise.complete.obs = TRUE, spearman = TRUE), spearman_cross_ref, ignore_attr = TRUE)
        }
    }
})
