# Helper: compute chi-squared statistic and p-value for a single row using
# R's stats::chisq.test(). Builds the full 2x2 contingency table from
# row values (a, c_val) and column sums.
chi2_reference <- function(a, c_val, col_sum1, col_sum2, yates = TRUE) {
    b <- col_sum1 - a
    d <- col_sum2 - c_val
    ct <- matrix(c(a, c_val, b, d), nrow = 2, byrow = FALSE)
    # If any marginal is zero, chisq.test will warn/error; handle gracefully
    tryCatch(
        {
            res <- stats::chisq.test(ct, correct = yates)
            c(chi2 = unname(res$statistic), pval = unname(res$p.value))
        },
        error = function(e) c(chi2 = 0, pval = 1),
        warning = function(w) {
            res <- suppressWarnings(stats::chisq.test(ct, correct = yates))
            c(chi2 = unname(res$statistic), pval = unname(res$p.value))
        }
    )
}

# Fully manual reference implementation (no chisq.test) to double-check
chi2_manual <- function(a, c_val, col_sum1, col_sum2, yates = TRUE) {
    a <- as.double(a)
    c_val <- as.double(c_val)
    col_sum1 <- as.double(col_sum1)
    col_sum2 <- as.double(col_sum2)
    b <- col_sum1 - a
    d <- col_sum2 - c_val
    N <- a + b + c_val + d
    r1 <- a + b
    r2 <- c_val + d
    c1 <- a + c_val
    c2 <- b + d
    denom <- r1 * r2 * c1 * c2
    if (denom == 0) {
        return(c(chi2 = 0, pval = 1))
    }
    ad_bc <- abs(a * d - b * c_val)
    if (yates) {
        numerator <- max(ad_bc - N / 2, 0)
    } else {
        numerator <- ad_bc
    }
    chi2_stat <- numerator^2 * N / denom
    pval <- stats::pchisq(chi2_stat, df = 1, lower.tail = FALSE)
    c(chi2 = chi2_stat, pval = pval)
}

# Compute expected results for an entire matrix using the manual reference
chi2_expected <- function(x, yates = TRUE) {
    col_sum1 <- sum(x[, 1])
    col_sum2 <- sum(x[, 2])
    result <- t(sapply(seq_len(nrow(x)), function(i) {
        chi2_manual(x[i, 1], x[i, 2], col_sum1, col_sum2, yates)
    }))
    colnames(result) <- c("chi2", "pval")
    if (!is.null(rownames(x))) {
        rownames(result) <- rownames(x)
    }
    result
}

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

test_that("tgs_chi2 with Yates correction matches chisq.test", {
    options(tgs_max.processes = 1)
    set.seed(42)
    x <- matrix(sample(0:50, 20, replace = TRUE), ncol = 2)
    rownames(x) <- paste0("gene", seq_len(nrow(x)))

    res <- tgs_chi2(x, yates = TRUE)

    col_sum1 <- sum(x[, 1])
    col_sum2 <- sum(x[, 2])
    for (i in seq_len(nrow(x))) {
        ref <- chi2_reference(x[i, 1], x[i, 2], col_sum1, col_sum2, yates = TRUE)
        expect_equal(unname(res[i, "chi2"]), unname(ref["chi2"]),
            tolerance = 1e-10,
            info = paste("Row", i, "chi2 with Yates")
        )
        expect_equal(unname(res[i, "pval"]), unname(ref["pval"]),
            tolerance = 1e-10,
            info = paste("Row", i, "pval with Yates")
        )
    }
})

test_that("tgs_chi2 without Yates correction matches chisq.test", {
    options(tgs_max.processes = 1)
    set.seed(42)
    x <- matrix(sample(0:50, 20, replace = TRUE), ncol = 2)
    rownames(x) <- paste0("gene", seq_len(nrow(x)))

    res <- tgs_chi2(x, yates = FALSE)

    col_sum1 <- sum(x[, 1])
    col_sum2 <- sum(x[, 2])
    for (i in seq_len(nrow(x))) {
        ref <- chi2_reference(x[i, 1], x[i, 2], col_sum1, col_sum2, yates = FALSE)
        expect_equal(unname(res[i, "chi2"]), unname(ref["chi2"]),
            tolerance = 1e-10,
            info = paste("Row", i, "chi2 without Yates")
        )
        expect_equal(unname(res[i, "pval"]), unname(ref["pval"]),
            tolerance = 1e-10,
            info = paste("Row", i, "pval without Yates")
        )
    }
})

test_that("tgs_chi2 works with dense integer matrix", {
    options(tgs_max.processes = 1)
    set.seed(123)
    x <- matrix(as.integer(sample(0:100, 40, replace = TRUE)), ncol = 2)
    expect_true(is.integer(x))

    res <- tgs_chi2(x, yates = TRUE)
    expected <- chi2_expected(x, yates = TRUE)
    expect_equal(res, expected, tolerance = 1e-10)
})

test_that("tgs_chi2 works with dense numeric (double) matrix", {
    options(tgs_max.processes = 1)
    set.seed(123)
    x <- matrix(as.double(sample(0:100, 40, replace = TRUE)), ncol = 2)
    expect_true(is.double(x))

    res <- tgs_chi2(x, yates = TRUE)
    expected <- chi2_expected(x, yates = TRUE)
    expect_equal(res, expected, tolerance = 1e-10)
})

test_that("tgs_chi2 works with dgCMatrix sparse matrix", {
    options(tgs_max.processes = 1)
    skip_if_not_installed("Matrix")
    set.seed(99)
    # Create a matrix with some zeros to make sparsity meaningful
    x_dense <- matrix(sample(c(0, 0, 0, 1:20), 40, replace = TRUE), ncol = 2)
    rownames(x_dense) <- paste0("row", seq_len(nrow(x_dense)))
    x_sparse <- Matrix::Matrix(x_dense, sparse = TRUE)
    expect_s4_class(x_sparse, "dgCMatrix")

    res_sparse <- tgs_chi2(x_sparse, yates = TRUE)
    expected <- chi2_expected(x_dense, yates = TRUE)
    expect_equal(res_sparse, expected, tolerance = 1e-10)

    # Also test without Yates
    res_sparse_no_yates <- tgs_chi2(x_sparse, yates = FALSE)
    expected_no_yates <- chi2_expected(x_dense, yates = FALSE)
    expect_equal(res_sparse_no_yates, expected_no_yates, tolerance = 1e-10)
})

test_that("tgs_chi2 edge case: row with zero counts in both columns", {
    options(tgs_max.processes = 1)
    # Row with (0, 0) should give chi2=0, pval=1 since c1=0 => denom=0
    x <- matrix(c(
        10, 5,
        0, 0,
        3, 7
    ), ncol = 2, byrow = TRUE)
    res <- tgs_chi2(x, yates = TRUE)

    # The zero-row should produce chi2=0, pval=1
    expect_equal(unname(res[2, "chi2"]), 0)
    expect_equal(unname(res[2, "pval"]), 1)

    # Other rows should match reference
    expected <- chi2_expected(x, yates = TRUE)
    expect_equal(res, expected, tolerance = 1e-10)
})

test_that("tgs_chi2 edge case: gene has all UMIs in one condition", {
    options(tgs_max.processes = 1)
    # A row where all counts are in column 1 (c_val=0), or all in column 2 (a=0)
    x <- matrix(c(
        10, 0,
        0, 15,
        5, 5,
        20, 0
    ), ncol = 2, byrow = TRUE)
    res <- tgs_chi2(x, yates = TRUE)
    expected <- chi2_expected(x, yates = TRUE)
    expect_equal(res, expected, tolerance = 1e-10)

    # Also check specific property: rows with extreme allocation should have chi2 >= 0
    expect_true(all(res[, "chi2"] >= 0))
    expect_true(all(res[, "pval"] >= 0 & res[, "pval"] <= 1))
})

test_that("tgs_chi2 edge case: single row matrix", {
    options(tgs_max.processes = 1)
    x <- matrix(c(10, 20), ncol = 2)
    rownames(x) <- "only_gene"

    res <- tgs_chi2(x, yates = TRUE)

    # Single row: a=10, c=20, col_sum1=10, col_sum2=20
    # b=0, d=0 => r2=0 => denom=0 => chi2=0, pval=1
    # (because there is only one row, all column sums equal that row)
    expect_equal(nrow(res), 1)
    expect_equal(ncol(res), 2)
    expect_equal(unname(res[1, "chi2"]), 0)
    expect_equal(unname(res[1, "pval"]), 1)
    expect_equal(rownames(res), "only_gene")
})

test_that("tgs_chi2 edge case: very large counts", {
    options(tgs_max.processes = 1)
    # Use large counts to check numerical stability
    x <- matrix(c(
        1e6, 1e6 + 100,
        500000, 500100,
        1e7, 1e7 - 50
    ), ncol = 2, byrow = TRUE)

    res <- tgs_chi2(x, yates = TRUE)
    expected <- chi2_expected(x, yates = TRUE)
    expect_equal(res, expected, tolerance = 1e-8)

    res_no_yates <- tgs_chi2(x, yates = FALSE)
    expected_no_yates <- chi2_expected(x, yates = FALSE)
    expect_equal(res_no_yates, expected_no_yates, tolerance = 1e-8)
})

test_that("tgs_chi2 preserves row names in output", {
    options(tgs_max.processes = 1)
    x <- matrix(c(10, 20, 30, 40, 50, 60), ncol = 2, byrow = TRUE)
    rownames(x) <- c("geneA", "geneB", "geneC")

    res <- tgs_chi2(x, yates = TRUE)
    expect_equal(rownames(res), c("geneA", "geneB", "geneC"))
})

test_that("tgs_chi2 output has correct dimensions and column names", {
    options(tgs_max.processes = 1)
    set.seed(7)
    x <- matrix(sample(0:30, 30, replace = TRUE), ncol = 2)

    res <- tgs_chi2(x, yates = TRUE)

    expect_true(is.matrix(res))
    expect_equal(nrow(res), nrow(x))
    expect_equal(ncol(res), 2)
    expect_equal(colnames(res), c("chi2", "pval"))
})

test_that("tgs_chi2 errors on non-2-column matrix", {
    options(tgs_max.processes = 1)

    # 3 columns
    x3 <- matrix(1:12, ncol = 3)
    expect_error(tgs_chi2(x3))

    # 1 column
    x1 <- matrix(1:5, ncol = 1)
    expect_error(tgs_chi2(x1))
})

test_that("tgs_chi2 dense and sparse give identical results", {
    options(tgs_max.processes = 1)
    skip_if_not_installed("Matrix")
    set.seed(2024)
    x_dense <- matrix(sample(c(0, 0, 1:10), 60, replace = TRUE), ncol = 2)
    rownames(x_dense) <- paste0("g", seq_len(nrow(x_dense)))
    x_sparse <- Matrix::Matrix(x_dense, sparse = TRUE)

    res_dense <- tgs_chi2(x_dense, yates = TRUE)
    res_sparse <- tgs_chi2(x_sparse, yates = TRUE)
    expect_equal(res_dense, res_sparse, tolerance = 1e-12)

    res_dense_noyates <- tgs_chi2(x_dense, yates = FALSE)
    res_sparse_noyates <- tgs_chi2(x_sparse, yates = FALSE)
    expect_equal(res_dense_noyates, res_sparse_noyates, tolerance = 1e-12)
})

test_that("tgs_chi2 Yates correction reduces chi2 compared to no correction", {
    options(tgs_max.processes = 1)
    set.seed(55)
    x <- matrix(sample(1:50, 20, replace = TRUE), ncol = 2)

    res_yates <- tgs_chi2(x, yates = TRUE)
    res_no_yates <- tgs_chi2(x, yates = FALSE)

    # Yates correction should produce chi2 <= the uncorrected version for every row
    expect_true(all(res_yates[, "chi2"] <= res_no_yates[, "chi2"] + 1e-12))
    # And pval >= the uncorrected version (since chi2 is smaller)
    expect_true(all(res_yates[, "pval"] >= res_no_yates[, "pval"] - 1e-12))
})

test_that("tgs_chi2 handles all-zero column (degenerate case)", {
    options(tgs_max.processes = 1)
    # Column 2 is all zeros => colSum2 = 0 => denom = 0 for all rows
    x <- matrix(c(10, 0, 20, 0, 30, 0), ncol = 2, byrow = TRUE)
    res <- tgs_chi2(x, yates = TRUE)

    expect_true(all(res[, "chi2"] == 0))
    expect_true(all(res[, "pval"] == 1))

    # Column 1 is all zeros
    x2 <- matrix(c(0, 10, 0, 20, 0, 30), ncol = 2, byrow = TRUE)
    res2 <- tgs_chi2(x2, yates = TRUE)

    expect_true(all(res2[, "chi2"] == 0))
    expect_true(all(res2[, "pval"] == 1))
})

test_that("tgs_chi2 results have pval in [0, 1] and chi2 >= 0", {
    options(tgs_max.processes = 1)
    set.seed(314)
    x <- matrix(sample(0:100, 100, replace = TRUE), ncol = 2)

    res <- tgs_chi2(x, yates = TRUE)
    expect_true(all(res[, "chi2"] >= 0))
    expect_true(all(res[, "pval"] >= 0 & res[, "pval"] <= 1))

    res2 <- tgs_chi2(x, yates = FALSE)
    expect_true(all(res2[, "chi2"] >= 0))
    expect_true(all(res2[, "pval"] >= 0 & res2[, "pval"] <= 1))
})

test_that("tgs_chi2 with matrix where row sums equal column sums (denom=0)", {
    options(tgs_max.processes = 1)
    # If a row consumes ALL of both column sums, c1 = N => c2 = 0 => denom = 0
    # This happens with a single-row matrix (covered above) or when one row
    # dominates.  Two-row case: row1 = (colSum1, colSum2) => row2 = (0,0)
    x <- matrix(c(100, 200, 0, 0), ncol = 2, byrow = TRUE)
    res <- tgs_chi2(x, yates = TRUE)
    expected <- chi2_expected(x, yates = TRUE)
    expect_equal(res, expected, tolerance = 1e-10)
})

test_that("tgs_chi2 no row names when input has none", {
    options(tgs_max.processes = 1)
    x <- matrix(c(5, 10, 15, 20), ncol = 2)
    # No rownames set
    expect_null(rownames(x))

    res <- tgs_chi2(x, yates = TRUE)
    expect_null(rownames(res))
})

test_that("tgs_chi2 handles NA in dense integer matrix", {
    options(tgs_max.processes = 1)
    x <- matrix(c(10L, NA_integer_, 30L, 20L, 25L, 35L), ncol = 2)
    rownames(x) <- c("gene1", "gene2", "gene3")
    res <- tgs_chi2(x)
    # Row 2 has NA -> should produce NA for chi2 and pval
    expect_true(is.na(res["gene2", "chi2"]))
    expect_true(is.na(res["gene2", "pval"]))
    # Rows 1 and 3 should have valid results
    expect_false(is.na(res["gene1", "chi2"]))
    expect_false(is.na(res["gene3", "chi2"]))
    expect_false(is.na(res["gene1", "pval"]))
    expect_false(is.na(res["gene3", "pval"]))
})

test_that("tgs_chi2 handles NaN in dense numeric matrix", {
    options(tgs_max.processes = 1)
    x <- matrix(c(10, NaN, 30, 20, 25, 35), ncol = 2)
    rownames(x) <- c("gene1", "gene2", "gene3")
    res <- tgs_chi2(x)
    # Row 2 has NaN -> should produce NA for chi2 and pval
    expect_true(is.na(res["gene2", "chi2"]))
    expect_true(is.na(res["gene2", "pval"]))
    # Rows 1 and 3 should have valid results
    expect_false(is.na(res["gene1", "chi2"]))
    expect_false(is.na(res["gene3", "chi2"]))
})

test_that("tgs_chi2 handles NA in second column", {
    options(tgs_max.processes = 1)
    x <- matrix(c(10L, 20L, 30L, NA_integer_, 25L, 35L), ncol = 2)
    res <- tgs_chi2(x)
    # Row 1 has NA in col 2
    expect_true(is.na(res[1, "chi2"]))
    expect_true(is.na(res[1, "pval"]))
    # Rows 2 and 3 are fine
    expect_false(is.na(res[2, "chi2"]))
    expect_false(is.na(res[3, "chi2"]))
})

test_that("tgs_chi2 errors on yates = NA", {
    options(tgs_max.processes = 1)
    x <- matrix(c(10, 20, 30, 40), ncol = 2)
    expect_error(tgs_chi2(x, yates = NA), "must not be NA")
})

test_that("tgs_chi2 handles all-zero sparse matrix", {
    options(tgs_max.processes = 1)
    # A dgCMatrix where every entry is 0 has an empty x slot
    x <- matrix(0L, nrow = 5, ncol = 2)
    rownames(x) <- paste0("gene", 1:5)
    smat <- Matrix::Matrix(x, sparse = TRUE)
    expect_equal(length(smat@x), 0L)  # confirm x slot is empty

    res <- tgs_chi2(smat)
    expect_equal(nrow(res), 5)
    expect_equal(ncol(res), 2)
    expect_equal(colnames(res), c("chi2", "pval"))
    # All rows should have chi2=0, pval=1 (denom is 0)
    expect_true(all(res[, "chi2"] == 0))
    expect_true(all(res[, "pval"] == 1))
    expect_equal(rownames(res), paste0("gene", 1:5))
})

test_that("tgs_chi2 errors on negative counts in dense matrix", {
    options(tgs_max.processes = 1)
    x <- matrix(c(-1, 3, 2, 4), ncol = 2)
    expect_error(tgs_chi2(x), "non-negative")
})

test_that("tgs_chi2 errors on negative counts in integer matrix", {
    options(tgs_max.processes = 1)
    x <- matrix(c(1L, -3L, 2L, 4L), ncol = 2)
    expect_error(tgs_chi2(x), "non-negative")
})
