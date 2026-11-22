test_that("tgs_cor_knn returns correct top-k correlations", {
    set.seed(60427)
    x <- matrix(rnorm(200), nrow = 20, ncol = 10)
    knn <- 3

    # Get knn result
    res <- tgs_cor_knn(x, NULL, knn)
    expect_s3_class(res, "data.frame")
    expect_named(res, c("col1", "col2", "cor", "rank"))

    # Verify: for each column, exactly knn neighbors are returned
    counts <- table(res$col1)
    expect_true(all(counts == knn))

    # Verify against full correlation + manual knn selection
    full_cor <- tgs_cor(x, tidy = TRUE)
    for (col in unique(res$col1)) {
        knn_res <- res[res$col1 == col, ]
        # Get top-k from full correlation for this column
        col_cors <- full_cor[full_cor$col1 == col | full_cor$col2 == col, ]
        # Normalize to always have our target in col1
        col_cors_norm <- data.frame(
            other = ifelse(col_cors$col1 == col, col_cors$col2, col_cors$col1),
            cor = col_cors$cor
        )
        col_cors_sorted <- col_cors_norm[order(-col_cors_norm$cor), ]
        top_k <- head(col_cors_sorted, knn)

        # The correlations should match (order may differ within same rank)
        expect_equal(sort(knn_res$cor, decreasing = TRUE), sort(top_k$cor, decreasing = TRUE), tolerance = 1e-10)
    }
})

test_that("tgs_cor_knn ranks are correct", {
    set.seed(60427)
    x <- matrix(rnorm(200), nrow = 20, ncol = 10)
    knn <- 3

    res <- tgs_cor_knn(x, NULL, knn)

    # Ranks within each column should be 1:knn
    for (col in unique(res$col1)) {
        ranks <- sort(res$rank[res$col1 == col])
        expect_equal(ranks, 1:knn)
    }

    # Higher correlation should have lower rank (rank 1 = highest)
    for (col in unique(res$col1)) {
        col_res <- res[res$col1 == col, ]
        col_res <- col_res[order(col_res$rank), ]
        expect_true(all(diff(col_res$cor) <= 0))
    }
})

test_that("tgs_cor_knn works with cross-correlation", {
    set.seed(60427)
    x <- matrix(rnorm(200), nrow = 20, ncol = 10)
    y <- matrix(rnorm(100), nrow = 20, ncol = 5)
    knn <- 2

    res <- tgs_cor_knn(x, y, knn)
    expect_s3_class(res, "data.frame")
    expect_named(res, c("col1", "col2", "cor", "rank"))

    # col1 should reference columns in x, col2 columns in y
    expect_true(all(res$col1 %in% seq_len(ncol(x))))
    expect_true(all(res$col2 %in% seq_len(ncol(y))))

    # Each column in x should have knn neighbors
    counts <- table(res$col1)
    expect_true(all(counts == knn))
})

test_that("tgs_cor_knn Spearman mode works", {
    set.seed(60427)
    x <- matrix(rnorm(200), nrow = 20, ncol = 10)
    knn <- 3

    res_pearson <- tgs_cor_knn(x, NULL, knn, spearman = FALSE)
    res_spearman <- tgs_cor_knn(x, NULL, knn, spearman = TRUE)

    # Results should differ (different correlation methods)
    expect_false(identical(res_pearson$cor, res_spearman$cor))

    # Both should be valid data frames
    expect_s3_class(res_spearman, "data.frame")
    expect_named(res_spearman, c("col1", "col2", "cor", "rank"))
})

test_that("tgs_cor_knn pairwise.complete.obs works with NAs", {
    set.seed(60427)
    x <- matrix(rnorm(200), nrow = 20, ncol = 10)
    x[1, 1] <- NA
    x[5, 3] <- NA
    x[10, 7] <- NA
    knn <- 3

    # Without pairwise.complete.obs, NAs propagate
    # (columns with any NA would have NA correlations)
    res_pco <- tgs_cor_knn(x, NULL, knn, pairwise.complete.obs = TRUE)
    expect_s3_class(res_pco, "data.frame")
    expect_true(all(!is.na(res_pco$cor)))
})

test_that("tgs_cor_knn threshold filters results", {
    set.seed(60427)
    x <- matrix(rnorm(200), nrow = 20, ncol = 10)
    knn <- 5
    threshold <- 0.3

    res <- tgs_cor_knn(x, NULL, knn, threshold = threshold)
    # All returned correlations should be above threshold
    expect_true(all(abs(res$cor) >= threshold))
})

test_that("tgs_cor_knn knn=1 returns single neighbor per column", {
    set.seed(60427)
    x <- matrix(rnorm(200), nrow = 20, ncol = 10)

    res <- tgs_cor_knn(x, NULL, knn = 1)
    counts <- table(res$col1)
    expect_true(all(counts == 1))
    expect_true(all(res$rank == 1))
})

test_that("tgs_cor_knn handles knn larger than available columns", {
    set.seed(60427)
    x <- matrix(rnorm(60), nrow = 20, ncol = 3)

    # knn = 2 (max possible for 3 columns: each column can have at most 2 neighbors)
    res <- tgs_cor_knn(x, NULL, knn = 2)
    expect_s3_class(res, "data.frame")
    counts <- table(res$col1)
    expect_true(all(counts == 2))
})

test_that("tgs_cor_knn works with single process", {
    set.seed(60427)
    x <- matrix(rnorm(200), nrow = 20, ncol = 10)
    knn <- 3

    old_opt <- getOption("tgs_max.processes")
    options(tgs_max.processes = 1)
    on.exit(options(tgs_max.processes = old_opt))

    res <- tgs_cor_knn(x, NULL, knn)
    expect_s3_class(res, "data.frame")
    expect_equal(nrow(res), ncol(x) * knn)
})

test_that("tgs_cor_knn works with multiple processes", {
    set.seed(60427)
    x <- matrix(rnorm(200), nrow = 20, ncol = 10)
    knn <- 3

    old_opt <- getOption("tgs_max.processes")
    options(tgs_max.processes = 4)
    on.exit(options(tgs_max.processes = old_opt))

    res <- tgs_cor_knn(x, NULL, knn)
    expect_s3_class(res, "data.frame")
    expect_equal(nrow(res), ncol(x) * knn)
})

test_that("tgs_cor_knn single vs multi process gives same results", {
    set.seed(60427)
    x <- matrix(rnorm(200), nrow = 20, ncol = 10)
    knn <- 3

    old_opt <- getOption("tgs_max.processes")

    options(tgs_max.processes = 1)
    res_single <- tgs_cor_knn(x, NULL, knn)

    options(tgs_max.processes = 4)
    res_multi <- tgs_cor_knn(x, NULL, knn)

    options(tgs_max.processes = old_opt)

    # Sort both for comparison
    res_single <- res_single[order(res_single$col1, res_single$rank), ]
    res_multi <- res_multi[order(res_multi$col1, res_multi$rank), ]
    rownames(res_single) <- NULL
    rownames(res_multi) <- NULL

    expect_equal(res_single, res_multi)
})

test_that("tgs_cor_knn rejects invalid inputs", {
    set.seed(60427)
    x <- matrix(rnorm(100), nrow = 10)

    # Non-numeric matrix
    expect_error(tgs_cor_knn(matrix(letters[1:20], nrow = 4), NULL, 2))

    # Invalid knn
    expect_error(tgs_cor_knn(x, NULL, -1))
    expect_error(tgs_cor_knn(x, NULL, 0))
    expect_error(tgs_cor_knn(x, NULL, 1.5))
})

test_that("tgs_cor_knn with integer matrix works", {
    set.seed(60427)
    x <- matrix(sample(1:100, 200, replace = TRUE), nrow = 20, ncol = 10)
    knn <- 3

    res <- tgs_cor_knn(x, NULL, knn)
    expect_s3_class(res, "data.frame")
    expect_equal(nrow(res), ncol(x) * knn)
})

test_that("tgs_cor_knn with named columns preserves names", {
    set.seed(60427)
    x <- matrix(rnorm(200), nrow = 20, ncol = 10)
    colnames(x) <- paste0("gene", 1:10)
    knn <- 3

    res <- tgs_cor_knn(x, NULL, knn)
    expect_s3_class(res, "data.frame")
    # col1 and col2 should be factor or character with column names
    expect_true(all(as.character(res$col1) %in% colnames(x)))
    expect_true(all(as.character(res$col2) %in% colnames(x)))
})

test_that("tgs_cor_knn does not depend on the number of processes", {
    set.seed(60427)
    x <- matrix(rnorm(10 * 120), nrow = 10)
    y <- matrix(rnorm(10 * 120), nrow = 10)
    old <- options(tgs_max.processes = 1)
    on.exit(options(old))
    expected <- tgs_cor_knn(x, NULL, 5)
    expected_cross <- tgs_cor_knn(x, y, 5)
    for (p in c(3, 7, 11, 13, 50)) {
        options(tgs_max.processes = p)
        expect_identical(tgs_cor_knn(x, NULL, 5), expected)
        expect_identical(tgs_cor_knn(x, y, 5), expected_cross)
    }
})

# reference: the highest correlations of each column of x with the columns of y, ties broken by the column number
cor_knn_ref <- function(x, y, knn, spearman = FALSE, threshold = 0) {
    # tgstat treats infinite values as missing
    x[!is.finite(x)] <- NA
    if (!is.null(y)) {
        y[!is.finite(y)] <- NA
    }
    cc <- suppressWarnings(stats::cor(x, if (is.null(y)) x else y, method = if (spearman) "spearman" else "pearson"))
    res <- lapply(seq_len(nrow(cc)), function(i) {
        v <- cc[i, ]
        ok <- !is.na(v) & abs(v) >= abs(threshold)
        if (is.null(y)) {
            ok[i] <- FALSE
        }
        cols <- which(ok)
        cols <- utils::head(cols[order(-v[cols], cols)], knn)
        data.frame(col1 = rep(i, length(cols)), col2 = cols, cor = unname(v[cols]), rank = seq_along(cols))
    })
    do.call(rbind, res)
}

# Spearman correlations of different pairs are often equal up to rounding, and then the neighbour that is chosen
# depends on the rounding: ties = TRUE does not compare the neighbours
expect_cor_knn_equal <- function(res, ref, ties = FALSE) {
    expect_equal(as.integer(res$col1), ref$col1)
    if (!ties) {
        expect_equal(as.integer(res$col2), ref$col2)
    }
    expect_equal(res$cor, ref$cor, tolerance = 1e-10)
    expect_equal(res$rank, ref$rank)
}

with_blas <- function(use_blas, expr) {
    old <- options(tgs_use.blas = use_blas)
    on.exit(options(old))
    expr
}

test_that("tgs_cor_knn with and without BLAS matches the reference", {
    set.seed(60427)
    x <- matrix(rnorm(40 * 60), 40, 60)
    y <- matrix(rnorm(40 * 25), 40, 25)
    y[, 1:10] <- y[, 1:10] + x[, 1:10]

    for (use_blas in c(TRUE, FALSE)) {
        for (spearman in c(FALSE, TRUE)) {
            expect_cor_knn_equal(
                with_blas(use_blas, tgs_cor_knn(x, NULL, 7, spearman = spearman)),
                cor_knn_ref(x, NULL, 7, spearman = spearman),
                ties = spearman
            )
            expect_cor_knn_equal(
                with_blas(use_blas, tgs_cor_knn(x, y, 7, spearman = spearman)),
                cor_knn_ref(x, y, 7, spearman = spearman),
                ties = spearman
            )
            expect_cor_knn_equal(
                with_blas(use_blas, tgs_cor_knn(x, NULL, 1000, spearman = spearman, threshold = 0.21)),
                cor_knn_ref(x, NULL, 1000, spearman = spearman, threshold = 0.21),
                ties = spearman
            )
            expect_cor_knn_equal(
                with_blas(use_blas, tgs_cor_knn(x, y, 1000, spearman = spearman, threshold = -0.21)),
                cor_knn_ref(x, y, 1000, spearman = spearman, threshold = 0.21),
                ties = spearman
            )
        }
    }
})

test_that("tgs_cor_knn skips constant columns and columns with missing values", {
    set.seed(60427)
    x <- matrix(rnorm(30 * 20), 30, 20)
    x[, 4] <- 3
    x[5, 9] <- NA
    x[7, 15] <- Inf
    storage.mode(x) <- "double"
    y <- x[, 20:1]

    for (use_blas in c(TRUE, FALSE)) {
        res <- with_blas(use_blas, tgs_cor_knn(x, NULL, 5))
        expect_cor_knn_equal(res, cor_knn_ref(x, NULL, 5))
        expect_false(any(c(4, 9, 15) %in% res$col1))
        expect_false(any(c(4, 9, 15) %in% res$col2))
        expect_true(all(table(res$col1) == 5))

        res <- with_blas(use_blas, tgs_cor_knn(x, y, 5, spearman = TRUE))
        expect_cor_knn_equal(res, cor_knn_ref(x, y, 5, spearman = TRUE), ties = TRUE)
    }
})

test_that("tgs_cor_knn with pairwise.complete.obs and missing values is the same with and without BLAS", {
    set.seed(60427)
    x <- matrix(rnorm(30 * 20), 30, 20)
    x[sample(length(x), 40)] <- NA
    res <- lapply(c(TRUE, FALSE), function(use_blas) {
        with_blas(use_blas, tgs_cor_knn(x, NULL, 5, pairwise.complete.obs = TRUE))
    })
    expect_identical(res[[1]], res[[2]])

    # without missing values pairwise.complete.obs changes nothing
    x <- matrix(rnorm(30 * 20), 30, 20)
    expect_cor_knn_equal(with_blas(TRUE, tgs_cor_knn(x, NULL, 5, pairwise.complete.obs = TRUE)), cor_knn_ref(x, NULL, 5))
})

test_that("tgs_cor_knn with BLAS handles integer matrices, column names and several blocks", {
    set.seed(60427)
    x <- matrix(sample(0:20, 15 * 12, replace = TRUE), 15, 12, dimnames = list(NULL, paste0("x", 1:12)))
    y <- matrix(sample(0:20, 15 * 8, replace = TRUE), 15, 8, dimnames = list(NULL, paste0("y", 1:8)))
    res <- with_blas(TRUE, tgs_cor_knn(x, y, 3))
    expect_cor_knn_equal(res, cor_knn_ref(x, y, 3))
    expect_equal(levels(res$col1), colnames(x))
    expect_equal(levels(res$col2), colnames(y))
    expect_identical(with_blas(TRUE, tgs_cor_knn(x, NULL, 3))$col2, with_blas(FALSE, tgs_cor_knn(x, NULL, 3))$col2)

    # 3000 columns need more than one block of 2^23 correlations
    x <- matrix(rnorm(12 * 3000), 12, 3000)
    for (spearman in c(FALSE, TRUE)) {
        expect_cor_knn_equal(
            with_blas(TRUE, tgs_cor_knn(x, NULL, 4, spearman = spearman)),
            cor_knn_ref(x, NULL, 4, spearman = spearman),
            ties = spearman
        )
    }
})

test_that("tgs_cor_knn with BLAS ranks in several processes", {
    set.seed(60427)
    x <- matrix(rnorm(2500 * 500), 2500, 500)
    old <- options(tgs_max.processes = 4)
    on.exit(options(old))
    expect_cor_knn_equal(with_blas(TRUE, tgs_cor_knn(x, NULL, 3, spearman = TRUE)), cor_knn_ref(x, NULL, 3, spearman = TRUE), ties = TRUE)
    expect_cor_knn_equal(with_blas(TRUE, tgs_cor_knn(x, x[, 1:100], 3, spearman = TRUE)), cor_knn_ref(x, x[, 1:100], 3, spearman = TRUE), ties = TRUE)
})

test_that("tgs_cor_knn accepts data frames and y is optional", {
    set.seed(60427)
    m <- matrix(rnorm(40 * 8), 40, 8, dimnames = list(NULL, letters[1:8]))
    df <- as.data.frame(m)
    for (use_blas in c(TRUE, FALSE)) {
        ref <- with_blas(use_blas, tgs_cor_knn(m, NULL, 3))
        expect_identical(with_blas(use_blas, tgs_cor_knn(m, knn = 3)), ref)
        expect_identical(with_blas(use_blas, tgs_cor_knn(df, knn = 3)), ref)
        expect_identical(with_blas(use_blas, tgs_cor_knn(df, df[, 1:4], 2)), with_blas(use_blas, tgs_cor_knn(m, m[, 1:4], 2)))
    }
})

test_that("tgs_cor_knn skips constant columns with pairwise.complete.obs and missing values", {
    set.seed(60427)
    x <- matrix(rnorm(30 * 10), 30, 10)
    x[, 4] <- 1 / 3
    x[c(2, 8), 4] <- NA
    x[5, 7] <- NA
    for (use_blas in c(TRUE, FALSE)) {
        res <- with_blas(use_blas, tgs_cor_knn(x, NULL, 4, pairwise.complete.obs = TRUE))
        expect_false(4 %in% res$col1)
        expect_false(4 %in% res$col2)
        expect_true(all(table(res$col1) == 4))
        expect_true(all(abs(res$cor) <= 1 + 1e-12))
        res <- with_blas(use_blas, tgs_cor_knn(x, x[, 10:1], 4, pairwise.complete.obs = TRUE))
        expect_false(4 %in% res$col1)
        expect_false(7 %in% res$col2)
    }
})

test_that("tgs_cor_knn skips a column that is constant over the rows it shares with another column", {
    set.seed(60427)
    for (const in c(3, 0.1, 1 / 3, 1e6 + 0.1, 0)) {
        for (num_shared in c(9, 20, 30)) {
            x <- matrix(rnorm(50 * 6), 50, 6)
            # column 2 is constant over the rows in which column 5 has values
            x[1:num_shared, 2] <- const
            x[(num_shared + 1):50, 5] <- NA
            y <- x[, 6:1]

            for (use_blas in c(TRUE, FALSE)) {
                for (spearman in c(FALSE, TRUE)) {
                    # all the other columns are returned: only the pair of columns 2 and 5 has no correlation
                    res <- with_blas(use_blas, tgs_cor_knn(x, NULL, 5, pairwise.complete.obs = TRUE, spearman = spearman))
                    expect_equal(as.vector(table(res$col1)), c(5, 4, 5, 5, 4, 5))
                    expect_false(any(res$col1 == 2 & res$col2 == 5 | res$col1 == 5 & res$col2 == 2))
                    expect_true(all(abs(res$cor) <= 1 + 1e-12))

                    res <- with_blas(use_blas, tgs_cor_knn(x, y, 6, pairwise.complete.obs = TRUE, spearman = spearman))
                    expect_equal(as.vector(table(res$col1)), c(6, 5, 6, 6, 5, 6))
                    expect_false(any(res$col1 == 2 & res$col2 == 2 | res$col1 == 5 & res$col2 == 5))
                }
            }
        }
    }
})

test_that("tgs_cor_knn returns neighbours with a correlation of exactly zero", {
    a <- c(1, -1, 1, -1, 0)
    b <- c(1, 1, -1, -1, 0)
    d <- c(1, -1, -1, 1, 0)
    x <- cbind(a, b, d, -a)
    dimnames(x) <- NULL
    x_na <- x
    x_na[5, 2] <- NA

    for (use_blas in c(TRUE, FALSE)) {
        for (spearman in c(FALSE, TRUE)) {
            for (pairwise in c(FALSE, TRUE)) {
                m <- if (pairwise) x_na else x
                res <- with_blas(use_blas, tgs_cor_knn(m, NULL, 3, pairwise.complete.obs = pairwise, spearman = spearman))
                expect_equal(as.vector(table(res$col1)), rep(3, 4))
                expect_equal(res$col2[res$col1 == 1], c(2, 3, 4))
                expect_equal(res$cor[res$col1 == 1], c(0, 0, -1))
                expect_equal(res$cor[res$col1 == 2], c(0, 0, 0))

                res <- with_blas(use_blas, tgs_cor_knn(m, m[, 4:1], 4, pairwise.complete.obs = pairwise, spearman = spearman))
                expect_equal(as.vector(table(res$col1)), rep(4, 4))
                expect_equal(res$col2[res$col1 == 1], c(4, 2, 3, 1))
                expect_equal(res$cor[res$col1 == 1], c(1, 0, 0, -1))
            }
        }
    }

    # a column that has no correlations still has no neighbours
    x[, 3] <- 2
    for (use_blas in c(TRUE, FALSE)) {
        res <- with_blas(use_blas, tgs_cor_knn(x, NULL, 3))
        expect_equal(as.vector(table(factor(res$col1, 1:4))), c(2, 2, 0, 2))
        expect_false(3 %in% res$col2)
        res <- with_blas(use_blas, tgs_cor_knn(x, x[, 4:1], 4))
        expect_equal(as.vector(table(factor(res$col1, 1:4))), c(3, 3, 0, 3))
        expect_false(2 %in% res$col2)
    }
})
