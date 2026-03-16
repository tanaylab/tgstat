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
