# reference implementations: ks.test and wilcox.test on each column, NA when one of the samples is empty
ks_ref <- function(x, labels, alternative = "two.sided") {
    x <- as.matrix(x)
    res <- t(apply(x, 2, function(v) {
        v1 <- v[labels %in% TRUE]
        v2 <- v[labels %in% FALSE]
        if (all(is.na(v1)) || all(is.na(v2))) {
            return(c(NA_real_, NA_real_))
        }
        r <- suppressWarnings(stats::ks.test(v1, v2, alternative = alternative, exact = FALSE))
        c(unname(r$statistic), r$p.value)
    }))
    res[is.nan(res)] <- NA_real_
    colnames(res) <- c("D", "pval")
    res
}

wilcox_ref <- function(x, labels, alternative = "two.sided") {
    x <- as.matrix(x)
    res <- t(apply(x, 2, function(v) {
        v1 <- v[labels %in% TRUE]
        v2 <- v[labels %in% FALSE]
        if (all(is.na(v1)) || all(is.na(v2))) {
            return(c(NA_real_, NA_real_, NA_real_))
        }
        r <- stats::wilcox.test(v1, v2, alternative = alternative, exact = FALSE, correct = TRUE)
        w <- unname(r$statistic)
        c(w, w / (sum(!is.na(v1)) * sum(!is.na(v2))), r$p.value)
    }))
    res[is.nan(res)] <- NA_real_
    colnames(res) <- c("W", "auc", "pval")
    res
}

expect_tests_equal <- function(x, labels) {
    for (alternative in c("two.sided", "less", "greater")) {
        expect_equal(tgs_ks(x, labels, alternative), ks_ref(x, labels, alternative), tolerance = 1e-10)
        expect_equal(tgs_wilcox(x, labels, alternative), wilcox_ref(x, labels, alternative), tolerance = 1e-10)
    }
}

test_that("tgs_ks and tgs_wilcox match ks.test and wilcox.test", {
    set.seed(60427)
    x <- matrix(rnorm(300 * 40), 300, 40)
    labels <- rep(c(TRUE, FALSE), c(120, 180))
    x[labels, 1:20] <- x[labels, 1:20] + 0.5
    x[!labels, 21:30] <- x[!labels, 21:30] + 0.5
    expect_tests_equal(x, labels)
    expect_tests_equal(x[1:7, ], labels[1:7])
})

test_that("tgs_ks and tgs_wilcox match base R with ties, missing values and infinite values", {
    set.seed(60427)
    x <- matrix(sample(0:5, 400 * 30, replace = TRUE), 400, 30)
    labels <- sample(c(TRUE, FALSE), 400, replace = TRUE, prob = c(0.3, 0.7))
    expect_true(is.integer(x))
    expect_tests_equal(x, labels)

    x[sample(length(x), 1000)] <- NA
    labels[sample(400, 40)] <- NA
    expect_tests_equal(x, labels)
    expect_equal(tgs_ks(x, labels), tgs_ks(x + 0, labels))
    expect_equal(tgs_wilcox(x, labels), tgs_wilcox(x + 0, labels))

    x <- x + 0
    x[sample(length(x), 300)] <- NaN
    x[sample(length(x), 300)] <- Inf
    x[sample(length(x), 300)] <- -Inf
    expect_tests_equal(x, labels)
})

test_that("tgs_ks and tgs_wilcox handle degenerate columns", {
    x <- cbind(
        constant = rep(3, 10),
        no_true = c(NA, NA, NA, 1:7),
        no_false = c(1:3, rep(NA, 7)),
        one_each = c(1, NA, NA, 2, rep(NA, 6)),
        tied_one_each = c(1, NA, NA, 1, rep(NA, 6)),
        regular = c(5, 1, 4, 2, 3, 6, 8, 7, 9, 10)
    )
    labels <- rep(c(TRUE, FALSE), c(3, 7))
    expect_tests_equal(x, labels)

    ks <- tgs_ks(x, labels)
    expect_equal(unname(ks["constant", ]), c(0, 1))
    expect_true(all(is.na(ks[c("no_true", "no_false"), ])))

    for (alternative in c("two.sided", "less", "greater")) {
        w <- tgs_wilcox(x, labels, alternative)
        expect_equal(unname(w["constant", c("W", "auc")]), c(10.5, 0.5))
        expect_equal(unname(w["constant", "pval"]), if (alternative == "two.sided") NA_real_ else 1)
        expect_true(all(is.na(w[c("no_true", "no_false"), ])))
    }

    expect_true(all(is.na(tgs_ks(x, rep(TRUE, 10)))))
    expect_true(all(is.na(tgs_wilcox(x, rep(NA, 10)))))
})

test_that("tgs_ks and tgs_wilcox accept a data frame and keep the column names", {
    set.seed(60427)
    df <- data.frame(a = rnorm(50), b = sample(1:4, 50, replace = TRUE), c = rnorm(50) + 1)
    labels <- rep(c(TRUE, FALSE), 25)
    expect_tests_equal(df, labels)
    expect_equal(rownames(tgs_ks(df, labels)), c("a", "b", "c"))
    expect_equal(rownames(tgs_wilcox(as.matrix(df), labels)), c("a", "b", "c"))
    expect_null(rownames(tgs_wilcox(unname(as.matrix(df)), labels)))
    expect_equal(colnames(tgs_ks(df, labels)), c("D", "pval"))
    expect_equal(colnames(tgs_wilcox(df, labels)), c("W", "auc", "pval"))
})

test_that("tgs_ks and tgs_wilcox reject invalid arguments", {
    x <- matrix(rnorm(20), 10, 2)
    labels <- rep(c(TRUE, FALSE), 5)
    for (f in list(tgs_ks, tgs_wilcox)) {
        expect_error(f(x))
        expect_error(f(x, as.numeric(labels)), "labels")
        expect_error(f(x, labels[-1]), "labels")
        expect_error(f(x, labels, "both"))
        expect_error(f(as.vector(x), rep(labels, 2)), "matrix")
        expect_error(f(matrix(letters[1:20], 10, 2), labels), "matrix")
        expect_error(f(data.frame(a = letters[1:10]), labels), "matrix")
        expect_error(f(x[, 0], labels), "matrix")
    }
})

test_that("tgs_ks and tgs_wilcox do not depend on the number of processes", {
    set.seed(60427)
    # 2000 x 1000 values: 4 processes when allowed
    x <- matrix(rnorm(2000 * 1000), 2000, 1000)
    x[sample(length(x), 20000)] <- NA
    labels <- sample(c(TRUE, FALSE, NA), 2000, replace = TRUE)
    res <- lapply(c(1, 4), function(np) {
        old <- options(tgs_max.processes = np)
        on.exit(options(old))
        list(tgs_ks(x, labels, "less"), tgs_wilcox(x, labels))
    })
    expect_identical(res[[1]], res[[2]])
    cols <- c(1:5, 996:1000)
    expect_equal(res[[2]][[1]][cols, ], ks_ref(x[, cols], labels, "less"), tolerance = 1e-10)
    expect_equal(res[[2]][[2]][cols, ], wilcox_ref(x[, cols], labels), tolerance = 1e-10)
})
