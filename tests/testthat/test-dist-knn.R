# reference: the knn nearest rows of y to each row of x, ties broken by the row number. The distances are those of
# dist(): missing dimensions are dropped and the sum of squares is scaled up, NA if there are none; non-finite values
# are missing, as in tgstat.
dist_knn_ref <- function(x, y, knn, threshold = Inf) {
    x[!is.finite(x)] <- NA
    ty <- t(if (is.null(y)) x else y)
    ty[!is.finite(ty)] <- NA
    res <- lapply(seq_len(nrow(x)), function(i) {
        sq <- (ty - x[i, ])^2
        count <- colSums(!is.na(sq))
        d <- sqrt(colSums(sq, na.rm = TRUE) * nrow(ty) / count)
        d[count == 0] <- NA
        if (is.null(y)) {
            d[i] <- NA
        }
        rows <- which(!is.na(d) & d <= threshold)
        rows <- utils::head(rows[order(d[rows], rows)], knn)
        list(rows = rows, dist = unname(d[rows]))
    })
    num_neighbours <- vapply(res, function(r) length(r$rows), integer(1))
    data.frame(
        row1 = rep(seq_len(nrow(x)), num_neighbours),
        row2 = as.integer(unlist(lapply(res, `[[`, "rows"))),
        dist = as.numeric(unlist(lapply(res, `[[`, "dist"))),
        rank = sequence(num_neighbours)
    )
}

expect_dist_knn_equal <- function(res, ref) {
    expect_named(res, c("row1", "row2", "dist", "rank"))
    expect_equal(as.integer(res$row1), ref$row1)
    expect_equal(as.integer(res$row2), ref$row2)
    expect_equal(res$dist, ref$dist, tolerance = 1e-10)
    expect_equal(res$rank, ref$rank)
}

with_blas <- function(use_blas, expr) {
    old <- options(tgs_use.blas = use_blas)
    on.exit(options(old))
    expr
}

test_that("tgs_dist_knn with and without BLAS matches the reference", {
    set.seed(60427)
    x <- matrix(rnorm(60 * 8), 60)
    y <- matrix(rnorm(25 * 8), 25)
    xi <- matrix(sample(0:1000, 60 * 8, replace = TRUE), 60)
    yi <- matrix(sample(0:1000, 25 * 8, replace = TRUE), 25)

    for (use_blas in c(TRUE, FALSE)) {
        expect_dist_knn_equal(with_blas(use_blas, tgs_dist_knn(x, knn = 7)), dist_knn_ref(x, NULL, 7))
        expect_dist_knn_equal(with_blas(use_blas, tgs_dist_knn(x, y, 7)), dist_knn_ref(x, y, 7))
        expect_dist_knn_equal(with_blas(use_blas, tgs_dist_knn(xi, knn = 7)), dist_knn_ref(xi, NULL, 7))
        expect_dist_knn_equal(with_blas(use_blas, tgs_dist_knn(xi, yi, 7)), dist_knn_ref(xi, yi, 7))
        expect_dist_knn_equal(with_blas(use_blas, tgs_dist_knn(x, y[1, , drop = FALSE], 3)), dist_knn_ref(x, y[1, , drop = FALSE], 3))
        expect_dist_knn_equal(with_blas(use_blas, tgs_dist_knn(x[1, , drop = FALSE], y, 3)), dist_knn_ref(x[1, , drop = FALSE], y, 3))

        # the nearest rows by tgs_dist
        d <- with_blas(use_blas, tgs_dist(x, y))
        res <- with_blas(use_blas, tgs_dist_knn(x, y, 1))
        expect_equal(as.integer(res$row2), unname(apply(d, 1, which.min)))
    }
})

test_that("tgs_dist_knn handles missing values", {
    set.seed(60427)
    x <- matrix(rnorm(50 * 6), 50)
    y <- matrix(rnorm(20 * 6), 20)
    x[sample(length(x), 40)] <- NA
    x[7, 3] <- Inf
    y[sample(length(y), 10)] <- NaN
    y[4, ] <- NA
    # rows that share no dimension with some of the others: their distance is NA
    x[10, ] <- c(NA, NA, NA, 1, 2, 3)
    x[11, ] <- c(1, 2, 3, NA, NA, NA)
    y[8, ] <- c(1, 2, 3, NA, NA, NA)

    for (use_blas in c(TRUE, FALSE)) {
        for (knn in c(5, 1000)) {
            res <- with_blas(use_blas, tgs_dist_knn(x, NULL, knn))
            expect_dist_knn_equal(res, dist_knn_ref(x, NULL, knn))
            expect_false(any(is.na(res$dist)))
            expect_false(11 %in% res$row2[res$row1 == 10])
            res <- with_blas(use_blas, tgs_dist_knn(x, y, knn))
            expect_dist_knn_equal(res, dist_knn_ref(x, y, knn))
            expect_false(4 %in% res$row2)
            expect_false(8 %in% res$row2[res$row1 == 10])
        }
    }
})

test_that("tgs_dist_knn returns factors of the row names", {
    set.seed(60427)
    x <- matrix(rnorm(30 * 4), 30, dimnames = list(paste0("x", 1:30), NULL))
    y <- matrix(rnorm(10 * 4), 10, dimnames = list(paste0("y", 1:10), NULL))
    for (use_blas in c(TRUE, FALSE)) {
        res <- with_blas(use_blas, tgs_dist_knn(x, y, 3))
        expect_dist_knn_equal(res, dist_knn_ref(x, y, 3))
        expect_equal(levels(res$row1), rownames(x))
        expect_equal(levels(res$row2), rownames(y))
        res <- with_blas(use_blas, tgs_dist_knn(x, NULL, 3))
        expect_equal(levels(res$row1), rownames(x))
        expect_equal(levels(res$row2), rownames(x))
        res <- with_blas(use_blas, tgs_dist_knn(x, unname(y), 3))
        expect_s3_class(res$row1, "factor")
        expect_type(res$row2, "integer")
    }
})

test_that("tgs_dist_knn applies the threshold and caps knn at the number of rows", {
    set.seed(60427)
    x <- matrix(rnorm(40 * 3), 40)
    y <- matrix(rnorm(15 * 3), 15)
    for (use_blas in c(TRUE, FALSE)) {
        for (threshold in c(0.5, 1, -1)) {
            res <- with_blas(use_blas, tgs_dist_knn(x, NULL, 5, threshold = threshold))
            expect_dist_knn_equal(res, dist_knn_ref(x, NULL, 5, abs(threshold)))
            expect_true(all(res$dist <= abs(threshold)))
            expect_dist_knn_equal(with_blas(use_blas, tgs_dist_knn(x, y, 5, threshold = threshold)), dist_knn_ref(x, y, 5, abs(threshold)))
        }
        expect_equal(nrow(with_blas(use_blas, tgs_dist_knn(x, knn = 3, threshold = 0))), 0)

        res <- with_blas(use_blas, tgs_dist_knn(x, NULL, 1000))
        expect_dist_knn_equal(res, dist_knn_ref(x, NULL, 39))
        expect_true(all(table(res$row1) == 39))
        res <- with_blas(use_blas, tgs_dist_knn(x, y, 1000))
        expect_dist_knn_equal(res, dist_knn_ref(x, y, 15))
        expect_true(all(table(res$row1) == 15))
        expect_dist_knn_equal(with_blas(use_blas, tgs_dist_knn(x[1:2, ], NULL, 5)), dist_knn_ref(x[1:2, ], NULL, 1))
    }
})

test_that("tgs_dist_knn gives distance 0 to identical rows and breaks ties by the row number", {
    for (seed in 1:3) {
        set.seed(seed)
        x <- matrix(runif(40 * 30) * 1000, 40)
        x[c(9, 20, 33), ] <- x[rep(25, 3), ]
        x[c(2, 3), ] <- x[c(1, 1), ]
        res <- lapply(c(TRUE, FALSE), function(use_blas) with_blas(use_blas, tgs_dist_knn(x, NULL, 4)))
        expect_identical(res[[1]][, c("row1", "row2", "rank")], res[[2]][, c("row1", "row2", "rank")])
        expect_equal(res[[1]]$dist, res[[2]]$dist)
        expect_false(any(is.nan(res[[1]]$dist)))
        for (r in res) {
            expect_dist_knn_equal(r, dist_knn_ref(x, NULL, 4))
            expect_identical(r$row2[r$row1 == 25][1:3], c(9L, 20L, 33L))
            expect_identical(r$row2[r$row1 == 20][1:3], c(9L, 25L, 33L))
            expect_identical(r$dist[r$row1 == 25][1:3], rep(0, 3))
            expect_identical(r$row2[r$row1 == 3][1:2], c(1L, 2L))
            expect_identical(r$dist[r$row1 == 3][1:2], rep(0, 2))
        }

        # with y = x each row is its own neighbour, at distance 0
        res_x <- lapply(c(TRUE, FALSE), function(use_blas) with_blas(use_blas, tgs_dist_knn(x, x, 5)))
        expect_identical(res_x[[1]][, c("row1", "row2", "rank")], res_x[[2]][, c("row1", "row2", "rank")])
        for (r in res_x) {
            expect_dist_knn_equal(r, dist_knn_ref(x, x, 5))
            expect_true(all(r$dist[r$row1 == r$row2] == 0))
            expect_equal(sum(r$row1 == r$row2), nrow(x))
        }
    }
})

test_that("tgs_dist_knn with y = x includes the row itself, unlike y = NULL", {
    set.seed(60427)
    x <- matrix(rnorm(50 * 5), 50)
    for (use_blas in c(TRUE, FALSE)) {
        res <- with_blas(use_blas, tgs_dist_knn(x, NULL, 4))
        res_x <- with_blas(use_blas, tgs_dist_knn(x, x, 5))
        expect_equal(res_x$row2[res_x$rank == 1], 1:50)
        expect_identical(res_x$dist[res_x$rank == 1], rep(0, 50))
        res_x <- res_x[res_x$rank > 1, ]
        res_x$rank <- res_x$rank - 1L
        rownames(res_x) <- NULL
        expect_equal(res_x, res)
    }
})

test_that("tgs_dist_knn with BLAS on more rows than one block matches the version without BLAS", {
    set.seed(60427)
    # 3000 x 3000 distances: blocks of 2^23 / 3000 rows
    x <- matrix(rnorm(3000 * 5), 3000)
    res <- with_blas(TRUE, tgs_dist_knn(x, NULL, 6))
    expect_dist_knn_equal(res, dist_knn_ref(x, NULL, 6))
    expect_identical(res[, c("row1", "row2", "rank")], with_blas(FALSE, tgs_dist_knn(x, NULL, 6))[, c("row1", "row2", "rank")])

    # 6000 x 1500 distances: blocks of 2^23 / 1500 rows
    x <- matrix(rnorm(6000 * 5), 6000)
    y <- matrix(rnorm(1500 * 5), 1500)
    res <- with_blas(TRUE, tgs_dist_knn(x, y, 6, threshold = 1))
    expect_dist_knn_equal(res, dist_knn_ref(x, y, 6, 1))
    expect_identical(res[, c("row1", "row2", "rank")], with_blas(FALSE, tgs_dist_knn(x, y, 6, threshold = 1))[, c("row1", "row2", "rank")])
})

test_that("tgs_dist_knn does not depend on the number of processes", {
    set.seed(60427)
    x <- matrix(rnorm(120 * 5), 120)
    y <- matrix(rnorm(30 * 5), 30)
    x[sample(length(x), 30)] <- NA
    old <- options(tgs_use.blas = FALSE, tgs_max.processes = 1)
    on.exit(options(old))
    expected <- tgs_dist_knn(x, NULL, 5)
    expected_cross <- tgs_dist_knn(x, y, 5)
    expect_dist_knn_equal(expected, dist_knn_ref(x, NULL, 5))
    expect_dist_knn_equal(expected_cross, dist_knn_ref(x, y, 5))
    for (p in c(3, 7, 11, 13, 50)) {
        options(tgs_max.processes = p)
        expect_identical(tgs_dist_knn(x, NULL, 5), expected)
        expect_identical(tgs_dist_knn(x, y, 5), expected_cross)
    }
})

test_that("tgs_dist_knn accepts data frames and rejects invalid input", {
    set.seed(60427)
    x <- matrix(rnorm(20 * 3), 20)
    y <- matrix(rnorm(8 * 3), 8)
    for (use_blas in c(TRUE, FALSE)) {
        expect_dist_knn_equal(with_blas(use_blas, tgs_dist_knn(as.data.frame(x), as.data.frame(y), 2)), dist_knn_ref(x, y, 2))
        expect_error(with_blas(use_blas, tgs_dist_knn(x, y[, 1:2], 2)), "different number of columns")
        expect_error(with_blas(use_blas, tgs_dist_knn(x, TRUE, 2)), "matrix of numeric values")
        expect_error(with_blas(use_blas, tgs_dist_knn(matrix(letters[1:6], 2), NULL, 2)), "matrix of numeric values")
        expect_error(with_blas(use_blas, tgs_dist_knn(x[1, , drop = FALSE], NULL, 2)), "at least two rows")
        for (knn in list(0, -1, 1.5, NA, 1:2)) {
            expect_error(with_blas(use_blas, tgs_dist_knn(x, y, knn)), "knn")
        }
    }
    expect_error(tgs_dist_knn(x), "Usage")
})
