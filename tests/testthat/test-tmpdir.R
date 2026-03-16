test_that("tgs_cor works with custom tgs_tmpdir option", {
    set.seed(60427)
    x <- matrix(rnorm(100), nrow = 10)

    tmpdir <- tempfile("tgstat_test_tmpdir")
    dir.create(tmpdir)
    on.exit(unlink(tmpdir, recursive = TRUE))

    old_opt <- getOption("tgs_tmpdir")
    options(tgs_tmpdir = tmpdir)
    on.exit(options(tgs_tmpdir = old_opt), add = TRUE)

    res <- tgs_cor(x)
    expect_equal(res, cor(x))
})

test_that("tgs_cor_knn works with custom tgs_tmpdir option", {
    set.seed(60427)
    x <- matrix(rnorm(200), nrow = 20, ncol = 10)
    knn <- 3

    tmpdir <- tempfile("tgstat_test_tmpdir")
    dir.create(tmpdir)
    on.exit(unlink(tmpdir, recursive = TRUE))

    old_opt <- getOption("tgs_tmpdir")
    options(tgs_tmpdir = tmpdir)
    on.exit(options(tgs_tmpdir = old_opt), add = TRUE)

    res <- tgs_cor_knn(x, NULL, knn)
    expect_s3_class(res, "data.frame")
    expect_equal(nrow(res), ncol(x) * knn)
})

test_that("tgs_dist works with custom tgs_tmpdir option", {
    set.seed(60427)
    x <- matrix(rnorm(100), nrow = 10)

    tmpdir <- tempfile("tgstat_test_tmpdir")
    dir.create(tmpdir)
    on.exit(unlink(tmpdir, recursive = TRUE))

    old_opt <- getOption("tgs_tmpdir")
    options(tgs_tmpdir = tmpdir)
    on.exit(options(tgs_tmpdir = old_opt), add = TRUE)

    res <- tgs_dist(x)
    expect_equal(as.matrix(res), as.matrix(dist(x)))
})

test_that("TMPDIR env var is respected when tgs_tmpdir option is not set", {
    set.seed(60427)
    x <- matrix(rnorm(100), nrow = 10)

    tmpdir <- tempfile("tgstat_test_tmpdir")
    dir.create(tmpdir)
    on.exit(unlink(tmpdir, recursive = TRUE))

    old_opt <- getOption("tgs_tmpdir")
    old_env <- Sys.getenv("TMPDIR", unset = NA)

    options(tgs_tmpdir = NULL)
    Sys.setenv(TMPDIR = tmpdir)
    on.exit(
        {
            options(tgs_tmpdir = old_opt)
            if (is.na(old_env)) Sys.unsetenv("TMPDIR") else Sys.setenv(TMPDIR = old_env)
        },
        add = TRUE
    )

    res <- tgs_cor(x)
    expect_equal(res, cor(x))
})

test_that("tgs_tmpdir option takes priority over TMPDIR env var", {
    set.seed(60427)
    x <- matrix(rnorm(200), nrow = 20, ncol = 10)
    knn <- 3

    tmpdir1 <- tempfile("tgstat_test_opt")
    tmpdir2 <- tempfile("tgstat_test_env")
    dir.create(tmpdir1)
    dir.create(tmpdir2)
    on.exit(unlink(c(tmpdir1, tmpdir2), recursive = TRUE))

    old_opt <- getOption("tgs_tmpdir")
    old_env <- Sys.getenv("TMPDIR", unset = NA)

    # Set both: option should win
    options(tgs_tmpdir = tmpdir1)
    Sys.setenv(TMPDIR = tmpdir2)
    on.exit(
        {
            options(tgs_tmpdir = old_opt)
            if (is.na(old_env)) Sys.unsetenv("TMPDIR") else Sys.setenv(TMPDIR = old_env)
        },
        add = TRUE
    )

    res <- tgs_cor_knn(x, NULL, knn)
    expect_s3_class(res, "data.frame")
    expect_equal(nrow(res), ncol(x) * knn)
})

test_that("multiprocess works with custom tgs_tmpdir", {
    set.seed(60427)
    x <- matrix(rnorm(200), nrow = 20, ncol = 10)
    knn <- 3

    tmpdir <- tempfile("tgstat_test_tmpdir")
    dir.create(tmpdir)
    on.exit(unlink(tmpdir, recursive = TRUE))

    old_tmpdir_opt <- getOption("tgs_tmpdir")
    old_proc_opt <- getOption("tgs_max.processes")
    options(tgs_tmpdir = tmpdir, tgs_max.processes = 4)
    on.exit(options(tgs_tmpdir = old_tmpdir_opt, tgs_max.processes = old_proc_opt), add = TRUE)

    res <- tgs_cor_knn(x, NULL, knn)
    expect_s3_class(res, "data.frame")
    expect_equal(nrow(res), ncol(x) * knn)
})
