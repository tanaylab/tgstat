test_that("an R error in a child process is reported and the child does not return to R", {
    pid <- Sys.getpid()
    x <- matrix(as.numeric(1:40), 4, 10)
    index <- factor(rep(1:2, 5))
    # the child builds the call of the function with the argument names, and R does not accept a name that long
    args <- structure(list(1), names = strrep("a", 10001))
    msg <- tryCatch(
        .tgs_call("tgs_matrix_tapply", x, index, function(v, ...) 1, "fun", args, 2L, new.env()),
        error = function(e) conditionMessage(e)
    )

    # a child that got here would go on running the tests as a second process
    if (Sys.getpid() != pid) {
        tools::pskill(Sys.getpid(), tools::SIGKILL)
    }

    expect_match(msg, "variable names are limited")
    expect_true(dir.exists(tempdir()))
})
