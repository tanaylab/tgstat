test_that("the BLAS guess looks at the library name and its directory name only", {
    guess <- tgstat:::.tgs_guess_use_blas
    expect_true(guess("/usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3"))
    expect_true(guess("/usr/lib/x86_64-linux-gnu/atlas/libblas.so.3"))
    expect_true(guess("/usr/lib64/libopenblasp.so.0"))
    expect_true(guess("/opt/intel/oneapi/mkl/latest/lib/libmkl_rt.so.2"))
    expect_false(guess("/usr/lib/R/lib/libRblas.so"))
    expect_false(guess("/usr/lib/x86_64-linux-gnu/blas/libblas.so.3"))
    expect_false(guess("/home/mklein/R/lib/libRblas.so"))
    expect_false(guess("/opt/conda/envs/atlas-rna/lib/libblas.so.3"))
    expect_type(guess(), "logical")
})
