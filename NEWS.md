# tgstat 2.5.0

* Added `tgs_cor_dist()`: the same as `tgs_dist(tgs_cor(x))`, but much faster and with less memory when `x` has many more columns than rows.
* `tgs_dist()` uses about half the memory on large inputs.
* Debian/Ubuntu installations with OpenBLAS (e.g. the rocker images) now use the faster BLAS code by default.
* Fixed `tgs_dist()` and `tgs_cor()` without BLAS, and `tgs_cor_knn()`, skipping one row or column for some input sizes when the number of processes is not a power of 2.
* Fixed `tgs_dist()` with BLAS returning NaN for identical points when the data has missing values.
* `tgs_cor()` with BLAS ranks and standardizes the columns in parallel and no longer holds the correlation matrix twice in memory: Spearman correlations are about 3 times faster and Pearson correlations up to 2 times.
* Added `tgs_dist_knn()`: the `knn` nearest rows of each row of `x`, among the rows of `y` or of `x` itself, by euclidean distance, without computing the full distance matrix.
* `tgs_dist()` gained a `y` argument: `tgs_dist(x, y)` returns the distances between the rows of `x` and the rows of `y`. `y` is now the second argument, so passing `diag` by position (`tgs_dist(x, TRUE)`) is an error.
* `tgs_dist()` with BLAS returns exactly 0 for identical rows; it used to return a small positive distance.
* `tgs_dist()`, `tgs_dist_knn()` and `tgs_cor_dist()` accept data frames.
* `tgs_cor()` and `tgs_cor_knn()` accept data frames. `tgs_cor()` accepts a vector as `y` (and then also as `x`), for the correlations of one vector with all the columns of a matrix. `y` of `tgs_cor_knn()` defaults to `NULL`.
* `tgs_cor_knn()` uses BLAS when `tgs_use.blas` is on, which makes it several times faster.
* `tgs_knn()` on a matrix is much faster and no longer needs memory for an index of the whole matrix.
* `tgs_matrix_tapply()` computes `max()`, `min()` and `median()` in C++, like `sum()` and `mean()`, which makes them much faster.
* Added a `margin` argument to `tgs_matrix_tapply()`: `margin = 1` groups the rows of the matrix instead of the columns, without transposing it (for `sum`, the same as `rowsum()`).
* `tgs_matrix_tapply()` with `sum()` or `mean()` on a sparse matrix takes time proportional to the number of non-zero values instead of the size of the matrix.
* Fixed `tgs_matrix_tapply()` returning only `NA` when arguments are passed to `fun` through `...` (including `mean()` with `trim`), and -2147483648 instead of `NA` when `fun` returns an integer `NA`.
* Fixed `tgs_matrix_tapply()` with `sum()` or `mean()` on an integer matrix with `NA`, which counted `NA` as -2147483648.
* A child process that ends without computing its part of the result is now an error in all the functions. It used to be ignored if the process exited (e.g. after an R error), leaving `NA` in the result.
* Fixed `tgs_cor()` returning `-Inf` or arbitrary small values instead of `NA` for constant columns.
* Fixed `tgs_cor_knn()` dropping valid neighbours when the data has a constant column.

# tgstat 2.4.0

* Added `tgs_chi2()`: vectorized chi-squared test for 2x2 contingency tables, with optional Yates' correction and support for sparse (`dgCMatrix`) input.
* Fixed compilation with clang 23 and the `PROTECT` issues reported by rchk.

# tgstat 2.3.32

* Replaced non-API C entry point `Rf_findVar` with `R_getVar` for R 4.6.0 compatibility.
* Added `.claude` to `.Rbuildignore`.

# tgstat 2.3.30

* FIFO temp directory is now configurable via `options(tgs_tmpdir = "/path")` or the `TMPDIR` environment variable, instead of being hardcoded to `/tmp`. Fixes permission errors on systems where `/tmp` is restricted (#17).
* Added a child-process liveness check in `wait_for_kid()`/`wait_for_kids()` to prevent potential hangs when child processes terminate unexpectedly without being reaped (#18).
* Added comprehensive tests for `tgs_cor_knn` and temp directory configuration.

# tgstat 2.3.29

* Removed non-API calls to R: `Rf_GetOption`.

# tgstat 2.3.28 

* Fixed *noRemap* additional issue by adding `Rf_` prefix in the c++ code.
* Removed non-API calls to R: `R_curErrorBuf`, `SET_TYPEOF` in the c++ code.


# tgstat 2.3.27

* Fixed URL of the package in the DESCRIPTION file.

# tgstat 2.3.26

* Guess `tgs.use_blas` from environment if not set.

# tgstat 2.3.25

* Fixed clang warnings.

# tgstat 2.3.24

* Fixed compilation errors on clang17.

# tgstat 2.3.23

* Added missing <memory> header. 

# tgstat 2.3.22

* Fixed compilation errors on clang16.

# tgstat 2.3.20

* added missing <cstdint> include in order to compile with gcc 13

# tgstat 2.3.19 

* Fixed rchk warning with the patch generously provided by Tomas Kalibera (thanks!).
* Fixed a compilation warning on M1 macs.

# tgstat 2.3.18 

* Fixed a compilation warning on M1 macs. 

# tgstat 2.3.17

- Fix compilation issues with R 4.2.0
- Use roxygen2 and markdown for documentation

# tgstat 2.3.16

- Fix compilation issues on Debian.

# tgstat 2.3.15

- Fix compilation errors on 32-bit platforms.
- Fix compilation and run-time issues on Solaris.

# tgstat 2.3.14

- Stop using non-portable bswap_32 and bswap_64.

# tgstat 2.3.13

- tgs_matrix_tapply: set correctly the column names of the resulted matrix.
- tgs_cor_knn: fix alignment issues on some platforms.

# tgstat 2.3.12

- Added authors with DOI to package description.

# tgstat 2.3.11

- Fixed compilation issues on some platforms.

# tgstat 2.3.10

- Bug fix in tgs_matrix_tapply: "fn_name is not a string" error when 'fun' is a function
defined inline.
- Bug fix in tgs_matrix_tapply: occasional "stack imbalance" warning.
