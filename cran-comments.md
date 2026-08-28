## R CMD check results

0 errors | 0 warnings | 0 notes

## Changes

* Fixed the installation failure with clang 23 / libc++, which no longer
  provides `<algorithm>` transitively. Missing standard headers were added
  explicitly across the C++ sources.
* Fixed the `PROTECT` issues reported by rchk in `tgs_knn()` and
  `tgs_cor_graph()`.
* Added `tgs_chi2()`.
