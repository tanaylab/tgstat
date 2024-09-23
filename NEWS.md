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
