# tgstat 2.3.13

- tgs_matrix_tapply: set correctly the column names of the resulted matrix.
- tgs_cor_knn: fix alignment issues on some platforms.
- Stop using non-portable bswap_32 and bswap_64.

# tgstat 2.3.12

- Added authors with DOI to package description.

# tgstat 2.3.11

- Fixed compilation issues on some platforms.

# tgstat 2.3.10

- Bug fix in tgs_matrix_tapply: "fn_name is not a string" error when 'fun' is a function
defined inline.
- Bug fix in tgs_matrix_tapply: occasional "stack imbalance" warning.
