# tgstat 2.0.2

* Major run-time optimization when spearman=T and pairwise.complete.obs=T and NaNs are presented in the data
* Bug fix: tgs_cor might hang when spearman=T and pairwise.complete.obs=T
* Bug fix: tgs_cor / tgs_dist might hang when tgs_debug option is TRUE


# tgstat 1.0.2

* New function: tgs_dist
* tgs_cor: in tidy format return column names (if exist) instead of their indices
* tgs_cor: minor run-time optimizations
* Bug fix in tgs_cor: results > threshold are returned instead of >= threshold

# tgstat 1.0.1

* tgs_cor: preserve dimnames in the result
