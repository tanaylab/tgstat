tgstat package - Tanay's group statistical utilities
====================================================

The goal of tgstat is to provide fast and efficient statistical tools.

Installation
------------

``` r
devtools::install_bitbucket('tanaylab/tgstat@default')
```

Correlation
-----------

``` r
library(tgstat)
rows <- 100
cols <- 1000
vals <- sample(1 : (rows * cols / 2), rows * cols, replace = T)
m <- matrix(vals, nrow = rows, ncol = cols)
m[sample(1 : (rows * cols), rows * cols / 1000)] <- NA

r1 <- tgs_cor(m, spearman = F)

r1[1:5, 1:5]
#>      [,1]       [,2]        [,3]        [,4]       [,5]
#> [1,]    1         NA          NA          NA         NA
#> [2,]   NA 1.00000000  0.03614276  0.12416333 0.10248252
#> [3,]   NA 0.03614276  1.00000000 -0.08007816 0.04201266
#> [4,]   NA 0.12416333 -0.08007816  1.00000000 0.05695502
#> [5,]   NA 0.10248252  0.04201266  0.05695502 1.00000000
```

``` r
r2 <- tgs_cor(m, pairwise.complete.obs = T, spearman = T)
r2[1:5, 1:5]
#>             [,1]       [,2]        [,3]        [,4]       [,5]
#> [1,]  1.00000000 0.07287901 -0.03257273  0.01080657 0.10550912
#> [2,]  0.07287901 1.00000000  0.05068107  0.12763276 0.11979598
#> [3,] -0.03257273 0.05068107  1.00000000 -0.08123612 0.05261326
#> [4,]  0.01080657 0.12763276 -0.08123612  1.00000000 0.03492349
#> [5,]  0.10550912 0.11979598  0.05261326  0.03492349 1.00000000
```
