<!-- badges: start -->
[![CRAN
status](https://www.r-pkg.org/badges/version/tgstat)](https://CRAN.R-project.org/package=tgstat)
[![Travis build
status](https://travis-ci.com/tanaylab/tgstat.svg?branch=master)](https://travis-ci.org/tanaylab/tgstat)
<!-- badges: end -->

tgstat
======

The goal of `tgstat` is to provide fast and efficient implementation of
certain R functions such as ‘cor’ and ‘dist’, along with specific
statistical tools.

Various approaches are used to boost the performance, including
multi-processing and use of optimized functions provided by the Basic
Linear Algebra Subprograms (BLAS) library.

Installation
------------

In order to install `tgstat`:

``` r
remotes::install_github("tanaylab/tgstat")
```

Examples
--------

``` r
set.seed(seed=1)
rows = 3000
cols = 3000
vals<-sample(1:(rows*cols/2), rows*cols, replace=T)
m<-matrix(vals, nrow=rows, ncol=cols)
m_with_NAs <- m
m_with_NAs[sample(1:(rows*cols), rows*cols / 10)] <- NA
dim(m)
#> [1] 3000 3000
```

### Fast computation of correlation matrices

Pearson correlation without BLAS, no NAs:

``` r
options(tgs_use.blas=F)
system.time(tgs_cor(m))
#>    user  system elapsed 
#>  19.817   1.573   1.099
```

Same with BLAS:

``` r
# tgs_cor, with BLAS, no NAs, pearson
options(tgs_use.blas=T)
system.time(tgs_cor(m))
#>    user  system elapsed 
#>   2.059   0.292   0.455
```

Base R version:

``` r
system.time(cor(m))
#>    user  system elapsed 
#>  18.120   0.096  18.260
```

Pearson correlation without BLAS, with NAs:

``` r
options(tgs_use.blas=F)
system.time(tgs_cor(m_with_NAs, pairwise.complete.obs=T))
#>    user  system elapsed 
#>  62.326   1.519   1.674
```

Same with BLAS:

``` r
options(tgs_use.blas=T)
system.time(tgs_cor(m_with_NAs, pairwise.complete.obs=T))
#>    user  system elapsed 
#>   6.145   0.784   0.651
```

Base R version:

``` r
system.time(cor(m_with_NAs, use="pairwise.complete.obs"))
#>    user  system elapsed 
#> 256.972   0.222 257.798
```

### Fast computation of distance matrices

Distance without BLAS, no NAs:

``` r
options(tgs_use.blas=F)
system.time(tgs_dist(m))
#> 84%...100%
#>    user  system elapsed 
#> 267.356   1.437   4.552
```

Same with BLAS:

``` r
options(tgs_use.blas=T)
system.time(tgs_dist(m))
#>    user  system elapsed 
#>   1.995   0.288   0.288
```

Base R:

``` r
system.time(dist(m, method="euclidean"))
#>    user  system elapsed 
#> 149.372   0.086 149.811
```

Notes regarding the usage of `BLAS`
-----------------------------------

`tgstat` runs best when R is linked with an optimized BLAS
implementation.

Many optimized BLAS implementations are available, both proprietary
(e.g. Intel’s MKL, Apple’s vecLib) and opensource (e.g. OpenBLAS,
ATLAS). Unfortunately, R often uses by default the reference BLAS
implementation, which is known to have poor performance.

Having `tgstat` rely on the reference BLAS will result in very poor
performance and is strongly discouraged. If your R implementation uses
an optimized BLAS, set `options(tgs_use.blas=TRUE)` to allow `tgstat` to
make BLAS calls. Otherwise, set `options(tgs_use.blas=FALSE)` (default)
which instructs `tgstat` to avoid BLAS and instead rely only on its own
optimization methods. If in doubt, it is possible to run one of `tgstat`
CPU intensive functions (e.g. `tgs_cor`) and compare its run time under
both `options(tgs_use.blas=FALSE)`.

Exact instructions for linking R with an optimized BLAS library are
system dependent and are out of scope of this document.
