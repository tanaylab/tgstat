tgstat package - Tanay's group statistical utilities
====================================================

The goal of tgstat is to provide fast and efficient statistical tools.

Installation
------------
``` r
install.packages('tgstat', repos='tanaylab.bitbucket.io/repo')
```

or:

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
#>             [,1]        [,2]        [,3]        [,4]       [,5]
#> [1,]  1.00000000  0.14245962 -0.09470294  0.01486385 0.22129726
#> [2,]  0.14245962  1.00000000 -0.06816740  0.06086571 0.29643287
#> [3,] -0.09470294 -0.06816740  1.00000000 -0.03645632 0.01065091
#> [4,]  0.01486385  0.06086571 -0.03645632  1.00000000 0.08422610
#> [5,]  0.22129726  0.29643287  0.01065091  0.08422610 1.00000000
```

``` r
r2 <- tgs_cor(m, pairwise.complete.obs = T, spearman = T)
r2[1:5, 1:5]
#>               [,1]        [,2]        [,3]          [,4]       [,5]
#> [1,]  1.0000000000  0.16414041 -0.09527753  0.0002160216 0.23579958
#> [2,]  0.1641404140  1.00000000 -0.06214221  0.0594179418 0.30271827
#> [3,] -0.0952775278 -0.06214221  1.00000000 -0.0371557156 0.01680168
#> [4,]  0.0002160216  0.05941794 -0.03715572  1.0000000000 0.07181518
#> [5,]  0.2357995800  0.30271827  0.01680168  0.0718151815 1.00000000
```

Especially useful for large matrices:

``` r
rows <- 1000
cols <- 1000
vals <- sample(1 : (rows * cols / 2), rows * cols, replace = T)
m <- matrix(vals, nrow = rows, ncol = cols)
m[sample(1 : (rows * cols), rows * cols / 1000)] <- NA

system.time(r1 <- tgs_cor(m, pairwise.complete.obs=TRUE, spearman = F))
#>    user  system elapsed 
#>   4.713   0.743   0.281
```

Compared to base R:

``` r
system.time(r2 <- cor(m, use='pairwise.complete.obs'))
#>    user  system elapsed 
#>  10.600   0.006  10.603
```

Euclidian distance
------------------

``` r
system.time(d <- tgs_dist(m))
#>    user  system elapsed 
#>   1.683   0.681   0.195
```

Compared to base R:

``` r
system.time(d <- dist(m, method='euclidian'))
#>    user  system elapsed 
#>   3.934   0.005   3.939
```

``` r
sessionInfo()
#> R version 3.3.2 (2016-10-31)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: CentOS Linux 7 (Core)
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> attached base packages:
#> [1] graphics  datasets  grid      stats     grDevices utils     methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] tgstat_1.0.2    cowplot_0.8.0   stringr_1.2.0   dplyr_0.7.2    
#>  [5] purrr_0.2.3     readr_1.1.1     tidyr_0.6.3     tibble_1.3.3   
#>  [9] ggplot2_2.2.1   tidyverse_1.1.1 reshape2_1.4.2  plyr_1.8.4     
#> [13] scales_0.4.1    optparse_1.3.2  pacman_0.4.6    misha_3.5.6    
#> 
#> loaded via a namespace (and not attached):
#>  [1] Rcpp_0.12.12         lubridate_1.6.0      lattice_0.20-35     
#>  [4] assertthat_0.2.0     rprojroot_1.2        digest_0.6.12       
#>  [7] psych_1.7.5          R6_2.2.2             cellranger_1.1.0    
#> [10] backports_1.1.0      evaluate_0.10.1      httr_1.2.1          
#> [13] GlobalOptions_0.0.12 rlang_0.1.2.9000     lazyeval_0.2.0      
#> [16] readxl_1.0.0         data.table_1.10.4    whisker_0.3-2       
#> [19] GetoptLong_0.1.6     pkgdown_0.1.0.9000   rmarkdown_1.6       
#> [22] desc_1.1.1           devtools_1.13.3      foreign_0.8-69      
#> [25] munsell_0.4.3        broom_0.4.2          modelr_0.1.0        
#> [28] pkgconfig_2.0.1      mnormt_1.5-5         htmltools_0.3.6     
#> [31] roxygen2_6.0.1       crayon_1.3.2         withr_2.0.0         
#> [34] commonmark_1.2       nlme_3.1-131         jsonlite_1.5        
#> [37] gtable_0.2.0         git2r_0.19.0         magrittr_1.5        
#> [40] stringi_1.1.5        bindrcpp_0.2         getopt_1.20.0       
#> [43] xml2_1.1.1           rjson_0.2.15         tools_3.3.2         
#> [46] forcats_0.2.0        glue_1.1.1           hms_0.3             
#> [49] pkgload_0.0.0.9000   parallel_3.3.2       yaml_2.1.14         
#> [52] colorspace_1.3-2     rvest_0.3.2          memoise_1.1.0       
#> [55] knitr_1.17           bindr_0.1            haven_1.0.0
```
