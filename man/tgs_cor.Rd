\name{tgs_cor}
\alias{tgs_cor}
\alias{tgs_cor_knn}
\title{
Calculates correlation or auto-correlation
}
\description{
Calculates correlation between two matrices columns or auto-correlation between a matrix columns.
}
\usage{
tgs_cor(x, y = NULL, pairwise.complete.obs = F, spearman = F,
        tidy = F, threshold = 0)

tgs_cor_knn(x, y, knn, pairwise.complete.obs = F, spearman = F,
            threshold = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{numeric matrix}
  \item{y}{numeric matrix}
  \item{pairwise.complete.obs}{see below}
  \item{spearman}{if 'TRUE' Spearman correlation is computed, otherwise Pearson}
  \item{tidy}{if 'TRUE' data is outputed in tidy format}
  \item{threshold}{absolute threshold above which values are outputed in tidy format}
  \item{knn}{the number of highest correlations returned per column}
}
\details{
'tgs_cor' is very similar to 'package:stats::cor'. Unlike the latter it uses all available CPU cores to compute the correlation in a much faster way. The basic implementation of 'pairwise.complete.obs' is also more efficient giving overall great run-time advantage.

Unlike 'package:stats::cor' 'tgs_cor' implements only two modes of treating data containing NA, which are equivalent to 'use="everything"' and 'use="pairwise.complete.obs". Please refer the documentation of this function for more details.

'tgs_cor(x, y, spearman = F)' is equivalent to 'cor(x, y, method = "pearson")'
'tgs_cor(x, y, spearman = T)' is equivalent to 'cor(x, y, method = "spearman")'
'tgs_cor(x, y, pairwise.complete.obs = T, spearman = T)' is equivalent to 'cor(x, y, use = "pairwise.complete.obs", method = "spearman")'
'tgs_cor(x, y, pairwise.complete.obs = T, spearman = F)' is equivalent to 'cor(x, y, use = "pairwise.complete.obs", method = "pearson")'

'tgs_cor' can output its result in "tidy" format: a data frame with three columns named 'col1', 'col2' and 'cor'. Only the correlation values which abs are equal or above the 'threshold' are reported. For auto-correlation (i.e. when 'y=NULL') a pair of columns numbered X and Y is reported only if X < Y.

'tgs_cor_knn' works similarly to 'tgs_cor'. Unlike the latter it returns only the highest 'knn' correlations for each column in 'x'. The result of 'tgs_cor_knn' is always outputed in "tidy" format.

One of the reasons to opt 'tgs_cor_knn' over a pair of calls to 'tgs_cor' and 'tgs_knn' is the reduced memory consumption of the former. For auto-correlation case (i.e. 'y=NULL') given that the number of columns NC exceeds the number of rows NR, then 'tgs_cor' memory consumption becomes a factor of NCxNC. In contrast 'tgs_cor_knn' would consume in the similar scenario a factor of max(NCxNR,NCxKNN). Similarly 'tgs_cor(x,y)' would consume memory as a factor of NCXxNCY, wherever 'tgs_cor_knn(x,y,knn)' would reduce that to max((NCX+NCY)xNR,NCXxKNN).
}
\value{
'tgs_cor_knn' or 'tgs_cor' with 'tidy=TRUE' return a data frame, where each row represents correlation between two pairs of columns from 'x' and 'y' (or two columns of 'x' itself if 'y==NULL'). 'tgs_cor' with the 'tidy=FALSE' returns a matrix of correlation values, where val[X,Y] represents the correlation between columns X and Y of the input matrices (if 'y' is not 'NULL') or the correlation between columns X and Y of 'x' (if 'y' is 'NULL').
}

\examples{
\donttest{
# Note: all the available CPU cores might be used

set.seed(seed = 0)
rows <- 100
cols <- 1000
vals <- sample(1 : (rows * cols / 2), rows * cols, replace = TRUE)
m <- matrix(vals, nrow = rows, ncol = cols)
m[sample(1 : (rows * cols), rows * cols / 1000)] <- NA

r1 <- tgs_cor(m, spearman = FALSE)
r2 <- tgs_cor(m, pairwise.complete.obs = TRUE, spearman = TRUE)
r3 <- tgs_cor_knn(m, NULL, 5, spearman = FALSE)
}

\dontshow{
options(tgs_use.blas=FALSE)
options(tgs_max.processes=1)

set.seed(seed = 0)
rows <- 100
cols <- 100
vals <- sample(1 : (rows * cols / 2), rows * cols, replace = TRUE)
m <- matrix(vals, nrow = rows, ncol = cols)
m[sample(1 : (rows * cols), rows * cols / 1000)] <- NA

r1 <- tgs_cor(m, spearman = FALSE)
r2 <- tgs_cor(m, pairwise.complete.obs = TRUE, spearman = TRUE)
r3 <- tgs_cor_knn(m, NULL, 5, spearman = FALSE)
}
}
\keyword{ ~correlation }

