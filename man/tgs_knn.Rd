\name{tgs_knn}
\alias{tgs_knn}
\title{
Returns k highest values of each column
}
\description{
Returns k highest values of each column.
}
\usage{
tgs_knn(x, knn, diag = F, threshold = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{numeric matrix or data frame (see below)}
  \item{knn}{the number of highest values returned per column}
  \item{diag}{if 'F' values of row 'i' and col 'j' are skipped for each i == j}
  \item{threshold}{filter out values lower than threshold}
}
\details{
'tgs_knn' returns the highest 'knn' values of each column of 'x' (if 'x' is a matrix). 'x' can be also a sparse matrix given in a data frame of 'col', 'row', 'value' format.

'NA' and 'Inf' values are skipped as well as the values below 'threshold'. If 'diag' is 'F' values of the diagonal (row == col) are skipped too.
}
\value{
A sparse matrix in a data frame format with 'col1', 'col2', 'val' and 'rank' columns. 'col1' and 'col2' represent the column and the row number of 'x'.
}

\examples{
\donttest{
# Note: all the available CPU cores might be used

set.seed(seed = 1)
rows <- 100
cols <- 1000
vals <- sample(1 : (rows * cols / 2), rows * cols, replace = TRUE)
m <- matrix(vals, nrow = rows, ncol = cols)
m[sample(1 : (rows * cols), rows * cols / 1000)] <- NA
r <- tgs_knn(m, 3)
}

\dontshow{
options(tgs_use.blas=FALSE)
options(tgs_max.processes=1)

set.seed(seed = 1)
rows <- 100
cols <- 100
vals <- sample(1 : (rows * cols / 2), rows * cols, replace = TRUE)
m <- matrix(vals, nrow = rows, ncol = cols)
m[sample(1 : (rows * cols), rows * cols / 1000)] <- NA
r <- tgs_knn(m, 3)
}
}
\keyword{ ~knn }

