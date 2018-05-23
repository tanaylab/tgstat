\name{tgs_finite}
\alias{tgs_finite}
\title{
Checks whether all the elements of the vector are finite
}
\description{
Checks whether all the elements of the vector are finite.
}
\usage{
tgs_finite(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{numeric or integer vector or matrix}
}
\details{
'tgs_finite' returns 'TRUE' if all the elements of 'x' are finite numbers. (See: 'is.finite'.)
}
\value{
'TRUE' if all the elements of 'x' are finite, otherwise 'FALSE'.
}

\examples{
tgs_finite(1:10)
tgs_finite(c(1:10, NaN))
tgs_finite(c(1:10, Inf))
}
\keyword{ ~finite }
