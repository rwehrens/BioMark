\name{selection}
\alias{selection}
\title{Accessor function to the selected variables of a BioMark object}
\description{Convenience function to get the indices of the selection in
  a BioMark object.}
\usage{selection(object, ...)}
\arguments{
  \item{object}{An object of class \code{BioMark}.}
  \item{\dots}{Further arguments, currently ignored.}
}
\value{A vector containing the indices of the selected variables.}
\seealso{\code{\link{get.biom}}}
\examples{
data(spikedApples)
X <- sqrt(spikedApples$dataMatrix)
Y <- factor(rep(0:1, each = 10))

## stability-based selection
set.seed(17)
apple.stab <- get.biom(X = X, Y = Y,, ncomp = 2:3, type = "stab")
selected.variables <- selection(apple.stab)
}
\keyword{models}
\keyword{multivariate}
