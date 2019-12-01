\name{SC-package}
\alias{SC-package}
\alias{SC}
\docType{package}
\title{
\packageTitle{SC}
}
\description{
\packageDescription{SC}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{SC}
\packageIndices{SC}
~~ An overview of how to use the package, including the most ~~
~~ important functions ~~
}
\author{
\packageAuthor{SC}

Maintainer: \packageMaintainer{SC}
}
\references{
Dallakyan A. and Pourahmadi M. - Learning Smooth CHolesky Factors of Large Covariance Matrices (2019)
}
\keyword{ package }
\examples{
n = 100
# random sample
X <- sample_gen(L = true, n = n)
# sample covariance matrix
S <- crossprod(scale(X, center = TRUE, scale = FALSE)) / n
# use identity matrix as initial estimate
L_fused = smoothchol(S, lambda1 = 0, lambda2 = 0.2, type = "fused")
}
