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

}
\author{
\packageAuthor{SC}

Maintainer: \packageMaintainer{SC}
}
\references{
Dallakyan A. and Pourahmadi M. - Learning Smooth Cholesky Factors of Large Covariance Matrices (2019)
}
\keyword{ package }
\examples{
n = 100
p = 50
band = 5
L_true = generateL(p = p, band = band, type = "b")$L
# random sample
X <- sample_gen(L = L_true, n = n)
# sample covariance matrix
S <- crossprod(scale(X, center = TRUE, scale = FALSE)) / n
# use identity matrix as an initial estimate
L_fused = smoothchol(S, lambda1 = 0, lambda2 = 0.2, type = "fused", band = band)
}
