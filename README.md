# SCPackage
## Introduction

The `SC` is a package that estimates a Cholesky factor of inverse covariance matrix via penalized maximum likelihood under local stationary assumption. In particular, given a n-by-p data matrix \eqn{X} , with each row an observation of a p dimensional random vector \eqn{X \sim N(0, \Omega^{-1}) = (L^T L)^{-1})}, this package implements a penalized likelihood-based approach of estimating \eqn{L} by smoothing its subdiagonals.
This document serves as an introduction of using the package.

The main function is `smoothchol`, which takes a sample covariance matrix of the observations and returns the estimate of \eqn{L}. 

## Installation

To install the latest version from Github, use

```s
library(devtools)
devtools::install_github("adallak/SCPackage")
```

## Usage
The package contains function `generateL` for generating true standard and modified Cholesky factor $L$. It takes as an input number of variables and number of bands and returns the Cholesky Factor. 
```s
library(SC)
set.seed(12)
p <- 50
L_true <- generateL(p = p, case = "b")$L
```

Having true Cholesky factor, we can then generate a (n x p)data matrix X with each row a random sample drawn independently from a Gaussian distribution of mean zero and covariance \eqn{\Sigma = (L^T L)^{-1}}. We use function `sample_gen` from the package `varband` to generate the data. 

```s
library(varband)
n = 100
# random sample
X <- sample_gen(L = L_true, n = n)
# sample covariance matrix
S <- crossprod(scale(X, center = TRUE, scale = FALSE)) / n
```


## Estimating L with a fixed tuning parameter

The `sc`function takes the following parameters:

* `S` - sample covariance matrix
* `lambda1` - controls the sparsity level
* `lambda_2` - controls the smoothness level
* `max_iter` - Number of iteration for block coordinate algorithm
* `init.x`   - Initial vectorized Cholesky factor \eqn{L}. If `NULL`, the default is \eqn{vec(\sqrt{diag(S)})}.
* `type` - type of the smoothing penalty
* `band` - if specified, algorithm iterates only over specifed number of subdiagonals and forces the rest of entries to be equal to zero.
* `ABSTOL` - Tolerance for algorithm convergence.

```s
L_fused = smoothchol(S, lambda1 = 0, lambda2 = 0.2, type = "fused")$L
L_l1trend = smoothchol(S, lambda1 = 0, lambda2 = 0.2, type = "l1trend")$L
L_HP = smoothchol(S, lambda1 = 0, lambda2 = 0.2, type = "HP")$L
```