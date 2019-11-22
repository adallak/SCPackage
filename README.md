# SCPackage
## Introduction
The `SC` is a package that estimates a Cholesky factor of inverse covariance matrix via penalized maximum likelihood under local stationary assumption. In particular, given a data matrix $X \in \mathbb{R}^{n \times p}$, with each row an observation of a $p$ dimensional random vector $X \sim N(0, \Omega^{-1} = (L^T L)^{-1})$, this package implements a penalized likelihood-based approach of estimating $L$ by smoothing its subdiagonals.
This document serves as an introduction of using the package.

The main function is `sc`, which takes a sample covariance matrix of the observations and returns the estimate of $L$. 

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
band <- 5
L_true <- generateL(p = p, band = band)
```

Having true Cholesky factor, we can then generate a data matrix $X \in \mathbb{R}^{n \times p}$ with each row a random sample drawn independently from a Gaussian distribution of mean zero and covariance $\Sigma = (L^T L)^{-1}$. We use function `sample_gen` from the package `varband` to generate the data. 

```s
library(varband)
n = 100
# random sample
X <- sample_gen(L = true, n = n)
# sample covariance matrix
S <- crossprod(scale(X, center = TRUE, scale = FALSE)) / n
```


## Estimating L with a fixed tuning parameter

The `sc`function takes the following parameters:
- `S` - sample covariance matrix
- `lambda1` - controls the sparsity level
- `lambda_2` - controls the smoothness
- `max_iter` - Number of iteration for block coordinate algorithm
- `init.x`   - Initial vectorized Cholesky vactor $L$. If `NULL`, the default is $vec(\sqrt{diag(S)})$.
- `type` - type of the smoothing penalty
- `band` - if specified, algorithm forces the rest of entries zero and iterates only over specifed subdiagonals.
- `ABSTOL` - Tolerance for algorithm convergence.

```s
# use identity matrix as initial estimate
L_fused = sc(S, lambda1 = 0, lambda2 = 0.2, type = "fused")
```