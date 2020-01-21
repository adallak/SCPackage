---
title: "sc_vignette"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{sc_vignette}
  %\VignetteEngine{knitr::knitr}
  %\VignetteEncoding{UTF-8}
---

<a id="top"></a>

# SC Vignette
### Aramayis Dallakyan

[Introduction](#intro)

<!-- [Installation](#install) -->


[Quick Start](#qs)

[Using Cross Validation](#cv)

<a id="intro"></a>

## Introduction

The `SC` is a package that estimates a Cholesky factor of inverse covariance matrix via penalized maximum likelihood under local stationary assumption. In particular, given a data matrix $X \in \mathbb{R}^{n \times p}$, with each row an observation of a $p$ dimensional random vector $X \sim N(0, \Omega^{-1} = (L^T L)^{-1})$, this package implements a penalized likelihood-based approach of estimating $L$ by smoothing its subdiagonals.
This document serves as an introduction of using the package.

The main function is `smoothchol`, which takes a sample covariance matrix of the observations and returns the estimate of $L$. 

<a id="qs"></a>

## Quick Start



```r
library(SC)
```

## Usage
The package contains function `generateL` for generating true standard and modified Cholesky factor $L$. The  It takes as an input number of variables and number of bands and returns the Cholesky Factor. The type function 


```r
set.seed(12)
p <- 30
L_true <- generateL(p = p, case = "c")$L
```

Having true Cholesky factor, we can then generate a data matrix $X \in \mathbb{R}^{n \times p}$ with each row a random sample drawn independently from a Gaussian distribution of mean zero and covariance $\Sigma = (L^T L)^{-1}$. We use function `sample_gen` from the package `varband` to generate the data. After centering the dat, the sample covariance is esitmated by $S = \frac{X^tX}{n}$


```r
library(varband)
n = 100
# random sample
X <- sample_gen(L = L_true, n = n)
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


```r
L_fused = smoothchol(S, lambda1 = 0, lambda2 = 0.5 , type = "fused", max_iter = 150)
```

<a id="cv"></a>

## Using Cross Validation

In this section we use `smoothcholCV` to select the tuning parameter for the proposed method using cross validation and plot the true and esitmated first subdiagonals.


```r
lambda2_seq = seq(0.01, 2, length.out = 60)
L_fused_cv = smoothcholCV(k = 5, X, lambda2_seq = lambda2_seq, pen.type = "fused", max_iter = 150, stand = TRUE )
L_trend_cv = smoothcholCV(k = 5, X, lambda2_seq = lambda2_seq, pen.type = "l1trend", max_iter = 150, stand = TRUE )
L_hp_cv = smoothcholCV(k = 5, X, lambda2_seq = lambda2_seq, pen.type = "HP", max_iter = 150, stan = TRUE )
## Converting L into T as described in Dallakyan and Pourahmadi (2019)
T_true = - L_true*(1/diag(L_true))
T_fused = - L_fused_cv$L_fit * (1 / diag(L_fused_cv$L_fit))
T_trend = - L_trend_cv$L_fit *(1 / diag(L_trend_cv$L_fit))
T_hp = - L_hp_cv$L_fit *(1/diag(L_fused_cv$L_fit))
## Extracting first subdiagonal
true_first = diag(T_true[-1, -p])
fused_first = diag(T_fused[-1, -p])
trend_first = diag(T_trend[-1, -p])
hp_first = diag(T_hp[-1, -p])
```

Now, to evaluate performance of our method we plot the true first subdiagonal and the estimated subdiagonal for all three estimators.

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png)

We can plot the true model and the estimated matrix by using `matimage` function from the package `varband`. 


```r
  matimage(L_true, main = "True L")
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png)

```r
  matimage(L_fused_cv$L_fit, main = "Smoothed Cholesky (Fused)")
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-2.png)

```r
  matimage(L_trend_cv$L_fit, main = "Smoothed Cholesky (l1trent)")
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-3.png)

```r
  matimage(L_hp_cv$L_fit, main = "Smoothed Cholesky (HP)")
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-4.png)
