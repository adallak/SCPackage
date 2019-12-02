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

[Installation](#install)

[Quick Start](#qs)

<a id="intro"></a>

## Introduction

The `SC` is a package that estimates a Cholesky factor of inverse covariance matrix via penalized maximum likelihood under local stationary assumption. In particular, given a data matrix $X \in \mathbb{R}^{n \times p}$, with each row an observation of a $p$ dimensional random vector $X \sim N(0, \Omega^{-1} = (L^T L)^{-1})$, this package implements a penalized likelihood-based approach of estimating $L$ by smoothing its subdiagonals.
This document serves as an introduction of using the package.

The main function is `sc`, which takes a sample covariance matrix of the observations and returns the estimate of $L$. 

<a id="install"></a>

## Installation

To install the latest version from Github, use


```r
#library(devtools)
#devtools::install_github("adallak/SCPackage")
```


## Quick Start



```r
library(SC)
```

## Usage
The package contains function `generateL` for generating true standard and modified Cholesky factor $L$. It takes as an input number of variables and number of bands and returns the Cholesky Factor. 


```r
set.seed(12)
p <- 50
band <- 5
L_true <- generateL(p = p, band = band)
```

Having true Cholesky factor, we can then generate a data matrix $X \in \mathbb{R}^{n \times p}$ with each row a random sample drawn independently from a Gaussian distribution of mean zero and covariance $\Sigma = (L^T L)^{-1}$. We use function `sample_gen` from the package `varband` to generate the data. After centering the dat, the sample covariance is esitmated by $S = \frac{X^tX}{n}$






