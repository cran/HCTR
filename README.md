<!-- README.md is generated from README.Rmd. Please edit that file -->

HCTR
====

The goal of HCTR is to create a new searching scheme for the
regularization parameter in penalized regression, such as Lasso,
adaptive Lasso, SCAD, and MCP.

Example
-------

This is a basic example which shows you how to (1) estimate false null
hypothesis proportion, and (2) create a new tuning region for the
regularization parameter.

``` r
## basic example code
library('HCTR')
# 1. Estimate proportion
set.seed(10)
X <- matrix(runif(n = 10000, min = 0, max = 1), nrow = 100)
result <- bounding.seq(p.value = X)
Y <- matrix(runif(n = 100, min = 0, max = 1), nrow = 100)
test <- est.prop(p.value = Y, cn = result)
# 2. Estimate a new tuning region
set.seed(10)
X <- matrix(rnorm(20000), nrow = 100)
beta <- rep(0, 200)
beta[1:100] <- 5
Y <- MASS::mvrnorm(n = 1, mu = X%*%beta, Sigma = diag(100))
fit <- glmnet::cv.glmnet(x = X, y = Y)
pihat <- 0.01
result <- est.lambda(cv.fit = fit, pihat = pihat, p = ncol(X))
```
