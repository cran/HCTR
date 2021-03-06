% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/final.selection.R
\name{final.selection}
\alias{final.selection}
\title{Final Selection}
\usage{
final.selection(cv.fit, pihat, p, cov.num = 0)
}
\arguments{
\item{cv.fit}{An object of either class "cv.glmnet" from glmnet::cv.glmnet() or class "cv.ncvreg" from ncvreg::cv.ncvreg(), which is a list generated by a cross-validation fit.}

\item{pihat}{eatimated proprtion from HCTR::est.prop().}

\item{p}{Total number of variables, except for covariates.}

\item{cov.num}{Number of covariates in model, default is 0. Covariate matrix, W, is assumed on the left side of variable matrix, X. The column index of covariates are before those of variables.}
}
\value{
A sequence of index of final selected variables in the final chosen model.
}
\description{
Returns the index of final selected variables in the final chosen model.
}
\examples{
set.seed(10)
X <- matrix(rnorm(20000), nrow = 100)
beta <- rep(0, 200)
beta[1:100] <- 5
Y <- MASS::mvrnorm(n = 1, mu = X\%*\%beta, Sigma = diag(100))
fit <- glmnet::cv.glmnet(x = X, y = Y)
pihat <- 0.01
result <- est.lambda(cv.fit = fit, pihat = pihat, p = ncol(X))
lambda.seq <- seq(from = result$lambda.min, to = result$lambda.max, length.out = 100)
# Note: The lambda sequences in glmnet and ncvreg are diffrent.
fit2 <- glmnet::cv.glmnet(x = X, y = Y, lambda = lambda.seq)
result2 <- final.selection(cv.fit = fit2, pihat = 0.01, p = ncol(X))
}
