#' @title p-values in high-dimensional linear model
#' @description Calculates p-values in high-dimentional linear models using multi-split method
#'
#' @param Y A numeric response vector, containing nobs variables.
#' @param X An input matrix, of dimension nobs x nvars.
#' @param W A covariate matrix, of dimension nobs x ncors, default is NULL.
#' @param type Penalized regression type, valid parameters include "Lasso", "AdaLasso", "SCAD", and "MCP".
#' @param B Multi-split times, default is 100.
#' @param fold.num The number of cross validation folds.
#'
#' @return A list of objects containing: (1) harmonic mean p-values; (2) original p-values; (3) index of selected samples; (4) index of selected variables
#'
#' @export
#' @import glmnet
#' @import ncvreg
#' @import harmonicmeanp
#' @import stats
#' @import MASS
#'
#' @examples
#' set.seed(10)
#' X <- matrix(rnorm(20000), nrow = 100)
#' beta <- rep(0, 200)
#' beta[1:100] <- 5
#' Y <- MASS::mvrnorm(n = 1, mu = X%*%beta, Sigma = diag(100))
#' result <- highdim.p(Y=Y, X=X, type = "Lasso", B = 2, fold.num = 10)
highdim.p <- function(Y, X, W = NULL, type, B = 100, fold.num) {
  X <- cbind(W, X)
  n <- nrow(X)
  p <- ncol(X)
  covar.num <- ncol(W)
  p.fac <- rep(1, p)
  if (!is.null(covar.num)) {
    p.fac[c(1:covar.num)] <- 0 # Set penalty term for covariates
  }

  p.values <- matrix(0, nrow = B, ncol = p)
  sample.list <- matrix(0, nrow = B, ncol = floor(n/2))
  selected.list <- list()

  i <- 1
  while (i <= B){
    index <- sort(sample(x = 1:n, size = floor(n/2)))
    sample.list[i,] <- index
    Y_1 <- Y[index]
    X_1 <- X[index,]
    Y_2 <- Y[-index]
    X_2 <- X[-index,]

    if (type == "Lasso") {
      final.index <- multi.lasso(X = X_1, Y = Y_1, p.fac = p.fac, fold.num = fold.num)
    } else if (type == "AdaLasso") {
      final.index <- multi.adlasso(X = X_1, Y = Y_1, covar.num = covar.num, fold.num = fold.num)
    } else if (type == "SCAD") {
      final.index <- multi.scad(X = X_1, Y = Y_1, p.fac = p.fac, fold.num = fold.num)
    } else if (type == "MCP") {
      final.index <- multi.mcp(X = X_1, Y = Y_1, p.fac = p.fac, fold.num = fold.num)
    } else {
      stop("Invalid regression type!")
    }

    selected.index <- final.index$selected.index
    unselected.index <- final.index$unselected.index

    selected.list[[i]] <- selected.index
    X_2 <- X_2[,selected.index]

    if ((!is.null(ncol(X_2)))&(NCOL(X_2)>0)) {
      fit_2 <- lm(Y_2 ~ X_2)
    } else {
      fit_2 <- lm(Y_2 ~ 1)
    }
    tmp.p.values <- summary(fit_2)$coefficients[-1,4]

    if (length(tmp.p.values)==length(selected.index)) {
      p.values[i,selected.index] <- tmp.p.values
      p.values[i,unselected.index] <- runif(n = (p-length(selected.index)), min = 0, max = 1)
      i <- (i+1)
    }
  }
  hmp <- function(p.values) { return(harmonicmeanp::p.hmp(p.values, L=B)) }
  p.hmp <- apply(p.values, 2, hmp)

  return (list("p.hmp"=p.hmp, "p.values"=p.values, "sample.list"=sample.list, "selected.list"=selected.list))
}
