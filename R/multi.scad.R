#' @title Multi-split SCAD
#' @description Multi-splitted variable selection using SCAD
#'
#' @param Y A numeric response vector, containing nobs variables.
#' @param X An input matrix, of dimension nobs x nvars.
#' @param p.fac A sequence of penalty factor applied on each variable.
#' @param fold.num The number of cross validation folds.
#'
#' @return A list of two numeric objects of index of (1) selected and (2) unselected variables.
#'
#' @import ncvreg
#' @import stats
#' @import MASS
#'
multi.scad <- function(X, Y, p.fac = NULL, fold.num) {
  if (is.null(p.fac)) {
    p.fac <- rep(1, ncol(X))
  }
  scad_cv <- ncvreg::cv.ncvreg(X = X, y = Y, penalty = "SCAD", penalty.factor = p.fac,
                               family = "gaussian", gamma = 3.7, nfolds = fold.num)

  n <- nrow(X)
  nzero <- (colSums(scad_cv$fit$beta != 0) - 1)
  lambda.index <- which(nzero < (n - floor(n/2)))
  lambda_hat <- scad_cv$lambda[lambda.index[which.min(scad_cv$cve[lambda.index])]]
  beta.est <- coef(scad_cv, lambda=lambda_hat)

  selected.index <- which(beta.est!=0,arr.ind = T)[-1]-1
  unselected.index <- which(beta.est==0,arr.ind = T)-1

  return (list("selected.index"=selected.index, "unselected.index"=unselected.index))
}
