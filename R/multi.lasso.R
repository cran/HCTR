#' @title Multi-split Lasso
#' @description Multi-splitted variable selection using Lasso
#'
#' @param Y A numeric response vector, containing nobs variables.
#' @param X An input matrix, of dimension nobs x nvars.
#' @param p.fac A sequence of penalty factor applied on each variable.
#' @param fold.num The number of cross validation folds.
#'
#' @return A list of two numeric objects of index of (1) selected and (2) unselected variables.
#'
#' @import glmnet
#' @import stats
#' @import MASS
multi.lasso <- function(X, Y, p.fac = NULL, fold.num) {
  fit <- glmnet::cv.glmnet(x=X, y=Y, family="gaussian", penalty.factor = p.fac,
                           type.measure = "mse", alpha =1, nfolds = fold.num)
  n <- nrow(X)
  p <- ncol(X)
  lambda.index <- which(fit$nzero < (n - floor(n/2)))
  lambda_hat <- fit$lambda[lambda.index[which.min(fit$cvm[lambda.index])]]
  beta.est <- glmnet::coef.cv.glmnet(fit, s = lambda_hat)

  selected.index <- beta.est@i[-1]
  if (length(selected.index) == 0) {
    unselected.index <- c(1:p)
  } else {
    unselected.index <- c(1:p)[-selected.index]
  }

  return (list("selected.index"=selected.index, "unselected.index"=unselected.index))
}
