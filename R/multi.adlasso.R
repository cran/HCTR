#' @title Multi-split Adaptive Lasso
#' @description Multi-splitted variable selection using Adaptive Lasso
#'
#' @param Y A numeric response vector, containing nobs variables.
#' @param X An input matrix, of dimension nobs x nvars.
#' @param covar.num Number of covariates in model, default is NULL.
#' @param fold.num The number of cross validation folds.
#'
#' @return A list of two numeric objects of index of (1) selected and (2) unselected variables.
#'
#' @import glmnet
#' @import stats
#' @import MASS
#'
multi.adlasso <- function(X, Y, covar.num = NULL, fold.num) {
  ridge_cv <- glmnet::cv.glmnet(x=X, y=Y, family="gaussian",
                                type.measure = "mse", alpha =0, nfolds = fold.num)
  best_ridge_coef <- as.numeric(coef(ridge_cv, s = ridge_cv$lambda.min))[-1]

  ad.p.fac <- 1 / abs(best_ridge_coef)
  if (!is.null(covar.num)) {
    ad.p.fac[c(1:covar.num)] <- 0
  }

  adlasso_cv <- glmnet::cv.glmnet(x=X, y=Y, family="gaussian",
                                  type.measure = "mse", alpha = 1, nfolds = fold.num,
                                  penalty.factor = ad.p.fac, keep = TRUE)

  n <- nrow(X)
  p <- ncol(X)
  lambda.index <- which(adlasso_cv$nzero < (n - floor(n/2)))
  lambda_hat <- adlasso_cv$lambda[lambda.index[which.min(adlasso_cv$cvm[lambda.index])]]
  beta.est <- stats::coef(adlasso_cv, s = lambda_hat)

  selected.index <- beta.est@i[-1]
  if (length(selected.index) == 0) {
    unselected.index <- c(1:p)
  } else {
    unselected.index <- c(1:p)[-selected.index]
  }

  return (list("selected.index"=selected.index, "unselected.index"=unselected.index))
}
