#' @title Permutation p-values
#' @description Calculates
#'
#' @param Y A numeric response vector, containing nobs variables.
#' @param X An input matrix, of dimension nobs x nvars.
#' @param W A covariate matrix, of dimension nobs x ncors, default is NULL.
#' @param type Penalized regression type, valid parameters include "Lasso", "AdaLasso", "SCAD", and "MCP".
#' @param B Multi-split times, default is 100.
#' @param fold.num The number of cross validation folds, default is 10.
#' @param perm.num Permutation times, default is 1000.
#'
#' @return A matrix containing harmonic mean p-values from permutation.
#'
#' @export
#'
#' @examples
#' set.seed(10)
#' X <- matrix(rnorm(20000), nrow = 100)
#' beta <- rep(0, 200)
#' beta[1:100] <- 5
#' Y <- MASS::mvrnorm(n = 1, mu = X%*%beta, Sigma = diag(100))
#' result <- pmpv(Y=Y, X=X, type = "Lasso", B = 2, fold.num = 10, perm.num = 10)
pmpv <- function(Y, X, W = NULL, type, B = 100, fold.num = 10, perm.num = 1000) {
  p <- ncol(X)
  p.noise <- matrix(rep(NA, perm.num*p), ncol = p)

  for (i in 1:perm.num) {
    Y_perm <- sample(Y)
    p.noise[i,] <- highdim.p(Y = Y_perm, X = X, W = W, type = type, B = B, fold.num = fold.num)$p.hmp
  }

  return(p.noise)
}

