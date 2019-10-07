#' @title Bounding Sequence
#' @description Calculates bounding sequence of higher crticism for proportion estimator using p-values
#'
#' @param p.value A matrix of p-values from permutation: row is from each permutation; column is from each variable.
#' @param alpha Probability of Type I error for bounding sequence, the default value is 1 / sqrt(log(p)), where p is number of p-values in each permutation.
#'
#' @return A bounding value of higher criticism with (1 - alpha) confidence.
#'
#' @export
#' @import stats
#' @importFrom Rdpack reprompt
#'
#' @examples
#' set.seed(10)
#' X <- matrix(runif(n = 10000, min = 0, max = 1), nrow = 100)
#' result <- bounding.seq(p.value = X)
#'
#' @references
#' \insertRef{jeng2019efficient}{HCTR}
bounding.seq <- function(p.value, alpha) {
  p <- ncol(p.value)
  perm.num <- nrow(p.value)

  if (missing(alpha)) {
    alpha <- (1 / sqrt(log(p)))
  }

  Vn <- rep(NA, perm.num)
  for(i in 1 : perm.num) {
    sortp <- sort(p.value[i, ], index.return = TRUE)$x
    ind <- 1 : p
    Un <- rep(0, p)
    Un <- (ind / p - sortp) / sqrt(sortp * (1 - sortp))
    Vn[i] <- max(Un[2 : floor(p / 2)])
  }

  c_n <- quantile(Vn, 1 - alpha)
  return (c_n)
}
