## Expectation propagation for logistic regression
## Assumes flat prior for simplicity
#' @param y vector of response, 1 or -1
#' @param X model matrix
ep_logistic <- function(y, X) {
  d <- ncol(X)
  stopifnot(isTRUE(all(y %in% c(-1, 1))))
  plogis(y * X %*% beta)
  for (i in seq_along(y)) {
    # 1) Form the cavity by removing obs. i and update natural parameters
    # 2) compute mean and variance of cavity from natural parameters
    # 3) compute mean and variance of linear projection x_i * beta
    # 4) Obtain the mean and variance of the marginal hybrid distribution
  }
}
