## Expectation propagation for logistic regression
## Assume flat prior for simplicity
#' @param y vector of response, 1 or -1
#' @param X model matrix
ep_logistic <- function(y, X) {
  d <- ncol(X)
  stopifnot(isTRUE(all(y %in% c(-1, 1))))
}
