library(cubature)


## Expectation propagation for logistic regression
## Assumes flat prior for simplicity
#' @param y vector of response
#' @param X model matrix
#' @param maxiter number of iterations
#' @param tol numerical tolerance
ep_logistic <- function(y, X, maxiter = 100, tol = 1e-2) {
  p <- ncol(X)
  n <- nrow(X)
  stopifnot(length(y) == n)
  # Make sure response is binary (bool)
  yb <- y == 1
  y <- ifelse(y == 1, 1L, -1L)
  # Starting values from maximum likelihood estimator / Gaussian approx
  logist_mle <- glm(yb ~ X - 1, family = binomial)
  coefs <- as.numeric(coef(logist_mle))
  # Initialize values for parameters
  Q <- Q_cur <- solve(vcov(logist_mle))
  r <- r_cur <- c(Q %*% coefs)
  # Keep track of factors for approximations
  a <- b <- rep(0, n)

  ep_update <- function(y, mu_lc, va_lc) {
    sd_lc <- sqrt(va_lc)
    # Calculate outside of the loop the cavity
    fn <- function(x) {
      dnorm(x, mean = mu_lc, sd = sd_lc) * plogis(y * x)
    }
    browser()
    # Compute normalizing constant
    cst <- cubature::hcubature(
      f = fn,
      lowerLimit = -Inf,
      upperLimit = Inf
    )$integral
    mu <- cubature::hcubature(
      f = function(x) {
        fn(x) * x
      },
      lowerLimit = -Inf,
      upperLimit = Inf
    )$integral /
      cst
    va <- cubature::hcubature(
      f = function(x) {
        fn(x) * (x - mu)^2
      },
      -Inf,
      Inf
    )$integral /
      cst
    # Return update
    c(a = mu / va - mu_lc / va_lc, b = 1 / va - 1 / va_lc)
  }
  stopifnot(tol > 0)
  # EP updates
  for (N in seq_len(as.integer(maxiter))) {
    for (i in seq_along(y)) {
      # 1) Form the cavity by removing obs. i and update natural parameters
      r_cav <- r - c(X[i, ] * a[i])
      Q_cav <- Q - tcrossprod(X[i, ]) * b[i]
      # 2) compute mean and variance of cavity from natural parameters
      C_m <- solve(Q_cav)
      mu_m <- c(C_m %*% r_cav)
      # 3) compute mean and variance of linear projection x_i * beta
      C_p <- c(X[i, , drop = FALSE] %*% C_m %*% t(X[i, , drop = FALSE]))
      m_p <- sum(X[i, ] * mu_m)
      # 4) Obtain the mean and variance of the marginal hybrid distribution
      new <- ep_update(y = y[i], mu_lc = m_p, va_lc = C_p)
      a[i] <- new['a']
      b[i] <- new['b']
      r <- r_cav + c(a[i] * X[i, ])
      Q <- Q_cav + b[i] * tcrossprod(X[i, ])
    }
    # Check convergence by tracking parameter changes
    if (norm(Q_cur - Q, type = "F") < tol & sqrt(sum(r - r_cur)^2) < tol) {
      break
    } else {
      Q_cur <- Q
      r_cur <- r
    }
  }
  vcov <- solve(Q)
  mu <- c(vcov %*% r)
  return(list(mean = mu, vcov = vcov, niter = N))
}


# Generate some fake data from a logistic model
set.seed(202504)
X <- cbind(1, rnorm(100), rpois(100, lambda = 4))
beta_t <- rnorm(3)
yb <- rbinom(n = 100, size = 1, prob = c(plogis(X %*% beta_t)))
y <- ifelse(yb == 0, -1, 1)

# Expectation propagation approximation
ep_logistic(y = y, X = X)
