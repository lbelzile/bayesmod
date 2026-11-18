library(numDeriv)

## Example of ADVI with logistic regression
## with an improper prior on betas
# Generate some fake data data
set.seed(202504)
X <- cbind(1, rnorm(100), rpois(100, lambda = 4))
beta_t <- rnorm(3)
yb <- rbinom(n = 100, size = 1, prob = c(plogis(X %*% beta_t)))
y <- ifelse(yb == 0, -1, 1)

# Logistic regression log likelihood
fn <- function(par, X, y) {
  # -sum(log1p(exp(-y * X %*% par)))
  sum(plogis(q = y * X %*% par, log.p = TRUE))
}
grad <- function(par, X, y) {
  colSums(y * plogis(q = c(-y * X %*% par)) * X)
}
# Check gradients
numDeriv::grad(func = fn, X = X, y = y, x = beta_t)
grad(par = beta_t, X = X, y = y)

hessian <- function(par, X, y) {
  f <- plogis(q = c(y * X %*% par))
  c <- sqrt(f * (1 - f))
  -crossprod(c * X)
}
# Check hessian matrix
numDeriv::hessian(func = fn, X = X, y = y, x = beta_t)
numDeriv::jacobian(func = grad, x = beta_t, X = X, y = y)
hessian(beta_t, y = y, X = X)

logistic_reg <- glm(
  formula = yb ~ X - 1,
  family = binomial
)
start <- list(
  mu = coef(logistic_reg),
  sigma = vcov(logistic_reg)
)

# Automatic differentiation variational inference
advi <- function(
  d = 1,
  start = list(mu = rep(0, d), sigma = diag(d)),
  fn,
  grad,
  hessian,
  meanfield = FALSE,
  step = 0.5,
  B = 1e3L,
  maxiter = 1e3L,
  delay = 2,
  eta_step = 0.1,
  alpha_step = 0.1,
  tol = 1e-3,
  random = TRUE,
  ...
) {
  if (missing(hessian)) {
    hessGrad <- TRUE
  } else {
    hessGrad <- FALSE
  }
  symM <- function(M) {
    0.5 * (M + t(M))
  }
  # Initialize algorithm
  mu <- start$mu
  if (isTRUE(meanfield)) {
    L <- 0.5 * log(diag(start$sigma))
    expL <- exp(L)
  } else {
    svdS <- La.svd(start$sigma)
    expL <- with(svdS, u %*% diag(sqrt(d)) %*% vt)
    # Matrix-log of covariance matrix
    L <- with(svdS, u %*% diag(0.5 * log(d)) %*% vt)
  }
  # Container for ELBO
  elbo <- rep(0, length.out = maxiter)
  cst <- 0.5 * d * (1 + log(2 * pi))
  # Stochastic gradient descent
  for (i in seq_len(maxiter)) {
    # Reinitialize containers
    grad_elbo <- rep(0, length.out = d)
    if (meanfield) {
      hessian_elbo <- rep(0, length.out = d)
    } else {
      hessian_elbo <- matrix(0, nrow = d, ncol = d)
    }
    for (b in seq_len(B)) {
      if (!isTRUE(random)) {
        set.seed(b) # make approximation deterministic
      }
      if (meanfield) {
        eta <- c(mu + expL * rnorm(d))
      } else {
        eta <- c(mu + expL %*% rnorm(d))
      }
      elbo[i] <- elbo[i] + fn(par = eta, ...)
      g <- grad(par = eta, ...)
      grad_elbo <- grad_elbo + g
      if (meanfield) {
        if (hessGrad) {
          hessian_elbo <- hessian_elbo + diag(tcrossprod(g, eta))
        } else {
          hessian_elbo <- hessian_elbo + diag(hessian(par = eta, ...))
        }
      } else {
        if (hessGrad) {
          hessian_elbo <- hessian_elbo + tcrossprod(g, eta)
        } else {
          hessian_elbo <- hessian_elbo + hessian(par = eta, ...)
        }
      }
    }
    # Renormalize entries
    g_mu <- grad_elbo / B
    if (meanfield) {
      elbo[i] <- elbo[i] / B + sum(L)
      if (hessGrad) {
        g_L <- expL * (hessian_elbo / B) + 1
      } else {
        g_L <- expL * expL * (hessian_elbo / B) + 1
      }
    } else {
      elbo[i] <- elbo[i] / B + sum(diag(L))
      if (hessGrad) {
        g_L <- symM(expL %*% (hessian_elbo / B)) + diag(d)
      } else {
        g_L <- symM(crossprod(expL) %*% (hessian_elbo / B)) + diag(d)
      }
    }
    # Adaptive step-size sequence from ADVI paper
    # Robbins-Monro sequence
    if (i == 1L) {
      s_mu <- alpha_step * g_mu^2
      s_L <- alpha_step * c(g_L)^2
    } else {
      s_mu <- alpha_step * g_mu^2 + (1 - alpha_step) * s_mu
      s_L <- alpha_step * c(g_L)^2 + (1 - alpha_step) * s_L
    }
    step <- eta_step * exp(-(0.5 + 1e-16) * log(i))
    step_mu <- step / (sqrt(s_mu) + delay)
    # The updates are done on the log scale!
    step_L <- step / (sqrt(s_L) + delay)
    # Adapt size of update
    update_mu <- step_mu * g_mu
    update_L <- step_L * c(g_L)
    # Update parameters
    mu <- mu + update_mu
    if (meanfield) {
      L <- L + update_L
      expL <- exp(L)
    } else {
      L <- L + matrix(update_L, ncol = d, nrow = d)
      ## Update exp(L)
      svdL <- La.svd(L)
      expL <- -with(svdL, u %*% diag(exp(-d)) %*% vt)
    }
    if (i >= 10L) {
      if (abs(elbo[i] - elbo[i - 1]) < tol) {
        break
      }
    }
  }
  if (meanfield) {
    expL <- diag(exp(L))
    L <- diag(L)
  }
  list(
    mean = mu,
    L = L,
    vcov = crossprod(expL),
    grad_mu = g_mu,
    grad_L = g_L,
    elbo = cst + elbo[1:i]
  )
}
# Fit the model
fit_advi <- advi(
  d = ncol(X),
  start = start,
  X = X,
  y = y,
  grad = grad,
  # meanfield = TRUE,
  hessian = hessian,
  fn = fn,
  B = 1e2L,
  tol = 1e-3,
  random = TRUE
)
# Evidence lower bound convergence
plot(
  fit_advi$elbo,
  ylab = "elbo",
  type = "l",
  xlab = "number of iterations"
)
fit_advi$mean
fit_advi$vcov


logist_stan <- cmdstanr::cmdstan_model(
  stan_file = "logistic.stan"
)
# Automatic-differentiation variational inference
mod <- logist_stan$variational(
  data = list(y = yb, X = X, p = ncol(X), N = nrow(X)),
  algorithm = "fullrank",
  # grad_samples = 1e3L,
  draws = 1e4L
)
mod$summary()
# Figure out the posterior by doing Monte Carlo with the draws
d <- mod$draws("beta")


fit_advi2 <- advi(
  d = ncol(X),
  start = list(mu = colMeans(d), sigma = cov(d)),
  X = X,
  y = y,
  grad = grad,
  # hessian = hessian,
  fn = fn,
  meanfield = TRUE,
  B = 1e2L,
  tol = 1e-3,
  random = FALSE
)
