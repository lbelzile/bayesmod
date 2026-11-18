library(ggplot2)
library(ggdensity)
library(patchwork)
theme_set(theme_classic())
#############################################################
## CAVI algorithm for Gaussian data with conjugate priors
#############################################################
## Simulate fake data
set.seed(1234)
y <- rnorm(n = 20, mean = 150, sd = sqrt(20))
n <- length(y)
a0 <- 0.01
b0 <- 0.01
mu0 <- 0
tau0 <- 1e-4

CAVI_gauss <- function(
  y,
  init = c(mean(y), length(y) / var(y), 1 / var(y)),
  tol = 1e-4
) {
  a0 <- 0.01
  b0 <- 0.01
  mu0 <- 0
  tau0 <- 1e-4
  n <- length(y)
  sum_y <- sum(y)
  E_mu <- init[1]
  var_mu <- init[2]
  E_tau <- init[3]
  B <- 20
  an <- (a0 + n / 2 - 1)
  elbo <- numeric(B)
  lcst <- a0 * log(b0) - lgamma(a0) + 0.5 * (1 + log(tau0) - n * log(2 * pi))
  for (i in 1:B) {
    var_mu <- 1 / (E_tau * (tau0 + n))
    E_mu <- (tau0 * mu0 + sum_y) / (tau0 + n)
    bn <- b0 +
      0.5 * (sum((y - E_mu)^2) + n * var_mu) +
      0.5 * tau0 * ((E_mu - mu0)^2 + var_mu)
    E_tau <- an / bn
    elbo[i] <- lcst - an * log(bn) + lgamma(an) + 0.5 * log(var_mu)
  }
  list(
    elbo = elbo,
    mu_mean = E_mu,
    mu_var = var_mu,
    tau_shape = an,
    tau_rate = bn
  )
}

CAVI_approx <- CAVI_gauss(y = y)
ggplot(
  data = data.frame(elbo = CAVI_approx$elbo),
  mapping = aes(x = seq_along(elbo), y = elbo)
) +
  geom_line() +
  geom_point() +
  labs(y = "elbo", x = "iteration number")


gt <- function(x, y, pars) {
  dnorm(x, mean = pars$mu_mean, sd = sqrt(pars$mu_var)) *
    dgamma(y, rate = pars$tau_rate, shape = pars$tau_shape)
}
pars_f <- list(
  mu_mean = (sum(y) + tau0 * mu0) / (n + tau0),
  mu_var = 1 / (n + tau0),
  tau_rate = b0 +
    0.5 * (sum(y^2) + tau0 * mu0^2 - (sum(y) + tau0 * mu0)^2 / (n + tau0)),
  tau_shape = 0.5 * (n - 1) + a0
)
ft <- function(x, y, pars) {
  dnorm(x, mean = pars$mu_mean, sd = sqrt(pars$mu_var) / sqrt(y)) *
    dgamma(y, rate = pars$tau_rate, shape = pars$tau_shape)
}

nsim <- 1e4L
tau_s <- rgamma(n = nsim, rate = pars_f$tau_rate, shape = pars_f$tau_shape)
mu_s <- rnorm(
  n = nsim,
  mean = pars_f$mu_mean,
  sd = sqrt(pars_f$mu_var) / sqrt(tau_s)
)

g1 <- ggplot() +
  ggdensity::geom_hdr_fun(
    fun = ft,
    probs = c(0.1, 0.25, 0.5, 0.75, 0.9),
    xlim = range(mu_s),
    ylim = range(tau_s),
    n = 1001L,
    args = list(pars = pars_f)
  ) +
  labs(x = expression(mu), y = expression(tau))
g2 <- ggplot() +
  geom_point(
    data = data.frame(tau = 1 / tau_s, mu = mu_s),
    mapping = aes(x = mu, y = tau)
  ) +
  ggdensity::geom_hdr_fun(
    fun = gt,
    probs = c(0.1, 0.25, 0.5, 0.75, 0.9),
    xlim = range(mu_s),
    ylim = range(tau_s),
    n = 1001L,
    args = list(pars = CAVI_approx)
  ) +
  labs(x = expression(mu), y = expression(tau))
g1 + g2 + plot_layout(guides = 'collect') & theme(legend.position = "bottom")

# Marginal densities
mu_loc <- (tau0 * mu0 + sum(y)) / (tau0 + n)
mu_scale <- sqrt(
  (2 * b0 + (n - 1) * var(y) + tau0 * n * (mean(y) - mu0)^2 / (tau0 + n)) /
    ((tau0 + n) * (n + 2 * a0))
)
mu_df <- 2 * a0 + n

g3 <- ggplot() +
  stat_function(
    fun = function(x) {
      dt((x - mu_loc) / mu_scale, df = mu_df)
    },
    xlim = range(mu_s)
  ) +
  stat_function(
    fun = dnorm,
    args = list(
      mean = CAVI_approx$mu_mean,
      sd = sqrt(CAVI_approx$mu_var)
    ),
    xlim = range(mu_s),
    linetype = "dashed"
  ) +
  labs(x = expression(mu), y = "density")
g4 <- ggplot() +
  stat_function(
    fun = dgamma,
    args = list(rate = pars_f$tau_rate, shape = pars_f$tau_shape),
    xlim = range(tau_s)
  ) +
  stat_function(
    fun = dgamma,
    args = list(
      rate = CAVI_approx$tau_rate,
      shape = CAVI_approx$tau_shape
    ),
    xlim = range(tau_s),
    linetype = "dashed"
  ) +
  labs(x = expression(tau), y = "density")
g3 + g4


#############################################################
## CAVI algorithm for probit regression
## example with data augmentation
#############################################################

# Data augmentation: "a" is latent Gaussian, truncated to be negative for failures and positive for successes
cavi_probit <- function(
  y, # response vector (0/1)
  X, # model matrix
  prior_beta_prec = diag(rep(0.01, ncol(X))),
  prior_beta_mean = rep(0, ncol(X)),
  maxiter = 1000L,
  tol = 1e-4
) {
  # Precompute fixed quantity
  sc <- solve(crossprod(X) + prior_beta_prec)
  pmu_prior <- prior_beta_prec %*% prior_beta_mean
  n <- length(y) # number of observations
  stopifnot(nrow(X) == n)
  mu_a <- rep(0, n)
  y <- as.logical(y)
  ELBO <- numeric(maxiter)
  lcst <- -0.5 *
    log(det(solve(prior_beta_prec) %*% crossprod(X) + diag(ncol(X))))
  for (b in seq_len(maxiter)) {
    mu_b <- c(sc %*% (t(X) %*% mu_a + pmu_prior))
    lp <- c(X %*% mu_b)
    mu_a <- lp + dnorm(lp) / (pnorm(lp) - ifelse(y, 0, 1))
    ELBO[b] <- sum(pnorm(lp, lower.tail = y, log.p = TRUE)) -
      0.5 *
        c(
          t(mu_b - prior_beta_mean) %*%
            prior_beta_prec %*%
            (mu_b - prior_beta_mean)
        ) -
      lcst
    if (b > 2 && (ELBO[b] - ELBO[b - 1]) < tol) {
      break
    }
  }
  list(mu_a = mu_a, mu_beta = mu_b, elbo = ELBO[1:b])
}

# Example with data from Experiment 2 of Duke and Amir (2023)
# on the effect of sequential decisions and purchasing formats
data(DA23_E2, package = "hecedsm")
X <- model.matrix(
  ~ scale(age) + format,
  data = DA23_E2,
  contrasts.arg = list(format = "contr.sum")
)
y <- DA23_E2$purchased
# Fit the probit model via coordinate-ascent variational inference
cavi_probit_DA <- cavi_probit(y = y, X = X)
# Compare posterior mean with frequentist estimates
cbind(cavi_probit_DA$mu_beta, coef(glm(y ~ X - 1, family = binomial("probit"))))
plot(
  cavi_probit_DA$elbo,
  type = "b",
  ylab = "ELBO",
  xlab = "number of iterations"
)


#############################################################
## CAVI algorithm for K-components Gaussian mixture models
#############################################################
# Parameters are
# 1) weights w (probability of components)
# 2) cluster means mu
# 3) cluster variance sigma_sq
# 4) binary indicators of clusters "a" (data augmentation)
CAVI_gauss_mixt <- function(
  K = 2L,
  y,
  prior_mean_mu = rep(0, K),
  prior_var_mu = rep(1e8, K),
  prior_shape_var = rep(0.01, K),
  prior_rate_var = rep(0.01, K),
  prior_shape_weights = 0.01,
  niter = 1000L
) {
  stopifnot(K >= 1)
  K <- as.integer(K)
  n <- length(y)
  # ELBO normalizing constant (only depends on hyperparameters)
  lcst <- 0.5 *
    K *
    (1 - n * log(2 * pi)) +
    lgamma(K * prior_shape_weights) -
    K * lgamma(prior_shape_weights) -
    lgamma(n + K * prior_shape_weights) +
    sum(prior_shape_var * log(prior_rate_var)) -
    sum(lgamma(prior_shape_var))
  # Initialization
  mu <- runif(K) * diff(range(y))
  sigma_sq <- alpha <- rep(1, K)
  A <- B <- rep(1, K) # K vector
  ELBO <- numeric(niter)
  nu <- matrix(NA, nrow = n, ncol = K)
  # CAVI runs
  for (b in seq_len(niter)) {
    for (k in seq_len(K)) {
      nu[, k] <- exp(
        digamma(alpha[k]) +
          0.5 * digamma(A[k]) -
          0.5 * log(B[k]) -
          0.5 * A[k] / B[k] * ((y - mu[k])^2 + sigma_sq[k])
      )
    }
    omega <- nu / rowSums(nu) # Probability of components for each obs
    om <- colSums(omega) # sample size in each cluster
    sigma_sq <- 1 / (1 / prior_var_mu + A * om / B) # posterior variance of cluster
    mu <- sigma_sq * (prior_mean_mu / prior_var_mu + A / B * colSums(omega * y)) # posterior mean of cluster
    alpha <- prior_shape_weights + om
    A <- prior_shape_var + 0.5 * om
    B <- prior_rate_var +
      0.5 * (colSums(omega * (outer(y, mu, FUN = "-")^2)) + om * sigma_sq)
    # Compute ELBO
    ELBO[b] <- lcst +
      sum(
        lgamma(A) -
          A * log(B) +
          lgamma(alpha) +
          0.5 * (log(sigma_sq) - log(prior_var_mu)) -
          0.5 * ((mu - prior_mean_mu)^2 + sigma_sq) / prior_var_mu
      ) -
      sum(omega * log(omega + 1e-80))
  }
  list(
    elbo = ELBO,
    mu_mean = mu,
    mu_var = sigma_sq,
    sigmasq_shape = A,
    sigmasq_rate = B,
    weights_alpha = alpha,
    probs = omega,
    mean_probs = alpha / sum(alpha),
    cluster_probs = colSums(omega) / sum(omega)
  )
}
## Data application
data(geyser, package = "MASS")
mixt <- CAVI_gauss_mixt(K = 2, y = geyser$duration)
plot(
  mixt$elbo,
  xlab = "number of iterations",
  ylab = "evidence lower bound",
  type = "b"
)
# Compute posterior mean of mean, variance and probability of each cluster
post_prob <- mixt$mean_probs
post_mean <- mixt$mu_mean
post_sd <- sqrt(mixt$sigmasq_rate / (mixt$sigmasq_shape - 1))
# Plot the mixture density at posterior mean
curve(
  post_prob[1] *
    dnorm(x, mean = post_mean[1], sd = post_sd[1]) +
    post_prob[2] * dnorm(x, mean = post_mean[2], sd = post_sd[2]),
  from = 1,
  to = 6,
  n = 1001,
  xlab = "duration of geyser eruption (in minutes)",
  ylab = "density"
)
rug(geyser$duration)


#############################################################
## CAVI algorithm for linear mixed effect model
#############################################################
## Considers a data augmentation with random effects u's
## and a factorizations with variance vs (mean + random effects)
## Assuming a diagonal structure (independent) on random effects
## - multivariate Gaussian for (beta, u)
## - independent inverse. gamma for each variance component
# Y ~ Gauss(X beta + Z u, sigma_sq_eps I)
# u ~ Gauss(0, blockdiag(sigma_sq_u))

cavi_varcomp_linmixmod <- function(
  y,
  Xform,
  Zform,
  data,
  maxiter = 1e3L,
  tol = 1e-8,
  prior_beta_mean,
  prior_beta_prec,
  prior_var_shape,
  prior_var_rate
) {
  stopifnot(is.matrix(X))
  y <- as.numeric(y)
  n <- length(y)
  X <- model.matrix(as.formula(Xform), data = data)
  mm <- model.frame(as.formula(Zform), data = data)
  ym <- model.response(mm)
  r <- ncol(mm) - ifelse(is.null(ym), 0, 1)
  Z <- Matrix::sparse.model.matrix(
    as.formula(paste(c(Zform, "-1"), collapse = "")),
    data = data
  )
  K <- apply(mm, 2, function(x) {
    length(unique(x))
  })
  stopifnot(ncol(Z) == sum(K))
  C <- cbind(X, Z)
  stopifnot(nrow(X) == n, nrow(Z) == n)
  p <- ncol(X)
  if (missing(prior_beta_mean)) {
    prior_beta_mean <- rep(0, p)
  } else {
    stopifnot(length(prior_beta_mean) == p)
  }
  if (missing(prior_beta_prec)) {
    prior_beta_prec <- Matrix::Diagonal(n = p, x = 0.01)
  } else {
    stopifnot(
      is.matrix(prior_beta_prec),
      ncol(prior_beta_prec) == p,
      nrow(prior_beta_prec) == p,
      isSymmetric(prior_beta_prec)
    )
  }
  if (missing(prior_var_shape)) {
    prior_var_shape <- rep(0.01, r + 1)
  } else {
    stopifnot(
      is.numeric(prior_var_shape),
      isTRUE(all(prior_var_shape > 0)),
      length(prior_var_shape) == (r + 1)
    )
  }
  if (missing(prior_var_rate)) {
    prior_var_rate <- rep(0.01, r + 1)
  } else {
    stopifnot(
      is.numeric(prior_var_rate),
      isTRUE(all(prior_var_rate > 0)),
      length(prior_var_rate) == (r + 1)
    )
  }
  # Check for variance components
  stopifnot(Matrix::rowSums(Z) == rep(1, n))
  lcst <- 0.5 *
    (p + sum(K)) -
    0.5 * n * log(2 * pi) +
    0.5 * c(Matrix::determinant(prior_beta_prec)$modulus) +
    sum(prior_var_shape * log(prior_var_rate)) -
    sum(log(prior_var_shape)) +
    sum(lgamma(prior_var_shape + 0.5 * c(n, K)))
  maxiter <- as.integer(maxiter)
  # Initial values
  B_u <- rexp(r)
  B_e <- rexp(1)
  # Fixed values (do not change)
  A_e <- prior_var_shape[1] + 0.5 * n
  A_u <- prior_var_shape[-1] + 0.5 * K
  elbo <- numeric(maxiter) # container for ELBO
  for (b in seq_len(maxiter)) {
    # Update coefficients (means)
    prec <- A_e /
      B_e *
      Matrix::crossprod(C) +
      Matrix::bdiag(c(
        prior_beta_prec,
        sapply(1:r, function(i) {
          Matrix::Diagonal(K[i], x = A_u[i] / B_u[i])
        })
      ))
    V <- solve(prec)
    mu <- A_e / B_e * as.numeric(V %*% crossprod(C, y))
    B_e <- prior_var_rate[1] +
      0.5 * as.numeric(crossprod(y - as.numeric(C %*% mu))) +
      sum(Matrix::diag(crossprod(C) %*% V))
    B_u <- prior_var_rate[2] +
      0.5 * as.numeric(crossprod(mu[-(1:p)])) +
      sum(diag(V)[-(1:p)])
    elbo[b] <- lcst -
      0.5 * as.numeric(Matrix::determinant(prec)$modulus) -
      0.5 *
        as.numeric(
          t(mu[1:p]) %*%
            prior_beta_prec %*%
            mu[1:p] +
            Matrix::determinant(
              prior_beta_prec %*% V[1:p, 1:p]
            )$modulus
        ) -
      A_e * log(B_e) -
      sum(A_u * log(B_u))
    if (b > 2 && (elbo[b] - elbo[b - 1]) < tol) {
      break
    }
  }
  list(
    elbo = elbo[seq_len(b)],
    X = X,
    Z = Z,
    sigmasq_shape = c(A_e, A_u),
    sigmasq_rate = c(B_e, B_u),
    coef_mean = mu,
    coef_var = V
  )
}

data(Orthodont, package = "nlme")
# Fix factors (unordered)
Orthodont$Subject <- factor(Orthodont$Subject, ordered = FALSE)
sd_y <- sd(Orthodont$distance)
y <- Orthodont$distance / sd_y # standardize response
fit_lmm <- cavi_varcomp_linmixmod(
  y = y,
  Xform = ~ age + Sex,
  Zform = ~Subject,
  data = Orthodont
)
lme <- lme4::lmer(y ~ age + Sex + (1 | Subject), data = Orthodont)
# Compare coefficients for mean + random effects
coefs <- c(lme4::fixef(lme), unlist(lme4::ranef(lme)))
max(abs(coefs - fit_lmm$coef_mean))
# Posterior mean of variance
lme4::VarCorr(lme)
sqrt(fit_lmm$sigmasq_rate / (fit_lmm$sigmasq_shape - 1))
