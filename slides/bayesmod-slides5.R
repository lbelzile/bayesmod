###########################################
######    Metropolis-Hastings    ##########
###########################################

# Load data
data(upworthy_question, package = "hecbayes")
# Compute sufficient statistics
data <- upworthy_question |>
  dplyr::group_by(question) |>
  dplyr::summarize(ntot = sum(impressions), y = sum(clicks))
# Create containers for MCMC
niter <- 1e4L
chain <- matrix(0, nrow = niter, ncol = 2L)
colnames(chain) <- c("beta", "kappa")
logpost_v <- numeric(niter)


# Code log posterior as sum of log likelihood and log prior
loglik <- function(par, counts = data$y, offset = data$ntot, ...) {
  lambda <- exp(c(par[1] + log(offset[1]), par[1] + par[2] + log(offset[2])))
  sum(dpois(x = counts, lambda = lambda, log = TRUE))
}
# Note common signature of function
logprior <- function(par, ...) {
  dnorm(x = par[1], mean = log(0.01), sd = 1.5, log = TRUE) +
    dnorm(x = par[2], log = TRUE)
}
logpost <- function(par, ...) {
  loglik(par, ...) + logprior(par, ...)
}


# Compute maximum a posteriori (MAP)
map <- optim(
  par = c(-4, 0.07),
  fn = logpost,
  control = list(fnscale = -1),
  offset = data$ntot,
  counts = data$y,
  hessian = TRUE
)
# Use MAP as starting value
cur <- map$par
# Keep track of log posterior to reduce calculations
logpost_cur <- logpost(cur)
# Proposal covariance
cov_map <- -2 * solve(map$hessian)
chol <- chol(cov_map)


set.seed(80601)
naccept <- 0L
for (i in seq_len(niter)) {
  # Multivariate normal proposal - symmetric random walk
  prop <- c(rnorm(n = 2) %*% chol + cur)
  logpost_prop <- logpost(prop)
  logR <- logpost_prop - logpost_cur
  if (logR > -rexp(1)) {
    # accept move
    cur <- prop
    logpost_cur <- logpost_v[i] <- logpost_prop
    naccept <- naccept + 1L
  }
  logpost_v[i] <- logpost_cur
  chain[i, ] <- cur
}

# Posterior summaries
summary(coda::as.mcmc(chain))

###########################################
#############     MALA     ################
###########################################

# Select data for a single question
qdata <- upworthy_question |>
  dplyr::filter(question == "yes") |>
  dplyr::mutate(y = clicks / impressions, no = impressions)


# Create functions with the same signature (...) for the algorithm
logpost <- function(par, data, ...) {
  mu <- par[1]
  sigma <- par[2]
  no <- data$no
  y <- data$y
  if (isTRUE(any(sigma <= 0, mu < 0, mu > 1))) {
    return(-Inf)
  }
  dnorm(x = mu, mean = 0.01, sd = 0.1, log = TRUE) +
    dexp(sigma, rate = 0.7, log = TRUE) +
    sum(dnorm(x = y, mean = mu, sd = sigma / sqrt(no), log = TRUE))
}


logpost_grad <- function(par, data, ...) {
  no <- data$no
  y <- data$y
  mu <- par[1]
  sigma <- par[2]
  c(
    sum(no * (y - mu)) / sigma^2 - (mu - 0.01) / 0.01,
    -length(y) / sigma + sum(no * (y - mu)^2) / sigma^3 - 0.7
  )
}


# Starting values - MAP
map <- optim(
  par = c(mean(qdata$y), 0.5),
  fn = function(x) {
    -logpost(x, data = qdata)
  },
  gr = function(x) {
    -logpost_grad(x, data = qdata)
  },
  hessian = TRUE,
  method = "BFGS"
)
# Check convergence
logpost_grad(map$par, data = qdata)


# Set initial parameter values
curr <- map$par
# Compute a mass matrix
Amat <- solve(map$hessian)
# Cholesky root - for random number generation
cholA <- chol(Amat)


# Create containers for MCMC
B <- 1e4L # number of iterations
warmup <- 1e3L # adaptation period
npar <- 2L
prop_sd <- rep(1, npar) # tuning parameter
chains <- matrix(nrow = B, ncol = npar)
damping <- 0.8
acceptance <- attempts <- 0
colnames(chains) <- names(curr) <- c("mu", "sigma")
# Proposal variance proportional to inverse hessian at MAP
prop_var <- diag(prop_sd) %*% Amat %*% diag(prop_sd)

for (i in seq_len(B + warmup)) {
  ind <- pmax(1, i - warmup)
  # Compute the proposal mean for the Newton step
  prop_mean <- c(
    curr +
      damping *
        Amat %*% logpost_grad(curr, data = qdata)
  )
  # prop <- prop_sd * c(rnorm(npar) %*% cholA) + prop_mean
  prop <- c(mvtnorm::rmvnorm(
    n = 1,
    mean = prop_mean,
    sigma = prop_var
  ))
  # Compute the reverse step
  curr_mean <- c(
    prop +
      damping *
        Amat %*% logpost_grad(prop, data = qdata)
  )
  # log of ratio of bivariate Gaussian densities
  logmh <- mvtnorm::dmvnorm(
    x = curr,
    mean = prop_mean,
    sigma = prop_var,
    log = TRUE
  ) -
    mvtnorm::dmvnorm(
      x = prop,
      mean = curr_mean,
      sigma = prop_var,
      log = TRUE
    ) +
    logpost(prop, data = qdata) -
    logpost(curr, data = qdata)
  if (logmh > log(runif(1))) {
    curr <- prop
    acceptance <- acceptance + 1L
  }
  attempts <- attempts + 1L
  # Save current value
  chains[ind, ] <- curr

  # MCMC loop
  if (i %% 100 & i < warmup) {
    # Check acceptance rate and increase/decrease variance
    out <- hecbayes::adaptive(
      attempts = attempts, # counter for number of attempts
      acceptance = acceptance,
      sd.p = prop_sd, #current proposal standard deviation
      target = 0.574
    ) # target acceptance rate
    prop_sd <- out$sd # overwrite current std.dev
    acceptance <- out$acc # if we change std. dev, this is set to zero
    attempts <- out$att # idem, otherwise unchanged
    prop_var <- diag(prop_sd) %*% Amat %*% diag(prop_sd)
  }
} # End of MCMC for loop

#' WAIC
#' @param loglik_pt B by n matrix of pointwise log likelihood
WAIC <- function(loglik_pt) {
  -mean(apply(loglik_pt, 2, mean)) + mean(apply(loglik_pt, 2, var))
}
