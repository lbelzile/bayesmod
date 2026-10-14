library(hecbayes)
# Load data
data(climatechange)
climatechange <- climatechange |>
  dplyr::mutate(RCP = factor(RCP))
# Frequentist analysis of variance with sum-to-zero constraints
lm(
  temp ~ RCP,
  data = climatechange,
  contrasts = list(RCP = contr.sum)
)

y <- climatechange$temp # response variable
group <- as.integer(climatechange$RCP) # dummies for groups
n <- length(y) # total number of observations
nj <- table(group) # number per group
K <- length(nj) # number of groups
B <- 1e4L # number of samples
draws <- matrix(nrow = B, ncol = K + 3)
# Gibbs sampler

# Parameter initialization
mu <- mean(y)
alpha <- rep(0, K)
sigma_a <- 0.5
sigma_y <- 1
for (b in seq_len(B)) {
  # Update the global mean mu
  mu <- rnorm(n = 1, mean = mean(y - alpha[group]), sd = sigma_y / sqrt(n))
  for (k in seq_len(K)) {
    prec_k <- 1 / sigma_a^2 + nj[k] / sigma_y^2
    alpha[k] <- rnorm(
      n = 1,
      mean = sum(y[group == k] - mu) / prec_k,
      sd = sqrt(1 / prec_k)
    )
  }
  sigma_a <- sqrt(
    1 / rgamma(n = 1, shape = 0.01 + 0.5 * sum(alpha^2), scale = K + 0.01)
  )
  sigma_y <- sqrt(
    1 /
      rgamma(
        n = 1,
        shape = 0.01 + 0.5 * sum((y - mu - alpha[group])^2),
        scale = n + 0.01
      )
  )
  draws[b, ] <- c(mu, alpha, sigma_a, sigma_y)
}


# Probabilistic programming with NIMBLE
library(nimble)
oneway_code <- nimbleCode({
  for (i in 1:n) {
    y[i] ~ dnorm(mu + alpha[group[i]], prec_y)
  }
  for (k in 1:K) {
    alpha[k] ~ dnorm(0, prec_alpha)
  }
  # Infamous priors - for comparability only
  mu ~ dnorm(0, sd = 100)
  prec_y ~ dgamma(0.01, 0.01)
  prec_alpha ~ dgamma(0.01, 0.01)
})
# Model is in terms of precision (reciprocal variance)
# here with conjugate priors (for comparison purposes)
oneway_model <- nimbleModel(
  code = oneway_code,
  constants = list(n = length(y), K = K, group = group),
  data = list(y = y),
  inits = list(
    mu = mean(y),
    alpha = rep(0, K),
    prec_y = 0.5,
    prec_alpha = 1
  )
)

params <- c("mu", "alpha", "prec_y", "prec_alpha")
niter <- 1e4L
nburnin <- 1e3L
nchains <- 4L

mcmc_out <- nimbleMCMC(
  code = oneway_model,
  monitors = params,
  niter = niter,
  nburnin = nburnin,
  nchains = nchains
)

plot(
  mcmc_out$chain1[, "mu"],
  type = "l",
  ylab = expression(mu),
  xlab = "iteration"
)
plot(
  mcmc_out$chain1[, "alpha[1]"],
  type = "l",
  ylab = expression(alpha[1]),
  xlab = "iteration"
)
plot(
  mcmc_out$chain1[, "mu"] + mcmc_out$chain1[, "alpha[1]"],
  type = "l",
  ylab = expression(mu + alpha[1]),
  xlab = "iteration"
)
# Neal's funnel (not quite visible due to poor exploration
plot(
  mcmc_out$chain1[, "alpha[1]"],
  log(1 / mcmc_out$chain1[, "prec_alpha"]),
  ylab = expression(sigma[alpha]^2),
  xlab = expression(alpha[1])
)
