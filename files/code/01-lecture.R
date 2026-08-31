## MATH 80601A - Bayesian modelling
## Code for lecture 1

## Hierarchical model

##,"Bar plot of 1K sample draws (density scale) with true probability
##from the negative binomial marginal model."
# Sample the means of the Poisson
n <- 1000L
k <- 0.5 # dispersion
mu <- 4 # mean
#  The gamma distribution has a mean of 4 and variance of 8
Lambda <- rgamma(n = n, shape = k * mu, rate = k)
# Functions r* for random number generation are vectorized wrt arguments
Y <- rpois(n = n, lambda = Lambda)
# Calculate the empirical mean and variance of the counts
# This is a form of Monte Carlo integration
rbind(
  empirical = c("mean" = mean(Y), "variance" = var(Y)),
  theoretical = c("mean" = mu, "variance" = mu + mu / k)
)
# Create a bar plot manually and compare with marginal CDF of neg. binom
plot(
  x = sort(unique(Y)),
  y = table(Y) / length(Y),
  type = "h",
  xlab = "x",
  ylab = "probability"
)
segments(
  x0 = 0:max(Y) - 0.25,
  x1 = 0:max(Y) + 0.25,
  y0 = dnbinom(0:max(Y), prob = 1 - 1 / (k + 1), size = k * mu),
  col = 2
)

## Change of variable formula and Jacobian

## Example 1.5 from the course notes (gamma to inverse gamma)

rate <- 0.5
shape <- 2
# Jacobian of transformation 1/x is 1/x^2
dinvgamma <- function(x, rate, shape) {
  dgamma(1 / x, rate = rate, shape = shape) / x^2
}
integrate(f = dinvgamma, lower = 0, upper = Inf, rate = rate, shape = shape)

par(mfrow = c(1, 2))
curve(
  expr = dgamma(x, rate = rate, shape = shape),
  from = 0,
  to = 10,
  ylab = "density",
  sub = "gamma"
)
curve(
  expr = dinvgamma(x, rate = rate, shape = shape),
  from = 0,
  to = 10,
  ylab = "density",
  sub = "gamma"
)

## First-order autoregressive process

## Time series of first-order autoregressive process draws with
## standard Gaussian innovations and different autocorrelation parameters
simulate_ar1 <- function(n, phi, mu = 0, sigma = 1) {
  y <- numeric(n) # container of size n
  # simulate from marginal if process is stationary
  y[1] <- ifelse(abs(phi) < 1, rnorm(n = 1, sd = sigma / sqrt(1 - phi^2)), 0)
  for (i in 2:n) {
    y[i] <- mu + phi * (y[i - 1] - mu) + rnorm(n = 1, sd = sigma)
  }
  return(y)
}
set.seed(2026)
par(mfrow = c(1, 3), bty = "l")
# Strong positive correlation
plot(
  simulate_ar1(n = 100, phi = 0.75),
  type = "l",
  ylab = "observation",
  xlab = "time"
)
# Negative correlation (oscillates around the mean)
plot(
  simulate_ar1(n = 100, phi = -0.5),
  type = "l",
  ylab = "observation",
  xlab = "time"
)
# Non-stationary (increasing variance)
plot(
  simulate_ar1(n = 100, phi = 1),
  type = "l",
  ylab = "observation",
  xlab = "time"
)

## Monte Carlo estimation
# "Monte Carlo estimates (running means) as a function of the sample size for
# Pr(Y>1), E(Y) and E(1/Y) for Y ~ gamma(0.5, 2)
set.seed(80601)
B <- 1e5L # number of simulation
Bseq <- seq_len(B)
# Parameters for shape and rate
alpha <- 0.5
beta <- 2
pv <- pgamma(1, shape = alpha, rate = beta)
samp <- rgamma(n = B, shape = alpha, rate = beta)
# Running means of the parameters
int1 <- cumsum(samp < 1) / Bseq
int2 <- cumsum(samp) / Bseq
int3 <- cumsum(1 / samp) / Bseq
par(mfrow = c(1, 3), bty = "l")
plot(
  x = Bseq[-(1:20)],
  y = int1[-(1:20)],
  xlab = "number of draws",
  ylab = "Monte Carlo estimate",
  type = "l"
)
plot(
  x = Bseq[-(1:20)],
  y = int2[-(1:20)],
  xlab = "number of draws",
  ylab = "Monte Carlo estimate",
  type = "l"
)
plot(
  x = Bseq[-(1:20)],
  y = int3[-(1:20)],
  xlab = "number of draws",
  ylab = "Monte Carlo estimate",
  type = "l"
)
# Compute std. error of estimator of total size (assuming it exists)
se1 <- sqrt(var(samp < 1) / B)
# 95% Wald-based confidence interval
confint <- mean(samp < 1) + se1 * qnorm(c(0.025, 0.975))
# Compare with truth
pv
