## Examples of Laplace approximation with exponential-gamma model
data(waiting, package = "hecbayes")
# Hyperparameters for prior
a <- 0.01
b <- 0.01
n <- length(waiting) # sample size
s <- sum(waiting) # sufficient statistic
# Exact formula for the log marginal likelihood
log_marg_lik <- lgamma(n + a) - lgamma(a) + a * log(b) - (n + a) * log(b + s)

# Laplace approximation
map <- (n + a - 1) / (s + b) # posterior mode/MAP
# Unnormalized log posterior
logpost <- function(x) {
  sum(dexp(waiting, rate = x, log = TRUE)) +
    dgamma(x, a, b, log = TRUE)
}
# Hessian evaluated at MAP
H <- -numDeriv::hessian(logpost, x = map)
# Approximation to log of marginal likelihood
log_marg_laplace <- 1 /
  2 *
  log(2 * pi) -
  c(determinant(H)$modulus) +
  logpost(map)


# Gaussian approximation to the posterior
# True posterior (obtained via conjugacy)
post <- function(x) {
  dgamma(x, shape = n + a, rate = sum(waiting) + b)
}
hessian <- function(x) {
  -(n + a) / x^2
}
# Plot true posterior
curve(
  post,
  from = 0.01,
  to = 1 / 15,
  bty = "l",
  ylab = "posterior density",
  xlab = expression("exponential rate" ~ lambda),
  n = 1001
)
#Gaussian approximation
curve(
  dnorm(x, mean = map, sd = sqrt(-1 / hessian(map))),
  add = TRUE,
  lty = 2,
  n = 1001
)


# Posterior expectation

laplace_mean <- exp(-1) /
  (s + b) *
  exp((n + a + 0.5) * log(n + a) - (n + a - 0.5) * log(n + a - 1))
post_mean <- (n + a) / (s + b)
