library(INLA)
library(ggplot2)
library(patchwork)

# Rainfall model binomial time series
# with circular random walk Gaussian prior
data(Tokyo)
n <- dim(Tokyo)[1]
# Formula (removing intercept)
formula <- y ~ f(time, model = "rw2", cyclic = TRUE) - 1
# INLA includes default priors
result <- inla(
   formula = formula,
   family = "binomial",
   Ntrials = n, data = Tokyo)
# The resulting object includes marginal densities
# posterior mean, etc. for each hyperparameter
# derived from the marginal approximations
summary(result)
plot(result)


# Stochastic volatility model
data(exchangerate, package = "hecbayes")
# Compute response from raw spot exchange rates at noon
y <- 100*diff(log(exchangerate$dexrate))
# 'y' is now a series of percentage of log daily differences
time <- seq_along(y)
data <- data.frame(y = y, time = time)
# Stochastic volatility model
# https://inla.r-inla-download.org/r-inla.org/doc/likelihood/stochvolgaussian.pdf
# The model uses a log link, and a (log)-gamma prior for the precision
f_stochvol <- formula(y ~ f(time, model = "ar1", param = list(prec = c(1, 0.001))))
mod_stochvol <- inla(f_stochvol, family = "stochvol", data = data)
summary(mod_stochvol)
plot(mod_stochvol)
marg_prec <- mod_stochvol$marginals.hyperpar[[1]]
marg_phi <- mod_stochvol$marginals.hyperpar[[2]]

g1 <- ggplot(data.frame(inla.smarginal(marg_prec))) +
  geom_line(aes(x, y)) +
  labs(x = expression(tau), y = "", subtitle = "posterior of precision") +
  scale_y_continuous(limits = c(0, NA), expand = expansion()) +
  theme_classic()

g2 <- ggplot(data.frame(inla.smarginal(marg_phi))) +
  geom_line(aes(x, y)) +
  labs(x = expression(rho), y = "", subtitle = "posterior correlation") +
  scale_y_continuous(limits = c(0, NA), expand = expansion()) +
  theme_classic()
g1 + g2


# Compute density, quantiles, etc. via inla.*marginal
## approximate 95% credible interval and marginal post median
INLA::inla.qmarginal(marg_phi, p = c(0.025, 0.5, 0.975))
# Change of variable to get variance from precision
marg_var <- inla.tmarginal(
  fun = function(x) { 1 / x },
  marginal = marg_prec)
INLA::inla.qmarginal(marg_var, p = c(0.025, 0.975))

# Posterior marginal mean and variance of phi
mom1 <- INLA::inla.emarginal(
    fun = function(x){x},
    marginal = marg_phi)
mom2 <- INLA::inla.emarginal(
    fun = function(x){x^2},
    marginal = marg_phi)
c(mean = mom1, sd = sqrt(mom2 - mom1^2))
