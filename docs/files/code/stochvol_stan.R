library(cmdstanr)
library(hecbayes)
library(posterior)
library(bayesplot)

data(exchangerate, package = "hecbayes")
# Compute response from raw spot exchange rates at noon
y <- 100 * diff(log(exchangerate$dexrate))
# 'y' is now a series of percentage of log daily differences
y <- y - mean(y) # mean center

# Download and compile the Stan model
url <- "https://lbelzile.github.io/bayesmod/files/code/stochvol1.stan"
file <- "stochvol1.stan"
download.file(url, destfile = file)
stochvol1_model <- cmdstanr::cmdstan_model(stan_file = file)

url <- "https://lbelzile.github.io/bayesmod/files/code/stochvol2.stan"
file <- "stochvol2.stan"
download.file(url, destfile = file)
stochvol2_model <- cmdstanr::cmdstan_model(stan_file = file)
# Sample from HMC-NUTS
fit <- stochvol2_model$sample(
  data = list(y = y, T = length(y)),
  iter_warmup = 1000,
  iter_sampling = 2500,
  chains = 4L
)
# Compute posterior summary statistics
fit$summary()
# Extract posterior samples
post_samp <- fit$draws()
# Plot a histogram
bayesplot::mcmc_hist(pars = "phi")
bayesplot::mcmc_trace_highlight(pars = "sigma")
