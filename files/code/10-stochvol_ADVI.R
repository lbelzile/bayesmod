library(cmdstanr)
library(hecbayes)
library(ggplot2)
data(exchangerate, package = "hecbayes")
# Compute response from raw spot exchange rates at noon
y <- 100*diff(log(exchangerate$dexrate))
# 'y' is now a series of percentage of log daily differences
y <- y - mean(y) # mean center

# Download and compile the Stan model
# We consider the model with the parametrization
# that incorporates support constraints
url <- "https://lbelzile.github.io/bayesmod/files/code/stochvol2.stan"
file <- "stochvol2.stan"
download.file(url, destfile = file)
stochvol_model <- cmdstanr::cmdstan_model(stan_file = file)

# Fit the model via ADVI
advi_mf <- stochvol_model$variational(
  data = list(y = y, T = length(y)),
  algorithm = "meanfield", # alternative "fullrank"
  seed = 80601)
# Note that there is no option to run a factorized
# density, and "fullrank" is too slow because it includes all
# latent parameters which scales O(n)
# Compute posterior summary statistics
summary <- advi_mf$summary()
postmean <- which(substr(summary$variable, 1,2) == "h[")

# Plot posterior mean of variability
ggplot(data = data.frame(x = exchangerate$date[-1],
                  y = exp(summary$mean[postmean]/2)),
       mapping = aes(x = x, y = y)) +
  geom_line() +
  labs(x = "date", y = "variability") +
  theme_classic()

# Plot marginal density of autocorrelation phi
ggplot(data = advi_mf$draws(variables = "phi",format = "df"),
       mapping = aes(x = phi)) +
  geom_density() +
  labs(x = expression(phi)) +
  theme_classic()
