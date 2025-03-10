# Fit Gaussian regression model via 'brms' interface 
# using STAN, sampling posterior draws via HMC
library(brms)
library(posterior)
library(ggplot2)
library(bayesplot)
library(patchwork)
data(diabetes, package = "lars")
fitmod <- brms::brm(
  formula = log(y) ~ x, 
  iter = 5000L, 
  warmup = 1000L,
  family = gaussian(), 
  data = diabetes, 
  prior = set_prior(horseshoe(df = 3, par_ratio = 0.1, main = TRUE), class = "b"),
  backend = "cmdstanr", # fitting engine
  control = list(adapt_delta = 0.9), #modify acceptance rate of HMC
  save_model = "fhs.stan"
)

bayesplot::color_scheme_set(scheme = "gray")
ggplot2::theme_set(ggplot2::theme_classic())
# Parameter estimates
summary(fitmod)
# Posterior predictive checks
pp_check(fitmod) + pp_check(fitmod, type = "loo_pit")

# Posterior density (quartiles + 89% intervals) via boxplots
mcmc_plot(fitmod)
# Traceplots
mcmc_plot(fitmod, type = "trace")
# ACFs
mcmc_plot(fitmod, type = "acf_bar")
