library(cmdstanr)
library(hecbayes)
library(posterior)
library(bayesplot)
library(ggplot2)

data(smartwatch, package = "hecbayes")
# Compile the Stan model
url <- "https://lbelzile.github.io/bayesmod/files/code/poisson_smartwatch.stan"
file <- "poisson_smartwatch.stan"
download.file(url, destfile = file)
smartwatch_model <- cmdstanr::cmdstan_model(stan_file = file)
# Set data and constants
data <-  with(
  smartwatch, 
  list(N = length(nviolation),
       J = length(levels(id)),
       K = length(levels(task)),
       y = as.integer(nviolation),
       id = as.integer(id),
       fixed = as.integer(task)))

# Sample from HMC-NUTS 
postsamp <- smartwatch_model$sample(
  data = data,
  iter_warmup = 1000, 
  iter_sampling = 2500, 
  chains = 4L, 
  refresh = 0L,
  # output_dir = "models", # output directory to save
  # output_basename = "Poisson_mixed", # nqame of files
  show_messages = FALSE)

# Summary of posterior samples
summary(postsamp)

# Extract posterior draws
posterior <- as.array(postsamp$post_warmup_draws)

bayesplot::color_scheme_set(scheme = "gray")
bayesplot::mcmc_areas(
  posterior, 
  pars = c("beta[1]", "beta[2]", "beta[3]", "beta[4]", "kappa"),
  prob = 0.50, # 50% intervals
  prob_outer = 0.99, # 99%
  point_est = "median"
) +
  scale_y_discrete(labels=c(
    "beta[1]" = expression(paste(beta[1],": phone")), 
    "beta[2]" = expression(paste(beta[2],": watch")), 
    "beta[3]" = expression(paste(beta[3],": speaker")), 
    "beta[4]" = expression(paste(beta[4],": texting")), 
    "kappa" = expression(kappa)
  ))
