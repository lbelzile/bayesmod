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
data <- with(
  smartwatch,
  list(
    N = length(nviolation),
    J = length(levels(id)),
    K = length(levels(task)),
    y = as.integer(nviolation),
    id = as.integer(id),
    fixed = as.integer(task)
  )
)

# Sample from HMC-NUTS
postsamp <- smartwatch_model$sample(
  data = data,
  iter_warmup = 1000,
  iter_sampling = 2500,
  chains = 4L,
  refresh = 0L,
  # output_dir = "models", # output directory to save
  # output_basename = "Poisson_mixed", # nqame of files
  show_messages = FALSE
)

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
  scale_y_discrete(
    labels = c(
      "beta[1]" = expression(paste(beta[1], ": phone")),
      "beta[2]" = expression(paste(beta[2], ": watch")),
      "beta[3]" = expression(paste(beta[3], ": speaker")),
      "beta[4]" = expression(paste(beta[4], ": texting")),
      "kappa" = expression(kappa)
    )
  )

# Load some posterior samples from the model
data(posterior_smartwatch, package = "hecbayes")
# Match coefficients to column of posterior
coef_beta <- paste0("beta[", as.integer(smartwatch$task), "]")
id_beta <- match(coef_beta, colnames(posterior_smartwatch))
coef_alpha <- paste0("alpha[", as.integer(smartwatch$id), "]")
id_alpha <- match(coef_alpha, colnames(posterior_smartwatch))
# Create containers for pointwise log likelihood and posterior predictive
B <- nrow(posterior_smartwatch)
y <- smartwatch$nviolation
n <- length(y)
loglik_pt <- matrix(nrow = B, ncol = n)
postpred <- matrix(nrow = B, ncol = n)
# Sample from posterior predictive / evaluate log likelihood
for (b in seq_len(B)) {
  loglik_pt[b, ] <- dpois(
    x = smartwatch$nviolation,
    lambda = exp(
      posterior_smartwatch[b, id_alpha] + posterior_smartwatch[b, id_beta]
    ),
    log = TRUE
  )
  postpred[b, ] <- rpois(
    n = n,
    lambda = exp(
      posterior_smartwatch[b, id_alpha] + posterior_smartwatch[b, id_beta]
    )
  )
}
# Watanabe's widely available information criterion
WAIC <- function(loglik_pt) {
  -mean(apply(loglik_pt, 2, mean)) + mean(apply(loglik_pt, 2, var))
}
WAIC(loglik_pt)
# LOO-CV P-P plot with importance sampling estimator
ploo <- numeric(n)
for (i in seq_along(y)) {
  ploo[i] <- weighted.mean(x = I(postpred[, i] <= y[i]), w = is_wgt[, i])
}
ggplot(
  data = data.frame(ploo = sort(ploo), punif = ppoints(n)),
  mapping = aes(x = punif, y = ploo)
) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  labs(x = "uniform", y = "LOO-PIT")

# Automatically with PSIS instead
loo_test <- loo::loo(loglik_pt, save_psis = TRUE)
# log unstandardized weights for leave-one-out
lw <- weights(loo_test$psis_object)
bayesplot::ppc_loo_pit_qq(
  y = y, # n vector of response
  yrep = postpred, # B x n matrix of posterior predictive values
  lw = lw
) + # B x n matrix of log weights
  theme_classic()
