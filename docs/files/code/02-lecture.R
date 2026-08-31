library(revdbayes)
library(ggplot2)
library(patchwork)
remotes::install_github("lbelzile/hecbayes")

###########################################################
#######  Example 2.4 - importance of selling format #######
###########################################################
# Load data set
data("sellingformat", package = "hecbayes")
# Create contingency table
(cont <- with(sellingformat, table(purchased, format)))
#  Put uniform prior on 0/1 for the binomial probability of success
#  the conjugate posterior is beta(y + 1, n - y + 1)

# Calculation of the marginal likelihood
log_marg_post_bern <- function(n, y) {
  lbeta(1 + y, 1 + n - y)
}

# Extract summary statistics (sufficient statistics)
n <- sum(cont)
n0 <- colSums(cont)['quantity-integrated']
n1 <- colSums(cont)['quantity-sequential']
y0 <- cont["1", "quantity-integrated"]
y1 <- cont["1", "quantity-sequential"]
y <- y0 + y1

# Alternative model (one proportion for each subgroup)
BF1 <- log_marg_post_bern(n = n0, y = y0) + # integrated
  log_marg_post_bern(n = n1, y = y1) # sequential decision
# Null model (same proportion regardless of experimental variable)
BF0 <- log_marg_post_bern(n = n, y = y) # pooled
# Bayes factor
exp(BF0 - BF1)

# Normalizing constant is marginal likelihood
# whose log is about -110
# so too small for accurate numerical integration
# Integrate instead a binomial likelihood
# and remove this from the normalizing constant later
marg_binom <- integrate(
  f = function(x) {
    dbinom(x = y0, size = n0, prob = x)
  },
  lower = 0,
  upper = 1
)
# Compare the numerical approximation with the true
c(
  "numerical" = log(marg_binom$value) - lchoose(n0, y0),
  "exact" = log_marg_post_bern(n = n0, y = y0)
)

## Maximum a posteriori
# Unnormalized log-posterior
logpost <- function(theta, y, n) {
  y * log(theta) + (n - y) * log(1 - theta) # prior is uniform
}
# Maximum a posteriori
map <- optimize(
  f = logpost,
  maximum = TRUE,
  interval = c(0, 1),
  y = y1,
  n = n1
)

par(mfrow = c(1, 2), bty = "l")
# Plot the posterior density
curve(
  expr = dbeta(x, shape1 = 24, shape2 = 177),
  from = 0,
  to = 0.4,
  n = 1001L,
  ylab = "posterior density",
  xlab = expression(theta)
)
abline(v = 23 / 197)
abline(v = map$maximum, lty = 2)
# Plot the posterior on the unconstrained space
curve(
  expr = dbeta(plogis(x), shape1 = 24, shape2 = 177) *
    plogis(x) *
    (1 - plogis(x)),
  from = -3,
  to = 0,
  n = 1001L,
  ylab = "posterior density",
  xlab = expression("logit" ~ theta)
)
# Mode is not invariant (due to Jacobian adjustment)
plogis(map$maximum)


# Sample from the posterior
post_p_int <- rbeta(n = 1e4L, shape1 = y0 + 1, shape2 = n0 - y0 + 1)
post_p_seq <- rbeta(n = 1e4L, shape1 = y1 + 1, shape2 = n1 - y1 + 1)
# Probability of superiority (akin to one-sided test of mu2 > mu1)
mean(post_p_int > post_p_seq)
# Reparametrization in terms of odds
post_odds_int <- (post_p_int / (1 - post_p_int))
post_odds_seq <- (post_p_seq / (1 - post_p_seq))
# Posterior odds
post_oddsratio <- post_odds_int / post_odds_seq

# Pinball loss
pinball <- function(y, qlev = 0.5) {
  ifelse(y < 0, (qlev - 1) * y, qlev * y)
}

# Posterior summary via loss functions
# Computing the 80% by minimizing a loss function
q80 <- optimize(
  f = function(x) {
    mean(pinball(post_oddsratio - x, 0.8))
  },
  interval = c(0, 1e10)
)$minimum
# Compare answer with empirical quantile
quantile(post_oddsratio, 0.8)

# 80% Highest posterior density interval
hdiD <- HDInterval::hdi(
  density(post_oddsratio),
  credMass = 0.80
)
# Equitailed confidence intervals
quantile(post_oddsratio, probs = c(0.1, 0.9))

# Plot posterior densities for probability of buying the product
cols <- MetBrewer::met.brewer("Hiroshige", 2)
g1 <- ggplot() +
  stat_function(
    fun = dbeta,
    xlim = c(0, 0.5),
    n = 1001,
    args = list(shape1 = y0 + 1, shape2 = n0 - y0 + 1),
    mapping = aes(col = "integrated")
  ) +
  stat_function(
    fun = dbeta,
    xlim = c(0, 0.5),
    n = 1001,
    args = list(shape1 = y1 + 1, shape2 = n1 - y1 + 1),
    mapping = aes(col = "sequential")
  ) +
  scale_color_manual(
    name = 'Sales format',
    breaks = c('sequential', 'integrated'),
    values = c('sequential' = cols[1], 'integrated' = cols[2])
  ) +
  scale_x_continuous(
    breaks = seq(0, 0.5, by = 0.25),
    labels = c("0", "0.25", "0.5")
  ) +
  labs(
    y = "",
    subtitle = "Posterior density",
    x = "probability of buying"
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion()
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.9, 0.9)
  )
# Plot posterior odds
g2 <- ggplot(
  data = data.frame(ratio = post_oddsratio),
  mapping = aes(x = ratio)
) +
  geom_density() +
  labs(
    x = "odds ratio of integrated vs sequential decisions",
    subtitle = "posterior density",
    y = ""
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion()
  )
g1 + g2


############################################################
#######  Example 2.6 VaR for Danish insurance losses #######
############################################################

data(danish, package = "evir")
# Using ratio-of-uniform, generate posterior samples
# from the binomial - generalized Pareto model
post_samp <- revdbayes::rpost_rcpp(
  n = 1000L,
  model = "bingp",
  data = danish,
  prior = revdbayes::set_prior(prior = "mdi", model = "gp"),
  thresh = 10
)
# Generates in $sim_vals the scale and shape parameters
# and in '$bin_sim_vals' the binomial probability 'probexc'
# Combine samples into a matrix
post_samp <- cbind(post_samp$sim_vals, post_samp$bin_sim_vals)
colnames(post_samp) <- c('scale', 'shape', "probexc")
post_samp <- as.data.frame(post_samp)

# Plot the bivariate posterior of the generalized Pareto parameters
g2 <- ggplot(
  data = post_samp,
  mapping = aes(x = scale, y = shape)
) +
  geom_point() +
  labs(x = expression(tau), y = expression(xi)) +
  theme_classic() +
  theme(axis.title.y = element_text(angle = 0))
g1 + g2

# Compute value at risk from generalized Pareto distribution quantile fn
VaR_post <- with(
  post_samp, # data frame of posterior draws
  revdbayes::qgp(
    # with columns 'probexc', 'scale', 'shape'
    p = 0.01 / probexc,
    loc = 10,
    scale = scale,
    shape = shape,
    lower.tail = FALSE
  )
)
# Loss function
loss <- function(qhat, q) {
  mean(ifelse(q > qhat, 0.5 * (0.99 * q - qhat), 0.75 * (qhat - 1.01 * q)))
}
# Create a grid of values over which to estimate the loss for VaR
nvals <- 101L
VaR_grid <- seq(
  from = quantile(VaR_post, 0.01),
  to = quantile(VaR_post, 0.99),
  length.out = nvals
)
# Create a container to store results
risk <- numeric(length = nvals)
for (i in seq_len(nvals)) {
  # Compute integral (Monte Carlo average over draws)
  risk[i] <- loss(q = VaR_post, qhat = VaR_grid[i])
}

VaR_post <- with(
  post_samp, # data frame of posterior draws
  revdbayes::qgp(
    # with columns 'probexc', 'scale', 'shape'
    p = 0.01 / probexc,
    loc = 10,
    scale = scale,
    shape = shape,
    lower.tail = FALSE
  )
)
# Loss functions
loss1 <- function(qhat, q) {
  mean(ifelse(q > qhat, 0.5 * (0.99 * q - qhat), 0.75 * (qhat - 1.01 * q)))
}
# squared error loss
loss2 <- function(qhat, q) {
  0.3 * mean((q - qhat)^2)
}
# Create a grid of values over which to estimate the loss for VaR
nvals <- 101L
VaR_grid <- seq(
  from = quantile(VaR_post, 0.01),
  to = quantile(VaR_post, 0.99),
  length.out = nvals
)
# Create a container to store results
risk1 <- risk2 <- numeric(length = nvals)
for (i in seq_len(nvals)) {
  # Compute integral (Monte Carlo average over draws)
  risk1[i] <- loss1(q = VaR_post, qhat = VaR_grid[i])
  risk2[i] <- loss2(q = VaR_post, qhat = VaR_grid[i])
}
g1 <- ggplot(data = data.frame(x = VaR_post)) +
  geom_density(mapping = aes(x = x)) +
  labs(
    subtitle = "posterior density",
    x = "value-at-risk 0.99 (in millions krone)",
    y = ""
  ) +
  scale_y_continuous(limits = c(0, 0.15), expand = c(0, 0)) +
  theme_classic()
# Plot loss function
g2 <- ggplot(
  data = data.frame(
    loss1 = risk1 - min(risk1),
    loss2 = risk2 - min(risk2),
    quantile = VaR_grid
  )
) +
  geom_line(mapping = aes(x = quantile, y = loss1)) +
  geom_line(mapping = aes(x = quantile, y = loss2), linetype = "dashed") +
  geom_vline(xintercept = VaR_grid[which.min(risk1)], linewidth = 0.1) +
  geom_vline(
    xintercept = mean(VaR_post),
    linetype = "dashed",
    linewidth = 0.1
  ) +
  scale_y_continuous(limits = c(0, 9), expand = c(0, 0)) +
  labs(
    x = "value-at-risk 0.99 (in millions krone)",
    subtitle = "custom loss function (full)\nand squared error loss (dashed)"
  ) +
  theme_classic()
g1 + g2
