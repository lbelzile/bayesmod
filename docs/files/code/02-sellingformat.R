library(ggplot2)
library(patchwork)
remotes::install_github("lbelzile/hecbayes")
# Load data set
data("sellingformat", package = "hecbayes")
# Create contingency table
(cont <- with(sellingformat, table(purchased, format)))
#  Put uniform prior on 0/1 for the binomial probability of success
#  the conjugate posterior is beta(y + 1, n - y + 1)

# Sample from the posterior
post_p_int <- rbeta(n = 1e4L, shape1 = 47, shape2 = 153)
post_p_seq <- rbeta(n = 1e4L, shape1 = 24, shape2 = 177)
# Probability of superiority (akin to one-sided test of mu2 > mu1)
mean(post_p_int > post_p_seq)
# Reparametrization in terms of odds
post_odds_int <- (post_p_int / (1 - post_p_int))
post_odds_seq <- (post_p_seq / (1 - post_p_seq))
# Posterior odds
post_oddsratio <- post_odds_int / post_odds_seq

# Plot posterior densities for probability of buying the product
cols <- MetBrewer::met.brewer("Hiroshige", 2)
g1 <- ggplot() +
  stat_function(
    fun = dbeta,
    xlim = c(0, 0.5),
    n = 1001,
    args = list(shape1 = 47, shape2 = 153),
    mapping = aes(col = "integrated")
  ) +
  stat_function(
    fun = dbeta,
    xlim = c(0, 0.5),
    n = 1001,
    args = list(shape1 = 24, shape2 = 177),
    mapping = aes(col = "sequential")
  ) +
  scale_color_manual(
    name = 'Sales format',
    breaks = c('sequential', 'integrated'),
    values = c('sequential' = cols[1], 'integrated' = cols[2])
  ) +
  scale_x_continuous(breaks = seq(0, 0.5, by = 0.25),
                     labels = c("0", "0.25", "0.5")) +
  labs(y = "", subtitle = "Posterior density", x = "probability of buying") +
  scale_y_continuous(limits = c(0, NA), expand = expansion()) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.9, 0.9))
# Plot posterior odds
g2 <- ggplot(data = data.frame(ratio = post_oddsratio),
             mapping = aes(x = ratio)) +
  geom_density() +
  labs(x = "odds ratio of integrated vs sequential decisions", subtitle = "posterior density", y = "") +
  scale_y_continuous(limits = c(0, NA), expand = expansion())
g1 + g2
