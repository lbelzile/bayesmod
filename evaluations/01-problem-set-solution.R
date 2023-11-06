library(dplyr)
library(tidyr)
library(ggplot2)
library(hecbayes)
data(voting, package = "hecbayes")
dirichlet <-  voting |>
    group_by(educ, vote) |>
    summarize(alpha = sum(weight))
# Marginals are beta
prop <- dirichlet |>
  group_by(educ) |>
  summarize(prop = alpha[vote == "always"]/sum(alpha))
# Extract parameter
pars <- matrix(dirichlet$alpha, nrow = 3, ncol = 3)
nsim <- 1e3L
post_samp <- array(dim = c(nsim, 3, 3))
oddsratio <- array(dim = c(nsim, 6))
for(i in 1:3){
  post_samp[,,i] <-  mev::rdir(n = nsim, alpha = pars[,i])
  oddsratio[,1:2 + 2*(i-1)] <- post_samp[,-3,i]/post_samp[,3,i]
}
# Identify columns for odds ratio
colnames(oddsratio) <- paste0(
  rep(levels(dirichlet$educ), each = 2), "_",
  rep(c("rarely/never vs always", "sporadic vs always"), length.out = 6))
# Create data frame with factors for ggplot2
odds_data <- tidyr::pivot_longer(as.data.frame(oddsratio),
                    cols = 1:6,
                    names_to = c("educ","voting"),
                    names_sep = "_",
                    values_to = "odds")
# Plot densities
ggplot(data = odds_data,
       mapping = aes(x = odds)) +
  geom_density() +
  facet_wrap(~ voting + educ) +
  theme_classic()
# Compute equitailed credible intervals
alpha <- 0.11
apply(oddsratio, 2, quantile, probs = c(alpha/2, 1-alpha/2))
