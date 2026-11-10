library(hecbayes)
# Load data
data(climatechange)
climatechange <- climatechange |>
  dplyr::mutate(RCP = factor(RCP))
# Frequentist analysis of variance with sum-to-zero constraints
lm(
  temp ~ RCP,
  data = climatechange,
  contrasts = list(RCP = contr.sum)
)

y <- climatechange$temp # response variable
group <- as.integer(climatechange$RCP) # dummies for groups
n <- length(y) # total number of observations
nj <- table(group) # number per group
K <- length(nj) # number of groups
B <- 1e4L # number of samples
draws <- matrix(nrow = B, ncol = K + 3)
# Gibbs sampler

# Parameter initialization
mu <- mean(y)
alpha <- rep(0, K)
sigma_a <- 0.5
sigma_y <- 1
for (b in seq_len(B)) {
  # Update the global mean mu
  mu <- rnorm(n = 1, mean = mean(y - alpha[group]), sd = sigma_y / sqrt(n))
  for (k in seq_len(K)) {
    prec_k <- 1 / sigma_a^2 + nj[k] / sigma_y^2
    alpha[k] <- rnorm(
      n = 1,
      mean = sum(y[group == k] - mu) / prec_k,
      sd = sqrt(1 / prec_k)
    )
  }
  sigma_a <- sqrt(
    1 / rgamma(n = 1, shape = 0.01 + 0.5 * sum(alpha^2), scale = K + 0.01)
  )
  sigma_y <- sqrt(
    1 /
      rgamma(
        n = 1,
        shape = 0.01 + 0.5 * sum((y - mu - alpha[group])^2),
        scale = n + 0.01
      )
  )
  draws[b, ] <- c(mu, alpha, sigma_a, sigma_y)
}
