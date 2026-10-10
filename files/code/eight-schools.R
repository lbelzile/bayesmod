library(ggplot2)

# Eight school data
es <- data.frame(
  school = LETTERS[1:8],
  y = c(28, 8, -3, 7, -1, 1, 18, 12),
  se = c(15, 10, 16, 11, 9, 11, 10, 18)
)
n <- nrow(es)
B <- 1e3L
set.seed(80601)
# Starting values
alpha <- rnorm(0, sd = 0.01)
mu <- mean(es$y)
sigma_al <- sqrt(var(es$y) - mean(es$se))
pars1 <- matrix(NA, nrow = B, ncol = n + 2)
var_mu <- 1 / sum(1 / es$se^2)
for (b in seq_len(B)) {
  var_al <- 1 / (1 / es$se^2 + 1 / sigma_al^2)
  pars1[b, 1:n] <- alpha <- rnorm(
    n = n,
    mean = var_al * (es$y - mu) / es$se^2,
    sd = sqrt(var_al)
  )
  pars1[b, n + 1] <- mu <- rnorm(
    n = 1,
    mean = var_mu * sum((es$y - alpha) / es$se^2),
    sd = sqrt(var_mu)
  )
  pars1[b, n + 2] <- sigma_al <- sum(alpha^2) / rchisq(1, df = n - 1)
}

## Sampler 2
# Same model, with joint updates
set.seed(80601)
# Starting values
alpha <- rnorm(0, sd = 0.01)
mu <- mean(es$y)
sigma_al <- sqrt(var(es$y) - mean(es$se))
pars2 <- matrix(NA, nrow = B, ncol = n + 2)
for (b in seq_len(B)) {
  var_al <- 1 / (1 / es$se^2 + 1 / sigma_al^2)
  # Sample from joint of alpha, mu
  # by taking p(alpha | sigma_al, y) * p(mu | alpha, sigma_al, y)
  pars2[b, 1:n] <- alpha <- rnorm(
    n = n,
    mean = var_al * es$y / es$se^2,
    sd = sqrt(var_al)
  )
  pars2[b, n + 1] <- mu <- rnorm(
    n = 1,
    mean = var_mu * sum((es$y - alpha) / es$se^2),
    sd = sqrt(var_mu)
  )
  pars2[b, n + 2] <- sigma_al <- sum(alpha^2) / rchisq(1, df = n - 1)
}

# Parameter expansion - V+PX sampler
a <- 1
pars3 <- matrix(NA, nrow = B, ncol = n + 3)
for (b in seq_len(B)) {
  var_al <- 1 / (1 / es$se^2 + 1 / sigma_al^2)
  # Sample from joint of alpha, mu
  # by taking p(alpha | sigma_al, y) * p(mu | alpha, sigma_al, y)
  alpha_st <- rnorm(n = n, mean = var_al * es$y / es$se^2, sd = sqrt(var_al))
  pars3[b, n + 1] <- mu <- rnorm(
    n = 1,
    mean = var_mu * sum((es$y - alpha) / es$se^2),
    sd = sqrt(var_mu)
  )
  pars3[b, n + 2] <- sigma_al_st <- sum(alpha_st^2) / rchisq(1, df = n - 1)
  pars3[b, n + 3] <- a <- rnorm(
    n = 1,
    mean = sum(alpha_st * (es$y - mu) / es$se^2) / sum(alpha_st^2 / es$se^2),
    sd = sqrt(sum(alpha_st^2 / es$se^2))
  )
  pars3[b, 1:n] <- alpha_st * a
  pars3[b, n + 2] <- sigma_al <- abs(a) * sigma_al_st
}


# Funnel behaviour and poor sampling of first model
g1 <- ggplot(
  data = data.frame(
    beta2 = c(pars1[, 2]),
    sigma = c(pars1[, n + 2])
  ),
  mapping = aes(x = beta2, y = log(sigma))
) +
  geom_point() +
  scale_color_grey() +
  scale_y_continuous(limits = c(0, 13)) +
  labs(
    x = expression(alpha[2]),
    y = expression("log" ~ sigma[alpha]),
    caption = "block"
  ) +
  theme_classic()
g2 <- ggplot(
  data = data.frame(
    beta2 = c(pars3[, 2]),
    sigma = c(pars3[, n + 2])
  ),
  mapping = aes(x = beta2, y = log(sigma))
) +
  geom_point() +
  scale_color_grey() +
  scale_y_continuous(limits = c(0, 13)) +
  labs(
    x = expression(alpha[2]),
    y = expression("log" ~ sigma[alpha]),
    caption = "expansion"
  ) +
  theme_classic()
g1 + g2
