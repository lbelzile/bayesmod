data(waiting, package = "hecbayes")
a <- 0.01
b <- 0.01
n <- length(waiting)
s <- sum(waiting)
log_marg_lik <- lgamma(n + a) - lgamma(a) + a * log(b) - (n + a) * log(b + s)
# Laplace approximation
map <- (n + a - 1) / (s + b)
logpost <- function(x) {
  n <- length(waiting)
  ifelse(
    x <= 0,
    -Inf,
    -sum(waiting) *
      x +
      length(waiting) * log(x) +
      # sum(dexp(waiting, rate = x, log = TRUE)) +
      dgamma(x, a, b, log = TRUE)
  )
}
# Hessian evaluated at MAP
Hess <- (n + a - 1) / map^2
log_marg_laplace <- 0.5 * log(2 * pi) - Hess + logpost(map)
# Skew-symmetric approximation
# on log-scale for lambda
# logmap <- log(map)
# logpost_log <- function(loglambda){
#  logpost(exp(loglambda))
# }
# Hess_log <- Hess * map^2

xs <- map +
  seq(
    from = -sqrt(1 / Hess) * 4,
    to = sqrt(1 / Hess) * 4,
    length.out = 101
  )
wskew <- function(theta) {
  p1 <- exp(logpost(theta))
  p2 <- exp(logpost(2 * map - theta))
  p1 / (p1 + p2)
}

skewpost_approx <- function(theta) {
  2 * dnorm(x = theta, mean = map, sd = sqrt(1 / Hess)) * wskew(theta)
}


integrate(skewpost_approx, lower = 1e-03, upper = 0.1)
library(ggplot2)
ggplot() +
  geom_function(
    fun = dgamma,
    args = list(shape = n + a, rate = sum(waiting) + b),
    linewidth = 1.5
  ) +
  geom_function(
    fun = skewpost_approx,
    col = "#e76254",
    linetype = "dashed",
    linewidth = 1.5,
    n = 1001
  ) +
  geom_function(
    fun = dnorm,
    args = list(mean = map, sd = sqrt(1 / Hess)),
    col = "#376795",
    linewidth = 1.5,
    linetype = "dotted",
    n = 1001
  ) +
  MetBrewer::scale_color_met_d(name = "Hiroshige") +
  scale_x_continuous(
    limits = c(0.005, 0.08)
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion()
  ) +
  labs(y = "density", x = expression(lambda)) +
  theme_classic()
