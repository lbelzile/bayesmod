# Forward sampling for model with Y|P ~ geometric(P) and P ~ beta(alpha1, alpha2)
alpha2 <- 2
alpha1 <- 2.5
n <- 1e5L
# Simulate first from P, then Y|P
P <- rbeta(n, shape1 = alpha1, shape2 = alpha2)
Y <- rgeom(n = n, prob = P)
# We can then discard P and recover the marginal of Y
mean(Y)
# Compute the running mean to check the Monte Carlo accuracy
plot(
  x = log10(seq_len(n)),
  y = cumsum(Y) / seq_len(n),
  type = "l",
  xlab = "number of samples (log10 scale)",
  ylab = "running mean"
)
