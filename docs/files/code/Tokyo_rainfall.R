remotes::install_github("lbelzile/hecbayes")
remotes::install_github("lbelzile/TruncatedNormal")
library(Matrix)
library(TruncatedNormal)
# Load data
data(tokyorain, package = "hecbayes")
# Aggregate data
tokyo <- tokyorain |>
  dplyr::group_by(day) |>
  dplyr::summarize(y = sum(y), n = dplyr::n())
nt <- 366L
# Circulant random walk of order two precision matrix
Q <- hecbayes::crw_Q(d = nt, type = "rw1", sparse = TRUE)
# Sparse Cholesky root
cholQ <- Matrix::chol(Q)
N <- Matrix::Diagonal(n = nt, x = tokyo$n)
# Create containers
B <- 1e4L # number of draws
beta_s <- matrix(nrow = B, ncol = nt)
x_s <- matrix(nrow = B, ncol = nt)
tau_s <- numeric(B)
# Initial values
beta <- rep(0, nt)
tau <- 1000
# Hyperprior parameter values
tau_a <- 1
tau_b <- 0.0001
# Gibbs sampling
for (b in seq_len(B)) {
  # Step 1: data augmentation
  x <- TruncatedNormal::rtnorm(
    n = 1,
    mu = beta[tokyorain$day],
    sd = 1,
    lb = ifelse(tokyorain$y == 0, -Inf, 0),
    ub = ifelse(tokyorain$y == 0, 0, Inf)
  )
  tx <- aggregate(x = x, by = list(tokyorain$day), FUN = sum)$x
  x_s[b, ] <- tx
  # Step 2: Simulate random effects in block
  beta <- beta_s[b, ] <- c(hecbayes::rGaussQ(
    n = 1,
    b = tx,
    Q = tau * Q + N
  ))
  # Simulate precision
  tau <- tau_s[b] <- rgamma(
    n = 1,
    shape = (nt - 2) / 2 + tau_a,
    rate = 0.5 * as.numeric((t(beta) %*% Q %*% beta)) + tau_b
  )
  # if beta is VERY smooth, then precision is large
}
colnames(beta_s) <- paste0("beta", 1:nt)
post_tokyo1 <- data.frame(cbind(tau = tau_s, beta_s))
post_beta <- colMeans(beta_s[, -1])
plot(1:366, pnorm(post_beta), type = "l")
lines(1:366, pnorm(colMeans(post_tokyo1[, -1])), col = 2)
# Initial values
beta <- rep(0, nt)
tau <- 1000
# Gibbs sampling
for (b in seq_len(B)) {
  # if(b %% 100L == 0L){print(paste(b, 'iterations completed')); }
  # Step 1: data augmentation
  x <- TruncatedNormal::rtnorm(
    n = 1,
    mu = beta[tokyorain$day],
    sd = 1,
    lb = ifelse(tokyorain$y == 0, -Inf, 0),
    ub = ifelse(tokyorain$y == 0, 0, Inf)
  )
  tx <- aggregate(x = x, by = list(tokyorain$day), FUN = sum)$x
  x_s[b, ] <- tx
  ## Step 2 b: alternative
  ## Update individually looping over indices
  st <- sample.int(nt, 1)
  # Compute mean vector for betas
  mbeta <- Matrix::solve(a = tau * Q + N, b = tx)
  nw <- c(1, -4, -4, 1)
  for (i in (st + seq_len(nt)) %% nt + 1L) {
    # Sample an index at random
    st <- sample.int(nt, 1)
    # Indices of the non-zero entries for row Q[i,]
    nh <- c(i - 3, i - 2, i, i + 1) %% 366 + 1
    prec <- tau * 6 + tokyo$n[i]
    condmean <- mbeta[i] - tau * sum(nw * (beta[nh] - mbeta[nh])) / prec
    beta[i] <- rnorm(n = 1, mean = condmean, sd = 1 / sqrt(prec))
    # beta[i] <- hecbayes::rcondmvnorm(
    #   n = 1,
    #   value = beta,
    #   ind = i,
    #   mean = mbeta,
    #   precision = tau* Q + N)
  }
  beta_s[b, ] <- beta

  # Simulate precision
  tau <- tau_s[b] <- rgamma(
    n = 1,
    shape = (nt - 1) / 2 + tau_a,
    rate = 0.5 * as.numeric(crossprod(cholQ %*% beta)) + tau_b
  )
}
post_tokyo2 <- data.frame(cbind(tau_s, beta_s))


logpost <- function(beta, tau) {
  sum(
    tokyo$y *
      pnorm(q = beta, log.p = TRUE) +
      (tokyo$n - tokyo$y) * pnorm(q = beta, lower.tail = TRUE, log.p = TRUE)
  ) +
    365 * log(tau) -
    0.5 * tau * as.numeric(crossprod(cholQ %*% beta)) +
    (tau_a - 1) * log(tau) -
    tau_b * tau
}

logpost2 <- function(beta, tau) {
  sum(dnorm(x, mean = beta[tokyorain$day], log = TRUE)) +
    (tau_a - 1) * log(tau) -
    tau_b * tau
}
beta <- rep(0, nt)
tau <- 1000
# Gibbs sampling
for (b in seq_len(B)) {
  # Step 1: data augmentation
  x <- TruncatedNormal::rtnorm(
    n = 1,
    mu = beta[tokyorain$day],
    sd = 1,
    lb = ifelse(tokyorain$y == 0, -Inf, 0),
    ub = ifelse(tokyorain$y == 0, 0, Inf)
  )
  tx <- aggregate(x = x, by = list(tokyorain$day), FUN = sum)$x
  x_s[b, ] <- tx
  # Simulate precision
  tau_prop <- hecbayes::rnorm(1, mean = tau, sd = sd(post_tokyo1$tau) / 4)
  beta_prop <- beta_s[b, ] <- c(hecbayes::rGaussQ(
    n = 1,
    b = tx,
    Q = tau * Q + N
  ))
  logR <- logpost2(beta_prop, tau_prop) - logpost2(beta, tau)
  if (logR > -rexp(1)) {
    beta <- beta_prop
    tau <- tau_prop
  }
  tau_s[b] <- tau
  beta_s[b, ] <- beta
}
