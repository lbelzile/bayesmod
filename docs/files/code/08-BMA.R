BMA_linreg <- function(y, X, intercept = TRUE, B = 1e4L) {
  n <- nrow(X)
  p <- ncol(X)
  # standardize inputs
  X <- apply(X, 2, scale)
  if (intercept) {
    stopifnot(isTRUE(all.equal(as.numeric(X[, 1]), rep(0, n))))
    X[, 1] <- rep(1, n)
  } else {
    X <- cbind(1, X)
  }
  # Hyperpriors (vague)
  mu_0 <- rep(0, p)
  Omega_0 <- diag(c(0.0001, rep(1, p - 1)))
  nu_0 <- 0.01
  tau_0 <- 0.01
  # Compute update parameters
  post_beta <- function(Xm, y, Omega0, mu0, tau0) {
    model <- lm(y ~ -1 + Xm)
    Omega_n <- crossprod(Xm) + Omega0
    Omega_n_inv <- solve(Omega_n)
    mu_n <- c(Omega_n_inv %*% (crossprod(Xm, y) + crossprod(Omega0, mu0)))
    tau_n <- tau0 +
      sum(resid(model)^2) +
      c(crossprod(Xm %*% (mu_n - coef(model)))) +
      sum((mu_n - mu0)^2)
    list(prec = Omega_n, vcov = Omega_n_inv, mean = mu_n, tau = tau_n)
  }
  # Score for model averaging
  score_M <- function(pars) {
    Cmat <- chol(pars$prec)
    c(crossprod(Cmat %*% pars$mean)) - sum(log(diag(Cmat)))
  }
  # Indices of covariates in the model
  # Initial model
  M_curr <- 1
  pars_curr <- post_beta(
    Xm = X[, M_curr, drop = FALSE],
    y = y,
    Omega0 = Omega_0[M_curr, M_curr, drop = FALSE],
    mu0 = mu_0[M_curr],
    tau0 = tau_0
  )
  # Containers
  omega_s <- numeric(B)
  beta_s <- matrix(0, nrow = B, ncol = p)
  # For loop for sampling
  for (b in 1:B) {
    # Step 1: update the model matrix
    new_ind <- sample(2:p, 1)
    if (new_ind %in% M_curr) {
      # Remove component (death)
      M_prop <- M_curr[-which(M_curr == new_ind)]
    } else {
      M_prop <- sort(c(M_curr, new_ind))
    }
    # Jacobian adjustment
    if (length(M_prop) == 1L) {
      logJ <- log(3) - log(2)
    } else if (length(M_curr) == 1L) {
      logJ <- log(2) - log(3)
    } else if (length(M_prop) == p) {
      logJ <- log(3) - log(2)
    } else if (length(M_curr) == p) {
      logJ <- log(2) - log(3)
    } else {
      logJ <- 0
    }
    score_curr <- score_M(pars_curr)
    # Proposal model
    pars_prop <- post_beta(
      Xm = X[, M_prop, drop = FALSE],
      y = y,
      Omega0 = Omega_0[M_prop, M_prop, drop = FALSE],
      mu0 = mu_0[M_prop],
      tau0 = tau_0
    )
    score_prop <- score_M(pars_prop)
    # Metropolis-Hastings step for change of model matrix
    logR <- score_prop - score_curr + logJ
    if (logR > log(runif(1))) {
      M_curr <- M_prop
      pars_curr <- pars_prop
    }
    # Step 2: update global precision and mean coefs
    omega_s[b] <- rgamma(
      n = 1,
      shape = 0.5 * (nu_0 + n),
      rate = pars_curr$tau^2 / 2
    )
    beta_s[b, M_curr] <- c(mvtnorm::rmvnorm(
      n = 1,
      mean = pars_curr$mean,
      sigma = pars_curr$vcov / omega_s[b]
    ))
  }
  return(list(beta = beta_s, sigma_sq = 1 / omega_s))
}

# Example of BMA with diabetes data set
data(diabetes, package = "lars")
bma <- BMA_linreg(
  y = diabetes$y,
  X = diabetes$x2,
  intercept = FALSE,
  B = 1e4L
)
# Warning: relatively slow...
burn <- 1000
pinclusion <- colMeans(bma$beta[-seq_len(burn), ] == 0)

