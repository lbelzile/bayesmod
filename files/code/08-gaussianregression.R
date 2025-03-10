library(ggplot2)
remotes::install_github("lbelzile/hecedsm")
data(LKUK24_S4, package = "hecedsm")
# Sum-to-zero constraint on factors
options(contrasts = c("contr.sum", "contr.poly"))
# Fit linear regression model
model <- lm(appropriation ~ politideo * chefdax * brandaction,
          data = LKUK24_S4)
# Extract model matrix, response and dimensions
y <- LKUK24_S4$appropriation
X <- model.matrix(model)
n <- nrow(X)
p <- ncol(X)
# Priors
# We set zero for contrasts and a small value
# for the global mean (since response is on [1,7] Likert scale)
# with lower values indicating more cultural appropriation)
mu_0 <- c(2.5, rep(0, p-1))
Omega_0 <- diag(p)
# Prior precision of 0.25 (variance of 4)
nu_0 <- 0.25
tau_0 <- 1
# Use formulas to obtain posterior precision and mean
Omega_n_inv <- solve(crossprod(X) + Omega_0)
mu_n <- Omega_n_inv %*% (crossprod(X, y) + crossprod(Omega_0, mu_0))
tau_n <- tau_0 + sum(resid(model)^2) + c(crossprod(X %*% (mu_n-coef(model)))) + sum((mu_n-mu_0)^2)
# Posterior draws from the model
omega <- rgamma(n = 1e3L, 
                shape = (nu_0 + n)/2, 
                tau_n^2/2)
# Use forward sampling to obtain the betas
# First simulate multivariate Gaussians
beta <- matrix(rnorm(n = 1e3L*p), ncol = p) %*% 
  chol(Omega_n_inv)
# Divide by gamma scale (scale mixture)
beta <- sweep(beta, 1, sqrt(omega), "/")
# add location parameter
beta <- sweep(beta, 2, mu_n, "+")
# Posterior quartiles for beta
beta_qu <- t(apply(beta, 2, quantile, probs = c(0.25,0.5,0.75)))
# Standard dev. for beta (from Student-t distribution properties)
beta_se <- sqrt((nu_0 + n)/(nu_0 + n - 2) *diag(tau_n/(nu_0 + n) * Omega_n_inv))
# Alternatively, we could also directly generate from multivariate Student-t
beta <- TruncatedNormal::rtmvt(
  n = 1e3L,
  mu = mu_n, 
  sigma = tau_n/(nu_0 + n) * Omega_n_inv, 
  df = nu_0 + n)
## Application to data set: marginal effects 
# Create a grid for each of 12 combinations of factors
dfg <- expand.grid(politideo =c("conservative", "liberal"),
            chefdax = c("not black", "black"),
            brandaction = c("peeking","permission", "control"))
mm <- model.matrix( ~ politideo * chefdax * brandaction,
              data = dfg)
# Subgroup means for each of the 12 categories
mu <- tcrossprod(beta, mm)
# Contrast weights, averaging over brandaction
w1 <- rep(c(1, 0, -1, 0), length.out = nrow(dfg))/3
w2 <- rep(c(0, 1, 0, -1), length.out = nrow(dfg))/3
# Posterior distribution of contrasts
tc <- mu %*% cbind(w1, w2)

# Plot of density
ggplot(data = 
         data.frame(x = c(tc),
                    f = factor(rep(c("liberal","conservative"), 
                                   each = nrow(beta)))),
       mapping = aes(x = x, y = f, group = f)) +
  ggdist::stat_dist_halfeye(.width = 0.5) +
  labs(y = "",
       x = "difference in cultural appropriation score",
       subtitle = "Comparison when Chef Dax is black (vs not)") +
theme_classic()
