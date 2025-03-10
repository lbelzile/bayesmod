## Fit univariate mixture model to 'faithful' data via Gibbs
## Number of mixture components is specified apriori
library(ggplot2)
library(patchwork)
data(faithful)
n <- nrow(faithful)
y <- faithful$waiting
# Fix hyperpriors
a1 <- 2; a2 <- 2; c <- 60; d <- 1/40; b1 <- 1; b2 <- 0.01
# Assign observations at random to groups
set.seed(80601)
# Fix some cutoff point and assign observations left/right to different clusters
cut <- runif(1, 0.1, 0.9)*diff(range(y)) + min(y)
group <- as.integer(y > cut)
# Sample proportion of observation in group zero
p <- sum(group == 0L)/n
# Initial mean and precision parameters, based on empirical estimates
mu <- c(mean(y[group == 0]), mean(y[group == 1]))
prec <- 1/c(var(y[group == 0]), var(y[group == 1]))
# Storage and number of replications
B <- 1e4L
theta <- matrix(nrow = B, ncol = 5L)
# Step 1: assign variables to clusters
for(b in 1:B){
  d1 <- dnorm(y, mean = mu[1], sd = 1/sqrt(prec[1])) # group 0 
  d2 <- dnorm(y, mean = mu[2], sd = 1/sqrt(prec[2])) # group 1
  # Data augmentation: group labels
  group <- rbinom(n = n, size = rep(1, n), prob = (1-p)*d2/(p*d1 + (1-p)*d2))
  # Step 2: update probability of cluster
  p <- rbeta(n = 1, shape1 = n - sum(group) + a1, sum(group) + a2)
  for(j in 1:2){
    # Extract observations in group g
    yg <- y[group == (j-1L)]
    # Number of observations in group g
    ng <- length(yg)
    # Update mean and precision of the group mean
    prec_mu <- ng*prec[j] + d
    mean_mu <- (sum(yg)*prec[j] + c*d)/prec_mu
    # Draw posterior of mu and precision
    mu[j] <- rnorm(n = 1, mean = mean_mu, sd = 1/sqrt(prec_mu))
    prec[j] <- rgamma(n = 1, 
                      shape = b1 + ng/2, 
                      rate = b2 + 0.5*sum((yg-mu[j])^2))
  }
  # Store values
  theta[b, ] <- c(p, mu, prec)
}
# Discard initial observations (burn in)
theta <- theta[-(1:100),]

# Plot histogram of observations
g1 <- ggplot(data = faithful, 
             mapping = aes(x = waiting)) +
  geom_histogram() +
  scale_y_continuous(limits = c(0, NA), 
                     expand = expansion(add = c(0, 0.5))) + 
  labs(x = "waiting time (in seconds)") +
  theme_classic()

# Define mixture density
mixtdens <- function(theta, x){ 
  theta[1]*dnorm(x, theta[2], sd = 1/sqrt(theta[4])) +
    (1-theta[1]) * dnorm(x, theta[3], sd = 1/sqrt(theta[5]))
}
# Extract limits and create a fine grid over whish to evaluate the density
xlim <- range(faithful$waiting) + c(-2,2)
xseq <- seq(xlim[1], xlim[2], length.out = 101)
map <- apply(X = theta[seq(1, nrow(theta), length.out = 400),], 
             MARGIN = 1, 
             FUN = function(theta_b){ mixtdens(theta_b, x = xseq)})
# Plot posterior density estimates
g2 <- ggplot(data = data.frame(
  x = rep(xseq, length.out = prod(dim(map))), 
  y = c(map), 
  .draw = rep(1:ncol(map), each = nrow(map))),
  mapping = aes(x = x, y = y, group = .draw)) + 
  # ggdist::stat_lineribbon(mapping = aes(fill = after_stat(.width)),
  #                 .width = ppoints(50)) +
  #  scale_fill_distiller() +
  geom_line(alpha = 0.2) +
  labs(x = "waiting time (in seconds)", y = "density") + 
  scale_x_continuous(limits = xlim, expand = expansion()) +
  scale_y_continuous(limits = c(0, NA), 
                     expand = expansion(mult = c(0, 0.1))) + 
  theme_classic()
g1 + g2
