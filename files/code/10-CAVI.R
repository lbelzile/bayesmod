# Simulate fake data
set.seed(80601)
y <- rnorm(n = 20, mean = 20, sd = 3)
hist(y)
init <- c(0,1,1)
# Prior parameters
a0 <- 0.01
b0 <- 0.01
mu0 <- 0
tau0 <- 1e-4
# Precompute summary statistics
n <- length(y); 
sum_y <- sum(y)
# Initialize algorithm
E_mu <- init[1]
var_mu <- init[2]
E_tau <- init[3]
# Determine the number of iterations
B <- 20
an <- a0+0.5*n
# Monitor the ELBO
elbo <- numeric(B)
for(i in 1:B){
  # Update each parameter in turn
  var_mu <- 1/(tau0 + n*E_tau)
  E_mu <- var_mu*(tau0*mu0 + E_tau*sum_y)
  bn <- (b0+0.5*sum((y-E_mu)^2) + 0.5*n*var_mu)
  E_tau <- an/bn
  # Recompute the ELBO at the end of each cycle
  elbo[i] <- -an*log(bn) +0.5*log(var_mu)-0.5*tau0*((mu0-E_mu)^2 + var_mu)
}

# Monitor the convergence of the ELBO
plot(1:B, elbo, type = "b")
# Plot the marginal approximations
par(mfrow = c(1,2))
curve(dnorm(x, mean = E_mu, sd = sqrt(var_mu)), 
            from = E_mu - 3*sqrt(var_mu), 
            to = E_mu + 3*sqrt(var_mu),
      ylab = "density", xlab = expression(mu), bty = "n")

curve(dgamma(x, shape = an, rate = bn), from = 0, to = 0.5,
      n = 1001,
      ylab = "density", xlab = expression(tau), bty = "n")


