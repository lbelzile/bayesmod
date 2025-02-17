stochVolCode <- nimbleCode({
  x[1] ~ dnorm(mu, sd = sigma / sqrt(1-phi*phi))
  y[1] ~ dnorm(0, sd = beta * exp(0.5 * x[1]))
  for (t in 2:T){
    x[t] ~ dnorm(mu + phi * (x[t-1] - mu), sd = sigma)
    y[t] ~ dnorm(0, sd = beta * exp(0.5 * x[t]))
  }
  phi ~ dunif(-1, 1)
  sigma ~ dt(mu = 0, sd = 5, df = 1)
  mu ~ dt(mu = 0, sd = 10, df = 1)
  beta <- exp(0.5*mu)
})
