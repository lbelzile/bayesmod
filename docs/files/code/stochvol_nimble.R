library(nimble, warn.conflicts = FALSE)
library(MCMCvis)
# See https://oliviergimenez.github.io/banana-book/intronimble.html

y <- 100*diff(log(exchangerate$dexrate))
y <- y - mean(y)

stochVol_code <- nimbleCode({
  x[1] ~ dnorm(mu, sd = sigma / sqrt(1-phi*phi))
  y[1] ~ dnorm(0, sd = beta * exp(0.5 * x[1]))
  for (t in 2:N){
    x[t] ~ dnorm(mu + phi * (x[t-1] - mu), sd = sigma)
    y[t] ~ dnorm(0, sd = beta * exp(0.5 * x[t]))
  }
  phi ~ dunif(-1, 1)
  sigma ~ dt(mu = 0, sigma = 5, df = 1)
  mu ~ dt(mu = 0, sigma = 10, df = 1)
  beta <- exp(0.5*mu)
})

stochVol_model <- nimbleModel(
  code = stochVol_code,
  constants = list(N = length(y)), data = list(y = y),
  inits = list(mu = 0, phi = 0.9, sigma = -5))

params <- c("mu","sigma","phi")
niter <- 5e3L
nburnin <- 1e3L
nchains <- 4L

mcmc_out <- nimbleMCMC(code = stochVol_model,
                          monitors = params,
                          niter = niter,
                          nburnin = nburnin,
                          nchains = nchains)

CstochVol_model <- compileNimble(stochVol_model)
CstochVol_model$simulate()
nimbleM

mean(mcmc_out$chain1[,'phi'])
