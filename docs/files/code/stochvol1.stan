data {
  int<lower=0> T;   // # sample size
  vector[T] y;      // mean zero return at time t
}
parameters {
  real mu;                     // mean log volatility
  real<lower=-1, upper=1> phi; // persistence of volatility, force stationarity
  real<lower=0> sigma;         // white noise shock scale
  vector[T] h;                 // log volatility at time t
}
model {
  phi ~ uniform(-1, 1);
  sigma ~ cauchy(0, 5);
  mu ~ cauchy(0, 10);
  h[1] ~ normal(mu, sigma / sqrt(1 - phi * phi)); // unconditional distribution
  for (t in 2:T) {
    h[t] ~ normal(mu + phi * (h[t - 1] -  mu), sigma);
  }
  y ~ normal(0, exp(h / 2));
}
