data {
  int<lower=0> T;   // # sample size
  vector[T] y;      // mean zero return at time t
}
parameters {
  real mu;                     // mean log volatility
  real<lower=-1, upper=1> phi; // persistence of volatility, force stationarity
  real<lower=0> sigma;         // white noise shock scale
  vector[T] h_std;             // standard Gaussian
}
transformed parameters {
  vector[T] h = h_std * sigma;  
  h[1] /= sqrt(1 - phi * phi);  // unconditional distribution
  h += mu;
  for (t in 2:T) {
    h[t] += phi * (h[t - 1] - mu);
  }
}
model {
  phi ~ uniform(-1, 1);
  sigma ~ cauchy(0, 5);
  mu ~ cauchy(0, 10);
  h_std ~ normal(0, 1); 
  y ~ normal(0, exp(h / 2));
  
}
