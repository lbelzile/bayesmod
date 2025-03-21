data {
  int<lower=1> N;
  int<lower=1> K;
  int<lower=2> L;
  vector[N] y;
  array[N] int<lower=0, upper=K> RCP;
  array[N] int<lower=0, upper=L> GCM;
  ordered[K] mu;
}

parameters {
  vector[K] alpha;
  sum_to_zero_vector[L] beta;
  real<lower=0> sigma;
  real<lower=0> tau;
  real<lower=0> omega;
}

model {
  y ~ normal(alpha[RCP] + beta[GCM], sigma);
  alpha ~ normal(mu, omega);
  beta ~ normal(0, sqrt(L * inv(L - 1)) * tau);
  // TODO: Add priors for sigma, tau, omega here
}
