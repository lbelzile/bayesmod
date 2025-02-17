data {
  int<lower=0> N;
  int<lower=1> K; // number of fixed effect levels for distraction type
  int<lower=1> J; // number of random effect levels for individual
  array[N] int<lower=0> y;    // vector of data
  array[N] int<lower=1, upper=J> id;
  array[N] int<lower=1, upper=K> fixed;
}
parameters {
  vector[K] beta;
   sum_to_zero_vector[J] alpha;
  real<lower=0> kappa;
}

model {
  beta ~ normal(0, 10);
  alpha ~ normal(0, sqrt(J * inv(J - 1)) * kappa);
  kappa ~ exponential(0.6);
  y ~ poisson_log(alpha[id] + beta[fixed]);
}

