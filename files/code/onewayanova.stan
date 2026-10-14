data {
  int<lower=1> N;
  int<lower=1> K;
  vector[N] y;
  array[N] int<lower=0, upper=K> group;
}

parameters {
  vector[K] alpha;
  real<lower=0> sigma;
  real<lower=0> tau;
  real mu;
}

model {
  alpha ~ normal(0, tau);
  y ~ normal(mu + alpha[group], sigma);
  sigma ~ student(0, 1, 3);
  tau ~ student(0, 1, 3);
}
