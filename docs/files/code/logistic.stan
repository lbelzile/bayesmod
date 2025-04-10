data {
  int<lower=0> N;
  int<lower=0> p;
  matrix[N, p] X;
  array[N] int<lower=0, upper=1> y;
}
parameters {
  vector[p] beta;
}
model {
  y ~ bernoulli_logit(X * beta);
}
