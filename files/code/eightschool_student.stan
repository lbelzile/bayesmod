data {
  int<lower=0> N;
  vector[N] y;
  vector[N] se_y;
}
parameters {
  vector[N] eta;
  real mu;
  real<lower=0> tau;
}
model {
  mu ~ normal(0, 100);
  eta ~ normal(0, tau);
  tau ~ student_t(3, 0, 5); 
  y ~ normal(mu + eta, se_y);
}
