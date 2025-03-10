// generated with brms 2.22.0
functions {
  /* Efficient computation of the horseshoe scale parameters
   * see Appendix C.1 in https://projecteuclid.org/euclid.ejs/1513306866
   * Args:
   *   lambda: local shrinkage parameters
   *   tau: global shrinkage parameter
   *   c2: slap regularization parameter
   * Returns:
   *   scale parameter vector of the horseshoe prior
   */
  vector scales_horseshoe(vector lambda, real tau, real c2) {
    int K = rows(lambda);
    vector[K] lambda2 = square(lambda);
    vector[K] lambda_tilde = sqrt(c2 * lambda2 ./ (c2 + tau^2 * lambda2));
    return lambda_tilde * tau;
  }
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int<lower=1> Kc;  // number of population-level effects after centering
  int<lower=1> Kscales;  // number of local scale parameters
  // data for the horseshoe prior
  real<lower=0> hs_df;  // local degrees of freedom
  real<lower=0> hs_df_global;  // global degrees of freedom
  real<lower=0> hs_df_slab;  // slab degrees of freedom
  real<lower=0> hs_scale_global;  // global prior scale
  real<lower=0> hs_scale_slab;  // slab prior scale
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] zb;  // unscaled coefficients
  real Intercept;  // temporary intercept for centered predictors
  // horseshoe shrinkage parameters
  real<lower=0> hs_global;  // global shrinkage parameter
  real<lower=0> hs_slab;  // slab regularization parameter
  vector<lower=0>[Kscales] hs_local;  // local parameters for the horseshoe prior
  real<lower=0> sigma;  // dispersion parameter
}
transformed parameters {
  vector[Kc] b;  // scaled coefficients
  vector<lower=0>[Kc] sdb;  // SDs of the coefficients
  vector<lower=0>[Kscales] scales;  // local horseshoe scale parameters
  real lprior = 0;  // prior contributions to the log posterior
  // compute horseshoe scale parameters
  scales = scales_horseshoe(hs_local, hs_global, hs_scale_slab^2 * hs_slab);
  sdb = scales[(1):(Kc)];
  b = zb .* sdb;  // scale coefficients
  lprior += student_t_lpdf(Intercept | 3, 140.5, 88.2);
  lprior += student_t_lpdf(hs_global | hs_df_global, 0, hs_scale_global * sigma)
    - 1 * log(0.5);
  lprior += inv_gamma_lpdf(hs_slab | 0.5 * hs_df_slab, 0.5 * hs_df_slab);
  lprior += student_t_lpdf(sigma | 3, 0, 88.2)
    - 1 * student_t_lccdf(0 | 3, 0, 88.2);
}
model {
  // likelihood including constants
  if (!prior_only) {
    target += normal_id_glm_lpdf(Y | Xc, Intercept, b, sigma);
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(zb);
  target += student_t_lpdf(hs_local | hs_df, 0, 1)
    - rows(hs_local) * log(0.5);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
}
