data {
  int<lower=1> N;                      // number of observations
  int<lower=1> L;                      // number of lakes
  int<lower=1> Y;                      // number of years (unique year ids)
  int<lower=1> K;                      // number of spline basis functions
  int<lower=1,upper=L> lake[N];        // lake index for each obs (1..L)
  int<lower=1,upper=Y> year[N];        // year index for each obs (1..Y)
  matrix[N, K] B;                      // basis matrix evaluated at each observation's day
  vector[N] log_obs;                   // observed log(presmolt_est)
  vector<lower=0>[N] log_se;           // reported log-scale sd for each observation
  int<lower=0,upper=1> is_observed[N]; // 1 if row i has an observed log_obs, 0 for prediction row

}

parameters {
  // seasonal coefficients, one vector per lake
  matrix[L, K] beta_raw;          // unscaled coeffs
  real<lower=0> sigma_beta;       // prior scale for seasonal coefficients

  // yearly multivariate state
  matrix[Y, L] z_raw;                 // N(0, 1) innovations

  // MV covariance factorization for Sigma: Sigma = diag(sig_y) * R * diag(sig_y)
  vector<lower=0>[L] sig_y;                 // marginal SD per lake for year-effects
  cholesky_factor_corr[L] L_R;              // cholesky factor of correlation matrix

  // AR1 coefficient
  real<lower=-0.99, upper=0.99> phi;

  // extra observation noise (on log scale)
  real<lower=0> sigma_obs;
  
  // Student-t degrees of freedom for observation model
  real<lower=2> nu_obs;
}

transformed parameters {
  // seasonal coefficients with scale
  matrix[L, K] beta;
  matrix[Y, L] z;
  
  // Cholesky for Sigma
  matrix[L, L] L_Sigma;
  L_Sigma = diag_pre_multiply(sig_y, L_R); // L_Sigma * L_Sigma' = Sigma

  for (l in 1:L)
    for (k in 1:K)
      beta[l, k] = sigma_beta * beta_raw[l, k];
      
  // first year:
  for (l in 1:L) z[1,l] = 0 + (L_Sigma * to_vector(z_raw[1]))[l];
  // AR1 with non-centered innovations:
  for (t in 2:Y) {
    vector[L] innov = L_Sigma * to_vector(z_raw[t]);
    z[t] = phi * z[t-1] + innov';
  }

}

model {
  // ---------- PRIORS ----------
  // beta_raw ~ N(0,1); sigma_beta sets magnitude
  to_vector(beta_raw) ~ normal(0, 1);
  to_vector(z_raw) ~ normal(0, 1);
  sigma_beta ~ student_t(3, 0, 1.0);   // tuneable

  // correlation prior
  L_R ~ lkj_corr_cholesky(2.0);

  // marginal SDs for year effects
  sig_y ~ student_t(3, 0, 1.0);        // tuneable, maybe bigger if log range large

  // AR1 prior
  phi ~ normal(0, 0.5);

  // sigma_obs prior
  sigma_obs ~ student_t(3, 0, 0.5); // adjusted 3rd arg from 0.2 to 0.5...
  
  // Robust-likelihood prior
  nu_obs ~ gamma(2, 0.1);

  // ---------- STATE DYNAMICS ----------
  // prior for first year
  {
    vector[L] mu0 = rep_vector(0.0, L);
  }

  // ---------- LIKELIHOOD ----------
  for (i in 1:N) {
    if (is_observed[i] == 1) {
      int l = lake[i];
      int t = year[i];
      real mu_i = dot_product(row(B, i), beta[l]') + z[t, l];
      real total_sd = sqrt(square(log_se[i]) + square(sigma_obs));
      
      log_obs[i] ~ student_t(nu_obs, mu_i, total_sd);
    }
  }
}

generated quantities {
  vector[N] log_true_all;
  
  for (i in 1:N) {
    int l = lake[i];
    int t = year[i];
    real s = dot_product(row(B, i), beta[l]');
    log_true_all[i] = s + z[t, l];
  }
}
