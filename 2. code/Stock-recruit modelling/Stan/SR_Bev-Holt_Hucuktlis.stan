data {
  int<lower=1> Y;  // Number of years
  vector[Y] S_obs; // Observed spawners
  vector[Y] R_obs; // Observed recruits
  vector[Y] sigma_S_obs; // Annual uncertainty for spawners
  vector[Y] sigma_R_obs; // Annual uncertainty for recruits
  int<lower=0,upper=1> enhanced[Y]; // Binary indicator for enhancement (0 = no, 1 = yes)
  real<lower=0> alpha_prior; 
  real<lower=0> beta_prior;
  real<lower=0> sigma_alpha_prior;
  real<lower=0> sigma_beta_prior;
}

parameters {
  real<lower=0> alpha_noenh;   // Productivity parameter (no enhancement)
  real<lower=0> beta_noenh;    // Density dependence parameter (no enhancement)
  real<lower=0> alpha_enh;     // Productivity parameter (enhanced)
  real<lower=0> beta_enh;      // Density dependence parameter (enhanced)
  real<lower=0> sigma_proc;    // Process error standard deviation
  real mu_S;                   // Mean of log(S_true)
  real mu_R;                   // Mean of log(R_true)
  real<lower=0> sigma_S;       // Latent state variation (spawners)
  real<lower=0> sigma_R;       // Latent state variation (recruits)
  vector[Y] S_raw;             // Standard normal latent variable for spawners
  vector[Y] R_raw;             // Standard normal latent variable for recruits
}

transformed parameters {
  vector[Y] S_true;         // Latent true spawner values
  vector[Y] R_true;         // Latent true recruit values
  vector[Y] lnresid;        // Log-scale recruitment residuals

  // Non-centered parameterization for latent states
  S_true = exp(mu_S + sigma_S * S_raw);
  R_true = exp(mu_R + sigma_R * R_raw);

  // Calculate log residuals
  for (t in 1:Y) {
    real expected_R_t;
    if (enhanced[t] == 0) {
      expected_R_t = (alpha_noenh * S_true[t]) / (1 + (alpha_noenh / beta_noenh) * S_true[t]);
    } else {
      expected_R_t = (alpha_enh * S_true[t]) / (1 + (alpha_enh / beta_enh) * S_true[t]);
    }
    lnresid[t] = log(R_true[t]) - log(expected_R_t);
  }
}

model {
  // Priors
  alpha_noenh ~ lognormal(log(alpha_prior), sigma_alpha_prior);
  beta_noenh  ~ lognormal(log(beta_prior), sigma_beta_prior);
  alpha_enh   ~ lognormal(log(alpha_prior), sigma_alpha_prior);
  beta_enh    ~ lognormal(log(beta_prior), sigma_beta_prior);
  sigma_proc ~ cauchy(0, 1);
  mu_S ~ normal(log(mean(S_obs)), 1);
  mu_R ~ normal(log(mean(R_obs)), 1);
  sigma_S ~ cauchy(0, 1);
  sigma_R ~ cauchy(0, 1);
  S_raw ~ normal(0, 1);
  R_raw ~ normal(0, 1);

  // Observation model
  S_obs ~ normal(log(S_true), sigma_S_obs);
  R_obs ~ normal(log(R_true), sigma_R_obs * 0.25);

  // Process model: Beverton-Holt structure with enhancement effect
  if (enhanced[1] == 0) {
    R_true[1] ~ normal((alpha_noenh * S_true[1]) / 
                       (1 + (alpha_noenh / beta_noenh) * S_true[1]),
                       sigma_proc);
  } else {
    R_true[1] ~ normal((alpha_enh * S_true[1]) / 
                       (1 + (alpha_enh / beta_enh) * S_true[1]),
                       sigma_proc);
  }

  for (t in 2:Y) {
    real mu_R_t;
    if (enhanced[t] == 0) {
      mu_R_t = (alpha_noenh * S_true[t]) / (1 + (alpha_noenh / beta_noenh) * S_true[t]);
    } else {
      mu_R_t = (alpha_enh * S_true[t]) / (1 + (alpha_enh / beta_enh) * S_true[t]);
    }
    R_true[t] ~ normal(mu_R_t, sigma_proc);
  }
}
