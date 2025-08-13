data {
  int<lower=1> Y;                   // Number of years
  vector[Y] S_obs;                  // (log) observed spawners
  vector[Y] R_obs;                  // (log) observed recruits
  vector[Y] sigma_S_obs;            // (log) annual uncertainty for spawners
  vector[Y] sigma_R_obs;            // (log) annual uncertainty for recruits
  real<lower=0> alpha_prior;        // plausible maximum recruits per spawner
  real<lower=0> Rmax_noenh_prior;   // plausible maximum recruitment (not enhanced)
  real<lower=0> Rmax_enh_prior;     // plausible maximum recruitment (enhanced)
  int<lower=0,upper=1> enhanced[Y]; // Binary indicator for enhancement (0 = no, 1 = yes)
}

parameters {
  real<lower=0> alpha_noenh;    // Productivity parameter (no enhancement)
  real<lower=0> alpha_enh;      // Productivity parameter (enhanced)
  real<lower=0> Rmax_noenh;     // Maximum recruitment (no enhancement)
  real<lower=0> Rmax_enh;       // Maximum recruitment (enhanced)
  real<lower=-1,upper=1> phi;   // AR1 correlation coefficient
  real log_sigma_R_proc;        // Process error standard deviation (log space)
  real mu_S;                    // Mean of log(S_true)
  real<lower=0> sigma_S;        // Latent state variation (spawners)
  vector[Y] S_raw;              // Standard normal latent variable for spawners
  vector[Y] z_v;                // Non-centered parameterization of v
  real<lower=0> k_R;            // Multiplicative scaling for reported SD
  real<lower=0> tau_R;          // additive obs SD on log scale
}

transformed parameters {
  vector[Y] S_true;             // Latent true spawner values
  vector[Y] R_true;             // Latent true recruit values
  vector[Y] lnresid;            // Log-scale recruitment residuals
  vector[Y] v;                  // AR1 process error
  real<lower=0> sigma_R_proc;   // Natural scale process error standard deviation
  real<lower=0> beta_noenh;
  real<lower=0> beta_enh;

  // Explicitly encode Rmax as prior on beta
  beta_noenh = Rmax_noenh;
  beta_enh = Rmax_enh;

  // Non-centered parameterization of latent recruitment deviations
  sigma_R_proc = exp(log_sigma_R_proc);
  v[1] = z_v[1] * sigma_R_proc;
  for (t in 2:Y) {
    v[t] = phi * v[t-1] + z_v[t] * sigma_R_proc;
  }

  // Non-centered parameterization for spawner latent states
  S_true = exp(mu_S + sigma_S * S_raw);

  // Beverton-Holt model for expected recruitment with AR1 deviations
  // Uses third parameterization (eqn 7.5.3) from Hilborn & Walters 1992 (p. 258)
  for (t in 1:Y) {
    real mu_R_t;
    if(enhanced[t] == 0) {
      mu_R_t = (alpha_noenh * S_true[t]) / (1 + (alpha_noenh / beta_noenh) * S_true[t]);
    } else {
      mu_R_t = (alpha_enh * S_true[t]) / (1 + (alpha_enh / beta_enh) * S_true[t]);
    }
    R_true[t] = mu_R_t * exp(v[t]);
    lnresid[t] = log(R_true[t]) - log(mu_R_t);
  }
}

model {
  // Priors
  alpha_enh ~ lognormal(log(alpha_prior), 0.5);
  alpha_noenh ~ lognormal(log(alpha_prior), 0.5);
  Rmax_enh ~ lognormal(log(Rmax_enh_prior), 0.5);
  Rmax_noenh ~ lognormal(log(Rmax_noenh_prior), 0.5);
  phi ~ normal(0, 0.5);
  log_sigma_R_proc ~ normal(log(4), 0.3);
  mu_S ~ normal(mean(S_obs), 1);
  sigma_S ~ std_normal() T[0,];
  S_raw ~ std_normal();
  z_v ~ std_normal();
  k_R ~ lognormal(log(1), 0.6);
  tau_R ~ std_normal() T[0,]; 

  // Observation model with adjusted R_obs likelihood
  S_obs ~ normal(log(S_true), sigma_S_obs);
  for (t in 1:Y) {
    real obs_sd = sqrt( (k_R * sigma_R_obs[t])^2 + tau_R^2 );
    R_obs[t] ~ normal(log(R_true[t]), obs_sd);
  }
}
