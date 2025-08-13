data {
  int<lower=1> Y;              // Number of years
  vector[Y] S_obs;             // (log) observed spawners
  vector[Y] R_obs;             // (log) observed recruits
  vector[Y] sigma_S_obs;       // (log) annual uncertainty for spawners
  vector[Y] sigma_R_obs;       // (log) annual uncertainty for recruits
  real<lower=0> alpha_prior;   // plausible maximum recruits per spawner
  real<lower=0> Rmax_prior;    // plausible maximum total recruitment
}

parameters {
  real<lower=0> alpha;          // Productivity parameter
  real<lower=0> Rmax;           // Maximum recruitment
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
  real<lower=0> beta;           // Beverton-Holt beta parameter
  real<lower=0> sigma_R_proc;   // Natural scale process error standard deviation
  
  // Make Rmax prior information explicit
  beta = Rmax;
  
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
    real mu_R_t = (alpha * S_true[t]) / (1 + (alpha / beta) * S_true[t]);
    R_true[t] = mu_R_t * exp(v[t]);
    lnresid[t] = log(R_true[t]) - log(mu_R_t);
  }
}

model {
  // Priors
  alpha ~ lognormal(log(alpha_prior), 0.5);
  Rmax ~ lognormal(log(Rmax_prior), 0.5);
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
