data {
  int<lower=1> Y;  // Number of years
  vector[Y] S_obs; // Observed spawners
  vector[Y] R_obs; // Observed recruits
  vector[Y] sigma_S_obs; // Annual uncertainty for spawners
  vector[Y] sigma_R_obs; // Annual uncertainty for recruits
}

parameters {
  real<lower=0> alpha;        // Productivity parameter
  real<lower=0> beta;         // Density dependence parameter
  real<lower=0> sigma_proc;   // Process error standard deviation
  real<lower=-1,upper=1> phi; // Autoregression coefficient
  real mu_S;                  // Mean of log(S_true)
  real mu_R;                  // Mean of log(R_true)
  real<lower=0> sigma_S;       // Latent state variation (spawners)
  real<lower=0> sigma_R;       // Latent state variation (recruits)
  vector[Y] S_raw;            // Standard normal latent variable for spawners
  vector[Y] R_raw;            // Standard normal latent variable for recruits
}

transformed parameters {
  vector[Y] S_true;  // Latent true spawner values
  vector[Y] R_true;  // Latent true recruit values

  // Non-centered parameterization for latent states
  S_true = exp(mu_S + sigma_S * S_raw);
  R_true = exp(mu_R + sigma_R * R_raw);
}

model {
  // Priors
  alpha ~ normal(2, 1);
  beta ~ lognormal(log(5e-4), 0.5);
  sigma_proc ~ cauchy(0, 1);
  phi ~ normal(0, 0.5);
  mu_S ~ normal(log(mean(S_obs)), 1);
  mu_R ~ normal(log(mean(R_obs)), 1);
  sigma_S ~ cauchy(0, 1);
  sigma_R ~ cauchy(0, 1);
  S_raw ~ normal(0, 1);
  R_raw ~ normal(0, 1);

  // Likelihood for measurement error with annual uncertainty
  S_obs ~ normal(log(S_true), sigma_S_obs);
  R_obs ~ normal(log(R_true), sigma_R_obs);

  // Process model with AR(1) structure
  R_true[1] ~ normal((alpha * S_true[1]) / (1 + beta * S_true[1]), sigma_proc);
  for (t in 2:Y) {
    real mu_R_t = (alpha * S_true[t]) / (1 + beta * S_true[t]);
    R_true[t] ~ normal(mu_R_t + phi * (R_true[t-1] - mu_R_t), sigma_proc);
  }
}
