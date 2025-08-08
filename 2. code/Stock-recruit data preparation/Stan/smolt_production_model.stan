data {
  int<lower=1> Y;                 // Number of years
  array[Y] real<lower=0> N_obs;   // Observed total lake abundance
  array[Y] int<lower=0> A1_obs;   // Observed age-1 migrants in samples
  array[Y] int<lower=0> A_total;  // Total aged fish in outmigration samples
  array[Y] int<lower=0> a1_obs;   // Observed age-1 fry counts (subset of years)
  array[Y] int<lower=0> a_total;  // Total fry aged (subset of years)
  
  // Weight data for outmigrants:
  // For age-1 outmigrants (O1)
  array[Y] real WO1_ln_mean;   // Observed log mean body weight for O1 fish in each year
  array[Y] real<lower=0> WO1_ln_sd;   // Observed log SD of body weight for O1 fish
  array[Y] int<lower=0> WO1_n;      // Sample size for O1 weight data
  // For age-2 outmigrants (O2)
  array[Y] real WO2_ln_mean;   // Observed mean body weight for O2 fish in each year
  array[Y] real<lower=0> WO2_ln_sd;   // Observed SD of body weight for O2 fish
  array[Y] int<lower=0> WO2_n;      // Sample size for O2 weight data
  
  real<lower=0> obs_error_prior;          // obs error on total lake abundance
  real mu_theta_prior;                    // Mean prior for smolting rate
  real<lower=0> sigma_theta_prior;        // SD prior for smolting rate
  real mu_M_prior;                        // Mean prior for mortality of age 2s in year 2
  real<lower=0> sigma_M_prior;            // SD prior for mortality
  real<lower=0> sigma_theta_sd_prior;     // SD prior for sigma_theta
  real<lower=0> sigma_M_sd_prior;         // SD prior for sigma_M
  real<lower=0> mu_N2_init_prior;         // prior for mean of N2s in first year
  real<lower=0> sigma_N2_init_prior;      // prior for sd of N2s in first year
  real<lower=0> mu_N1_prior;              // prior for mean of age1 fry abundance
  real<lower=0> sigma_N1_prior;           // prior for sd of age1 fry abundance
  //real<lower=0> w1_theta_prior;           // prior for (log) mean of age1 smolt weight
  //real<lower=0> w1_sigma_prior;           // prior for (log) sd of age1 smolt weight
  //real<lower=0> w2_theta_prior;           // prior for (log) mean of age2 smolt weight
  //real<lower=0> w2_sigma_prior;           // prior for (log) sd of age2 smolt weight

  
  // Companion vectors for data presence/absence:
  array[Y] int is_observed_ats;   // = 1 if N_obs data are observed in that year
  array[Y] int is_observed_count; // = 1 if A_total and A1_obs are observed in that year
  array[Y] int is_observed_fry;   // = 1 if fry age data is observed in that year
  array[Y] int is_observed_WO1;   // = 1 if age1 smolt weight data are observed in that year
  array[Y] int is_observed_WO2;   // = 1 if age2 smolt weight data are observed in that year
                                  // Note: even if no age-2 smolts were observed, we still
                                  // need estimated weights because the model will likely estimate
                                  // that a small fraction of outmigrants were age-2s
}
 
parameters {
  real mu_theta;                     // Population-level mean of age-1 outmigration proportion in logit space
  real mu_M;                         // Population-level mean of mortality in logit space
  real<lower=0> mu_N1;               // Population-level mean of age-1 fry abundance
 
  real<lower=0> sigma_theta;         // Yearly variation in outmigration
  real<lower=0> sigma_M;             // Yearly variation in mortality
  real<lower=0> sigma_N1;            // Yearly variation in age-1 fry abundance
 
  array[Y] real theta_raw;           // Yearly outmigration deviation - z score
  array[Y] real M_raw;               // Yearly mortality deviation - z score
  array[Y] real log_N1_raw;          // Yearly age-1 fry abundance - z score
  real log_obs_error;                // Observation error 
 
  real<lower=0> N2_init;             // initialization for age 2 fish in year 1

  vector[Y] w1_ln;                   // True (log) weight for O1 fish in each year
  vector[Y] w2_ln;                   // True (log) weight for O2 fish in each year
  
  // Regression coefficients for predictive relationship between abundance and weight
  real alpha_w1;                     
  real beta_w1;
  real<lower=0> w_sigma1;
  real alpha_w2;
  real beta_w2;
  real<lower=0> w_sigma2;

}
 
// deterministic statements
transformed parameters {                 
  array[Y] real<lower=0, upper=1> theta; // proportion of age1s that smolt "smolting rate"
  array[Y] real<lower=0, upper=1> M;     // Second winter mortality for holdovers
  array[Y] real<lower=0> N1;             // Age-1 fry in lake
  array[Y] real<lower=0> N_lake;         // Total fry in lake
  array[Y] real<lower=0> N2;             // Age-2 fry in lake
  array[Y] real<lower=0> O1;             // Age-1 outmigrants
  array[Y] real<lower=0> O2;             // Age-2 outmigrants
  array[Y] real<lower=0, upper=1> p1;    // Age-1 smolt proportion
  array[Y] real<lower=0, upper=1> p2;    // Age-2 smolt proportion
  array[Y] real<lower=0> BO1;            // Biomass of age-1 outmigrants
  array[Y] real<lower=0> BO2;            // Biomass of age-2 outmigrants 
  real<lower=0> obs_error = exp(log_obs_error);

  
  // Hierarchical transformation
  for (y in 1:Y) {
    theta[y] = inv_logit(mu_theta + sigma_theta * theta_raw[y]); // mu_theta in logit space 
    M[y] = inv_logit(mu_M + sigma_M * M_raw[y]);
    N1[y] = exp(mu_N1 + sigma_N1 * log_N1_raw[y]);
  }
  

  // State dynamics
  for (y in 1:Y) {
    if (y == 1) {
      N2[y] = N2_init;  // No age-2 data in first year
    } else {
      N2[y] = N1[y-1] * (1 - theta[y-1]) * (1 - M[y-1]); // calculation of age2s in y+1
    }
    // Implicit calculation from total; *key assumption*
    N_lake[y] = N1[y] + N2[y];
    // Outmigrants
    O1[y] = theta[y] * N1[y];
    O2[y] = N2[y]; // all age2s outmigrate
    // Annual age proportions of smolts
    p1[y] = O1[y] / N_lake[y];
    p2[y] = O2[y] / N_lake[y];
    // Biomass: multiply abundance by true weight
    BO1[y] = O1[y] * exp(w1_ln[y]); // ensure these are inputted in log space
    BO2[y] = O2[y] * exp(w2_ln[y]);
  }
}

 
model {
  // Priors
  log_obs_error ~ normal(log(obs_error_prior), 0.5);
  mu_theta ~ normal(mu_theta_prior, sigma_theta_prior);
  mu_M ~ normal(mu_M_prior, sigma_M_prior);
  // Use half-normal distributions
  sigma_theta ~ normal(0, sigma_theta_sd_prior) T[0,];
  sigma_M ~ normal(0, sigma_M_sd_prior) T[0,];
  // total lake age1 fry abundance (in log space)
  mu_N1 ~ normal(mu_N1_prior, sigma_N1_prior);
  sigma_N1 ~ normal(0, 1) T[0,];
  // initialization of age 2s
  // Priors must be transformed to log space in R first, as per:
  // https://msalganik.wordpress.com/2017/01/21/making-sense-of-the-rlnorm-function-in-r
  N2_init ~ lognormal(mu_N2_init_prior, sigma_N2_init_prior);
  // Hierarchical deviations - z scores
  theta_raw ~ std_normal();
  M_raw ~ std_normal();
  log_N1_raw ~ std_normal();
  // Priors on regression parameters for abundance-weight relationship
  alpha_w1 ~ normal(0, 2);
  beta_w1  ~ normal(0, 1);  // Expect negative slope (weight declines as abundance increases)
  w_sigma1 ~ normal(0, 1) T[0,];
  alpha_w2 ~ normal(0, 2);
  beta_w2  ~ normal(0, 1);
  w_sigma2 ~ normal(0, 1) T[0,];
  // Predictive priors on w1 and w2 based on log abundance of O1 and O2
  for (y in 1:Y) {
    w1_ln[y] ~ normal(alpha_w1 + beta_w1 * log(O1[y] + 1e-6), w_sigma1);
    w2_ln[y] ~ normal(alpha_w2 + beta_w2 * log(O2[y] + 1e-6), w_sigma2);
  }

  // Observation model: Negative binomial for total lake abundance (ATS est)
  for (y in 1:Y) {
    if (is_observed_ats[y] == 1) {
      //N_obs[y] ~ lognormal(log(N_lake[y]) - obs_error^2/2, obs_error);
      N_obs[y] ~ lognormal(log(N_lake[y]), obs_error); // try without bias correction
    }
  }
  
  // Observation model: Binomial for age composition in outmigrants
  for (y in 1:Y) { 
    if (is_observed_count[y] == 1) {
      A1_obs[y] ~ binomial(A_total[y], O1[y] / (O1[y] + O2[y])); // total number of observed age1s in outmigration
    }
  }
  
  // Observation model: Fry age composition (in lake) where available
  for (y in 1:Y) {
    if (is_observed_fry[y] == 1) {
    a1_obs[y] ~ binomial(a_total[y], N1[y] / N_lake[y]);
    }
  }
  
  // Observation models for weight data for outmigrants
  // Ignores missing observations and uninformative (rare) cases where
  // sample size for weights is 1 and therefore SD = 0
  for (y in 1:Y) {
    if (is_observed_WO1[y] == 1 && WO1_n[y] > 1 && WO1_ln_sd[y] > 0) {
    WO1_ln_mean[y] ~ normal(w1_ln[y], WO1_ln_sd[y] / sqrt(WO1_n[y]));
    }
    if (is_observed_WO2[y] == 1 && WO2_n[y] > 1 && WO2_ln_sd[y] > 0) {
    WO2_ln_mean[y] ~ normal(w2_ln[y], WO2_ln_sd[y] / sqrt(WO2_n[y]));
    }
  }
}
