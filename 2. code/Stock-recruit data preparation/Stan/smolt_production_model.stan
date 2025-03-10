data {
  int<lower=1> Y;                // Number of years
  array[Y] int<lower=0> N_obs;    // Observed total lake abundance
  array[Y] int<lower=0> A1_obs;   // Observed age-1 migrants in samples
  array[Y] int<lower=0> A_total;  // Total aged fish in outmigration samples
  real<lower=0> obs_error_prior; // obs error on total lake abundance (sigma)
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
  // Companion vector for count data: 
  // = 1 if A_total and A1_obs are observed in that year, = 0 otherwise
  array[Y] int is_observed_count;
}
 
parameters {
  real<lower=0> obs_error;           // Observation error for total lake abundance
  real mu_theta;                     // Population-level mean of age 1 outmigration proportion in logit space
  real mu_M;                         // Population-level mean of mortality in logit space
 
  real<lower=0> sigma_theta;         // Yearly variation in outmigration
  real<lower=0> sigma_M;             // Yearly variation in mortality
 
  array[Y] real theta_raw;           // Yearly outmigration deviation - z score
  array[Y] real M_raw;               // Yearly mortality deviation - z score
 
  real<lower=0> N2_init;             // initialization for age 2 fish in year 1
  array[Y] real<lower=0> N1;         // Yearly total age1 fry abundance
 
}
 
// deterministic statements
transformed parameters {                 
  array[Y] real<lower=0, upper=1> theta; // proportion of age1s that smolt "smolting rate"
  array[Y] real<lower=0, upper=1> M;     // Second winter mortality for holdovers
 
  array[Y] real<lower=0> N_lake;         // Total fry in lake
  array[Y] real<lower=0> N2;             // Age-2 fry in lake
  array[Y] real<lower=0> O1;             // Age-1 outmigrants
  array[Y] real<lower=0> O2;             // Age-2 outmigrants
  array[Y] real<lower=0, upper=1> p1;    // Age-1 smolt proportion
  array[Y] real<lower=0, upper=1> p2;    // Age-2 smolt proportion
 
  
  // Hierarchical transformation
  for (y in 1:Y) {
    theta[y] = inv_logit(mu_theta + sigma_theta * theta_raw[y]); // mu_theta in logit space 
    M[y] = inv_logit(mu_M + sigma_M * M_raw[y]);
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
  }
}
 
model {
  // Priors
  obs_error ~ exponential(obs_error_prior);
 
  mu_theta ~ normal(mu_theta_prior, sigma_theta_prior);
  mu_M ~ normal(mu_M_prior, sigma_M_prior);
  // Use half-normal distributions
  sigma_theta ~ normal(0, sigma_theta_sd_prior) T[0,];
  sigma_M ~ normal(0, sigma_M_sd_prior) T[0,];
 
  // initialization of age 2s
  // Priors must be transformed to log space in R first, as per:
  // https://msalganik.wordpress.com/2017/01/21/making-sense-of-the-rlnorm-function-in-r
  N2_init ~ lognormal(mu_N2_init_prior, sigma_N2_init_prior);
  // total lake age1 fry abundance
  N1 ~ lognormal(mu_N1_prior, sigma_N1_prior);
  // Hierarchical deviations - z scores
  theta_raw ~ normal(0,1);
  M_raw ~ normal(0,1);
 
  // Observation model: Negative binomial for total lake abundance (ATS est)
  for (y in 1:Y) {
    N_obs[y] ~ lognormal(log(N_lake[y]) - obs_error^2/2, obs_error);
  }
  
  // Observation model: Binomial for age composition in outmigrants
  //for(y in 1:N_obs) y = index number of which years you have observations
  for (y in 1:Y) { 
    if (is_observed_count[y] == 1) {
      A1_obs[y] ~ binomial(A_total[y], O1[y] / (O1[y] + O2[y])); // total number of observed age1s in outmigration
    }
  }
}
