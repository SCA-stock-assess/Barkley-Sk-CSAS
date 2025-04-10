// Code generated by ChatGPT based on description of the problem and 
// uploading the raw data file

data {
  int<lower=1> T;                        // Number of years
  vector<lower=0>[T] rate_age1;            // Observed rate for Age1 fish in year t (counts per sampling day)
  vector<lower=0>[T] rate_age2;            // Observed rate for Age2 fish in year t (for t>=2, holdover only)
  vector<lower=1>[T] sample_days;          // Sampling effort (days) in each year
  int<lower=0, upper=1> is_observed[T];    // Indicator if data are observed in year t
  matrix<lower=0>[T, 2] alpha_prior;       // Informative prior for initial state (2 components: [Age1 leaving, holdover])
  real<lower=0> sigma_prior;               // Lake-specific prior for process noise sigma
}

parameters {
  simplex[2] theta[T];                     // Latent state for each year:
                                           // theta[t,1] = proportion leaving as Age1 (immediate smolts)
                                           // theta[t,2] = holdover proportion (those that remain in the lake)
  real<lower=0> sigma;                     // Process noise (affects state evolution and observation variability)
}

model {
  // Prior for the initial state (year 1) using adult returns as informative prior.
  theta[1] ~ dirichlet(alpha_prior[1]);
  
  // State evolution: Dirichlet random walk on the 2-component simplex.
  for (t in 2:T) {
    theta[t] ~ dirichlet(theta[t-1] * sigma + 1e-6);
  }
  
  // Observation model for Age1:
  // Observed rate_age1 in year t is modeled as lognormally distributed with mean proportional to theta[t,1] * sample_days[t].
  for (t in 1:T) {
    if (is_observed[t] == 1) {
      rate_age1[t] ~ lognormal(log(theta[t,1] * sample_days[t]), sigma);
    }
  }
  
  // Observation model for Age2:
  // For t >= 2, the observed rate_age2 comes solely from the holdover (theta[t-1,2]) of the previous year.
  for (t in 2:T) {
    if (is_observed[t] == 1) {
      rate_age2[t] ~ lognormal(log(theta[t-1,2] * sample_days[t]), sigma);
    }
  }
  
  // Prior for process noise sigma:
  sigma ~ normal(sigma_prior, 0.5) T[0,];
}

generated quantities {
  // The smolting rate for each year is defined as the fraction of juvenile fish that leave as Age1,
  // which is the first component of theta.
  vector[T] smolting_rate;
  for (t in 1:T) {
    smolting_rate[t] = theta[t,1];
  }
}

