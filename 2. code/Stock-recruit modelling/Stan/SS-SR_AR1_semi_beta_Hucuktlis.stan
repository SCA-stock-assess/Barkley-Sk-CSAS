// Stan version of age-structured state-space spawner-recruitment model with AR-1 process variation (adapted from Fleischman et al. CJFAS. 2013)

data{
  int nyrs;           // number of calender years
  int a_min;          // minimum age class
  int a_max;          // maximum age class
  int A;              // number of age classes
  int nRyrs;          // number of recruitment years
  int A_obs[nyrs, A]; // observed age composition in counts by age class
  vector[nyrs] S_obs; // observed spawners
  vector[nyrs] H_obs; // observed harvest
  vector[nyrs] S_cv;  // spawner observation error CV
  vector[nyrs] H_cv;  // harvest observation error CV
  int f[nyrs];     // binary variable indicating whether lake was fertilized in a given year
  real Smax_p;        // prior for Smax
  real Smax_p_sig;    // prior for Smax
  vector[nyrs] use;   // binary variable determining whether to use data for state space model

}

transformed data{
  real Smax_p_corr;
  real Smax_p_sig_corr;

  Smax_p_sig_corr = sqrt(log(1+(Smax_p_sig^2)/(Smax_p^2))); //this converts sigma on the untransformed scale to a log scale
  // ChatGPT thinks the below calculation of Smax_p_corr is incorrect and suggested the following alteration:
  // Smax_p_corr = log(Smax_p) - 0.5 * pow(sqrt(log(1 + pow(Smax_p_sig, 2) / pow(Smax_p, 2))), 2);
  Smax_p_corr = log(Smax_p)-0.5*(Smax_p_sig_corr)^2; //convert smax prior to per capita slope - transform to log scale with bias correction
}

parameters{
  // [Is it appropriate for the log parameters to have a lower bound of 0?]
  vector<lower=0>[nRyrs] lnR;             // log recruitment states
  real<lower=0> lnalpha_fert;             // lnalpha for fertilized state
  real<lower=0> lnalpha_unfert;           // lnalpha for unfertilized state
  real<lower=0> Smax_fert;                // Smax for fertilized state
  real<lower=0> Smax_unfert;              // Smax for unfertilized state  
  real<lower=0> sigma_R;                  // process error
  real<lower=0> sigma_R0;                 // process error for first a.max years with no spawner link
  real<lower=-1,upper=1> phi;             // lag-1 correlation in process error
  real lnresid_0;                         // first residual for lag-1 correlation in process error
  real<lower=0> mean_ln_R0;               // "true" mean log recruitment in first a.max years with no spawner link
  vector<lower=0.01,upper=0.99>[nyrs] U;  // harvest rate
  vector<lower=0,upper=1>[3] prob;        // maturity schedule probs
  real<lower=0,upper=1> D_scale;          // governs variability of age proportion vectors across cohorts
  matrix<lower=0.01>[nRyrs, A] g;         // individual year/age class gamma variates for generating age at maturity proportions
}

transformed parameters{
  vector<lower=0>[nyrs] N;              // run size states
  vector<lower=0>[nyrs] S;              // spawner states
  vector<lower=0>[nyrs] C;              // catch states
  vector[nyrs] lnS;                     // log spawner states
  vector[nyrs] lnC;                     // log catch states
  vector<lower=0>[nRyrs] R;             // recruitment states
  real<lower=0> sigma_R_corr;           // log-normal bias-corrected process error
  vector[nRyrs] lnresid;                // log recruitment residuals
  vector[nRyrs] lnRm_1;                 // log recruitment states in absence of lag-one correlation
  vector[nRyrs] lnRm_2;                 // log recruitment states after accounting for lag-one correlation
  matrix<lower=0>[nyrs, A] N_ta;        // returns by age matrix
  matrix<lower=0, upper=1>[nRyrs, A] p; // age at maturity proportions
  vector<lower=0,upper=1>[4] pi;        // maturity schedule probs
  real<lower=0> D_sum;                  // inverse of D_scale which governs variability of age proportion vectors across cohorts
  vector<lower=0>[A] Dir_alpha;         // Dirichlet shape parameter for gamma distribution used to generate vector of age-at-maturity proportions
  matrix<lower=0, upper=1>[nyrs, A] q;  // age composition by year/age classr matrix
  real beta_fert;                       // Ricker b for fertilized state
  real beta_unfert;                     // Ricker b for unfertilized state

  // Maturity schedule: use a common maturation schedule to draw the brood year specific schedules
  pi[1] = prob[1];
  pi[2] = prob[2] * (1 - pi[1]);
  pi[3] = prob[3] * (1 - pi[1] - pi[2]);
  pi[4] = 1 - pi[1] - pi[2] - pi[3];
  D_sum = 1/D_scale^2;

  for (a in 1:A) {
    Dir_alpha[a] = D_sum * pi[a];
    for (y in 1:nRyrs) {
      p[y,a] = g[y,a]/sum(g[y,]);
    }
  }

  // simple calculations
  R = exp(lnR);
  beta_fert = 1/Smax_fert;
  beta_unfert = 1/Smax_unfert;
  sigma_R_corr = (sigma_R*sigma_R)/2;

  // Calculate the numbers at age matrix as brood year recruits at age (proportion that matured that year)
  for (t in 1:nyrs) {
    for(a in 1:A){
      N_ta[t,a] = R[t+A-a] * p[t+A-a,a];
    }
  }

  // Calculate returns, spawners and catch by return year
  for(t in 1:nyrs) {
    N[t] = sum(N_ta[t,1:A]);
    S[t] = N[t] * (1 - U[t]);
    lnS[t] = log(S[t]);
    C[t] = N[t] * U[t];
    lnC[t] = log(C[t]);
  }

  // Calculate age proportions by return year
  for (t in 1:nyrs) {
    for(a in 1:A){
      q[t,a] = N_ta[t,a]/N[t];
    }
  }

  // Ricker SR with AR1 process on log recruitment residuals for years with brood year spawners
  for (i in 1:nRyrs) {
    lnresid[i] = 0.0;
    lnRm_1[i] = 0.0;
    lnRm_2[i] = 0.0;
  }


// allow different alpha/beta for each fertilization state
// could try above approach but allowing only alpha or only beta to vary according to fertilzation
// fertilization maybe has more an effect on beta than alpha?
  for (y in (A+a_min):nRyrs) {
    real beta = f[y - a_max] ? beta_fert : beta_unfert;
    real lnalpha = f[y - a_max] ? lnalpha_fert : lnalpha_unfert;
    lnRm_1[y] = lnS[y-a_max] + lnalpha - beta * S[y-a_max];
    lnresid[y] = lnR[y] - lnRm_1[y];
  }


  lnRm_2[A+a_min] =  lnRm_1[A+a_min] + phi * lnresid_0;

  for (y in (A+a_min+1):nRyrs) {
    lnRm_2[y] =  lnRm_1[y] + phi * lnresid[y-1];
  }
}

model{
  // Priors
  lnalpha_fert ~ normal(1, 2);
  lnalpha_unfert ~ normal(1, 2);
  Smax_fert ~ lognormal(Smax_p_corr, Smax_p_sig_corr); // 1/beta
  Smax_unfert ~ lognormal(Smax_p_corr, Smax_p_sig_corr); // 1/beta
  sigma_R ~ normal(0,2);
  lnresid_0 ~ normal(0,10);
  mean_ln_R0 ~ normal(0,10);
  sigma_R0 ~ inv_gamma(2,1); 
  prob[1] ~ beta(1,1);
  prob[2] ~ beta(1,1);
  prob[3] ~ beta(1,1);
  D_scale ~ beta(1,1);

  // Likelihoods
  // Gamma variates for each year and age class which are used to determine age at maturity proportions
  for (y in 1:nRyrs) {
    for (a in 1:A) {
      //g[y,a] ~ gamma(Dir_alpha[a],1);
      target += gamma_lpdf(g[y,a]|Dir_alpha[a], 1);
    }
  }

  // First `a.max` years of recruits, for which there is no spawner link
  lnR[1:a_max] ~ normal(mean_ln_R0, sigma_R0);

  // State model
  for (y in (A + a_min):nRyrs) {
    // Adjust index to align "use" years with recruitment years
    int adj_y = y - (A - 1); 
    
    if (use[adj_y]) {
      lnR[y] ~ normal(lnRm_2[y], sigma_R_corr);
    }
  }
  // could create a loop where this is applied only on chosen years (e.g. 
  // remove the 1993 data point)
  // wrap above in if statement, if(a+amin /= 1993) {fit normal model to obs}
  // no "else" required...

  // Observation model
  for(t in 1:nyrs){
  //A_obs[t,1:A]) ~ multinomial(q[t,1:A]);
    target += multinomial_lpmf(A_obs[t,1:A]|to_vector(q[t,1:A]));
    U[t] ~ beta(1,1);
    //exclude chosen years from the observation model for harvest and spawners
    if(use[t]) {
      H_obs[t] ~ lognormal(lnC[t], sqrt(log((H_cv[t]^2) + 1)));
      S_obs[t] ~ lognormal(lnS[t], sqrt(log((S_cv[t]^2) + 1)));
    }
  }
}

generated quantities {

}
