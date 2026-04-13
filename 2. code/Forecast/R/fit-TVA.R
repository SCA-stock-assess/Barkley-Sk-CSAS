# fit time-varying productivity models to same data used in rest of the analysis
library(here)
source(here("2. code/0. functions/common_functions.R"))

# run part of the SR fitting just to get the data 
source2(here("2. code/Stock-recruit modelling/R/Bayesian_state-space_alt-beta-prior.R"), 1, 250)
  #^ breaks. looks like I need to run spawner-smolt_analyesi.R first but that also breaks 
    #at L919

if(FALSE){
#fit models (1 by 1 for simplicity's sake)
SPR_fit <- stan(file = here("2. code/Forecast/Stan/SS-SR_TVA_semi_beta_Somass.stan"), 
                data = stocks_stan_data$SPR, 
                seed = 1)

GCL_fit <- stan(file = here("2. code/Forecast/Stan/SS-SR_TVA_semi_beta_Somass.stan"), 
                data = stocks_stan_data$GCL, 
                seed = 1)

HUC_fit <-  stan(file = here("2. code/Forecast/Stan/SS-SR_TVA_semi_beta_Somass.stan"), 
                 data = stocks_stan_data$HUC, 
                 seed = 1)
}

#save model within forecast subfolder

#brief model validation, chains, ESS, R-hats 

