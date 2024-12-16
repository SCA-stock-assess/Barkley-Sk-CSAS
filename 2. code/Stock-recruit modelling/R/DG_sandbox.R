#Model fitting sandbox - just a place to play with parameterizations
  #can work through stocks and models 1 by 1 to fine tune paramaterizations. 
GCL_AR1 <- stan(
  file = here(
    "2. code",
    "Stock-recruit modelling",
    "Stan",
    #"SS-SR_AR1.stan"
    "SS-SR_AR1_semi_beta.stan" #toggle which to fit 
  ),
  data = stocks_stan_data$GCL
)

AR1.model.summary <- as.data.frame(rstan::summary(GCL_AR1)$summary)
model.pars.AR1 <- rstan::extract(GCL_AR1)

mcmc_combo(GCL_AR1, pars = c("beta", "Smax", "lnalpha", "sigma_R", "phi"),
           combo = c("dens_overlay", "trace"),
           gg_theme = legend_none())

mcmc_combo(GCL_AR1, pars = c("D_scale", "D_sum"),
           combo = c("dens_overlay", "trace"),
           gg_theme = legend_none())
