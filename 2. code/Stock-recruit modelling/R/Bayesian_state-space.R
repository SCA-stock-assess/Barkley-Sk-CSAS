# Packages ----------------------------------------------------------------

pkgs <- c("here", "tidyverse", "readxl", "rstan", "bayesplot")
#install.packages(pkgs)

library(here)
library(tidyverse); theme_set(theme_bw())
library(rstan)
library(readxl)
library(bayesplot)


# Load and format data for Stan input -------------------------------------


# Load stock-recruit time series by return year
run_ts <- here("1. data", "return by age time series.xlsx") |> 
  read_xlsx() |> 
  mutate(
    run = catch + escapement,
    brood_year = return_year - ttl_age
  )


# Calculate brood year returns
brood_ts <- run_ts |> 
  summarize(
    .by = c(brood_year, stock),
    return = sum(catch, escapement)
  )


# Build stock-recruit table
sr <- run_ts |> 
  mutate(
    .by = c(return_year, stock),
    spawners = sum(escapement)
  ) |> 
  pivot_wider(
    id_cols = c(return_year, stock, spawners),
    names_from = ttl_age,
    names_prefix = "N.age.",
    values_from = run,
    values_fn = sum
  ) |> 
  rowwise() |> 
  mutate(
    run = sum(c_across(contains("age"))),
    # Convert age #s to proportions
    across(contains("age"), \(x) x/run),
  ) |> 
  ungroup() |> 
  left_join(
    brood_ts,
    by = join_by(
      stock, 
      return_year == brood_year
    )
  ) |> 
  rename(
    "year" = return_year,
    "S" = spawners,
    "N" = run,
    "R" = return
  )


# Declare a function that transforms data for 1 stock into correct 
# Stan input format. Need to do this because will be required for both
# SPR and GCL.
make_stan_data <- function(stock) {
  
  fn_data <- sr[sr$stock == stock, ]
  
  A_obs <- fn_data |> 
    select(contains("age")) |> 
    as.matrix() |>
    base::`*`(1000) |> 
    round() 
    # Using 1000/year as placeholder... would need to do considerable
    # data gathering from historic files to get the time series of 
    # number of samples going back to 1977. Does this actually 
    # matter enough to be worth doing so?
    
  S_obs <- fn_data$S
  H_obs <- fn_data$N-fn_data$S
  H_obs[H_obs<=0] <- 0.01 # Replace 0's with 0.01, otherwise you get ln(0) which breaks
  
  a_min <- min(run_ts$ttl_age) # youngest age at maturity
  a_max <- max(run_ts$ttl_age) # oldest age at maturity
  nyrs <- length(S_obs) # number of spawning years
  A <- a_max - a_min + 1 # total age classes
  nRyrs <- nyrs + A - 1 # number of recruitment years 
  
  stan.data <- list(
    "nyrs" = nyrs,
    "a_min" = a_min,
    "a_max" = a_max,
    "A" = A,
    "nRyrs" = nyrs + A - 1,
    "A_obs" = A_obs,
    "S_obs" = S_obs,
    "H_obs" = H_obs,
    "S_cv" = rep(0.05,length(S_obs)),
    "H_cv" = rep(0.05,length(H_obs)), 
    "Smax_p" = 0.75*max(fn_data$S), #what do we think Smax is? 
    "Smax_p_sig" = 1*max(fn_data$S) #can tinker with these values making the factor smaller mean's you've observed density dependence (i.e. the ricker "hump")
  )
  
  return(stan.data)
}

# Save Stan data for GCL and SPR (more useful as a list?)
stocks_stan_data <- unique(sr$stock) |> 
  purrr::set_names() |> 
  map(make_stan_data)

# Fit Stan models for both stocks -----------------------------------------

# Embrace Stan model in a function call to iterate over stocks
fit_stan_mod <- function(stan_data) {
  stan(
    file = here(
      "2. code",
      "Stock-recruit modelling",
      "Stan",
      "SS-SR_AR1.stan"
    ),
    model_name = "SS-SR_AR1",
    data = stan_data
  )
}


# Try stan model on GCL data
# `FALSE` disables the code from running
# Switch to `TRUE` to run
if(FALSE) {
  GCL_AR1 <- stan(
    file = here(
      "2. code",
      "Stock-recruit modelling",
      "Stan",
      "SS-SR_AR1.stan"
      #"SS-SR_AR1_semi_beta.stan #toggle which to fit 
     ),
    model_name = "SS-SR_AR1", #toggle which version of model you want 
    data = stocks_stan_data$GCL
  )
  
  # Save the fitted model object
  saveRDS(
    GCL_AR1,
    file = here(
      "3. outputs",
      "stock-recruit modelling",
      "GCL_AR1.rds"
    )
  )      
}


# Try stan model on SPR data
# `FALSE` disables the code from running
# Switch to `TRUE` to run
if(FALSE) {
  SPR_AR1 <- stan(
    file = here(
      "2. code",
      "Stock-recruit modelling",
      "Stan",
      "SS-SR_AR1.stan"
    ),
    model_name = "SS-SR_AR1",
    data = stocks_stan_data$SPR
  )
  
  # Save the fitted model object
  saveRDS(
    SPR_AR1,
    file = here(
      "3. outputs",
      "stock-recruit modelling",
      "SPR_AR1.rds"
    )
  )      
}


# Load fitted models from files (if the code above wasn't run)
if(!exists("GCL_AR1")) {
  GCL_AR1 <- readRDS(
    here(
      "3. outputs",
      "stock-recruit modelling",
      "GCL_AR1.rds"
    )
  )
}


# Load fitted models from files (if the code above wasn't run)
if(!exists("SPR_AR1")) {
  SPR_AR1 <- readRDS(
    here(
      "3. outputs",
      "stock-recruit modelling",
      "SPR_AR1.rds"
    )
  )
}


# Combine the models into a list for diagnostics
AR1_mods <- list(
  "GCL" = GCL_AR1,
 "SPR" = SPR_AR1
)


# Save summary data from the fitted models
AR1_models_summary <- AR1_mods |> 
  map(rstan::summary) |> 
  map(\(x) x$summary) |> 
  map(as.data.frame) |> 
  list_rbind(names_to = "stock") |> 
  nest(.by = stock, .key = "model_summary")


model_pars_AR1 <- AR1_mods |> 
  map(
    \(x) x |> 
      rstan::extract() |> 
      map(as.data.frame) |> 
      enframe()
  ) |> 
  list_rbind(names_to = "stock")


# some diagnostics

# check the chains directly
# leading pars
lead_pars <- c("beta", "lnalpha", "sigma_R", "lnresid_0", "phi")

AR1_mods |> 
  map(
    \(x) mcmc_combo(
      x, 
      pars = lead_pars,
      combo = c("dens_overlay", "trace"),
      gg_theme = legend_none()
    ) 
  )


# age pars
AR1_mods |> 
  map(
    \(x) mcmc_combo(
      x, 
      pars = c("D_scale", "D_sum"),
      combo = c("dens_overlay", "trace"),
      gg_theme = legend_none()
    )
  )

AR1_mods |> 
  map(
    \(x) mcmc_combo(
      x, 
      pars = paste0("Dir_alpha[", 1:4, "]"),
      combo = c("dens_overlay", "trace"),
      gg_theme = legend_none()
    )
  )


# how do correlations in leading parameters look?
AR1_mods |> 
  map(\(x) pairs(x, pars = lead_pars))




# Inference from fitted models --------------------------------------------


# Create dataframes of spawner abundances and predicted recruitment
make_srplot_data <- function(stock, stan_mod, stan_data) {
  
  model_pars <- stan_mod |> 
    rstan::extract()
  
  sr_years <- stan_data$nyrs # nRyrs doesn't work...
  max_samples <- dim(model_pars$lnalpha)
  
  spwn <- exp(model_pars$lnS)
  spwn.quant <- apply(spwn, 2, quantile, probs=c(0.05,0.5,0.95))[,1:(sr_years-4)]
  
  rec <-exp(model_pars$lnR)
  rec.quant <- apply(rec, 2, quantile, probs=c(0.05,0.5,0.95))[,8:dim(model_pars$R)[2]]
  
  brood_t <- as.data.frame(cbind(1:(sr_years-4),t(spwn.quant), t(rec.quant)))
  colnames(brood_t) <- c("BroodYear","S_lwr","S_med","S_upr","R_lwr","R_med","R_upr")
  
  brood_t <- as.data.frame(brood_t) |> 
    mutate(BroodYear = BroodYear + min(sr[sr$stock == stock,]$year) - 1)
  
  # SR relationship
  spw <- seq(0,max(brood_t[,4]),length.out=100)
  SR_pred <- matrix(NA,100,max_samples)
  
  for(i in 1:max_samples){
    r <- sample(seq(1,max_samples),1,replace=T)
    a <- model_pars$lnalpha[r]
    b <- model_pars$beta[r]
    SR_pred[,i] <- (exp(a)*spw*exp(-b*spw))
  }
  
  SR_pred <- cbind(spw,t(apply(SR_pred,c(1),quantile,probs=c(0.05,0.5,0.95),na.rm=T)))
  colnames(SR_pred) <- c("Spawn", "Rec_lwr","Rec_med","Rec_upr")
  SR_pred <- as.data.frame(SR_pred)
  
  output <- list(
    "SR_pred" = SR_pred,
    "brood_t" = brood_t
  )
  
  return(output)
  
}


# Make data for plotting
plot_data <- unique(sr$stock) |> 
  purrr::set_names() |> 
  map(
    \(x) make_srplot_data(
      stock = x,
      stan_mod = AR1_mods[[x]],
      stan_data = stocks_stan_data[[x]]
    )
  )


# Plot
plot_data |> 
  imap(
    \(x, idx) ggplot() +
      geom_errorbar(
        data = x$brood_t, 
        aes(x= S_med, y = R_med, ymin = R_lwr, ymax = R_upr),
        colour="grey", 
        width = 0,
        size = 0.3
      ) +
      geom_abline(intercept = 0, slope = 1,col="dark grey") +
      geom_ribbon(
        data = x$SR_pred, 
        aes(x = Spawn, ymin = Rec_lwr, ymax = Rec_upr),
        fill = "grey80", 
        alpha = 0.5, 
        linetype = 2, 
        colour = "gray46"
      ) +
      geom_line(
        data = x$SR_pred, 
        aes(x = Spawn, y = Rec_med), 
        color = "black", 
        size = 1
      ) +
      geom_errorbarh(
        data = x$brood_t, 
        aes(y = R_med, xmin = S_lwr, xmax = S_upr),
        height = 0, 
        colour = "grey", 
        size = 0.3
      ) +
      geom_point(
        data = x$brood_t, 
        aes(x = S_med, y = R_med, color=BroodYear), 
        size = 3
      ) +
      scale_x_continuous(
        limits = c(0, max(x$brood_t[,4])),
        labels = scales::label_number(),
        expand = expansion(c(0, 0.05))
      ) +
      scale_y_continuous(
        limits = c(0, max(x$brood_t[,6])),
        labels = scales::label_number(),
        expand = expansion(c(0, 0.05)),
        oob = scales::oob_keep
      ) +
      scale_colour_viridis_c()+
      labs(
        x = "Spawners",
        y = "Recruits",
        title = paste("Stock:", idx)
      ) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size=9),
        legend.text = element_text(size=8)
      )
  )
