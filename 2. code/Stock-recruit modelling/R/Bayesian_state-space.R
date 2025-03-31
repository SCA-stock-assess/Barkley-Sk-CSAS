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
sr <- here(
  "3. outputs", 
  "Stock-recruit data",
  'Barkley_Sockeye_stock-recruit_infilled.xlsx'
  ) |> 
  read_xlsx(sheet = "S-R data") |> 
  # Remove earlier years in the time series
  filter(!year < 1972)


# Create dataframe of historical average age compositions
stock_age_averages <- sr |> 
  summarize(
    .by = stock,
    across(contains("N.age"), \(x) sum(x, na.rm = TRUE))
  ) |> 
  rowwise() |> 
  mutate(total = sum(c_across(contains("N.age")))) |> 
  mutate(
    across(contains("N.age"), \(x) x/total),
    .keep = "unused"
  )


# Infill missing ages with historical averages 
sr_age_infill <- sr |> 
  filter(if_all(contains("N.age"), is.na)) |> 
  select(-contains("N.age")) |> 
  left_join(stock_age_averages) |> 
  mutate(age.samples = 1) |> 
  bind_rows(sr) |> 
  filter(!if_all(contains("N.age"), is.na)) |> 
  rowwise() |> 
  mutate(
    # Calculate adult spawners for years where ages were infilled
    adult_S = if_else(is.na(adult_S), S*sum(c_across(matches("N.*(4|5|6)"))), adult_S),
    # Convert all age columns to proportions
    across(
      contains("N.age"), 
      \(x) x/sum(c_across(contains("N.age")))
    )
  ) |> 
  ungroup() |> 
  mutate(
    .by = stock,
    max.age.samples = max(age.samples),
    # Assign effective sample size (ESS) values for age data by stock and year, 
    # ESS = 100 implies high confidence
    age.ess = case_when(
      stock %in% c("GCL", "SPL") & year > 2009 ~ 100,
      age.samples < 15 ~ 10,
      stock == "HUC" & year > 2020 ~ 100,
      stock == "HUC" & age.samples/max.age.samples < 0.25 ~ 50,
      stock %in% c("GCL", "SPL") & age.samples < 400 ~ 60,
      age.samples/max.age.samples < 0.25 ~ 75,
      .default = 75
    )
  )


# Ensure all age columns sum to 1
sr_age_infill |> 
  rowwise() |> 
  summarize(age_comp = sum(c_across(contains("N.age")))) |> 
  count(age_comp)

  
# Declare a function that transforms data for 1 stock into correct 
# Stan input format. Need to do this because will be required for both
# SPR and GCL.
make_stan_data <- function(stock_name) {
  
  fn_data <- sr_age_infill |> 
    filter(stock == stock_name) |> 
    arrange(year) # ensure years are correctly ordered
  
  A_obs <- fn_data |> 
    # Multiply age compositions by assigned effective sample sizes
    mutate(across(contains("N.age"), ~.x*age.ess)) |> 
    select(contains("N.age")) |> 
    as.matrix() |> 
    round() 
    
  brood_years <- fn_data$year
  S_obs <- fn_data$adult_S
  H_obs <- fn_data$H
  H_obs[H_obs<=0] <- 0.01 # Replace 0's with 0.01, otherwise you get ln(0) which breaks
  H_cv <- fn_data$H_cv
  S_cv <- fn_data$S_cv
  
  # Vector describing whether the lake was fertilized each year
  f <- fn_data |> 
    mutate(
      # Simplify the information for Great Central and Sproat
      f = case_when(
        stock == "GCL" ~ 1, # changes nothing: all GCL years already = 1
        stock == "SPR" ~ 0, # overwrites 1985, the one year Sproat was fertilized
        stock == "HUC" ~ fertilized
      )
    ) |> 
    pull(f)
  
  # Extract a vector of the observed total age classes
  a_years <- fn_data |> 
    select(matches("N\\.age\\.\\d")) |> 
    colnames() |> 
    str_extract("\\d") |> 
    as.integer()
  
  a_min <- min(a_years) # youngest age at maturity
  a_max <- max(a_years) # oldest age at maturity
  nyrs <- length(brood_years) # number of spawning years/length of data time series
  A <- a_max - a_min + 1 # total age classes
  nRyrs <- nyrs + A - 1 # recruitment years
  # note that additional years are added at the beginning of the time series
  # (i.e. nRyrs does not extend into the future)
  
  # Declare which years to include in state-space SR model
  use <- fn_data |> 
    # Exclude 1993 for Hucuktlis model
    mutate(use = if_else(stock == "HUC" & year == 1993, 0, 1)) |> 
    pull(use)
  
  stan.data <- list(
    "nyrs" = nyrs,
    "brood_years" = brood_years,
    "a_min" = a_min,
    "a_max" = a_max,
    "A" = A,
    "nRyrs" = nRyrs,
    "A_obs" = A_obs,
    "S_obs" = S_obs,
    "H_obs" = H_obs,
    "S_cv" = S_cv,
    "H_cv" = H_cv, 
    "use" = use,
    "f" = f,
    "Smax_p" = 0.75*max(fn_data$adult_S), #what do we think Smax is? 
    "Smax_p_sig" = 1*max(fn_data$adult_S) #can tinker with these values making the factor smaller mean's you've observed density dependence (i.e. the Ricker "hump")
  )
  
  return(stan.data)
}


# Save Stan data for GCL and SPR (more useful as a list?)
stocks_stan_data <- unique(sr_age_infill$stock) |> 
  purrr::set_names() |> 
  map(make_stan_data)


# Fit Stan models for Somass stocks and validate fits -------------------


# Embrace Stan model in a function call to iterate over stocks
fit_stan_mod <- function(stan_data, stock, model_filename) {
  stan(
    file = here(
      "2. code",
      "Stock-recruit modelling",
      "Stan",
      # toggle which model to fit
      #"SS-SR_AR1.stan"
      model_filename
    ),
    iter = 5000, # minor issues achieving ESS for Sproat; increasing iterations
    control = list(max_treedepth = 15),
    model_name = stock,
    data = stan_data
  )
}


#  Stan model on GCL data
# `FALSE` disables the code from running
# Switch to `TRUE` to run
if(FALSE) {
  GCL_AR1 <- fit_stan_mod(
    stocks_stan_data$GCL, 
    stock = "GCL",
    model_filename = "SS-SR_AR1_semi_beta_Somass.stan"
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
  SPR_AR1 <- fit_stan_mod(
    stocks_stan_data$SPR, 
    stock = "SPR",
    model_filename = "SS-SR_AR1_semi_beta_Somass.stan"
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
Somass_mods <- list(
  "GCL" = GCL_AR1,
  "SPR" = SPR_AR1
)


# some diagnostics

# n_eff versus Rhat
Somass_mods |> 
  map(\(x) summary(x)$summary) |>  
  map(as.data.frame) |> 
  list_rbind(names_to = "stock") |> 
  ggplot(aes(x = n_eff, y = Rhat))+
  facet_wrap(~stock) +
  geom_point() +
  stat_density_2d(
    alpha = 0.6,
    geom = "polygon", 
    contour = TRUE,
    aes(fill = after_stat(level)), 
    bins = 6
  ) +
  geom_hline(yintercept = 1.01, lty = 2)+
  geom_vline(xintercept = 400, lty = 2) +
  scale_fill_distiller(palette = "RdPu", direction = 1) +
  guides(fill = "none")


# check the chains directly
# leading pars
lead_pars <- c("beta", "lnalpha", "sigma_R", "lnresid_0", "phi")

(lead_pars_p <- Somass_mods |> 
  map(
    \(x) mcmc_combo(
      x, 
      pars = lead_pars,
      combo = c("dens_overlay", "trace"),
      gg_theme = legend_none()
    ) 
  )
)


# Save the plots
lead_pars_p |> 
  iwalk(
    \(x, idx) ggsave(
      plot = x,
      filename = here(
        "3. outputs",
        "Stock-recruit modelling",
        "Bayesian state-space diagnostics",
        paste0(
          "Leading_parameters_",
          idx,
          ".png"
        )
      )
    )
  )
  

# age pars
(age_pars_p <- Somass_mods |> 
  map(
    \(x) mcmc_combo(
      x, 
      pars = c("D_scale", "D_sum"),
      combo = c("dens_overlay", "trace"),
      gg_theme = legend_none()
    )
  )
)


# Save the plots
age_pars_p |> 
  iwalk(
    \(x, idx) ggsave(
      plot = x,
      filename = here(
        "3. outputs",
        "Stock-recruit modelling",
        "Bayesian state-space diagnostics",
        paste0(
          "Age_parameters_",
          idx,
          ".png"
        )
      )
    )
  )

  
# Dirichlet parameters  
(Dir_pars_p <- Somass_mods |> 
  map(
    \(x) mcmc_combo(
      x, 
      pars = paste0("Dir_alpha[", 1:4, "]"),
      combo = c("dens_overlay", "trace"),
      gg_theme = legend_none()
    )
  )
)


# Save the plots
Dir_pars_p |> 
  iwalk(
    \(x, idx) ggsave(
      plot = x,
      filename = here(
        "3. outputs",
        "Stock-recruit modelling",
        "Bayesian state-space diagnostics",
        paste0(
          "Dirichlet_parameters_",
          idx,
          ".png"
        )
      )
    )
  )
  

# how do correlations in leading parameters look?
(corr_p <- Somass_mods |> 
  map(\(x) pairs(x, pars = c("beta", "lnalpha", "phi")))
)


# Save the correlation plots
corr_p |> 
  iwalk(
    \(x, idx) ggsave(
      plot = x,
      filename = here(
        "3. outputs",
        "Stock-recruit modelling",
        "Bayesian state-space diagnostics",
        paste0(
          "Correlation_leading-parameters_",
          idx,
          ".png"
        )
      ),
      width = 9,
      height = 9,
      units = "in",
      dpi = "print"
    )
  )



# Plot Somass stock-recruit curves with latent states ----------------------


# Create dataframes of spawner abundances and predicted recruitment
make_srplot_data <- function(stock, stan_mod, stan_data) {
  
  model_pars <- rstan::extract(stan_mod)
  
  a_min <- stan_data$a_min
  A <- stan_data$A
  nyrs <- stan_data$nyrs 
  nRyrs <- stan_data$nRyrs
  max_samples <- dim(model_pars$lnalpha)
  
  spwn <- model_pars$S
  spwn.quant <- apply(spwn, 2, quantile, probs=c(0.05,0.5,0.95))[,1:(nyrs-a_min)]
  
  rec <- model_pars$R
  rec.quant <- apply(rec, 2, quantile, probs=c(0.05,0.5,0.95))[,(A+a_min):nRyrs]
  
  brood_years <- stan_data$brood_years[1:(nyrs-a_min)]
  brood_t <- as.data.frame(cbind(brood_years, t(spwn.quant), t(rec.quant)))
  colnames(brood_t) <- c("BroodYear","S_lwr","S_med","S_upr","R_lwr","R_med","R_upr")
  
  
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
Somass_sr_plot_data <- names(Somass_mods) |> 
  purrr::set_names() |> 
  map(
    \(x) make_srplot_data(
      stock = x,
      stan_mod = Somass_mods[[x]],
      stan_data = stocks_stan_data[[x]]
    )
  )


# Function to make SR plot
make_sr_plots <- function(sr_plot_data, stock) {
  ggplot() +
    geom_errorbar(
      data = sr_plot_data$brood_t, 
      aes(x= S_med, y = R_med, ymin = R_lwr, ymax = R_upr),
      colour="grey", 
      width = 0,
      size = 0.3
    ) +
    geom_abline(intercept = 0, slope = 1,col="dark grey") +
    geom_ribbon(
      data = sr_plot_data$SR_pred, 
      aes(x = Spawn, ymin = Rec_lwr, ymax = Rec_upr),
      fill = "grey80", 
      alpha = 0.5, 
      linetype = 2, 
      colour = "grey46"
    ) +
    geom_line(
      data = sr_plot_data$SR_pred, 
      aes(x = Spawn, y = Rec_med), 
      color = "black", 
      size = 1
    ) +
    geom_errorbarh(
      data = sr_plot_data$brood_t, 
      aes(y = R_med, xmin = S_lwr, xmax = S_upr),
      height = 0, 
      colour = "grey", 
      size = 0.3
    ) +
    geom_point(
      data = sr_plot_data$brood_t, 
      aes(x = S_med, y = R_med, color=BroodYear), 
      size = 3
    ) +
    scale_x_continuous(
      limits = c(0, max(sr_plot_data$brood_t[,4])),
      labels = scales::label_number(),
      breaks = scales::pretty_breaks(n = 3),
      expand = expansion(c(0, 0.05))
    ) +
    scale_y_continuous(
      limits = c(0, max(sr_plot_data$brood_t[,6])),
      labels = scales::label_number(),
      expand = expansion(c(0, 0.05)),
      oob = scales::oob_keep
    ) +
    scale_colour_viridis_c() +
    #coord_equal(ratio = 1) + # interesting way to show the data, but ultimately just looks crowded 
    labs(
      x = "Spawners",
      y = "Recruits",
      title = paste("Stock:", stock)
    ) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key.size = unit(0.4, "cm"),
      legend.title = element_text(size=9),
      legend.text = element_text(size=8),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}


# Plot
(sr_plots <- Somass_sr_plot_data |> 
    imap(make_sr_plots)
)


# Save plots
sr_plots |> 
  iwalk(
    \(x, idx) ggsave(
      plot = x,
      filename = here(
        "3. outputs",
        "Plots",
        paste0(
          "State-space_stock-recruit_predictions_",
          idx,
          ".png"
        )
      ),
      dpi = "print"
    )
  )



# Test method of data exclusion on Hucuktlis data -------------------------


# Lists of stan data for Hucuktlis with and without the 1993 observation
HUC_stan_full <- make_stan_data("HUC")
HUC_stan_full$use <- rep(1, length(HUC_stan_full$use))
HUC_stan_trim <- make_stan_data("HUC") # default 'use' in make_stan_data excludes 1993

HUC_stan_data <- list(
  "HUC_full" = HUC_stan_full,
  "HUC_trim" = HUC_stan_trim
)  


# Fit two models: one with and one without 1993 data in the state space model  
# Round 1 because round 2 (following) will include the fertilization component
if(FALSE) {
  
  HUC_round1_mods <- HUC_stan_data |> 
    set_names(\(x) paste0(x, "_nofert")) |> # Rename so model names have "_nofert" suffix
    imap(
      \(x, idx) fit_stan_mod(
        stan_data = x,
        stock = idx, 
        # Using the Somass model for simplicity in this section
        # The difference is that only 1 value each for lnalpha and beta are estimated
        # (i.e. irrespective of fertilization state, which comes in the next section)
        model_filename = "SS-SR_AR1_semi_beta_Somass.stan"
      )
    )
  
  # Save the fitted models
  HUC_round1_mods |> 
    iwalk(
      \(x, idx) saveRDS(
        x,
        file = here(
          "3. outputs",
          "stock-recruit modelling",
          paste0(idx, "_AR1.rds")
        )
      )     
    )
}


# Read fitted model objects from RDS if code above is not run
if(!exists("HUC_round1_mods")) {
  HUC_round1_mods <- names(HUC_stan_data) |> 
    paste0("_nofert") |> 
    purrr::set_names() |> 
    map(
      \(x) readRDS(
        here(
          "3. outputs",
          "stock-recruit modelling",
          paste0(x, "_AR1.rds")
        )
      )
    )
}


# Skip straight to SR plots to see if/how much the line shifts
HUC_sr_plots_data <- pmap(
  list(
    purrr::set_names(names(HUC_round1_mods)),
    HUC_round1_mods,
    HUC_stan_data
  ),
  make_srplot_data
)


HUC_sr_plots_data |> 
  imap(make_sr_plots)
# Doesn't look like the relationship is affected much by exclusion of 1993...



# Fit Hucuktlis model with fertilization-specific a and b vals --------


# This time use the Hucuktlis-specific Stan model
if(FALSE) {

  HUC_round2_mods <- HUC_stan_data |> 
    set_names(\(x) paste0(x, "_fert")) |> # Rename so model names have "_fert" suffix
    imap(
      \(x, idx) fit_stan_mod(
        stan_data = x,
        stock = idx, 
        model_filename = "SS-SR_AR1_semi_beta_Hucuktlis.stan"
      )
    )

  # Save the fitted models
  HUC_round2_mods |> 
    iwalk(
      \(x, idx) saveRDS(
        x,
        file = here(
          "3. outputs",
          "stock-recruit modelling",
          paste0(idx, "_AR1.rds")
        )
      )     
    )
}


# Read fitted model objects from RDS if code above is not run
if(!exists("HUC_round2_mods")) {
  HUC_round2_mods <- names(HUC_stan_data) |> 
    paste0("_fert") |> 
    purrr::set_names() |> 
    map(
      \(x) readRDS(
        here(
          "3. outputs",
          "stock-recruit modelling",
          paste0(x, "_AR1.rds")
        )
      )
    )
}


# Make a list of all 4 fitted Hucuktlis models
HUC_mods <- c(HUC_round1_mods, HUC_round2_mods)


# Run the diagnostic plots on all 4 Hucuktlis models
# n_eff versus Rhat
HUC_mods |> 
  map(\(x) summary(x)$summary) |>  
  map(as.data.frame) |> 
  list_rbind(names_to = "model") |> 
  ggplot(aes(x = n_eff, y = Rhat))+
  facet_wrap(~model) +
  geom_point() +
  stat_density_2d(
    alpha = 0.6,
    geom = "polygon", 
    contour = TRUE,
    aes(fill = after_stat(level)), 
    bins = 6
  ) +
  geom_hline(yintercept = 1.01, lty = 2)+
  geom_vline(xintercept = 400, lty = 2) +
  scale_fill_distiller(palette = "RdPu", direction = 1) +
  guides(fill = "none")


# check the chains directly
(lead_pars_HUC_p <- HUC_mods |> 
    map(
      function(x) {
        
        lead_pars_mod <- str_subset(
          names(x),
          "lnalpha|beta|phi|sigma_R$|lnresid_0|phi"
        )
        
        mcmc_combo(
          x, 
          pars = lead_pars_mod,
          combo = c("dens_overlay", "trace"),
          gg_theme = legend_none()
        )
      }
    )
)


# Save the plots
lead_pars_HUC_p |> 
  iwalk(
    \(x, idx) ggsave(
      plot = x,
      filename = here(
        "3. outputs",
        "Stock-recruit modelling",
        "Bayesian state-space diagnostics",
        paste0(
          "Leading_parameters_",
          idx,
          ".png"
        )
      )
    )
  )


# age pars
(age_pars_HUC_p <- HUC_mods |> 
    map(
      \(x) mcmc_combo(
        x, 
        pars = c("D_scale", "D_sum"),
        combo = c("dens_overlay", "trace"),
        gg_theme = legend_none()
      )
    )
)


# Save the plots
age_pars_HUC_p |> 
  iwalk(
    \(x, idx) ggsave(
      plot = x,
      filename = here(
        "3. outputs",
        "Stock-recruit modelling",
        "Bayesian state-space diagnostics",
        paste0(
          "Age_parameters_",
          idx,
          ".png"
        )
      )
    )
  )


# Dirichlet parameters  
(Dir_pars_HUC_p <- HUC_mods |> 
    map(
      \(x) mcmc_combo(
        x, 
        pars = paste0("Dir_alpha[", 1:4, "]"),
        combo = c("dens_overlay", "trace"),
        gg_theme = legend_none()
      )
    )
)


# Save the plots
Dir_pars_HUC_p |> 
  iwalk(
    \(x, idx) ggsave(
      plot = x,
      filename = here(
        "3. outputs",
        "Stock-recruit modelling",
        "Bayesian state-space diagnostics",
        paste0(
          "Dirichlet_parameters_",
          idx,
          ".png"
        )
      )
    )
  )


# how do correlations in leading parameters look?
(corr_HUC_p <- HUC_mods |> 
    map(
      function(x) {
        corr_pars <- str_subset(
          names(x),
          "lnalpha|beta|phi"
        )
        
        pairs(
          x, 
          pars = corr_pars
        )
      }
    )
)


# Save the correlation plots
corr_HUC_p |> 
  iwalk(
    \(x, idx) ggsave(
      plot = x,
      filename = here(
        "3. outputs",
        "Stock-recruit modelling",
        "Bayesian state-space diagnostics",
        paste0(
          "Correlation_leading-parameters_",
          idx,
          ".png"
        )
      ),
      width = 9,
      height = 9,
      units = "in",
      dpi = "print"
    )
  )


# Time series plots of model residuals ------------------------------------


# Residuals from fitted models
resids <- c(Somass_mods, HUC_mods) |> 
  map(\(x) rstan::extract(x, "lnresid")) |> 
  map(\(x) map(x, \(y) as_tibble(y, rownames = "draw", .name_repair = NULL))) |> 
  list_flatten(name_spec = "{outer}") |> 
  map(\(x) pivot_longer(x, !draw, names_to = "Ryr")) |> 
  list_rbind(names_to = "model") |> 
  mutate(
    .by = model,
    Ryr = as.integer(str_extract(Ryr, "\\d+")),
    max_Ryr = max(Ryr),
    max_yr = max(sr$year),
    brood_year = max_yr - max_Ryr + Ryr
  ) |> 
  filter(brood_year <= max(sr$year - 6)) |> 
  summarize(
    .by = c(model, brood_year),
    probs = list(quantile(value, probs = c(0.1, 0.25, 0.5, 0.75, 0.9)))
  ) |> 
  unnest_wider(probs) |> 
  rename(
    "lwr" = 3,
    "midlwr" = 4,
    "mid" = 5,
    "midupr" = 6,
    "upr" = 7
  )


# Plot
ggplot(resids, aes(x=brood_year, y = mid)) +
  facet_wrap(
    ~model, 
    ncol = 1, 
    strip.position = "right", 
    scales = "free_y"
  ) +
  geom_abline(
    intercept = 0, 
    slope = 0, 
    col = "dark grey", 
    lty = 2
  ) +
  geom_ribbon(
    aes(ymin = lwr, ymax = upr),  
    fill = "darkgrey", 
    alpha = 0.5
  ) +
  geom_ribbon(
    aes(ymin = midlwr, ymax = midupr),  
    fill = "black", 
    alpha=0.2
  ) + 
  geom_line(lwd = 1.1) +
  labs(
    x = "Return year",
    y = "Recruitment residuals", 
  ) +
  theme(
    legend.position = "none",
    panel.grid = element_blank()
  ) 


# Combined SR plots for all CUs -------------------------------------------


# Posterior draws of alpha and beta for all models
ab_posterior <- c(HUC_mods, Somass_mods) |> 
  map(rstan::extract) |> 
  map(\(x) keep_at(x, \(y) str_detect(y, "lnalpha|beta"))) |> 
  map(\(x) map(x, \(y) as_tibble(y, rownames = "draw"))) |> 
  map(\(x) list_rbind(x, names_to = "parameter")) |> 
  list_rbind(names_to = "model") |> 
  mutate(
    stock = str_extract(model, "GCL|SPR|HUC"),
    fert = case_when(
      str_detect(parameter, "_fert") ~ 1, 
      stock == "GCL" ~ 1,
      .default = 0
    ),
    data_scope = if_else(str_detect(model, "trim"), "trim", "full"),
    parameter = str_remove_all(parameter, "_.*")
  ) |> 
  nest(ab_draws = c(parameter, value, draw)) |> 
  rowwise() |> 
  mutate(ab_draws = list(pivot_wider(ab_draws, names_from = parameter)))


# Posterior draws of spawners and recruits
sr_draws <- c(HUC_mods, Somass_mods) |> 
  map(\(x) rstan::extract(x, c("S", "R"))) |> 
  map(\(x) map(x, \(y) as_tibble(y, rownames = "draw", .name_repair = NULL))) |> 
  map(\(x) map(x, \(y) pivot_longer(y, !draw, names_to = "yr"))) |> 
  map(\(x) list_rbind(x, names_to = "parameter")) |> 
  list_rbind(names_to = "model") |> 
  mutate(
    .by = c(model, parameter),
    yr = as.integer(str_extract(yr, "\\d+")),
    max_yr = max(yr),
    max_sr_yr = max(sr$year),
    brood_year = max_sr_yr - max_yr + yr,
    stock = str_extract(model, "GCL|SPR|HUC"),
    data_scope = if_else(str_detect(model, "trim"), "trim", "full")
  ) |> 
  # Constrain range of brood years to only those that are observed in the time series
  filter(
    .by = stock,
    between(brood_year, min(brood_year) + 3, max(brood_year) - 6)
  ) |> 
  select(-contains("yr")) |> 
  pivot_wider(names_from = parameter, values_from = value) |> 
  # Add fertilization metadata from input dataframe
  left_join(
    select(.data = sr, brood_year = year, stock, fert = fertilized),
    by = c("brood_year", "stock")
  ) |> 
  summarize(
    .by = c(model, stock, fert, brood_year, data_scope),
    across(
      c(S, R), 
      .fns = list(
        `5` = ~quantile(.x, 0.05),
        `50` = ~quantile(.x, 0.5),
        `95` = ~quantile(.x, 0.95)
      ),
      .names = "{.col}_{.fn}"
    )
  ) |> 
  nest(brood_t = c(brood_year, fert, matches("_\\d+")))
  

# Merge data for SR plotting
all_mods_data <- left_join(
  ab_posterior,
  sr_draws,
  by = c("stock", "model", "data_scope")
) |> 
  rowwise() |> 
  mutate(
    long_name = factor(
      model,
      levels = c(
        "GCL", 
        "SPR", 
        "HUC_full_nofert", 
        "HUC_trim_nofert", 
        "HUC_full_fert", 
        "HUC_trim_fert"
      ),
      labels = c(
        "Great Central Lake",
        "Sproat Lake",
        "Hucuktlis Lake (all data)",
        "Hucuktlis Lake (excl. 1993)",
        "Hucuktlis Lake (w/fertilization)",
        "Hucuktlis Lake (w/fertilization; excl. 1993)"
      ) 
    ),
    # Calculate predicted spawner values from posterior lnalpha and beta draws
    pred_frame = list(
      data.frame("S" = seq(0, max(brood_t$S_95), length.out = 100)) |> 
        expand_grid(ab_draws) |> 
        mutate(
          R = exp(lnalpha)*S*exp(-beta*S),
          long_name = long_name,
          fert = factor(fert)
        ) |> 
        summarize(
          .by = c(S, long_name, fert),
          quantiles = list(quantile(R, c(0.05, 0.5, 0.95)))
        ) |> 
        unnest_wider(quantiles)
    ),
    brood_t = list(mutate(brood_t, long_name = long_name))
  ) |> 
  ungroup() 


# Data frame with latent states of spawners and recruits
all_mods_brood_t <- list_rbind(all_mods_data$brood_t) |> 
  mutate(across(matches("_\\d+"), \(x) x/1000))


# Data frame with predicted recruits per spawner by model
all_mods_pred_frame <- list_rbind(all_mods_data$pred_frame) |> 
  mutate(
    across(c(S, matches("%")), \(x) x/1000),
    fert = if_else(str_detect(long_name, "fertilization"), fert, NA)
  )



# Make tidy SR plots
(all_mods_sr_plots <- ggplot(
  data = all_mods_pred_frame,
  aes(x = S, y = `50%`)
) +
    facet_wrap(
      ~long_name,
      scales = "free"
    ) +
    geom_abline(intercept = 0, slope = 1,col="dark grey") +
    geom_line(aes(colour = fert)) +
    geom_ribbon(
      aes(
        ymin = `5%`, 
        ymax = `95%`,
        fill = fert,
        colour = fert
      ),
      alpha = 0.3,
      lty = 2
    ) +
    geom_errorbar(
      data = all_mods_brood_t, 
      aes(x= S_50, y = R_50, ymin = R_5, ymax = R_95),
      colour="grey", 
      width = 0,
      size = 0.3
    ) +
    geom_errorbarh(
      data = all_mods_brood_t, 
      aes(y = R_50, x = S_50, xmin = S_5, xmax = S_95),
      height = 0, 
      colour = "grey", 
      size = 0.3
    ) +
    geom_point(
      data = all_mods_brood_t, 
      aes(x = S_50, y = R_50, colour = factor(fert)), 
      size = 2,
      alpha = 0.6
    ) +
    scale_colour_manual(
      name = "Fertilization\nstate",
      values = c("red", "blue"),
      aesthetics = c("colour", "fill")
    ) +
    scale_x_continuous(
      limits = c(0, NA),
      labels = scales::label_number(),
      breaks = scales::pretty_breaks(n = 3),
      expand = expansion(c(0, 0))
    ) +
    scale_y_continuous(
      limits = c(0, NA),
      labels = scales::label_number(),
      expand = expansion(c(0, 0.05)),
      oob = scales::oob_keep
    ) +
    labs(
      x = "Spawners (1000s)",
      y = "Recruits (1000s)"
    ) +
    theme(
      strip.background = element_blank()
    )
)


# Save the plotted model fits
ggsave(
  plot = all_mods_sr_plots,
  here(
    "3. outputs", 
    "Plots",
    "All_Barkley_Sk_AR1_fits.png"
  ),
  height = 7, 
  width = 10,
  units = "in",
  dpi = "print"
)

