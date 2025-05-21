# Packages ----------------------------------------------------------------


pkgs <- c("here", "tidyverse", "readxl", "rstan", "ggrepel", "MASS")
#install.packages(pkgs)


library(here)
library(tidyverse); theme_set(theme_bw())
library(readxl)
library(ggrepel)


# Load in data and fitted models ------------------------------------------


# Spawner-recruit data time series
sr_data <- here(
  "3. outputs",
  "Stock-recruit data",
  "Barkley_Sockeye_stock-recruit_infilled.xlsx"
) |> 
  read_xlsx(sheet = "S-R data") |> 
  # Calculate harvest rates
  mutate(hr = H/N)


# Fitted state-space spawner-recruit models
AR1_mods <- list.files(
  here(
    "3. outputs",
    "Stock-recruit modelling"
  ),
  pattern = "AR1.rds",
  full.names = TRUE
) |> 
  set_names(nm = ~str_extract(.x, "(SPR|GCL|HUC).*_AR1")) |> 
  map(readRDS) 

AR1_frame <- AR1_mods |> 
  enframe(name = "spec", value = "model") 
  
  
# Historic management forecast error
fcst_err <- here(
  "1. data",
  "Somass_Sockeye_forecast_error.xlsx"
) |> 
  read_xlsx() |> 
  rename_with(tolower) |> 
  filter(forecast == "MGT") 


# Show distribution of forecast errors
ggplot(fcst_err, aes(x = err_pct)) +
  geom_density(fill = "grey") +
  scale_y_continuous(expand = c(0, 0.05))


# Bootstrap mean & SD of forecast errors
fcst_err_boot <- fcst_err |> 
  select(err_pct) |> 
  expand_grid(sample = seq(1:1000)) |> 
  nest(.by = sample) |> 
  rowwise() |> 
  mutate(boot = list(sample(data, length(data), replace = TRUE))) |> 
  select(boot) |> 
  unnest(boot) |> 
  summarize(
    mean = mean(err_pct),
    sd = sd(err_pct)
  )


# Harvest rate implementation error
hr_err <- sr_data |> 
  # Calculate annual harvest rates on Somass and Hucuktlis
  mutate(group = if_else(stock == "HUC", "Hucuktlis", "Somass")) |> 
  summarize(
    .by = c(group, year),
    across(c(N, H), sum)
  ) |> 
  mutate(hr = H/N) |> 
  pivot_wider(
    id_cols = year,
    names_from = group,
    values_from = hr,
    names_glue = "{group}_hr"
  ) |> 
  right_join(fcst_err) |> 
  # Calculate harvest rate error
  mutate(hr_err = Somass_hr-target_hr) 


# Bootstrap mean & SD of HR implementation error
hr_err_boot <- hr_err |> 
  select(hr_err) |> 
  expand_grid(sample = seq(1:1000)) |> 
  nest(.by = sample) |> 
  rowwise() |> 
  mutate(boot = list(sample(data, length(data), replace = TRUE))) |> 
  select(boot) |> 
  unnest(boot) |> 
  summarize(
    mean = mean(hr_err),
    sd = sd(hr_err)
  )


# Fit simple LM to relate Hucuktlis HR to Somass HR
Hucuktlis_hr_lm <- hr_err |> 
  distinct(year, Somass_hr, Hucuktlis_hr) |> 
  filter(year > 2011) %>% # Constrain to post-Maa-nulth Treaty time period 
  lm(Hucuktlis_hr ~ Somass_hr, data = .)


summary(Hucuktlis_hr_lm) # Fit looks pretty strong, 
# should be good for simulation purposes


# Plot Hucuktlis versus Somass HR in recent years
(hr_corr_p <- sr_data |> 
    # Calculate annual harvest rates on Somass and Hucuktlis
    mutate(group = if_else(stock == "HUC", "Hucuktlis", "Somass")) |> 
    filter(year > 2011) |> 
    summarize(
      .by = c(group, year),
      across(c(N, H), sum)
    ) |> 
    mutate(hr = H/N) |> 
    pivot_wider(
      id_cols = year,
      names_from = group,
      values_from = hr,
      names_glue = "{group}_hr"
    ) |> 
    ggplot(aes(x = Somass_hr, y = Hucuktlis_hr)) +
    geom_abline(slope = 1, lty = 2) +
    geom_smooth(method = "lm") +
    geom_point(aes(colour = year), size = 2) +
    geom_text_repel(aes(label = year)) +
    scale_x_continuous(
      limits = c(0, 0.75),
      expand = c(0, 0),
      labels = scales::percent,
      oob = scales::oob_keep
    ) +
    scale_y_continuous(
      limits = c(0, 0.75),
      expand = c(0, 0),
      labels = scales::percent,
      oob = scales::oob_keep
    ) +
    scale_color_viridis_c() + 
    coord_fixed(1) +
    labs(
      x = "Somass harvest rate",
      y = "Hucuktlis harvest rate"
    )
)


# Save the plot
ggsave(
  plot = hr_corr_p,
  filename = here(
    "3. outputs",
    "Plots",
    "Hucuktlis_vs_Somass_HRs_2012-2024.png"
  ),
  width = 6,
  units = "in",
  dpi = "print"
)
  



# Set up the simulation ----------------------------------------


# Subsample posterior draws from Sproat to match GCL and Hucuktlis
SPR_AR1_sub <- AR1_mods[["SPR_AR1"]] |> 
  rstan::extract() |> 
  # Subset the various posterior components to 12000 rows according to class
  map(
    function(x) {
      if(is.array(x) & length(dim(x)) == 3) {
        x[sample(nrow(x), size = 12000, replace = FALSE),,]
      } 
      else if(is.matrix(x)) {
        x[sample(nrow(x), size = 12000, replace = FALSE),]
      }
      else {
        x[sample(nrow(x), size = 12000, replace = FALSE)]
      }
    }
  )


# Select only 1 model for Hucuktlis
AR1_subset <- list(
  "SPR" = SPR_AR1_sub,
  "GCL" = rstan::extract(AR1_mods[["GCL_AR1"]]),
  "HUC" = rstan::extract(AR1_mods[["HUC_full_enh_AR1"]])
)


# Select which lnalpha and beta values to use for Hucuktlis
names(AR1_subset[["HUC"]]) <- AR1_subset[["HUC"]] |> 
  names() |> 
  str_remove_all("_unenh") # erase the "_unfert" suffix
# Renaming lnalpha_unfert and beta_unfert identifies these as the values
# of choice for this model, since all code indexes into "lnalpha" and "beta"

# Initialize for loop products
pi_samps <- array(NA, dim = c(nrow(AR1_subset[[1]]$beta), ncol(AR1_subset[[1]]$pi), length(AR1_subset)))
p_samps <- array(NA, dim = c(nrow(AR1_subset[[1]]$beta), ncol(AR1_subset[[1]]$pi), length(AR1_subset), 3))
sig_R_samps <- NULL
model_samps_list <- list()


# Iterate over each model in AR1_subset
for (i in seq_along(AR1_subset)) {
  
  # Extract the current model
  model <- AR1_subset[[i]]
  
  # Number of age classes
  A <- ncol(model$pi)
  
  # CU name
  cu <- names(AR1_subset)[i]
  
  # Number of spawner years 
  nyrs <- ncol(model$S)
  
  # Number of recruitment years
  nRyrs <- ncol(model$R)
  
  # Extract most recemt spawner and recruitment states
  S_states <- model$S[, (nyrs - A + 1):nyrs]
  R_states <- model$R[, (nRyrs - A + 2):nRyrs]
  
  # Extract the last year's recruitment residual
  last_resid <- model$lnresid[, nRyrs]
  
  # Combine parameters into a single matrix
  sub_samps <- cbind(
    exp(model$lnalpha),
    model$beta,
    S_states,
    R_states,
    last_resid
  )
  
  # Assign column names to sub_samps
  colnames(sub_samps) <- c(
    paste0("alpha_", cu),
    paste0("beta_", cu),
    paste0("S_yr", (nyrs - A + 1):nyrs, "_", cu),
    paste0("R_yr", (nRyrs - A + 2):nRyrs, "_", cu),
    paste0("last_resid_", cu)
  )
  
  # Extract pi samples and assign column names
  pi_samps <- model$pi
  colnames(pi_samps) <- paste("pi_a", 1:A, cu, sep = "_")
  
  # Extract p samples for the last 3 recruitment years and assign column names
  p_samps_list <- list()
  for (j in 1:3) {
    p_samps <- model$p[, nyrs + j, ]
    colnames(p_samps) <- paste0("p_yr", nyrs + j, "_a", 1:A, "_", cu)
    p_samps_list[[j]] <- p_samps
  }
  
  # Stack the list of 3 p_samps matrices into a 3D array: [draws x ages x 3 years]
  p_3yr_array <- simplify2array(p_samps_list)  # dims: [n_draws, A, 3]
  
  # Apply median across the 3rd dimension (years) for each draw and age class
  p_avg <- apply(p_3yr_array, c(1, 2), median)
  
  # Normalize each row to sum to 1
  p_avg <- p_avg / rowSums(p_avg)
  
  # Assign column names
  colnames(p_avg) <- paste0("p_avg_a", 1:A, cu)
  
  # Combine into final model_samps matrix
  model_samps <- cbind(sub_samps, pi_samps, p_avg)
  
  # Store the model_samps matrix in the list
  model_samps_list[[cu]] <- model_samps
  
  # Time series of recruitment residuals (most recent 43 years)
  sig_R_samps <- cbind(sig_R_samps, apply(AR1_subset[[i]]$lnresid, 2, median)[(nRyrs-43):nRyrs])
  # using 43 is a bit of a hacky solution to match up the uneven time series 
  # lengths at the end year 
}


# Variance covariance matrix for recruitment residuals
epi <- cov(sig_R_samps)
colnames(epi) <- names(model_samps_list)

# Define function used to slice posterior states from the models ----------------------------


# below fun adapted from BC's kusko code (https://github.com/brendanmichaelconnors/Kusko-harvest-diversity-tradeoffs/blob/master/functions.R#L237)
# NB: used ChatGPT to adjust to match model_samps structure
process_iteration <- function(model_samps_row) {
  
  # Extract parameter names
  param_names <- names(model_samps_row)
  
    # Extract parameters
  alpha <- unname(model_samps_row[grep("^alpha_", param_names)])
  beta <- unname(model_samps_row[grep("^beta_", param_names)])
  last_resid <- unname(model_samps_row[grep("^last_resid_", param_names)])
  
  # Extract p values
  ps <- unname(model_samps_row[grep("^p_avg_a", param_names)])
  
  # Extract state variables S and R
  S_indices <- grep("^S_", param_names)
  R_indices <- grep("^R_", param_names)

  S <- matrix(model_samps_row[S_indices], nrow = length(S_indices), ncol = 1)
  R <- matrix(model_samps_row[R_indices], nrow = length(R_indices), ncol = 1)
  
  # Create output list
  output <- list(
    alpha = alpha,
    beta = beta,
    last_resid = last_resid,
    S = S,
    R = R,
    ps = ps
  )
  
  return(output)
}


# Have a look at an example run of the process_iteration function
lapply(model_samps_list, function(samps) {
  sampled_row <- samps[sample(nrow(samps), 1), ]
  process_iteration(sampled_row)
})



# Stipulate Harvest Control Rules -----------------------------------------


# Current management plan
status_quo_HCR <- function(fcst_run) {
  
  if (fcst_run < 200000)  {Somass_HR <- 0}
  else if (fcst_run < 250000)  {Somass_HR <- 0.15}
  else if (fcst_run < 300000)  {Somass_HR <- 0.2}
  else if (fcst_run < 350000)  {Somass_HR <- 0.229166667}
  else if (fcst_run < 400000)  {Somass_HR <- 0.25}
  else if (fcst_run < 450000)  {Somass_HR <- 0.291666667}
  else if (fcst_run < 500000)  {Somass_HR <- 0.324074074}
  else if (fcst_run < 550000)  {Somass_HR <- 0.35}
  else if (fcst_run < 600000)  {Somass_HR <- 0.397727273}
  else if (fcst_run < 650000)  {Somass_HR <- 0.4375}
  else if (fcst_run < 700000)  {Somass_HR <- 0.471153846}
  else if (fcst_run < 750000)  {Somass_HR <- 0.5}
  else if (fcst_run < 800000)  {Somass_HR <- 0.522222222}
  else if (fcst_run < 850000)  {Somass_HR <- 0.541666667}
  else if (fcst_run < 900000)  {Somass_HR <- 0.558823529}
  else if (fcst_run < 950000)  {Somass_HR <- 0.574074074}
  else if (fcst_run < 1000000)  {Somass_HR <- 0.587719298}
  else if (fcst_run < 1050000)  {Somass_HR <- 0.6}
  else if (fcst_run < 1100000)  {Somass_HR <- 0.618253968}
  else if (fcst_run < 1150000)  {Somass_HR <- 0.634848485}
  else if (fcst_run < 1200000)  {Somass_HR <- 0.65}
  else if (fcst_run < 1250000)  {Somass_HR <- 0.658928571}
  else if (fcst_run < 1300000)  {Somass_HR <- 0.667142857}
  else if (fcst_run < 1350000)  {Somass_HR <- 0.674725275}
  else if (fcst_run < 1400000)  {Somass_HR <- 0.681746032}
  else if (fcst_run < 1450000)  {Somass_HR <- 0.688265306}
  else if (fcst_run < 1500000)  {Somass_HR <- 0.694334975}
  else {Somass_HR <- 0.7}
 
 Somass_HR_realized <- Somass_HR + Somass_HR * rnorm(1, mean = hr_err_boot$mean, sd = hr_err_boot$sd)
 if(Somass_HR_realized < 0) {Somass_HR_realized <- 0}
 if(Somass_HR_realized >= 1) {Somass_HR_realized <- 0.99}
  
 # Predict the Hucuktlis Harvest Rate
 huc_pred <- predict(
   Hucuktlis_hr_lm, 
   newdata = data.frame(Somass_hr = Somass_HR_realized), 
   se.fit = TRUE
 )
 
 Hucuktlis_HR <- if(Somass_HR_realized == 0) {0} else {
   # Can go below 0 but this is fixed in the simulation
   rnorm(1, mean = huc_pred$fit, sd = huc_pred$residual.scale)
 }
 
  return(
    list(
      fcst_run = fcst_run,
      Somass = Somass_HR_realized,
      Hucuktlis = Hucuktlis_HR
    )
  )
}


# See how the HCR function performs on some random data
rlnorm(n = 1000, meanlog = log(500000), sdlog = 0.5) |>
  set_names(1:1000) |> 
  map(status_quo_HCR) |> 
  map(enframe) |> 
  list_rbind(names_to = "iter") |> 
  unnest(value) |> 
  pivot_wider() |> 
  pivot_longer(c(Somass, Hucuktlis)) |> 
  ggplot(aes(x = fcst_run, y = value)) +
  geom_line(aes(colour = name)) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = c(0, 0),
    labels = scales::percent
  ) +
  scale_x_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_colour_viridis_d(end = 0.8, direction = -1) +
  labs(x = "Forecast Somass run", y = "Harvest rate", colour = "Stock")


# Simulation function -----------------------------------------------------



simulate_forward <- function(
    model_samps_list, 
    phi = 0.75, 
    cov_matrix, 
    n_years = 25,
    HCR_fn
) {
  
  ns <- length(model_samps_list)      # Number of CUs
  A <- 4                              # Age classes (3-6)
  total_years <- A + n_years          # Buffer for age lag
  
  # Prep storage
  R <- matrix(0, total_years, ns)
  predR <- matrix(0, total_years, ns)
  v <- matrix(0, total_years, ns)
  N <- array(0, dim = c(total_years, A, ns))
  S <- matrix(0, total_years, ns)
  HR <- matrix(0, total_years, ns)
  C <- matrix(0, total_years, ns)
  
  # ----- 1. Initialize from posterior samples -----
  init_vals <- lapply(model_samps_list, function(samps) {
    sampled_row <- samps[sample(nrow(samps), 1), ]
    process_iteration(sampled_row)
  })
  
  alpha <- sapply(init_vals, function(x) x$alpha)
  beta <- sapply(init_vals, function(x) x$beta)
  last_resid <- sapply(init_vals, function(x) x$last_resid)
  ps <- lapply(init_vals, function(x) x$ps)
  
  for (i in 1:ns) {
    R[1:3, i] <- init_vals[[i]]$R
    S[4:7, i] <- init_vals[[i]]$S
    v[4, i] <- last_resid[i]
  }
  
  # ----- 2. Simulate recruitment errors -----
  epi <- MASS::mvrnorm(n = total_years, mu = rep(0, ns), Sigma = cov_matrix)
  
  # ----- 3. Fill years 4–7 -----
  for (i in 1:ns) {
    mu <- log(alpha[i]) + log(S[4, i]) - beta[i] * S[4, i]
    predR[4, i] <- exp(mu)
    R[4, i] <- predR[4, i] * exp(phi * v[4, i] + epi[4, i])
    v[4, i] <- log(R[4, i]) - log(predR[4, i])
  }
  
  for (t in 5:7) {
    for (i in 1:ns) {
      mu <- log(alpha[i]) + log(S[t, i]) - beta[i] * S[t, i]
      predR[t, i] <- exp(mu)
      R[t, i] <- predR[t, i] * exp(phi * v[t - 1, i] + epi[t, i])
      v[t, i] <- log(R[t, i]) - log(predR[t, i])
    }
  }
  
  for (t in 4:7) {
    for (a in 1:A) {
      true_age <- a + 2
      for (i in 1:ns) {
        brood_year <- t - true_age
        if (brood_year > 0) {
          N[t, a, i] <- R[brood_year, i] * ps[[i]][a]
        }
      }
    }
  }
  
  v[is.nan(v)] <- 0
  
  # ----- 4. Forward simulation loop -----
  for (t in 8:total_years) {
    
    ## I. Calculate N[t, , ]
    for (a in 1:A) {
      true_age <- a + 2
      for (i in 1:ns) {
        brood_year <- t - true_age
        if (brood_year > 0) {
          N[t, a, i] <- R[brood_year, i] * ps[[i]][a]
        }
      }
    }
    
    ## II. Harvest Rule
    somass_N <- sum(N[t, , 1:2])
    #if (is.na(somass_N)) { somass_N <- 0 }
    somass_fcst <- somass_N + somass_N * rnorm(1, fcst_err_boot$mean, fcst_err_boot$sd) # Add forecast error
    # Note, it is ok if somass_fcst occasionally dips below 0, the HCR_fn will set HR = 0 for negative inputs
    HR_t <- HCR_fn(somass_fcst) # Includes implementation error
    HR[t, 1:2] <- HR_t$Somass
    HR[t, 3] <- HR_t$Hucuktlis
    if (HR[t, 3] < 0) { HR[t, 3] <- 0 } # Ensure Hucuktlis HR doesn't drop below 0
    
    ## III. Spawners and Catch
    for (i in 1:ns) {
      S[t, i] <- sum(N[t, , i] * (1 - HR[t, i]))
      C[t, i] <- sum(N[t, , i] * HR[t, i])
    }
    
    ## IV. Recruitment
    for (i in 1:ns) {
      mu <- log(alpha[i]) + log(S[t, i]) - beta[i] * S[t, i]
      predR[t, i] <- exp(mu)
      R[t, i] <- predR[t, i] * exp(phi * v[t - 1, i] + epi[t, i])
      v[t, i] <- log(R[t, i]) - log(predR[t, i])
    }
  }
  
  # ----- 5. Combine into a dataframe -----
  output <- data.frame()
  for (i in 1:ns) {
    df_cu <- data.frame(
      year = 1:total_years,
      CU = paste0("CU", i),
      S = S[, i],
      R = R[, i],
      C = C[, i],
      HR = HR[, i],
      N_age3 = N[, 1, i],
      N_age4 = N[, 2, i],
      N_age5 = N[, 3, i],
      N_age6 = N[, 4, i]
    )
    output <- rbind(output, df_cu)
  }
  
  rownames(output) <- NULL
  return(output)
}


# Test out one iteration of the simulation
sim_result <- simulate_forward(
  model_samps_list, 
  phi = 0.75,
  cov_matrix = epi,
  HCR_fn = status_quo_HCR
)


# Plot the resulting time series
sim_result |> 
  rowwise() |> 
  mutate(N = sum(c_across(contains("N_age")))) |> 
  pivot_longer(cols = c(S, C, R, HR, N)) |> 
  ggplot(aes(x = year, y = value)) +
  facet_grid(name ~ CU, scales = "free_y") +
  geom_line()


# Run the simulation 1000 times -------------------------------------------

sims <- c(1:1e3) |> 
  set_names() |> 
  map(
    \(x) 
    simulate_forward(
      model_samps_list, 
      phi = 0.75, 
      cov_matrix = epi,
      HCR_fn = status_quo_HCR
    )
  ) |> 
  list_rbind(names_to = "sim")


# Plot the resulting time series
sims |> 
  rowwise() |> 
  mutate(N = sum(c_across(contains("N_age")))) |> 
  pivot_longer(cols = c(S, C, R, HR, N)) |> 
  ggplot(aes(x = year, y = value)) +
  facet_grid(name ~ CU, scales = "free_y") +
  geom_line(aes(group = sim), linewidth = 0.2, alpha = 0.2)
