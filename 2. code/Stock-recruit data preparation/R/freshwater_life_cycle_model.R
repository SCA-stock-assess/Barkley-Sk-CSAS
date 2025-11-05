# Packages and setup --------------------------------------------------------


pkgs <- c(
  "tidyverse", "rstan", "here", "bayesplot", "tidybayes", 
  "geomtextpath", "writexl", "readxl", "ggridges", "Hmisc"
)
#install.packages(pkgs)


library(here)
library(tidyverse); theme_set(theme_bw())
library(rstan)
library(readxl)
library(bayesplot)
library(tidybayes)
library(geomtextpath)
library(ggridges)
library(writexl)


# Logit transformation
logit <- function(p) {log(p/(1-p))}

# Lake names 
lake_names <- c("Great Central", "Sproat", "Hucuktlis")


# Use parallel processing when running Stan models
options(mc.cores = parallel::detectCores())


# Load and summarize smolt size data -----------------------------------------


# Raw data
smolt_sizes <- c("GCL", "SPR", "HEN") |> 
  purrr::set_names() |> 
  map(
    \(x) here(
      "1. data",
      "Hyatt, Stiff, Rankin smolt data",
      paste0(x, "_Smolt_Export.xlsx")
    )
  ) |> 
  imap(
    \(x, idx) read_xlsx(
      path = x,
      sheet = paste(idx, "Smolt Data (all)")
    )
  ) |> 
  list_rbind(names_to = "lake") |> 
  janitor::clean_names() |> 
  mutate(
    lake = factor(
      lake,
      levels = c("GCL", "SPR", "HEN"),
      labels = lake_names
    ),
    gear = tolower(gear),
    ln_wt = log(fresh_std_weight)
  )


# Does size-at-age change over dates? If so, need to correct for CSRs
smolt_sizes |> 
  mutate(day = as.numeric(format(sample_date, "%j"))) |> 
  summarize(
    .by = c(year, day, lake, fnlage),
    mean_wt = mean(fresh_std_weight),
    n = n()
  ) |> 
  filter(fnlage %in% c(1, 2)) |> 
  ggplot(aes(x = day, y = mean_wt, colour = year)) +
  facet_grid(
    fnlage ~ lake,
    scales = "free"
  ) +
  geom_point(aes(size = n)) +
  geom_line(aes(group = year)) +
  geom_smooth(method = "lm", aes(weight = n)) +
  scale_color_viridis_c()
# Sizes probably decline slightly over time. Should correct for CSR. 


# Annual smolt survey sample sizes
smolt_catch <- here(
  "1. data",
  "Hyatt, Stiff, Rankin smolt data",
  "Smolt Sample Metadata 25.01.23.xls"
) |> 
  map2(
    c("GCL", "SPR", "HEN"),
    \(x, y) read_xls(x, sheet = y)
  ) |> 
  list_rbind() |> 
  janitor::clean_names() |> 
  mutate(
    lake = factor(
      lake,
      levels = c("GCL", "SPR", "HEN"),
      labels = lake_names
    ),
    gear = tolower(gear_type),
    catch = if_else(
      is.na(total_catch) | (total_catch < total_retained),
      total_retained, 
      total_catch
    )
  ) |> 
  select(lake, sample_date, gear, catch)
  

# Join the catch data by date and calculate CSRs
smolt_sizes_csr <- smolt_sizes |> 
  add_count(sample_date, lake, gear, name = "samples") |> 
  left_join(smolt_catch) |> 
  mutate(csr = if_else(samples >= catch, 1, catch/samples))


# Calculate annual weight summaries corrected for CSRs
smolt_sizes_summary <- smolt_sizes_csr |> 
  filter(
    fnlage %in% c(1, 2),
    !is.na(ln_wt)
  ) |> 
  summarize(
    .by = c(lake, year, fnlage),
    ln_mean = sum(csr * ln_wt) / sum(csr),
    ln_sd = sqrt(sum(csr * (ln_wt - ln_mean)^2) / sum(csr)),
    n = n()
  ) |> 
  complete(lake, year, fnlage) |> 
  pivot_wider(
    names_from = fnlage,
    values_from = c(ln_mean, ln_sd, n),
    names_glue = "WO{fnlage}_{.value}"
  )


# Load and prepare raw fry and smolt abundance data  ----------------


# Load the smolt outmigration data
smolt_data <- here(
  "3. outputs",
  "Stock-recruit data",
  "Inferred_annual_smolt_age_composition.csv"
) |> 
  read.csv() |> 
  rename_with(\(x) str_replace_all(x, "\\.", "x")) |> 
  # Trim time series to include only the studied years
  filter(!(year > 2013 & if_all(contains("rate"), is.na)))  |> 
  # Convert annual rates to counts
  mutate(
    across(
      c(rate_age1, rate_age2), 
      \(x) round(x*sample_days, 0),
      .names = "{str_replace(col, 'rate', 'count')}"
    ),
    smolting_rate = 1-(rate_age1x/(rate_age1x + rate_age1)),
    # Remove exact 0s and 1s from age compositions and smolting rate
    across(
      c(matches("^prop_age"), smolting_rate),
      \(x) case_when(
        x == 1 ~ 0.999999,
        x == 0 ~ 1e-6,
        .default = x
      )
    ),
    # Ensure age compositions sum to 1
    across(
      c("prop_age1", "prop_age2", "prop_age1x"),
      \(x) x / (prop_age1 + prop_age2 + prop_age1x)
    )
  )


# Plot calculated age proportions
smolt_data |> 
  pivot_longer(
    matches("(rate|prop)_"),
    names_sep = "_",
    names_to = c(".value", "age")
  ) |> 
  ggplot(aes(x = year, y = prop)) +
  facet_wrap(
    ~ lake,
    ncol = 1,
    strip.position = "right"
  ) +
  geom_point(
    aes(fill = age, size = rate),
    colour = "black",
    shape = 21
  ) +
  scale_y_continuous(labels = scales::percent) +
  coord_cartesian(expand = FALSE) +
  labs(y = "Percent composition") +
  theme(panel.spacing.y = unit(1, "lines"))


# fry ages (from results received in 2025)
fry_age <- here(
  "1. data",
  "smolt size data.xlsx"
) |> 
  read_xlsx(sheet = "fry_ages") |> 
  mutate(
    lake = case_when(
      lake == "gcl" ~ "Great Central",
      lake == "huc" ~ "Hucuktlis",
      lake == "spr" ~ "Sproat"
    )
  ) |> 
  summarize(
    .by = c(smolt_year, lake),
    across(matches("count_fry\\d"), sum)
  )


# Load the annual ATS estimates (corrected for variation in survey date)
ats_data <- here(
  "3. outputs",
  "Stock-recruit data",
  "GAMM-estimated_pre-smolt_time_series.csv"
) |> 
  read.csv() |> 
  rename("year" = smolt_year) |> 
  rename_with(\(x) paste0("ats_", x), .cols = c(est, lwr_80, upr_80, cv)) |> 
  # Remove NAs at start and end of time series
  filter(!is.na(ats_est)) |> 
  group_by(lake) |> 
  complete(year = full_seq(min(year):max(year), 1)) |> 
  ungroup() |> 
  mutate(ats_se = ats_est*ats_cv)


# View which years are missing estimates
ats_data |>
  ggplot(aes(x = year, y = ats_est)) +
  facet_wrap(
    ~lake,
    ncol = 1,
    scales = "free_y"
  ) +
  geom_pointrange(
    aes(
      ymin = ats_lwr_80, 
      ymax = ats_upr_80,
      colour = interpolated
    )
  ) +
  scale_colour_manual(values = c("black", "grey"))


# Join ATS, smolt, and fry data
lakes_data <- ats_data |> 
  left_join(
    select(
      .data = smolt_data, 
      lake, year, contains("count"), smolting_rate, prop_age2
    )
  ) |> 
  # Add infilled smolt weight data
  left_join(smolt_sizes_summary) |> 
  # Add fry ages
  full_join(
    fry_age,
    by = c("lake", "year" = "smolt_year")
  ) |> 
  mutate(
    is_observed_count = !if_any(matches("count_age"), is.na),
    is_observed_fry = !if_any(matches("count_fry"), is.na),
    is_observed_ats = interpolated == "No",
    is_observed_WO1 = !if_any(matches("WO1"), is.na),
    is_observed_WO2 = !if_any(matches("WO2"), is.na), 
    count_total = count_age1 + count_age2,
    lake = factor(lake, levels = lake_names)
  ) |> 
  arrange(lake, year) |> 
  # Remove two interpolated years at beginning of Sproat time series
  filter(!(lake == "Sproat" & year < 1980))



# Format data for Stan ----------------------------------------------------


# Function to create observed Stan data for each lake
# As a general note: priors should reflect what is *possible*
# not what we think is likely.
make_stan_data <- function(
    lake_name,
    # User-defined second year mortality and theta priors 
    # (with sensible default values)
    mu_M_prior = logit(0.5),
    sigma_M_prior = 0.7,
    mu_theta_prior = logit(0.9),
    sigma_theta_prior = 1
) {
  
  # Filter data for a single lake
  lake_data <- lakes_data %>%
    filter(lake == lake_name) |> 
    # Convert NAs to 0s
    mutate(across(matches("ats|count|W"), \(x) if_else(is.na(x), 0, x)))

  # Identify missing rows in observations
  is_observed_count <- lake_data$is_observed_count
  is_observed_ats <- lake_data$is_observed_ats
  is_observed_fry <- lake_data$is_observed_fry
  is_observed_WO1 <- lake_data$is_observed_WO1
  is_observed_WO2 <- lake_data$is_observed_WO2
  
  # Vector of age1 counts (note that these counts are EXPANDED to account for 
  # total smolt captures on the days samples were taken) and total counts
  A1_obs <- lake_data$count_age1
  A_total <- lake_data$count_total
  
  # Vector of observed fry ages from in-lake trawl samples
  a1_obs <- lake_data$count_fry1
  a_total <- lake_data$count_fry1 + lake_data$count_fry2
  
  # Vector of observed total lake abundance
  N_obs <- lake_data$ats_est
  
  # Age 1 and age 2 smolt weight priors
  w1_theta_prior <- mean(lake_data$WO1_ln_mean, na.rm = TRUE)
  w1_sigma_prior <- sd(lake_data$WO1_ln_mean, na.rm = TRUE)
  w2_theta_prior <- mean(lake_data$WO2_ln_mean, na.rm = TRUE)
  w2_sigma_prior <- sd(lake_data$WO2_ln_mean, na.rm = TRUE)
  
  # Smolting rate prior
  theta_prior <- lake_data |> 
    mutate(
      holdover_rate = 1-smolting_rate,
      # Adjust holdover rate according to assumed overwinter mortality
      # Reflects our understanding that the prior will always be an 
      # overestimate of the smolting rate
      holdover_rate_adj = holdover_rate/(1-binomial(link = "logit")$linkinv(mu_M_prior)),
      smolting_rate_adj = 1-holdover_rate_adj,
      theta = logit(smolting_rate_adj)
    ) |> 
    summarize(
      mu_theta = mean(theta, na.rm = TRUE),
      sigma_theta = sd(theta, na.rm = TRUE)
    )
  
  # Alternative specification of smolting rate prior
  # PT suggestion after reviewing the informative prior above and
  # observing that it seems to draw the model away from what we assume
  # is going on for GCL
  theta_prior_alt <- data.frame(
    mu_theta = mu_theta_prior,
    sigma_theta = sigma_theta_prior
  )
  
  # Age-1 fry proportion prior
  pN1_prior <- lake_data |> 
    mutate(pN1 = logit(1-prop_age2)) |> 
    summarize(
      mu_pN1 = mean(pN1, na.rm = TRUE),
      sigma_pN1 = sd(pN1, na.rm = TRUE)
    )
  
  # Initial age2 fry population prior
  N2_init_prior <- lake_data |> 
    filter(year <= min(year) + 4) |>
    # Use average proportion from earliest five years
    summarize(
      prop_age2 = median(prop_age2),
      across(c(ats_est, ats_se), \(x) x[which.min(year)])
    ) |> 
    mutate(
      N2_init = prop_age2 * ats_est,
      N2_init_sigma = prop_age2 * ats_se,
      # Need to transform N2_init mean and SD for lognormal distribution:
      # https://msalganik.wordpress.com/2017/01/21/making-sense-of-the-rlnorm-function-in-r/
      N2_init_log_mean = log(N2_init^2 / sqrt(N2_init_sigma^2 + N2_init^2)),
      N2_init_log_sigma = sqrt(log(1 + (N2_init_sigma^2 / N2_init^2))),
    )
  
  # Age 1 fry abundance prior
  N1_prior <- lake_data |> 
    mutate(N1 = (1-prop_age2)*ats_est) |> 
    summarize(
      mu_N1_ln = mean(log(N1), na.rm = TRUE),
      sigma_N1_ln = sd(log(N1), na.rm = TRUE)
    )
  
  # Observation error on lake ATS estimates
  obs_error <- lake_data |> 
    # summarize(
    #   ats_se = mean(ats_se),
    #   ats_est = mean(ats_est)
    # ) |> 
    mutate(obs_error = sqrt(log(1 + (ats_se^2 / ats_est^2)))) |> 
    pull(obs_error)
  
  # Stan data list
  stan_data <- list(
    Y = nrow(lake_data),
    year = lakes_data$year,
    A1_obs = A1_obs,
    A_total = A_total,
    a1_obs = a1_obs,
    a_total = a_total,
    N_obs = N_obs,
    # Age-specific weight data are already log-transformed in the raw data
    # (i.e. values are mean and sd of the log-transformed individual weights)
    WO1_ln_mean = lake_data$WO1_ln_mean,
    WO1_ln_sd = lake_data$WO1_ln_sd,
    WO1_n = lake_data$WO1_n,
    WO2_ln_mean = lake_data$WO2_ln_mean,
    WO2_ln_sd = lake_data$WO2_ln_sd,
    WO2_n = lake_data$WO2_n,
    is_observed_count = as.integer(is_observed_count),
    is_observed_ats = as.integer(is_observed_ats),
    is_observed_fry = as.integer(is_observed_fry),
    is_observed_WO1 = as.integer(is_observed_WO1),
    is_observed_WO2 = as.integer(is_observed_WO2),
    obs_error = obs_error,
    #mu_theta_prior = theta_prior$mu_theta,
    #sigma_theta_prior = theta_prior$sigma_theta,
    mu_theta_prior = theta_prior_alt$mu_theta, # alternative user specification
    sigma_theta_prior = theta_prior_alt$sigma_theta, # alternative user specification
    mu_M_prior = mu_M_prior,
    sigma_M_prior = sigma_M_prior,
    mu_N2_init_prior = N2_init_prior$N2_init_log_mean,
    sigma_N2_init_prior = N2_init_prior$N2_init_log_sigma,
    #w1_theta_prior = w1_theta_prior,
    #w1_sigma_prior = w1_sigma_prior,
    #w2_theta_prior = w2_theta_prior,
    #w2_sigma_prior = w2_sigma_prior,
    #mu_pN1_prior = pN1_prior$mu_pN1,
    #mu_pN1_prior = logit(0.8),
    #sigma_pN1_prior = 2,
    #sigma_pN1_prior = pN1_prior$sigma_pN1 # might try wider prior?
    mu_N1_prior = N1_prior$mu_N1_ln,
    sigma_N1_prior = N1_prior$sigma_N1_ln
  )
  
  return(stan_data)
}


# Create nested list of Stan input data for each lake
lakes_stan_data <- unique(lakes_data$lake) |> 
  set_names() |> 
  map2(
    # Second winter mortality rate priors for GCL, HUC, SPR, respectively
    .y = logit(c(0.3, 0.6, 0.6)),
    \(x, y) make_stan_data(
      lake_name = x, 
      #mu_M_prior = y # Try with all default values first
    )
  )


# View prior distributions for theta, mu_M, and N2_init
some_priors <- lakes_stan_data |> 
  # Discard all parameters with year-specific estimates
  map(\(x) discard(.x = x, .p = \(y) length(y) > 1)) |> 
  map(\(x) imap(x, as_tibble_col)) |> 
  map(\(x) list_cbind(unname(x))) |> 
  list_rbind(names_to = "lake") |> 
  rowwise() |> 
  mutate(
    n = 1e4,
    mu_M = list(
      rnorm(
        n = n,
        mean = mu_M_prior,
        sd = sigma_M_prior
      )
    ),
    mu_theta = list(
      rnorm(
        n = n,
        mean = mu_theta_prior,
        sd = sigma_theta_prior
      )
    ),
    N2_init = list(
      rlnorm(
        n = n,
        meanlog = mu_N2_init_prior,
        sdlog = sigma_N2_init_prior
      )
    ),
    N1 = list(
      rlnorm(
        n = n,
        meanlog = mu_N1_prior,
        sdlog = sigma_N1_prior
      )
    )
  ) |>
  ungroup() |> 
  select(lake, mu_M, mu_theta, N2_init, N1) |> 
  unnest(everything()) |> 
  pivot_longer(!lake) %>%
  split(.$name) 

some_priors |> 
  imap(
    \(x, idx) x |> 
      ggplot(
        # Show distributions of M and theta priors in response (%) scale
        aes(
          if(idx %in% c("mu_M", "mu_theta")) {
            x = binomial()$linkinv(value)
          } else {
          x = value
          }
        )
      ) +
      facet_wrap(
        ~ lake,
        nrow = 1,
        scales = "free"
      ) +
      geom_density(fill = "lightgrey") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      scale_x_continuous(expand = c(0, 0), name = idx)
  )


# Fit Stan models and assess convergence ------------------------------------


# Compile and fit the Stan model
fit_stan_model <- function(stan_data, index, iter = 4000) {
  stan(
    file = here(
      "2. code",
      "Stock-recruit data preparation",
      "Stan",
      "smolt_production_model_AR1_mv-thetaM.stan"
    ),
    model_name = index,
    data = stan_data,
    iter = iter,
    chains = 4,
    control = list(
      adapt_delta = 0.995
      ,max_treedepth = 15
    )
  )
}


# Try fitting the model with just one lake first
# test_fit <- fit_stan_model(lakes_stan_data$`Great Central`, "Great Central")


# Fit models for each lake (toggle to TRUE to run)
if(
  FALSE
  #TRUE
) {
  lakes_stan_fits <- lakes_stan_data |> 
    imap(fit_stan_model, 6000) # Adjust such that iter is lake-specific?
}


# Save fitted models as RDS objects (toggle to TRUE to run)
if(
  #TRUE
  FALSE
) {
  lakes_stan_fits |> 
    iwalk(
      \(x, idx) x |> 
        saveRDS(
          file = here(
            "3. outputs",
            "Stock-recruit data",
            paste("Bayesian state-space", idx, "smolts_revised.rds")
          )
        )
    )
}


# Load fitted models from RDS files (if not already in global environment)
lakes_stan_fits <- if(exists("lakes_stan_fits")) {
  lakes_stan_fits
  } else {
  names(lakes_stan_data) |>
    purrr::set_names() |> 
    map(
      \(x) list.files(
        here(
          "3. outputs",
          "Stock-recruit data"
        ),
        pattern = paste("Bayesian state-space", x, "smolts_revised.rds"),
        full.names = TRUE
      )
    ) |> 
    map(readRDS)
}


#assess convergence###
worst_Rhat <- lakes_stan_fits |> 
  map(\(x) summary(x)$summary) |>  
  map(as.data.frame) |> 
  list_rbind(names_to = "lake") |> 
  mutate(Rhat = round(Rhat, 3)) |> 
  arrange(desc(Rhat))

worst_Rhat %>% 
  filter(n_eff>3) %>% 
  ggplot(aes(x = n_eff, y = Rhat))+
  facet_wrap(~lake) +
  geom_point()+
  geom_hline(yintercept = 1.01, lty = 2)+
  geom_vline(xintercept = 400, lty = 2)

head(worst_Rhat, n = 20)


# Trace plots
lakes_stan_fits |> 
  imap(
    \(x, idx)
    traceplot(
      x, 
      pars = rownames(filter(worst_Rhat, lake == idx))[1:20]
    ) +
      ggtitle(label = idx)
  )


# Pair plots
lakes_stan_fits |> 
  map(
    \(x) bayesplot::mcmc_pairs(
      x, 
      # AR1 params
      pars = c("obs_error_scale", "lp__", "phi", "sigma_proc")
    )
  )

lakes_stan_fits |> 
  map(
    \(x) bayesplot::mcmc_pairs(
      x, 
      # Holdover dynamics
      pars = c("mu_M", "mu_theta", "sigma_M", "sigma_theta")
    )
  )



# Plot predicted versus observed values -----------------------------------


# Dataframe with first year in time series for each lake
# used to properly code years in the Stan outputs
min_yrs <- lakes_data |> 
  summarize(
    .by = lake,
    min_yr = min(year)
  )


# Posterior values that are not year-specific
post_interannual <- lakes_stan_fits |> 
  map(extract) |> 
  # Discard all parameters with year-specific estimates
  map(\(x) discard(.x = x, .p = \(y) length(dim(y)) > 1)) |> 
  map(\(x) map(x, \(y) as_tibble(y, rownames = "draw"))) |> 
  map(\(x) list_rbind(x, names_to = "parameter")) |> 
  list_rbind(names_to = "lake") 


# Posterior values that are year-specific
post_annual <- lakes_stan_fits |> 
  map(extract) |> 
  # Keep only parameters with year-specific estimates
  map(\(x) keep(.x = x, .p = \(y) length(dim(y)) == 2 & dim(y)[2] > 10)) |> 
  map(\(x) map(x, \(y) as_tibble(y, .name_repair = NULL, rownames = "draw"))) |> 
  map(\(x) list_rbind(x, names_to = "parameter")) |> 
  list_rbind(names_to = "lake") |> 
  pivot_longer(
    cols = !c(lake, parameter, draw),
    names_to = "year",
    values_to = "value",
    values_drop_na = TRUE,
    names_transform = \(x) str_extract(x, "\\d+")
  ) |> 
  left_join(min_yrs) |> 
  mutate(year = as.numeric(year) - 1 + min_yr) |> 
  select(-min_yr)


# theta-M correlation coefficients
post_omega <- lakes_stan_fits |> 
  map(extract) |> 
  # Keep only parameters with year-specific estimates
  map(\(x) keep_at(x = x, at = "Omega")) |> 
  map(\(x) map(x, \(y) as_tibble(y[,1,2], rownames = "draw"))) |> 
  map(\(x) list_rbind(x, names_to = "parameter")) |> 
  list_rbind(names_to = "lake")


# Calculate brood year total outmigrants and biomass from individual draws
post_by_ttls <- post_annual |> 
  filter(str_detect(parameter, "(B|)O\\d")) |> 
  # Convert year (currently smolt year) to brood_year
  mutate(year = if_else(str_detect(parameter, "1"), year-2, year-3)) |> 
  pivot_wider(
    names_from = parameter,
    values_from = value
  ) |> 
  mutate(
    BYO = O1 + O2,
    BYB = BO1 + BO2
  ) |> 
  pivot_longer(
    cols = c(BYO, BYB),
    names_to = "parameter",
    values_to = "value"
  ) |> 
  select(lake, year, draw, parameter, value)


# Calculate smolt year total outmigrants and biomass from individual draws
post_sy_ttls <- post_annual |> 
  filter(str_detect(parameter, "(B|)O\\d")) |> 
  pivot_wider(
    names_from = parameter,
    values_from = value
  ) |> 
  mutate(
    SYO = O1 + O2,
    SYB = BO1 + BO2
  ) |> 
  pivot_longer(
    cols = c(SYO, SYB),
    names_to = "parameter",
    values_to = "value"
  ) |> 
  select(lake, year, draw, parameter, value)


# Bring both posterior dataframes together
posterior_df <- bind_rows(
  post_annual, 
  post_interannual,
  post_omega,
  post_by_ttls,
  post_sy_ttls
) |> 
  # Make parameter names more informative
  mutate(
    parameter_long = case_when(
      parameter == "O1" ~ "age1 smolts",
      parameter == "O2" ~ "age2 smolts",
      parameter == "BO1" ~ "age1 smolt biomass (g)",
      parameter == "BO2" ~ "age2 smolt biomass (g)",
      parameter == "w1" ~ "age1 smolt mean weight (g)",
      parameter == "w2" ~ "age2 smolt mean weight (g)",
      parameter == "p1" ~ "age1 smolt proportion",
      parameter == "p2" ~ "age2 smolt proportion",
      parameter == "N1" ~ "age1 fry",
      parameter == "N2" ~ "age2 fry",
      parameter == "N_lake" ~ "lake total fry abundance",
      parameter == "theta" ~ "smolting rate",
      parameter == "M" ~ "age1 to age2 mortality",
      parameter == "Omega" ~ "smolting rate-mortality correlation coefficient",
      parameter == "phi" ~ "AR1 correlation coefficient",
      parameter == "BYB" ~ "brood year smolt biomass production (g)",
      parameter == "BYO" ~ "brood year smolt production",
      parameter == "SYO" ~ "smolt year total outmigrant abundance",
      parameter == "SYB" ~ "smolt year total outmigrant biomass (g)",
      .default = "NO_ENTRY"
    )
  )


# Plot posterior versus prior distributions for mu_theta and mu_M
posterior_df |> 
  filter(parameter %in% c("mu_M", "mu_theta")) |> 
  ggplot(aes(x = binomial()$linkinv(value))) +
  facet_wrap(
    parameter ~ lake,
    scales = "free",
    labeller = labeller(.multi_line = FALSE)
  ) +
  geom_textdensity(
    label = "posterior",
    hjust = 0.3
  ) +
  geom_textdensity(
    data = some_priors |> 
      keep_at(c("mu_M", "mu_theta")) |> 
      list_rbind(names_to = "parameter"),
    colour = "red",
    label = "prior",
    hjust = 0.7
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "value")
  

# Plot posterior distributions for sigma_theta, sigma_M, Omega, and phi
posterior_df |> 
  filter(parameter %in% c("sigma_M", "sigma_theta", "Omega", "phi")) |> 
  ggplot(aes(x = value)) +
  facet_grid(
    parameter ~ lake,
    scales = "free",
    labeller = labeller(.multi_line = FALSE)
  ) +
  geom_vline(
    xintercept = 0,
    lty = 2
  ) +
  geom_density() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "value") +
  theme(panel.spacing.x = unit(1.5, "lines"))


# Annual parameters estimated
annual_params <- posterior_df |> 
  filter(!is.na(year)) |> 
  distinct(parameter) |> 
  pull() |> 
  str_subset("_raw", negate = TRUE)


# Dataframe with estimated parameters from raw data
obs_params <- lakes_data |> 
  arrange(lake, year) |> 
  mutate(
    .by = lake,
    O1 = (1-prop_age2)*smolting_rate*ats_est,
    O2 = prop_age2*ats_est,
    p1 = (1-prop_age2)*smolting_rate,
    p2 = prop_age2,
    theta = smolting_rate,
    M = case_when(
      lake == "Great Central" ~ lakes_stan_data[["Great Central"]]$mu_M_prior,
      lake == "Sproat" ~ lakes_stan_data[["Sproat"]]$mu_M_prior,
      lake == "Hucuktlis" ~ lakes_stan_data[["Hucuktlis"]]$mu_M_prior
    ) |> 
      binomial(link = "logit")$linkinv(),
    w1 = exp(WO1_ln_mean),
    w2 = exp(WO2_ln_mean),
    N_lake = ats_est,
    BO1 = O1*exp(WO1_ln_mean),
    BO2 = O2*exp(WO2_ln_mean),
    BYB = BO1 + lead(BO2),
    BYO = O1 + lead(O2),
    SYB = BO1+BO2,
    SYO = ats_est
  ) |> 
  select(lake, year, any_of(annual_params)) |> 
  pivot_longer(
    cols = !c(lake, year),
    names_to = "parameter",
    values_to = "obs"
  ) |> 
  # Ensure year is brood year for the BY total values
  mutate(year = if_else(parameter %in% c("BYB", "BYO"), year-2, year))


# Plot predicted versus estimated values for parameters of interest across years
# Ensure mean smolt weight is added to these plots
(out_plots <- posterior_df |> 
    filter(
      !is.na(value),
      parameter %in% annual_params,
    ) |> 
    summarize(
      .by = c(lake, parameter, parameter_long, year),
      quantiles = list(quantile(value, probs = c(0.025, 0.1, 0.5, 0.9, 0.975)))
    ) |> 
    unnest_wider(quantiles) |>
    # Add observed values
    left_join(obs_params) |> 
    filter(!parameter_long == "NO_ENTRY") %>%
    split(.$parameter_long) |> 
    imap(
      \(x, idx) x |> 
        ggplot(aes(x = year, y = `50%`)) +
        facet_wrap(
          ~lake,
          scales = "free_y",
          ncol = 1,
          strip.position = "right"
        ) +
        geom_linerange(
          aes(
            ymin = `2.5%`,
            ymax = `97.5%`
          ),
          linewidth = 0.3,
          alpha = 0.6
        ) +
        geom_pointrange(
          aes(
            ymin = `10%`,
            ymax = `90%`
          ),
          linewidth = 1
        ) +
        geom_point(aes(y = obs), color = "red") +
        labs(y = idx) +
        ggtitle(paste("Posterior versus estimated", idx))
    )
)


# Save output plots
out_plots |> 
  iwalk(
    \(x, idx) ggsave(
      plot = x,
      filename = here(
        "3. outputs",
        "Stock-recruit data",
        "Bayesian smolt model posterior",
        paste("Posterior vs estimated", idx, ".png")
      ),
      dpi = "print",
      width = 6.5, 
      height = 5.5,
      units = "in"
    )
  )


# Export posterior values of interest (with imputed years) -------------------


# Save full posterior_df as an RDS object
saveRDS(
  object = posterior_df,
  file = here(
    "3. outputs",
    "Stock-recruit data",
    "Freshwater_LifeCycle_model_full-posterior.RDS"
  )
)
# Large file, circa 130MB


# Summarize posterior for export
annual_estimates <- posterior_df |> 
  filter(
    parameter %in% str_subset(annual_params, "p\\d", negate = TRUE),
    !parameter %in% c("z_N1", "w1_ln", "w2_ln", "logN1", "obs_sigma"),
    !is.na(value)
  ) |> 
  mutate(
    family = case_when(
      str_detect(parameter, "N|O|BO|BY|SY") ~ "lognormal",
      parameter %in% c("theta", "M", "p1", "p2") ~ "logitnormal",
      parameter %in% c("w1", "w2") ~ "normal",
      .default = NA
    ),
    log = log(value),
    logit = logit(value)
  ) |> 
  summarize(
    .by = c(lake, parameter, year, family),
    quantiles = list(quantile(value, probs = c(0.025, 0.1, 0.5, 0.9, 0.975))),
    across(
      c(value, log, logit), 
      .fns = list("mean" = mean, "sd" = sd),
      .names = c("{.fn}_{.col}")
    )
  ) |> 
  unnest_wider(quantiles) |> 
  mutate(
    mu = case_when(
      family == "lognormal" ~ mean_log,
      family == "logitnormal" ~ mean_logit,
      family == "normal" ~ mean_value
    ),
    sigma = case_when(
      family == "lognormal" ~ sd_log,
      family == "logitnormal" ~ sd_logit,
      family == "normal" ~ sd_value
    )
  ) |> 
  select(-matches("(mean|sd)_(value|log|logit)"))


# Metadata sheet
metadata <- tibble(parameter = unique(annual_estimates$parameter)) |> 
  mutate(
    definition = case_when(
      parameter == "O1" ~ "Posterior estimates for age1 outmigrating smolts",
      parameter == "O2" ~ "Posterior estimates for age2 outmigrating smolts",
      parameter == "w1" ~ "posterior estimates for age1 smolt average weight (g)",
      parameter == "w2" ~ "Posterior estimates for age2 smolts average weight (g)",
      parameter == "BO1" ~ "Posterior estimates for age1 smolt biomass (g)",
      parameter == "BO2" ~ "Posterior estimates for age2 smolt biomass (g)",
      parameter == "BYO" ~ "Brood year sum of posterior estimates for O1 & O2 (see above)",
      parameter == "BYB" ~ "Brood year sum of posterior estimates for BO1 & BO2 (see above)",
      parameter == "SYO" ~ "Smolt year sum of posterior estimates for O1 & O2 (see above)",
      parameter == "SYB" ~ "Smolt year sum of posterior estimates for BO1 & BO2 (see above)",
      parameter == "N1" ~ "Posterior estimates for age1 fry abundance in the lake on 1 March",
      parameter == "N2" ~ "Posterior estimates for age2 fry abundance in the lake on 1 March",
      parameter == "N_lake" ~ "Posterior estimates for total fry abundance in the lake on 1 March",
      parameter == "theta" ~ "Posterior estimates for age1 smolting rate",
      parameter == "M" ~ "Posterior estimates for second year fry mortality (age1 to age2)",
      .default = "MISSING"
    ),
    year_type = case_when(
      str_detect(parameter, "BY") ~ "brood year",
      .default = "smolt year"
    )
  )


# Export posterior estimates and metadata to excel
list(
  "metadata" = metadata, 
  "model_estimates" = annual_estimates
) |> 
  write_xlsx(
    path = here(
      "3. outputs",
      "Stock-recruit data",
      "Bayesian_state-space_smolt-production_estimated_time_series.xlsx"
    )
  )


# Posterior summaries for ResDoc tables and writing -----------------------


# Median fry and smolt abundances
post_annual |> 
  filter(parameter %in% c("O1", "O2", "N1", "N2")) |> 
  pivot_wider(names_from = parameter) |> 
  mutate(
    SYN = N1 + N2,
    SYO = O1 + O2,
    .keep = "unused"
  ) |> 
  pivot_longer(c(SYN, SYO)) |> 
  mutate(value = value/1e6) |> 
  summarize(
    .by = c(lake, name),
    median = median(value, na.rm = TRUE),
    q25 = quantile(value, 0.25, na.rm = TRUE),
    q75 = quantile(value, 0.75, na.rm = TRUE)
  )


# Plot annual estimates of fry and smolt abundances
(ann_O_N_p <- post_annual |> 
  filter(parameter %in% c("O1", "O2", "N1", "N2")) |> 
  pivot_wider(names_from = parameter) |> 
  mutate(
    SYN = N1 + N2,
    SYO = O1 + O2,
    lake = factor(lake, levels = lake_names),
    .keep = "unused"
  ) |> 
  pivot_longer(c(SYN, SYO)) |> 
  mutate(value = value/1e6) |> 
  summarize(
    .by = c(lake, name, year),
    median = median(value, na.rm = TRUE),
    q10 = quantile(value, 0.10, na.rm = TRUE),
    q90 = quantile(value, 0.90, na.rm = TRUE)
  ) |> 
  mutate(name = if_else(name == "SYO", "Smolts", "Fry")) |> 
  ggplot(aes(x = year, y = median)) +
  facet_wrap(
    ~ lake,
    ncol = 1,
    strip.position = "right",
    scales = "free_y"
  ) +
  geom_pointrange(
    aes(
      ymin = q10,
      ymax = q90,
      colour = name
    ),
    position = position_dodge(width = 0.5)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(
    limits = c(0,NA),
    labels = scales::label_number(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_colour_manual(values = c("black", "grey")) +
  labs(
    y = "Abundance (millions)",
    x = "Outmigration year",
    colour = "Life\nstage"
  ) +
  theme(
    strip.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank(),
    panel.spacing.y = unit(1, "lines")
  )
)


ggsave(
  ann_O_N_p,
  filename = here(
    "3. outputs",
    "Plots",
    "Fry_smolt_abundance_Bayesian_estimates.png"
  ),
  width = 8,
  height = 5,
  units = "in",
  dpi = "print"
)


# Median age compositions of fry from the three lakes
post_annual |> 
  filter(parameter %in% c("N1", "N2")) |> 
  pivot_wider(names_from = parameter) |> 
  mutate(
    pN1 = N1/(N1 + N2),
    pN2 = 1 - pN1,
    lake = factor(lake, levels = lake_names),
    .keep = "unused"
  ) |> 
  pivot_longer(c(pN1, pN2)) |> 
  summarize(
    .by = c(lake, name),
    median = median(value, na.rm = TRUE),
    q25 = quantile(value, 0.25, na.rm = TRUE),
    q75 = quantile(value, 0.75, na.rm = TRUE)
  )


# Table of key parameters governing process dynamics
posterior_df |> 
  filter(
    parameter %in% c(
      "mu_theta", "mu_M", "sigma_theta", "sigma_M", "rho_theta_M",
      "mu_N1", "phi", "sigma_proc", 
      "obs_error_scale",
      "N2_init"
    )
  ) |> 
  # Back-transform values to natural/interpretation scale
  mutate(
    value = case_when(
      parameter %in% c("mu_theta", "mu_M") ~ plogis(value),
      parameter == "mu_N1" ~ exp(value),
      .default = value
    )
  ) |> 
  summarize(
    .by = c(lake, parameter),
    median = median(value, na.rm = TRUE),
    q10 = quantile(value, 0.1, na.rm = TRUE),
    q90 = quantile(value, 0.9, na.rm = TRUE)
  ) |> 
  mutate(
    across(median:q90, \(x) round(x, 3)),
    value = paste0(median, " (", q10, "-", q90, ")")
  ) |> 
  pivot_wider(
    id_cols = parameter,
    names_from = lake
  ) |> 
  write.table(
    "clipboard",
    row.names = FALSE
  )

# Smolt weight versus abundance relationships -----------------------------


# Lookup for enhanced years
enh_lu <- here(
  "3. outputs",
  "Total return time series",
  "Barkley Sockeye total annual returns by CU_collated.xlsx"
) |> 
  read_xlsx() |> 
  mutate(
    across(matches("hat.*release"), \(x) if_else(is.na(x), 0, x)),
    enhanced = if_else(fertilized == 1 | if_any(matches("hat.*release")) > 0, 1, 0),
    lake = factor(
      CU,
      levels = c("GCL", "SPR", "HUC"),
      labels = lake_names
    )
  ) |> 
  select(year, lake, enhanced)


# Create a dataframe of approximated annual average smolt weights from the 
# observed sample data 
BY_smolt_w <- smolt_sizes_csr |> 
  filter(
    fnlage %in% c(1, 2),
    !is.na(fresh_std_weight)
  ) |> 
  # Convert year to brood year
  mutate(year = if_else(fnlage == 1, year-2, year-3)) |> 
  rename("w" = fresh_std_weight) |> 
  summarize(
    .by = c(lake, year),
    median_w = Hmisc::wtd.quantile(w, weights = csr, probs = 0.5, na.rm = TRUE),
    l_80_w = Hmisc::wtd.quantile(w, weights = csr, probs = 0.1, na.rm = TRUE),
    u_80_w = Hmisc::wtd.quantile(w, weights = csr, probs = 0.9, na.rm = TRUE)
  ) |> 
  complete(lake, year)


# Join observed weight data to posterior smolt abundance estimates
obs_wt_posterior_BYO <- posterior_df |> 
  filter(parameter == "BYO") |>
  mutate(
    lake = factor(lake, levels = lake_names),
    BYO = value/1e6
  ) |> 
  summarize(
    .by = c(lake, year),
    across(
      BYO,
      .fns = list(
        median = ~median(.x, na.rm = TRUE),
        l_80 = ~quantile(.x, 0.1, na.rm = TRUE),
        u_80 = ~quantile(.x, 0.9, na.rm = TRUE)
      ),
      .names = "{.fn}_{.col}"
    ) 
  ) |> 
  left_join(enh_lu) |> 
  left_join(BY_smolt_w) |> 
  mutate(method = "Average smolt weight (g)\nfrom sample data")


# Calculate abundance and weights from model posteriors
posterior_BYO_wt <- posterior_df |> 
  filter(parameter %in% c("BYO", "BYB")) |> 
  pivot_wider(
    id_cols = c(lake, draw, year),
    names_from = parameter,
    values_from = value
  ) |> 
  mutate(
    w = BYB/BYO,
    lake = factor(lake, levels = lake_names),
    BYO = BYO/1e6
  ) |> 
  summarize(
    .by = c(lake, year),
    across(
      c(BYO, w),
      .fns = list(
        median = ~median(.x, na.rm = TRUE),
        l_80 = ~quantile(.x, 0.1, na.rm = TRUE),
        u_80 = ~quantile(.x, 0.9, na.rm = TRUE)
      ),
      .names = "{.fn}_{.col}"
    ) 
  ) |> 
  left_join(enh_lu) |> 
  left_join(select(lakes_data, lake, year, is_observed_WO1)) |> 
  #filter(is_observed_WO1) |> # Only include years where weights were measured?
  mutate(method = "Average smolt weight (g) from\nlife cycle model posterior")


(w_abun_p <- posterior_BYO_wt |> 
    bind_rows(obs_wt_posterior_BYO) |> 
    ggplot(aes(x = median_BYO, y = median_w)) +
    facet_grid(
      method~lake,
      scales = "free_x",
      switch = "y"
    ) +
    geom_linerange(
      aes(
        xmin = l_80_BYO,
        xmax = u_80_BYO
      ),
      colour = "grey50",
      alpha = 0.5
    ) +
    geom_pointrange(
      aes(
        ymin = l_80_w,
        ymax = u_80_w,
        fill = factor(enhanced),
      ),
      alpha = 0.5,
      colour = "grey50",
      shape = 21
    ) +
    scale_fill_manual(values = c("red", "blue")) +
    scale_x_continuous(
      limits = c(0, NA),
      labels = scales::label_number(accuracy = 1),
      breaks = scales::pretty_breaks(n = 4),
      expand = expansion(mult = c(0, 0.05))
    ) +
    scale_y_continuous(
      limits = c(0, 8),
      labels = scales::label_number(),
      breaks = scales::pretty_breaks(n = 4),
      expand = expansion(mult = c(0, 0.05)),
      oob = scales::oob_keep
    ) +
    labs(
      x = "Brood year smolt abundance (millions)\nfrom life cycle model posterior",
      y = NULL,
      fill = "Enhancement\nstate"
    ) +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(1, "lines"),
      strip.placement = "outside",
      panel.spacing.y = unit(1, "lines")
    )
)


# Save 
ggsave(
  w_abun_p,
  filename = here(
    "3. outputs",
    "Plots",
    "Smolt_weight_versus_abundance.png"
  ),
  width = 8,
  height = 5.2,
  units = "in",
  dpi = "print"
)
