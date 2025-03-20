# Packages ----------------------------------------------------------------


pkgs <- c(
  "tidyverse", "rstan", "here", "bayesplot", "tidybayes", 
  "geomtextpath", "writexl", "readxl", "ggridges"
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



# Load and infill smolt size data -----------------------------------------


# Raw data, with annual mean lengths and weights infilled for missing years
smolt_sizes <- here(
  "1. data",
  "smolt size data.xlsx"
) |> 
  read_xlsx(sheet = "smolt_morpho_age") |> 
  filter(smolt_year >= 1977) |> 
  pivot_longer(
    matches("age\\d"),
    names_pattern = "(age\\d)_(.*)",
    names_to = c("age", "measure")
  ) |> 
  pivot_wider(
    names_from = measure,
    values_from = value
  ) |> 
  pivot_longer(
    matches("^(len|w)"),
    names_sep = "_",
    names_to = c("measure", "statistic")
  ) |> 
  pivot_wider(
    names_from = statistic,
    values_from = value
  ) |> 
  mutate(cv = sd/mean)


# Does using average CV make sense for SD infilling?
smolt_sizes |> 
  ggplot(aes(x = smolt_year, y = cv)) +
  facet_grid(measure ~ lake) +
  geom_point(aes(colour = age)) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1),
    expand = c(0, 0)
  )
# Data look a bit noisy for weights


# Does CV scale with mean?
smolt_sizes |> 
  ggplot(aes(x = mean, y = cv, colour = age)) +
  facet_grid(lake ~ measure, scales = "free_x") +
  geom_point() +
  geom_smooth()
# Not really

# It looks like it's probably fine to use average CV to infill missing SD values


# Missing age-2 Hucuktlis weights
# see if there is any possibility to infill based on age1 weights or lagged
# age1 weights
smolt_sizes |> 
  filter(lake == "huc") |> 
  pivot_wider(
    id_cols = c(lake, smolt_year, measure),
    names_from = age,
    values_from = mean
  ) |> 
  arrange(lake, measure, smolt_year) |> 
  mutate(
    .by = c(lake, measure),
    age1_lag = lag(age1, 1L)
  ) |> 
  pivot_longer(
    matches("age1"),
    names_to = "predictor",
    values_to = "x"
  ) |> 
  ggplot(aes(x = log(x), y = log(age2))) +
  facet_wrap(
    predictor ~ measure, 
    scales = "free",
    labeller = labeller(.multi_line = FALSE)
  ) +
  geom_smooth(method = "lm") +
  geom_point()
# Age1 sizes seem reasonable as a predictor


# Use lm to predict age2 values for Hucuktlis
huc_infill <- smolt_sizes |> 
  filter(lake == "huc") |> 
  pivot_wider(
    id_cols = c(lake, smolt_year, measure),
    names_from = age,
    values_from = mean
  ) |> 
  nest(.by = measure) |> 
  rowwise() |> 
  mutate(
    model = list(lm(age2 ~ age1, data = data)),
    pred = list(predict(model, data))
  ) |> 
  unnest(cols = c(data, pred)) |> 
  mutate(
    infill = if_else(is.na(age2), pred, age2),
    age = "age2"
  ) |> 
  select(measure, lake, smolt_year, age, infill)


# Add infilled age2 values for Hucuktlis and infill SD and N
smolt_sizes_infilled <- smolt_sizes |> 
  left_join(huc_infill) |> 
  mutate(mean = if_else(is.na(mean), infill, mean)) |> 
  mutate(
    .by = c(lake, measure, age),
    mean_cv = mean(cv, na.rm = TRUE),
    median_n = round(median(n, na.rm = TRUE), 0),
    # Assume SD is similar to historic average CV
    sd = if_else(is.na(sd), mean*mean_cv, sd),
    # Assume N is equivalent to the historic median N
    n = if_else(is.na(n), median_n, n)
  ) |> 
  select(-infill, -mean_cv, -median_n)


# Dataframe with just summary statistics for weight
smolt_wt <- smolt_sizes_infilled |> 
  filter(measure == "w") |> 
  pivot_wider(
    id_cols = c(smolt_year, lake),
    names_from = age,
    values_from = c(mean, sd, n),
    names_glue = "{age}_{.value}"
  ) |> 
  rename_with(\(x) str_replace(x, "age", "WO")) |> 
  mutate(lake = factor(lake, labels = c("Great Central", "Hucuktlis", "Sproat")))


# Load and prepare raw smolt abundance data  ----------------------------------


# Load the smolt outmigration data
smolt_data <- here(
  "3. outputs",
  "Stock-recruit data",
  "Inferred_annual_smolt_age_composition.csv"
) |> 
  read.csv() |> 
  rename_with(\(x) str_replace_all(x, "\\.", "x")) |> 
  filter(year < 2017)  |> # post-2017 won't have data
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


# Look at distribution of observed smolting rates
smolt_data |> 
  ggplot(aes(x = logit(smolting_rate))) +
  facet_wrap(~lake) +
  geom_density(
    fill = "grey",
    bw = 1
  )
# Median values best? Appears to be some bimodality


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


# Load smolt size data
smolt_size <- here(
  "1. data",
  "smolt size data.xlsx"
) |> 
  read_xlsx(sheet = "smolt_morphometrics") |> 
  pivot_longer(
    cols = !smolt_year,
    names_pattern = "(.{3})_(.*)",
    names_to = c("lake", "metric")
  ) |> 
  pivot_wider(names_from = metric) |> 
  mutate(
    lake = case_when(
      lake == "gcl" ~ "Great Central",
      lake == "huc" ~ "Hucuktlis",
      lake == "spr" ~ "Sproat"
    ),
    K = 100*(w_g/len_mm^3)
  )


# fry size data
fry_size <- here(
  "1. data",
  "smolt size data.xlsx"
) |> 
  read_xlsx(sheet = "fry_morphometrics") |> 
  pivot_longer(
    cols = !smolt_year,
    names_pattern = "(.{3})_(.*)",
    names_to = c("lake", "metric")
  ) |> 
  pivot_wider(names_from = metric) |> 
  rename_with(\(x) paste0("fry_", x), .cols = !c(lake, smolt_year)) |> 
  mutate(
    lake = case_when(
      lake == "gcl" ~ "Great Central",
      lake == "huc" ~ "Hucuktlis",
      lake == "spr" ~ "Sproat"
    )
  )


# Plot smolt size versus smolting rate
smolt_size |> 
  left_join(fry_size) |> 
  mutate(fw_growth_idx = len_mm - fry_len_mm) |> 
  left_join(
    smolt_data, 
    by = c("smolt_year" = "year", "lake")
  ) |>
  pivot_longer(c(w_g, len_mm, K, fw_growth_idx)) |> 
  ggplot(aes(x = value, y = smolting_rate)) +
  facet_grid(
    lake~name, 
    scales = "free",
    switch = "x"
    ) +
  geom_point() +
  geom_smooth() +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1),
    expand = c(0, 0),
    oob = scales::oob_keep
  ) +
  theme(
    strip.placement = "outside",
    strip.background.x = element_blank(),
    axis.title.x = element_blank()
  )
# No strong patterns evident here. Both the dependent and 
# independent variable are subject to considerable uncertainty, but 
# if there was a strong link between growth and smolting rate, a 
# clearer pattern would likely be discernable.
# Possible conclusion: environmental conditions in the lakes do not
# have a strong influence on smolting rate?



# Load the annual ATS estimates
ats_data <- here(
  "3. outputs",
  "Stock-recruit data",
  "GAM-estimated_pre-smolt_time_series.csv"
) |> 
  read.csv() |> 
  rename("year" = smolt_year) |> 
  rename_with(\(x) paste0("ats_", x), .cols = c(est, lwr, upr, se, cv)) |> 
  # Remove NAs at start and end of time series
  filter(!is.na(ats_est)) |> 
  group_by(lake) |> 
  complete(year = full_seq(min(year):max(year), 1)) |> 
  ungroup()


# Are there any years with missing estimates?
ats_data |> 
  ggplot(aes(x = year, y = ats_est)) +
  facet_wrap(
    ~lake,
    ncol = 1,
    scales = "free_y"
  ) +
  geom_pointrange(aes(ymin = ats_lwr, ymax = ats_upr))


# Join ATS and smolt data
lakes_data <- ats_data |> 
  left_join(
    select(
      .data = smolt_data, 
      lake, year, contains("count"), smolting_rate, prop_age2
    )
  ) |> 
  # Add infilled smolt weight data
  left_join(
    smolt_wt,
    by = c("lake", "year" = "smolt_year")
  ) |> 
  mutate(
    is_observed_count = !if_any(matches("count"), is.na),
    is_observed_ats = !if_any(matches("ats"), is.na),
    count_total = count_age1 + count_age2,
    # Convert NAs to 0s
    across(matches("ats|count"), \(x) if_else(is.na(x), 0, x))
  ) |> 
  arrange(year) |> 
  # Post-2017 won't have data
  filter(
    year < 2017,
    !(lake == "Hucuktlis" & year == 2016) # No weight data available for Hucuktlis in 2016
  )



# Format data for Stan ----------------------------------------------------


# Function to create observed Stan data for each lake
# As a general note: priors should reflect what is *possible*
# not what we think is likely.
make_stan_data <- function(
    lake_name,
    # User-defined second year mortality and theta priors 
    # (with sensible default values)
    mu_M_prior = logit(0.5),
    sigma_M_prior = 0.5,
    mu_theta_prior = logit(0.8),
    sigma_theta_prior = 1,
    sigma_theta_sd_prior = 1,
    sigma_M_sd_prior = 2
) {
  
  # Filter data for a single lake
  lake_data <- lakes_data %>%
    filter(lake == lake_name)

  # Identify missing rows in observations
  is_observed_count <- lake_data$is_observed_count
  is_observed_ats <- lake_data$is_observed_ats
  
  # Vector of age1 counts
  A1_obs <- lake_data$count_age1
  A_total <- lake_data$count_total
  
  # Vector of observed total lake abundance
  N_obs <- round(lake_data$ats_est, 0)
  
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
  
  # Age-1 fry proportion prior
  pN1_prior <- lake_data |> 
    mutate(pN1 = logit(1-prop_age2)) |> 
    summarize(
      mu_pN1 = mean(pN1, na.rm = TRUE),
      sigma_pN1 = sd(pN1, na.rm = TRUE)
    )
  
  # Alternative specification of smolting rate prior
  # PT suggestion after reviewing the data-driven approach above and
  # observing that it seems to draw the model away from what we assume
  # is going on for GCL
  theta_prior_alt <- data.frame(
    mu_theta = mu_theta_prior,
    sigma_theta = sigma_theta_prior
  )
  
  # Initial age2 fry population prior
  N2_init_prior <- lake_data |> 
    filter(year == min(year)) |> 
    mutate(
      N2_init = prop_age2 * ats_est,
      N2_init_sigma = prop_age2 * ats_se,
      # Need to transform N2_init mean and SD for lognormal distribution:
      # https://msalganik.wordpress.com/2017/01/21/making-sense-of-the-rlnorm-function-in-r/
      N2_init_log_mean = log(N2_init^2 / sqrt(N2_init_sigma^2 + N2_init^2)),
      N2_init_log_sigma = sqrt(log(1 + (N2_init_sigma^2 / N2_init^2))),
    )
  
  # Total lake fry abundance prior
  N1_prior <- lake_data |> 
    mutate(N1 = (1-prop_age2)*ats_est) |> 
    summarize(
      mean_N1 = mean(N1, na.rm = TRUE),
      sd_N1 = sd(N1, na.rm = TRUE)
    ) |> 
    mutate(
      mu_N1 = log(mean_N1^2 / sqrt(sd_N1^2 + mean_N1^2)),
      sigma_N1 = sqrt(log(1 + (sd_N1^2 / mean_N1^2)))
    )
  
  # Observation error on lake ATS estimates
  obs_error_prior <- lake_data |> 
    summarize(
      ats_se = mean(ats_se),
      ats_est = mean(ats_est)
    ) |> 
    mutate(ats_se = ats_se/5) |> # Reduce the SE estimates
    mutate(obs_error_prior = sqrt(log(1 + (ats_se^2 / ats_est^2)))) |> 
    pull(obs_error_prior)
  
  # Stan data list
  stan_data <- list(
    Y = nrow(lake_data),
    A1_obs = A1_obs,
    A_total = A_total,
    N_obs = N_obs,
    WO1_mean = lake_data$WO1_mean,
    WO1_sd = lake_data$WO1_sd,
    WO1_n = lake_data$WO1_n,
    WO2_mean = lake_data$WO2_mean,
    WO2_sd = lake_data$WO2_sd,
    WO2_n = lake_data$WO2_n,
    is_observed_count = as.integer(is_observed_count),
    is_observed_ats = as.integer(is_observed_ats),
    obs_error_prior = obs_error_prior,
    #mu_theta_prior = theta_prior$mu_theta,
    #sigma_theta_prior = theta_prior$sigma_theta,
    mu_theta_prior = theta_prior_alt$mu_theta, # alternative user specification
    sigma_theta_prior = theta_prior_alt$sigma_theta, # alternative user specification
    mu_M_prior = mu_M_prior,
    sigma_M_prior = sigma_M_prior,
    sigma_theta_sd_prior = sigma_theta_sd_prior,
    sigma_M_sd_prior = sigma_M_sd_prior,
    mu_N2_init_prior = N2_init_prior$N2_init_log_mean,
    sigma_N2_init_prior = N2_init_prior$N2_init_log_sigma,
    #mu_pN1_prior = pN1_prior$mu_pN1,
    #mu_pN1_prior = logit(0.8),
    #sigma_pN1_prior = 2,
    #sigma_pN1_prior = pN1_prior$sigma_pN1 # might try wider prior?
    mu_N1_prior = N1_prior$mu_N1,
    sigma_N1_prior = N1_prior$sigma_N1
  )
  
  return(stan_data)
}


# Create nested list of Stan input data for each lake
lakes_stan_data <- unique(lakes_data$lake) |> 
  set_names() |> 
  map2(
    # Second winter mortality rate priors for GCL, HUC, SPR, respectively
    .y = logit(c(0.25, 0.65, 0.5)),
    \(x, y) make_stan_data(
      lake_name = x, 
      mu_M_prior = y # Try with all default values first
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
fit_stan_model <- function(stan_data, index) {
  stan(
    file = here(
      "2. code",
      "Stock-recruit data preparation",
      "Stan",
      "smolt_production_model.stan"
    ),
    model_name = index,
    data = stan_data,
    iter = 4000,
    chains = 4,
    control = list(
      adapt_delta = 0.99
      ,max_treedepth = 15
    )
  )
}


# Try fitting the model with just one lake first
# fit_stan_model(lakes_stan_data$`Great Central`, "Great Central")
# working, as of 6 March

# Fit models for each lake (toggle to TRUE to run)
if(
  FALSE
  #TRUE
) {
  lakes_stan_fits <- lakes_stan_data |> 
    imap(fit_stan_model)
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
            paste("Bayesian state-space", idx, "smolts.rds")
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
        pattern = paste("Bayesian state-space", x, "smolts.rds"),
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
  map(\(x) pairs(x, pars = c("mu_theta", "mu_M")))


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
  map(\(x) discard(.x = x, .p = \(y) is.matrix(y))) |> 
  map(\(x) map(x, \(y) as_tibble(y, rownames = "draw"))) |> 
  map(\(x) list_rbind(x, names_to = "parameter")) |> 
  list_rbind(names_to = "lake") 


# Posterior values that are year-specific
post_annual <- lakes_stan_fits |> 
  map(extract) |> 
  # Keep only parameters with year-specific estimates
  map(\(x) keep(.x = x, .p = \(y) is.matrix(y))) |> 
  map(\(x) map(x, \(y) as_tibble(y, .name_repair = NULL, rownames = "draw"))) |> 
  map(\(x) list_rbind(x, names_to = "parameter")) |> 
  list_rbind(names_to = "lake") |> 
  pivot_longer(
    cols = !c(lake, parameter, draw),
    names_to = "year",
    values_to = "value",
    names_transform = \(x) str_extract(x, "\\d+")
  ) |> 
  left_join(min_yrs) |> 
  mutate(year = as.numeric(year) - 1 + min_yr) |> 
  select(-min_yr)


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


# Bring both posterior dataframes together
posterior_df <- bind_rows(
  post_annual, 
  post_interannual,
  post_by_ttls
) |> 
  # Make parameter names more informative
  mutate(
    parameter_long = case_when(
      parameter == "O1" ~ "age1 smolts",
      parameter == "O2" ~ "age2 smolts",
      parameter == "BO1" ~ "age1 smolt biomass",
      parameter == "BO2" ~ "age2 smolt biomass",
      parameter == "w1" ~ "age1 smolt mean weight",
      parameter == "w2" ~ "age2 smolt mean weight",
      parameter == "p1" ~ "age1 smolt proportion",
      parameter == "p2" ~ "age2 smolt proportion",
      parameter == "N1" ~ "age1 fry",
      parameter == "N2" ~ "age2 fry",
      parameter == "N_lake" ~ "lake total fry abundance",
      parameter == "theta" ~ "smolting rate",
      parameter == "M" ~ "age1 to age2 mortality",
      parameter == "BYB" ~ "brood year smolt biomass production",
      parameter == "BYO" ~ "brood year smolt production",
      .default = "NO_ENTRY"
    )
  )


# Plot posterior versus prior distributions for mu_theta and mu_M
posterior_df |> 
  filter(parameter %in% c("mu_M", "mu_theta")) %>% 
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
  

# Plot posterior distributions for sigma_theta & sigma_M
posterior_df |> 
  filter(parameter %in% c("sigma_M", "sigma_theta")) |> 
  ggplot(aes(x = value)) +
  facet_wrap(
    parameter ~ lake,
    scales = "free",
    labeller = labeller(.multi_line = FALSE)
  ) +
  geom_textdensity(
    label = "posterior",
    hjust = 0.8
  ) +
  geom_textvline(
    data = lakes_stan_data |> 
      map(\(x) keep_at(x, paste0(c("sigma_M", "sigma_theta"), "_sd_prior"))) |> 
      map(as_tibble) |> 
      list_rbind(names_to = "lake") |> 
      pivot_longer(
        !lake,
        names_to = "label",
        values_to = "value"
      ) |> 
      mutate(parameter = str_remove(label, "_sd_prior")),
    aes(
      xintercept = value, 
      label = label
    ),
    colour = "red",
    hjust = 0.7
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "value")


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
    w1 = WO1_mean,
    w2 = WO2_mean,
    N_lake = ats_est,
    BO1 = O1*WO1_mean,
    BO2 = O2*WO2_mean,
    BYB = BO1 + lead(BO2),
    BYO = O1 + lead(O2),
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
        ggtitle(paste("Predicted versus estimated", idx))
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


# Make the same plots but with ridges
(out_plots <- posterior_df |> 
    filter(
      !is.na(value),
      parameter %in% unique(obs_params$parameter)
    ) |> 
    # Add observed values
    left_join(obs_params) %>%
    split(.$parameter_long) |> 
    imap(
      \(x, idx) x |> 
        ggplot(aes(y = year, x = value, group = year)) +
        facet_wrap(
          ~lake,
          scales = "free_x",
          nrow = 1
        ) +
        stat_summary(
          aes(x = obs),
          geom = "point",
          fun = \(y) median(y, na.rm = TRUE),
          colour = "red",
          shape = "|",
          vjust = 1.5
        ) +
        geom_density_ridges(
          rel_min_height = 0.001,
          alpha = 0.5
        ) +
        scale_y_reverse() +
        labs(x = idx) +
        ggtitle(paste("Predicted versus estimated", idx)) +
        theme_ridges() +
        theme(strip.background = element_blank())
    )
)


# Export posterior values of interest (with imputed years) -------------------


# Summarize posterior for export
annual_estimates <- posterior_df |> 
  filter(
    parameter %in% str_subset(annual_params, "p\\d", negate = TRUE),
    !is.na(value)
  ) |> 
  mutate(
    family = case_when(
      str_detect(parameter, "N|O|BO|BY") ~ "lognormal",
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
metadata <- tibble(
  parameter = unique(annual_estimates$parameter)
) |> 
  mutate(
    definition = case_when(
      parameter == "O1" ~ "Posterior estimates for age1 outmigrating smolts",
      parameter == "O2" ~ "Posterior estimates for age2 outmigrating smolts",
      parameter == "w1" ~ "Posterior estimates for age1 smolt average weight (g)",
      parameter == "w2" ~ "Posterior estimates for age2 smolts average weight (g)",
      parameter == "BO1" ~ "Posterior estimates for age1 smolt biomass (g)",
      parameter == "BO2" ~ "Posterior estimates for age2 smolt biomass (g)",
      parameter == "BYO" ~ "Brood year sum of posterior estimates for O1 & O2 (see above)",
      parameter == "BYB" ~ "Brood year sum of posterior estimates for BO1 & BO2 (see above)",
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
  
  