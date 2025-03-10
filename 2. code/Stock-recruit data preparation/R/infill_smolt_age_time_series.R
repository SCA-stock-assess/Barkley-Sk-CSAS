# Packages ----------------------------------------------------------------


pkgs <- c(
  "tidyverse", "rstan", "here", "bayesplot", "tidybayes", 
  "geomtextpath", "writexl"
)
#install.packages(pkgs)


library(here)
library(tidyverse); theme_set(theme_bw())
library(rstan)
library(bayesplot)
library(tidybayes)
library(geomtextpath)
library(writexl)


# Logit transformation
logit <- function(p) {log(p/(1-p))}


# Load and prepare raw data for Stan --------------------------------------


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
  mutate(
    is_observed_count = !if_any(matches("count"), is.na),
    is_observed_ats = !if_any(matches("ats"), is.na),
    count_total = count_age1 + count_age2,
    # Convert NAs to 0s
    across(matches("ats|count"), \(x) if_else(is.na(x), 0, x))
  ) |> 
  arrange(year) |> 
  # Post-2017 won't have data
  filter(year < 2017)


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
    sigma_M_sd_prior = 1
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
    mutate(ats_se = ats_se/2) |> # Reduce the SE estimates
    mutate(obs_error_prior = sqrt(log(1 + (ats_se^2 / ats_est^2)))) |> 
    pull(obs_error_prior)
  
  # Stan data list
  stan_data <- list(
    Y = nrow(lake_data),
    A1_obs = A1_obs,
    A_total = A_total,
    N_obs = N_obs,
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

# Fit models for each lake
lakes_stan_fits <- lakes_stan_data |> 
  imap(fit_stan_model)


# Save summary data from the fitted models
models_summary <- lakes_stan_fits |> 
  map(rstan::summary) |> 
  map(\(x) x$summary) |> 
  map(\(x) as_tibble(x, rownames = "parameter")) |> 
  list_rbind(names_to = "stock") |> 
  nest(.by = stock, .key = "model_summary")


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


# Posterior values
post <- lakes_stan_fits |> 
  map(extract) |> 
  # Discard all parameters with year-specific estimates
  map(\(x) discard(.x = x, .p = \(y) is.matrix(y))) |> 
  map(\(x) imap(x, as_tibble_col)) |> 
  map(\(x) list_cbind(unname(x))) |> 
  list_rbind(names_to = "lake")


# Plot posterior versus prior distributions for mu_theta and mu_M
post |> 
  select(lake, mu_M, mu_theta) |> 
  pivot_longer(!lake) %>% 
  split(.$name) |> 
  imap(
    \(x, idx) x |> 
      ggplot(
        aes(
          #x = value
          x = binomial()$linkinv(value)
        )
      ) +
      facet_wrap(
        ~lake,
        scales = "free"
      ) +
      geom_textdensity(
        label = "posterior",
        hjust = 0.3
      ) +
      geom_textdensity(
        data = some_priors[[idx]],
        colour = "red",
        label = "prior",
        hjust = 0.7
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      labs(x = idx)
  )


# Plot posterior distributions for sigma_theta & sigma_M
post |> 
  select(lake, sigma_M, sigma_theta) |> 
  pivot_longer(!lake) %>% 
  split(.$name) |> 
  imap(
    function(x, idx) {
      prior <- lakes_stan_data |> 
        map(\(x) keep_at(x, \(y) str_detect(y, idx))) |> 
        map(as_tibble) |> 
        list_rbind()
      
      x |> 
      ggplot(aes(x = value)) +
      facet_wrap(
        ~lake,
        scales = "free"
      ) +
      geom_textdensity(
        label = "posterior",
        hjust = 0.8
      ) +
      geom_textvline(
        data = prior,
        aes(xintercept = .data[[paste0(idx, "_sd_prior")]]),
        colour = "red",
        label = paste0(idx, "_sd_prior"),
        hjust = 0.7
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      labs(x = idx)
    }
  )


# Predicted versus observed Age1 out-migrants
lakes_stan_fits |> 
  imap(
    function(x, idx) {
      O1 <- lakes_data |> 
        filter(lake == idx) |> 
        mutate(age1s = (1-prop_age2)*smolting_rate*ats_est) |> 
        pull(age1s)
      Y <- lakes_stan_data[[idx]]$Y
      
      spread_draws(x, O1[year]) %>% 
        ggplot(aes(x = year, y = O1)) +
        stat_pointinterval(.width = c(0.5, 0.8, 0.95))+
        geom_point(data = data.frame(year = 1:Y, O1 = O1), color = "red") +
        scale_y_continuous(labels = scales::label_number()) +
        ggtitle(paste(idx, "predicted versus estimated age1 out-migrants"))
    }
  )


# Predicted versus observed Age2 out-migrants
lakes_stan_fits |> 
  imap(
    function(x, idx) {
      O2 <- lakes_data |> 
        filter(lake == idx) |> 
        mutate(age2s = prop_age2*ats_est) |> 
        pull(age2s)
      Y <- lakes_stan_data[[idx]]$Y
      
      spread_draws(x, O2[year]) %>% 
        ggplot(aes(x = year, y = O2)) +
        stat_pointinterval(.width = c(0.5, 0.8, 0.95))+
        geom_point(data = data.frame(year = 1:Y, O2 = O2), color = "red") +
        scale_y_continuous(labels = scales::label_number()) +
        ggtitle(paste(idx, "predicted versus estimated age2 out-migrants"))
    }
  )


# Predicted versus observed age-1 smolt proportion
lakes_stan_fits |> 
  imap(
    function(x, idx) {
      p1 <- lakes_data |> 
        filter(lake == idx) |> 
        mutate(prop_age1s = (1-prop_age2)*smolting_rate) |> 
        pull(prop_age1s)
      Y <- lakes_stan_data[[idx]]$Y
      
      spread_draws(x, p1[year]) %>% 
        ggplot(aes(x = year, y = p1)) +
        stat_pointinterval(.width = c(0.5, 0.8, 0.95))+
        geom_point(data = data.frame(year = 1:Y, p1 = p1), color = "red") +
        scale_y_continuous(
          expand = c(0, 0),
          limits = c(0, 1),
          labels = scales::percent
        ) +
        ggtitle(paste(idx, "predicted versus estimated age1 proportion"))
    }
  )


# Predicted Age2 proportions
lakes_stan_fits |> 
  imap(
    function(x, idx) {
      p2 <- lakes_data |> 
        filter(lake == idx) |> 
        pull(prop_age2)
      Y <- lakes_stan_data[[idx]]$Y
      
      spread_draws(x, p2[year]) %>% 
        ggplot(aes(x = year, y = p2)) +
        stat_pointinterval(.width = c(0.5, 0.8, 0.95))+
        geom_point(data = data.frame(year = 1:Y, p2 = p2), color = "red") +
        scale_y_continuous(
          expand = c(0, 0),
          limits = c(0, 1),
          labels = scales::percent
        ) +
        ggtitle(paste(idx, "predicted versus estimated age2 proportion"))
    }
  )


# Predicted versus estimated smolting rate
lakes_stan_fits |> 
  imap(
    function(x, idx) {
      theta <- lakes_data |> 
        filter(lake == idx) |> 
        pull(smolting_rate)
      Y <- lakes_stan_data[[idx]]$Y
      
      spread_draws(x, theta[year]) %>% 
        ggplot(aes(x = year, y = theta)) +
        stat_pointinterval(.width = c(0.5, 0.8, 0.95))+
        geom_point(data = data.frame(year = 1:Y, theta = theta), color = "red") +
        scale_y_continuous(
          expand = c(0, 0),
          limits = c(0, 1),
          labels = scales::percent
        ) +
        ggtitle(paste(idx, "predicted versus estimated smolting rate"))
    }
  )


# Predicted versus observed lake abundance
lakes_stan_fits |> 
  imap(
    function(x, idx) {
      N_lake <- lakes_data |> 
        filter(lake == idx) |> 
        pull(ats_est)
      Y <- lakes_stan_data[[idx]]$Y
      
      spread_draws(x, N_lake[year]) %>% 
        ggplot(aes(x = year, y = N_lake)) +
        stat_pointinterval(.width = c(0.5, 0.8, 0.95))+
        geom_point(data = data.frame(year = 1:Y, N_lake = N_lake), color = "red") +
        ggtitle(paste(idx, "predicted versus estimated lake abundance"))
    }
  )


# Predicted mortality rate
lakes_stan_fits |> 
  imap(
    \(x, idx) spread_draws(x, M[year]) %>% 
        ggplot(aes(x = year, y = M)) +
        stat_pointinterval(.width = c(0.5, 0.8, 0.95))+
      scale_y_continuous(
        expand = c(0, 0),
        limits = c(0, 1),
        labels = scales::percent
      ) +
        ggtitle(paste(idx, "predicted mortality rate"))
  )



  

# Extract posterior (with imputed) values and export ----------------------


# Dataframe with first year in time series for each lake
min_yrs <- lakes_data |> 
  summarize(
    .by = lake,
    min_yr = min(year)
  )


# Big dataframe with all draws for all lakes
posterior_df <- lakes_stan_fits |> 
  map(extract) |> 
  map(
    \(x) keep(
      .x = x, 
      .p = \(y) is.matrix(y))
  ) |> 
  map(\(x) map(x, \(y) as_tibble(y, .name_repair = NULL))) |> 
  map(\(x) list_rbind(x, names_to = "parameter")) |> 
  list_rbind(names_to = "lake") |> 
  pivot_longer(
    cols = !c(lake, parameter),
    names_to = "year",
    values_to = "value",
    names_transform = \(x) str_extract(x, "\\d+")
  ) |> 
  left_join(min_yrs) |> 
  mutate(year = as.numeric(year) - 1 + min_yr) |> 
  select(-min_yr) |> 
  filter(parameter %in% c("O1", "O2", "N1", "N2", "N_lake", "theta", "M"))
  

# Summarize posterior for export
annual_estimates <- posterior_df |> 
  filter(!is.na(value)) |> 
  summarize(
    .by = c(lake, parameter, year),
    quantiles = list(quantile(value, probs = c(0.025, 0.1, 0.5, 0.9, 0.975)))
    # How about an SD-esque statistic? 
  ) |> 
  unnest_wider(quantiles)


# Metadata sheet
metadata <- tibble(
  parameter = unique(annual_estimates$parameter)
) |> 
  mutate(
    definition = case_when(
      parameter == "O1" ~ "Posterior estimates for age1 outmigrating smolts",
      parameter == "O2" ~ "Posterior estimates for age2 outmigrating smolts",
      parameter == "N1" ~ "Posterior estimates for age1 fry abundance in the lake on 1 March",
      parameter == "N2" ~ "Posterior estimates for age2 fry abundance in the lake on 1 March",
      parameter == "N_lake" ~ "Posterior estimates for total fry abundance in the lake on 1 March",
      parameter == "theta" ~ "Posterior estimates for age1 smolting rate",
      parameter == "M" ~ "Posterior estimates for second year fry mortality (age1 to age2)",
      .default = "MISSING"
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
  
  