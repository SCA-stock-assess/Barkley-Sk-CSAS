# Packages ----------------------------------------------------------------


pkgs <- c("tidyverse", "rstan", "here", "bayesplot", "tidybayes", "geomtextpath")
#install.packages(pkgs)


library(here)
library(tidyverse); theme_set(theme_bw())
library(rstan)
library(bayesplot)
library(tidybayes)
library(geomtextpath)


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
make_stan_data <- function(
    lake_name,
    # User-defined second year mortality and theta priors 
    # (with sensible default values)
    mu_M_prior = -1.5,
    sigma_M_prior = 0.5,
    mu_theta_prior = 4,
    sigma_theta_prior = 2
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
  
  # Observation error on lake ATS estimates
  obs_error_prior <- lake_data$ats_se
  
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
    N2_init_prior = N2_init_prior$N2_init_log_mean,
    N2_init_sigma_prior = N2_init_prior$N2_init_log_sigma
  )
  
  return(stan_data)
}


# Create nested list of priors for each lake
lakes_stan_data <- unique(lakes_data$lake) |> 
  set_names() |> 
  map2(
    # Second winter mortality rate priors for GCL, HUC, SPR, respectively
    .y = logit(c(0.25, 0.5, 0.4)),
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
        meanlog = N2_init_prior,
        sdlog = N2_init_sigma_prior
      )
    )
  ) |>
  ungroup() |> 
  select(lake, mu_M, mu_theta, N2_init) |> 
  unnest(everything()) |> 
  pivot_longer(!lake) %>%
  split(.$name) 

some_priors |> 
  imap(
    \(x, idx) x |> 
      ggplot(aes(x = value)) +
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
    iter = 2000,
    chains = 4,
    control = list(adapt_delta = 0.9)
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


# Plot posterior versus prior distributions for mu_theta and 
# mu_M
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
        ggtitle(paste(idx, "predicted versus estimated age1 out-migrants"))
    }
  )


# Predicted versus observed Age2 out-migrants
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
        ggtitle(paste(idx, "predicted versus estimated age1 proportion"))
    }
  )


# Predicted Age1 proportions
lakes_stan_fits |> 
  imap(
    function(x, idx) {
      2 <- lakes_data |> 
        filter(lake == idx) |> 
        mutate(age1s = prop_age2*ats_est) |> 
        pull(age1s)
      Y <- lakes_stan_data[[idx]]$Y
      
      spread_draws(x, O2[year]) %>% 
        ggplot(aes(x = year, y = O2)) +
        stat_pointinterval(.width = c(0.5, 0.8, 0.95))+
        geom_point(data = data.frame(year = 1:Y, O2 = O2), color = "red") +
        ggtitle(paste(idx, "predicted versus estimated age2 out-migrants"))
    }
  )


# Predicted versus true smolting rate
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
        ggtitle(paste(idx, "predicted versus estimated smolting rate"))
    }
  )


# Predicted mortality rate
lakes_stan_fits |> 
  imap(
    \(x, idx) spread_draws(x, theta[year]) %>% 
        ggplot(aes(x = year, y = theta)) +
        stat_pointinterval(.width = c(0.5, 0.8, 0.95))+
        ggtitle(paste(idx, "predicted mortality rate"))
  )


# # Make dataframe with predicted and observed values for each year
# make_pred_frame <- function(fit, index) {
#   
#   # Extract posterior samples for theta
#   theta_samples <- extract(fit)$theta
#   
#   # Posterior mean estimates
#   theta_summary <- apply(
#     theta_samples, 
#     c(2, 3), # Specifies margins of output to summarize
#     # For this Stan output, 1 is the iterations, 2 is the time steps,
#     # and 3 is the categories, i.e. age classes
#     function(x) c(
#       mean = mean(x),
#       lower80 = quantile(x, 0.1, names = FALSE),
#       upper80 = quantile(x, 0.9, names = FALSE),
#       lower95 = quantile(x, 0.025, names = FALSE),
#       upper95 = quantile(x, 0.975, names = FALSE)
#     )
#   )
#   
#   # Convert to dataframe
#   theta_df <- map_dfr(
#     1:dim(theta_summary)[2], 
#     function(t) {
#       tibble(
#         year = filter(data, lake == index)$year[t],
#         fitted_age1_mean = theta_summary["mean", t, 1],
#         fitted_age1_lower80 = theta_summary["lower80", t, 1],
#         fitted_age1_upper80 = theta_summary["upper80", t, 1],
#         fitted_age1_lower95 = theta_summary["lower95", t, 1],
#         fitted_age1_upper95 = theta_summary["upper95", t, 1],
#         fitted_age2_mean = theta_summary["mean", t, 2],
#         fitted_age2_lower80 = theta_summary["lower80", t, 2],
#         fitted_age2_upper80 = theta_summary["upper80", t, 2],
#         fitted_age2_lower95 = theta_summary["lower95", t, 2],
#         fitted_age2_upper95 = theta_summary["upper95", t, 2]
#       )
#     }
#   )
#   
#   # Add fitted values to data
#   pred_frame <- data |> 
#     filter(lake == index) |> 
#     select(year, lake, matches("rate|prop")) |> 
#     left_join(theta_df) |> 
#     rename_with(.cols = matches("prop|rate"), \(x) paste0(x, "_mean")) |> 
#     pivot_longer(
#       matches("prop|rate|fitted"),
#       names_sep = "_",
#       names_to = c("source", "age", "interval")
#     ) |> 
#     pivot_wider(
#       names_from = c(source, interval),
#       names_sep = "_",
#       values_from = value
#     )
# }
# 
# 
# # Save dataframe with predictions for each of the three lakes
# lakes_pred_frame <- lakes_stan_fits |> 
#   imap(make_pred_frame) |> 
#   list_rbind()
# 
# 
# # Plot predicted versus observed data
# lakes_pred_frame |> 
#   ggplot(aes(x = year, y = prop_mean)) +
#   facet_wrap(
#     ~ lake,
#     ncol = 1,
#     strip.position = "right"
#   ) +
#   geom_linerange(
#     aes(
#       ymin = fitted_lower95, 
#       ymax = fitted_upper95,
#       colour = age
#     ),
#     linewidth = 0.25,
#     alpha = 0.6
#   ) +
#   geom_pointrange(
#     aes(
#       y = fitted_mean,
#       ymin = fitted_lower80,
#       ymax = fitted_upper80,
#       colour = age
#     ),
#     shape = "–",
#     linewidth = 0.75,
#     #alpha = 0.6
#   ) +
#   geom_point(
#     aes(fill = age, size = rate_mean),
#     colour = "black",
#     shape = 21,
#     alpha = 0.6
#   ) +
#   scale_y_continuous(labels = scales::percent) +
#   coord_cartesian(expand = FALSE) +
#   labs(y = "Percent composition") +
#   theme(panel.spacing.y = unit(1, "lines"))



 