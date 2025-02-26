# Packages ----------------------------------------------------------------


pkgs <- c("tidyverse", "rstan", "here", "bayesplot")
#install.packages(pkgs)


library(here)
library(tidyverse); theme_set(theme_bw())
library(rstan)
library(bayesplot)


# Load and prepare raw data for Stan --------------------------------------


# Load the data
data <- here(
  "3. outputs",
  "Stock-recruit data",
  "Inferred_annual_smolt_age_composition.csv"
) |> 
  read.csv() |> 
  rename_with(\(x) str_replace_all(x, "\\.", "x")) |> 
  filter(year < 2017) # post-2017 won't have data


# Plot observed proportions
data |> 
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


# Specify sigma priors for each lake
lakes_sigma_priors <- c(1, 0.5, 0.3) |> 
  set_names(nm = unique(data$lake))


# Weakly informative prior from adult age compositions
adult_prior <- here(
  "3. outputs",
  "Stock-recruit data",
  "Annual_adult_fw-age_composition.csv"
) |> 
  read.csv() |>
  rename_with(\(x) str_replace_all(x, "\\.", "x")) |> 
  arrange(lake, year) |> 
  mutate(across(c(returns_age1, returns_age2), \(x) x * ttl)) |> 
  mutate(
    .by = lake,
    returns_age1x = lead(returns_age2),
    # Replace 0s with small positive number
    across(matches("returns_"), \(x) if_else(x == 0, 1e-6, x))
  ) |> 
  filter(year %in% unique(data$year))


# Function to create Stan data for each lake
make_stan_data <- function(lake_name) {
  
  # Filter data for a single lake
  lake_data <- data %>%
    filter(lake == lake_name) %>%
    select(year, sample_days, rate_age1, rate_age2) %>%
    arrange(year)
  
  # Identify missing rows
  is_observed <- ifelse(is.na(lake_data$rate_age1) | is.na(lake_data$rate_age2), 0, 1)
  
  # Replace NAs with placeholder values (needed for Stan input)
  rate <- lake_data %>%
    select(rate_age1, rate_age2) %>%
    mutate(
      across(everything(), ~replace_na(.x, 1e-6)),
      # Replace 0s with very small positive number
      across(everything(), \(x) if_else(x == 0, 1e-6, x))
    )
  
  rate_age1 <- pull(rate, rate_age1)
  rate_age2 <- pull(rate, rate_age2)
  
  # Define concentration parameter (weight/sample size proxy)
  conc <- lake_data$sample_days
  conc[is.na(conc)] <- mean(conc, na.rm = TRUE)  # Replace missing sample sizes with average
  
  # Lake-specific sigma prior
  sigma_prior <- lakes_sigma_priors[[lake_name]]
  
  # Lake-specific adult FW age comps prior
  alpha_prior <- adult_prior |> 
    filter(lake == lake_name) |> 
    select(returns_age1, returns_age2) |> 
    as.matrix()
  
  # Handle missing values in alpha_prior
  for (i in 1:nrow(alpha_prior)) {
    if (any(is.na(alpha_prior[i, ]))) {
      alpha_prior[i, ] <- colMeans(alpha_prior, na.rm = TRUE)  # Fill with mean values
    }
  }
  
  # Normalize alpha_prior to avoid zero issues
  alpha_prior <- alpha_prior + 1e-6
  alpha_prior <- alpha_prior / rowSums(alpha_prior) * rowSums(alpha_prior, na.rm = TRUE)
  
  # Convert adult return counts to proportions
  alpha_prior <- alpha_prior / rowSums(alpha_prior)
  

  # Stan data list
  stan_data <- list(
    T = nrow(lake_data),
    rate_age1 = rate_age1,
    rate_age2 = rate_age2,
    is_observed = as.integer(is_observed),
    sigma_prior = sigma_prior,
    alpha_prior = alpha_prior,
    sample_days = conc
  )
  
  return(stan_data)
}


# List of Stan data for each lake
lakes_stan_data <- unique(data$lake) |> 
  set_names() |> 
  map(make_stan_data)



# Fit Stan models and plot predictions ------------------------------------


# Compile and fit the Stan model
fit_stan_model <- function(stan_data, index) {
  stan(
    file = here(
      "2. code",
      "Stock-recruit data preparation",
      "Stan",
      "dirichlet_state_space_model.stan"
    ),
    model_name = index,
    data = stan_data,
    iter = 2000,
    chains = 3,
    control = list(adapt_delta = 0.9)
  )
}


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


model_pars <- lakes_stan_fits |> 
  map(
    \(x) x |> 
      rstan::extract() |> 
      map(as.data.frame) |> 
      enframe()
  ) |> 
  list_rbind(names_to = "stock")


# some diagnostics

# check the chains directly
# leading parameters
# lead_pars <- c("theta", "sigma", "epsilon_raw", "epsilon", "lp__")
# 
# (lead_pars_p <- lakes_stan_fits |> 
#     map(
#       \(x) mcmc_combo(
#         x, 
#         pars = c("sigma"),
#         combo = c("dens_overlay", "trace"),
#         gg_theme = legend_none()
#       ) 
#     )
# )
# 
# 
# # Dirichlet parameters  
# (Dir_pars_p <- lakes_stan_fits |> 
#     map(
#       \(x) mcmc_combo(
#         x, 
#         pars = paste0("Dir_alpha[", 1:4, "]"),
#         combo = c("dens_overlay", "trace"),
#         gg_theme = legend_none()
#       )
#     )
# )



# Plot predicted versus observed values -----------------------------------


# Make dataframe with predicted and observed values for each year
make_pred_frame <- function(fit, index) {
  
  # Extract posterior samples for theta
  theta_samples <- extract(fit)$theta
  
  # Posterior mean estimates
  theta_summary <- apply(
    theta_samples, 
    c(2, 3), # Specifies margins of output to summarize
    # For this Stan output, 1 is the iterations, 2 is the time steps,
    # and 3 is the categories, i.e. age classes
    function(x) c(
      mean = mean(x),
      lower80 = quantile(x, 0.1, names = FALSE),
      upper80 = quantile(x, 0.9, names = FALSE),
      lower95 = quantile(x, 0.025, names = FALSE),
      upper95 = quantile(x, 0.975, names = FALSE)
    )
  )
  
  # Convert to dataframe
  theta_df <- map_dfr(
    1:dim(theta_summary)[2], 
    function(t) {
      tibble(
        year = filter(data, lake == index)$year[t],
        fitted_age1_mean = theta_summary["mean", t, 1],
        fitted_age1_lower80 = theta_summary["lower80", t, 1],
        fitted_age1_upper80 = theta_summary["upper80", t, 1],
        fitted_age1_lower95 = theta_summary["lower95", t, 1],
        fitted_age1_upper95 = theta_summary["upper95", t, 1],
        fitted_age2_mean = theta_summary["mean", t, 2],
        fitted_age2_lower80 = theta_summary["lower80", t, 2],
        fitted_age2_upper80 = theta_summary["upper80", t, 2],
        fitted_age2_lower95 = theta_summary["lower95", t, 2],
        fitted_age2_upper95 = theta_summary["upper95", t, 2]
      )
    }
  )
  
  # Add fitted values to data
  pred_frame <- data |> 
    filter(lake == index) |> 
    select(year, lake, matches("rate|prop")) |> 
    left_join(theta_df) |> 
    rename_with(.cols = matches("prop|rate"), \(x) paste0(x, "_mean")) |> 
    pivot_longer(
      matches("prop|rate|fitted"),
      names_sep = "_",
      names_to = c("source", "age", "interval")
    ) |> 
    pivot_wider(
      names_from = c(source, interval),
      names_sep = "_",
      values_from = value
    )
}


# Save dataframe with predictions for each of the three lakes
lakes_pred_frame <- lakes_stan_fits |> 
  imap(make_pred_frame) |> 
  list_rbind()


# Plot predicted versus observed data
lakes_pred_frame |> 
  ggplot(aes(x = year, y = prop_mean)) +
  facet_wrap(
    ~ lake,
    ncol = 1,
    strip.position = "right"
  ) +
  geom_linerange(
    aes(
      ymin = fitted_lower95, 
      ymax = fitted_upper95,
      colour = age
    ),
    linewidth = 0.25,
    alpha = 0.6
  ) +
  geom_pointrange(
    aes(
      y = fitted_mean,
      ymin = fitted_lower80,
      ymax = fitted_upper80,
      colour = age
    ),
    shape = "–",
    linewidth = 0.75,
    #alpha = 0.6
  ) +
  geom_point(
    aes(fill = age, size = rate_mean),
    colour = "black",
    shape = 21,
    alpha = 0.6
  ) +
  scale_y_continuous(labels = scales::percent) +
  coord_cartesian(expand = FALSE) +
  labs(y = "Percent composition") +
  theme(panel.spacing.y = unit(1, "lines"))

