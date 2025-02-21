# Packages ----------------------------------------------------------------


pkgs <- c("tidyverse", "rstan", "here")
#install.packages(pkgs)


library(here)
library(tidyverse); theme_set(theme_bw())
library(rstan)



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


# Function to create Stan data for each lake
make_stan_data <- function(lake_name) {
  
  # Filter data for a single lake
  lake_data <- data %>%
    filter(lake == lake_name) %>%
    select(year, starts_with("rate_")) %>%
    arrange(year)
  
  # Ensure weights are conveyed into the model
  weights <- lake_data %>%
    rowwise() %>%
    mutate(
      weight = sum(c_across(starts_with("rate_"))),
      # Replace weights with a placeholder where data is missing
      # Placeholder value, 1, doesn't affect unobserved years
      weight = if_else(is.na(weight), 1, weight)
    ) %>%
    pull(weight)
  
  # Prepare data for Stan
  Y <- lake_data %>%
    select(starts_with("rate_")) %>%
    as.matrix()
  
  # Identify missing rows
  is_observed <- complete.cases(Y)
  
  # Replace NAs with placeholder values (needed for Stan input)
  Y[!is_observed, ] <- rep(1, ncol(Y))
  
  # Adjust 0 values to accomodate log-transformation in Dirichlet fit
  Y[Y == 0] <- 1e-6
  Y <- Y / rowSums(Y)
  
  # Calculate prior mean composition (ignoring NAs)
  prior_mean <- colMeans(Y[is_observed, ], na.rm = TRUE)
  prior_mean <- prior_mean / sum(prior_mean) # Ensure it is a simplex
  
  # Lake-specific sigma prior
  sigma_prior <- lakes_sigma_priors[[lake_name]]
  
  # Stan data list
  stan_data <- list(
    T = nrow(Y),
    K = ncol(Y),
    rate = Y,
    is_observed = as.integer(is_observed),
    prior_mean = prior_mean,
    sigma_prior = sigma_prior,
    conc = weights
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
    #warmup = 1500, # probably need lower ratio of warmup to total iterations
    chains = 4,
    control = list(adapt_delta = 0.9)
  )
}


# Fit models for each lake
lakes_stan_fits <- lakes_stan_data |> 
  imap(fit_stan_model)
# Currently yields numerous warnings, many of which imply poor fits
# High numbers of divergent transitions, estimated Bayesian Fractions of
# Missing Information are low in several chains
# Large R-hat values
# Low tail effective sample sizes
# Transitions exceeding maximum treedepth


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
        fitted_age2_upper95 = theta_summary["upper95", t, 2],
        fitted_age1x_mean = theta_summary["mean", t, 3],
        fitted_age1x_lower80 = theta_summary["lower80", t, 3],
        fitted_age1x_upper80 = theta_summary["upper80", t, 3],
        fitted_age1x_lower95 = theta_summary["lower95", t, 3],
        fitted_age1x_upper95 = theta_summary["upper95", t, 3]
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

