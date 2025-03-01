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



# Patrick Thompson's code -------------------------------------------------


library(tidybayes)

#simulate fake data##


# Nick messing around to understand sigma_theta
sigma_theta_p <- function(mu_theta) {
  data.frame(
  mu_theta = mu_theta,
  sigma_theta = seq(0.1, 1.5, by = 0.1)
) |> 
  rowwise() |> 
  mutate(sim_theta = list(plogis(rnorm(1e5, mu_theta, sigma_theta)))) |> 
  unnest(sim_theta) |> 
  ggplot(aes(x = sim_theta)) +
  facet_wrap(~sigma_theta, labeller = label_both) +
  geom_histogram(binwidth = 0.01) +
  geom_vline(xintercept = plogis(mu_theta), lty = 2) +
  coord_cartesian(expand = FALSE)
}
# Similar simulations likely useful in understanding the other parameters


# Print some plots showing sigma_theta and mu_theta combos
c(logit(seq(0.15, 0.95, by = 0.1))) |> 
  map(sigma_theta_p)


# Do simulations on realistic values for the three lakes

# Function to simulate some data
make_sim_data <- function(
    Y = 30, # years
    phi = 0.1, # dispersion parameter
    mu_theta, # logit-transformed smolting rate mean estimate
    sigma_theta, # interannual smolting rate variation following a normal distribution in logit space
    mu_M, # logit-transformed overwinter mortality rate
    sigma_M, # logit-transformed overwinter mortality variation
    N1_mu, # Mean number of age1 fry prior to smolting
    N1_phi = 10, # dispersion parameter around number smolting
    init_N2 # Initial number of age1 holdovers
) {
  
  N1 <- rnbinom(n = Y, mu = N1_mu, size = N1_phi)
  N2 <- rep(NA, Y)
  N2[1] <- init_N2
  
  theta <- plogis(rnorm(Y, mu_theta, sigma_theta))
  hist(theta)
  
  M <- plogis(rnorm(Y, mu_M, sigma_M)) # Mortality age1 to age2 in lake
  
  N_lake <- O1 <- O2 <- rep(NA, Y)
  for(y in 1:Y){
    if(y > 1){
      N2[y] = N1[y-1] * (1 - theta[y-1]) * M[y-1]
    }
    N_lake[y] = N1[y] + N2[y]
    
    # Outmigrants
    O1[y] = theta[y] * N1[y]
    O2[y] = N2[y]
  }
  
  obs_error <- 0.001
  N_obs <- rnbinom(Y, mu = N_lake, 1/obs_error)
  
  A_total <- rbinom(n = Y, size = N_obs, prob = 0.02) # Observed fish outmigrating
  A1_obs <- rbinom(n = Y, size = A_total, prob = O1 / (O1 + O2)) # Observed Age1 fish outmigrating
  
  dat_list <- list(
    Y = Y, 
    N_obs = N_obs, 
    N_lake = N_lake,
    A1_obs = A1_obs, 
    A_total = A_total, 
    obs_error_prior = 1/0.001,
    mu_theta_prior = 2,
    sigma_theta_prior = 0.5,
    mu_M_prior = -1.5,
    sigma_M_prior = 0.5,
    N2_init_prior = 300,
    N2_init_sigma_prior = 0.2
  )
  
  return(dat_list)
}


# Logit transformation
logit <- function(p) {log(p/(1-p))}


# Data for the three lakes
sim_data <- data.frame(
    lake = c("Great Central", "Sproat", "Hucuktlis"),
    mu_theta = logit(c(0.85, 0.9, 0.95)), 
    sigma_theta = c(0.5, 0.5, 0.9), 
    mu_M = logit(c(.2, .45, .7)), # haphazard guesses 
    sigma_M = c(0.4, 0.5, 0.7), 
    N1_mu = c(7e6, 6e6, 1.5e6), 
    N1_phi = c(8, 10, 5), # smaller phi = more variable
    init_N2 = c(0.15*7e6, 0.1*6e6, 0.05*1.5e6)
) %>%
  split(.$lake) |> 
  map(\(x) select(x, -lake)) |> 
  map(as.vector) |> 
  map(\(x) do.call(make_sim_data, x))
  

# Look at lake observed versus total fry
sim_data |> 
  map(\(x) plot(x$N_obs ~ x$N_lake))


#run stan model##
model <- stan_model(
  here(
    "2. code",
    "Stock-recruit data preparation",
    "Stan",
    "age_prop_beta_PT.stan"
  )
)
fit <- sampling(model, data = dat_list, chains = 4, cores = 4, iter = 1000)

#assess convergence###
worst_Rhat <- summary(fit)$summary %>% 
  as.data.frame() %>% 
  mutate(Rhat = round(Rhat, 3)) %>% 
  arrange(desc(Rhat))

worst_Rhat %>% 
  filter(n_eff>3) %>% 
  ggplot(aes(x = n_eff, y = Rhat))+
  geom_point()+
  geom_hline(yintercept = 1.01, lty = 2)+
  geom_vline(xintercept = 400, lty = 2)

head(worst_Rhat, n = 20)

traceplot(fit, pars = rownames(worst_Rhat)[1:20])
pairs(fit, pars = c("mu_theta", "mu_M"))

#explore posterior###
post <- extract(fit)

plot(density(post$mu_M))
abline(v = mu_M)

plot(density(post$mu_theta))
abline(v = mu_theta)

plot(density(post$init_N2))
abline(v = init_N2)

spread_draws(fit, O1[year]) %>% 
  ggplot(aes(x = year, y = O1)) +
  stat_pointinterval(.width = c(0.5, 0.8, 0.95))+
  geom_point(data = data.frame(year = 1:Y, O1 = O1), color = "red")

spread_draws(fit, O2[year]) %>% 
  ggplot(aes(x = year, y = O2)) +
  stat_pointinterval(.width = c(0.5, 0.8, 0.95))+
  geom_point(data = data.frame(year = 1:Y, O2 = O2), color = "red")

spread_draws(fit, theta[year]) %>% 
  ggplot(aes(x = year, y = theta)) +
  stat_pointinterval(.width = c(0.5, 0.8, 0.95))+
  geom_point(data = data.frame(year = 1:Y, theta = theta), color = "red") +
  scale_y_continuous(
    limits = c(0, 1),
    expand = c(0,0)
  )

spread_draws(fit, M[year]) %>% 
  ggplot(aes(x = year, y = M)) +
  stat_pointinterval(.width = c(0.5, 0.8, 0.95))+
  geom_point(data = data.frame(year = 1:Y, M = M), color = "red")
 