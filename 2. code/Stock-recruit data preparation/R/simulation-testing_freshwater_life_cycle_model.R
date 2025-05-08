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


# Simulation approach based on Patrick Thompson's code ------------------- 


# Simulate some plausible data to validate the Stan model
# and better understand how priors and natural variability
# can affect estimate quality
set.seed(10)


# Nick messing around to understand sigma_theta values
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
c(logit(seq(0.15, 0.95, by = 0.1))) |> # mu_theta values
  map(sigma_theta_p)


# Make simulated data based on realistic values for the three lakes

# Function to simulate some data and stipulate model inputs
make_sim_data <- function(
    
  Y = 30, # number of years
  
  # Input variables (i.e. values the model will try to estimate)
  mu_theta, # logit-transformed smolting rate mean estimate
  sigma_theta, # interannual smolting rate variation following a normal distribution in logit space
  mu_M, # logit-transformed overwinter mortality rate
  sigma_M, # logit-transformed overwinter mortality variation
  N1_mu, # Mean number of age1 fry prior to smolting
  N1_phi = 10, # dispersion parameter around number smolting
  N2_init, # Initial number of age1 holdovers
  
  # Priors (with sensible default values)
  obs_error_prior = 5e5, # Not sure about this
  mu_theta_prior = 2,
  sigma_theta_prior = 0.5,
  mu_M_prior = -1.5,
  sigma_M_prior = 0.5,
  N2_init_prior = N2_init/2,
  N2_init_sigma_prior = 0.2
  
) {
  
  # Simulate numbers of age-1 and age-2 fry
  N1 <- rnbinom(n = Y, mu = N1_mu, size = N1_phi)
  N2 <- rep(NA, Y)
  N2[1] <- N2_init
  
  # Simulation smolting rate
  theta <- plogis(rnorm(Y, mu_theta, sigma_theta))
  
  # Mortality age1 to age2 in lake
  M <- plogis(rnorm(Y, mu_M, sigma_M)) 
  
  # True numbers of fry in the lake
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
  
  # Error around ATS estimate
  obs_error <- 0.25*N1_mu
  
  # Observed ATS estimate
  # Assumed to reflect true number (i.e. N_lake) with some added variability
  # No idea why log()*10 seems appropriate for dispersion parameter...
  N_obs <- rnbinom(Y, mu = N_lake, size = log(obs_error)*10) 
  
  # Observed fish outmigrating
  A_total <- rbinom(n = Y, size = N_obs, prob = 2e-5) # total outmigrants in sample
  A1_obs <- rbinom(n = Y, size = A_total, prob = O1 / (O1 + O2)) # Age1 fish in sample
  
  # Compile all desired outputs into a list
  dat_list <- list(
    Y = Y, 
    theta = theta,
    M = M,
    N_obs = N_obs, 
    N_lake = N_lake,
    A1_obs = A1_obs, 
    A_total = A_total, 
    O1 = O1,
    O2 = O2,
    obs_error_prior = obs_error_prior, # Not sure about this
    mu_theta_prior = mu_theta_prior,
    sigma_theta_prior = sigma_theta_prior,
    mu_M_prior = mu_M_prior,
    sigma_M_prior = sigma_M_prior,
    N2_init_prior = N2_init_prior,
    N2_init_sigma_prior = N2_init_sigma_prior
  )
  
  return(dat_list)
}


# Declare input values for simulation
sim_inputs <- data.frame(
  lake = c("Great Central", "Sproat", "Hucuktlis"),
  mu_theta = logit(c(0.80, 0.85, 0.95)), # Estimates from the calculations
  sigma_theta = c(0.5, 0.5, 0.9), # Best guesses of the estimated distributions
  mu_M = logit(c(.2, .45, .7)), # haphazard guesses 
  sigma_M = c(0.4, 0.5, 0.7), # Looks plausible when plotted
  N1_mu = c(7e6, 6e6, 1.5e6), # Somewhat reflective of the observed data
  N1_phi = c(8, 10, 5), # smaller phi = more variable
  N2_init = c(0.15*7e6, 0.1*6e6, 0.05*1.5e6) # Based on estimated smolt compositions
)

# Declare priors for model fits 
priors <- data.frame(
  lake = c("Great Central", "Sproat", "Hucuktlis"),
  obs_error_prior = c(7e6*0.3, 6e6*0.3, 1.5e6*0.3), # Not sure about this
  mu_theta_prior = logit(c(0.85, 0.9, 0.99)), # Assume priors will be biased high
  sigma_theta_prior = rep(1, 3),
  mu_M_prior = logit(c(0.25, 0.25, 0.5)), # We will assume Hucuktlis will be higher
  sigma_M_prior = rep(0.5, 3),
  N2_init_prior = c(5e5, 4e5, 5e4),
  N2_init_sigma_prior = rep(0.2, 3) # No idea about this
)

# Data for the three lakes
sim_data <- sim_inputs |> 
  left_join(priors) %>%
  split(.$lake) |> 
  map(\(x) select(x, -lake)) |> 
  map(as.vector) |> 
  map(\(x) do.call(make_sim_data, x))


# Look at lake observed versus actual fry
sim_data |> 
  imap(\(x, name) plot(x$N_obs ~ x$N_lake, main = name))
# Looks plausible given observed variability


# Function to fit Stan model
fit_stan <- function(stan_data) {
  sampling(
    object = stan_model(
      here(
        "2. code",
        "Stock-recruit data preparation",
        "Stan",
        "age_prop_beta_PT.stan"
      )
    ), 
    data = stan_data, 
    chains = 4, 
    cores = 4, 
    iter = 1500
  )
}


# Save fitted stan models for the three simulated datasets
sim_stan_fits <- map(sim_data, fit_stan)


#assess convergence###
worst_Rhat <- sim_stan_fits |> 
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
sim_stan_fits |> 
  imap(
    \(x, idx)
    traceplot(
      x, 
      pars = rownames(filter(worst_Rhat, lake == idx))[1:20]
    ) +
      ggtitle(label = idx)
  )


# Pair plots
sim_stan_fits |> 
  map(\(x) pairs(x, pars = c("mu_theta", "mu_M")))


#explore posterior###
post <- sim_stan_fits |> 
  map(extract) |> 
  # Discard all parameters with year-specific estimates
  map(\(x) discard(.x = x, .p = \(y) is.matrix(y))) |> 
  map(\(x) imap(x, as_tibble_col)) |> 
  map(\(x) list_cbind(unname(x))) |> 
  list_rbind(names_to = "lake")


target <- sim_inputs |> 
  select(any_of(colnames(post)))


prior_dist <- priors |> 
  rowwise() |> 
  mutate(
    mu_M = list(
      rnorm(
        1e5,
        mean = mu_M_prior,
        sd = sigma_M_prior
      )
    ),
    mu_theta = list(
      rnorm(
        1e5,
        mean = mu_theta_prior,
        sd = sigma_theta_prior
      )
    ),
    N2_init = list(
      rlnorm(
        1e5,
        meanlog = log(N2_init_prior),
        sdlog = N2_init_sigma_prior
      )
    )
  ) |>
  ungroup() |> 
  select(lake, mu_M, mu_theta, N2_init) |> 
  unnest(everything()) |> 
  pivot_longer(!lake) %>%
  split(.$name)


# Show posterior estimates versus simulation input values and priors
# Complex function employed to customize plots according to which prior
# information is pertinent, and to which have targets
post |> 
  pivot_longer(!lake) %>% 
  split(.$name) |> 
  imap(
    function(x, idx) {
      
      if(idx %in% colnames(target)) {
        target_vals <- select(target, lake, .data[[idx]])
      }
      
      plot <- x |> 
        ggplot(aes(x = value)) +
        facet_wrap(
          ~lake,
          scales = "free_x"
        ) +
        geom_density() +
        scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
        labs(x = idx)
      
      if(idx %in% names(prior_dist)) {
        plot2 <- plot + 
          geom_density(
            data = prior_dist[[idx]],
            colour = "red"
          ) +
          geom_textvline(
            data = target_vals,
            aes(xintercept = .data[[idx]]),
            label = "Simulation input",
            lty = 2
          )
        
      } else if(idx %in% colnames(target)) {
        plot2 <- plot +
          geom_textvline(
            data = target_vals,
            aes(xintercept = .data[[idx]]),
            label = "Simulation input",
            lty = 2
          )
      } else {
        plot2 <- plot
      }
      
      return(plot2)
    }
  )


# Predicted versus true Age1 out-migrants
sim_stan_fits |> 
  imap(
    function(x, idx) {
      O1 <- sim_data[[idx]]$O1
      Y <- sim_data[[idx]]$Y
      
      spread_draws(x, O1[year]) %>% 
        ggplot(aes(x = year, y = O1)) +
        stat_pointinterval(.width = c(0.5, 0.8, 0.95))+
        geom_point(data = data.frame(year = 1:Y, O1 = O1), color = "red") +
        ggtitle(paste(idx, "predicted versus true age1 out-migrants"))
    }
  )


# Predicted versus true Age2 out-migrants
sim_stan_fits |> 
  imap(
    function(x, idx) {
      O2 <- sim_data[[idx]]$O2
      Y <- sim_data[[idx]]$Y
      
      spread_draws(x, O2[year]) %>% 
        ggplot(aes(x = year, y = O2)) +
        stat_pointinterval(.width = c(0.5, 0.8, 0.95))+
        geom_point(data = data.frame(year = 1:Y, O2 = O2), color = "red") +
        ggtitle(paste(idx, "predicted versus true age2 out-migrants"))
    }
  )


# Predicted versus true smolting rate
sim_stan_fits |> 
  imap(
    function(x, idx) {
      theta <- sim_data[[idx]]$theta
      Y <- sim_data[[idx]]$Y
      
      spread_draws(x, theta[year]) %>% 
        ggplot(aes(x = year, y = theta)) +
        stat_pointinterval(.width = c(0.5, 0.8, 0.95))+
        geom_point(data = data.frame(year = 1:Y, theta = theta), color = "red") +
        ggtitle(paste(idx, "predicted versus true smolting rate"))
    }
  )


# Predicted versus true second winter mortality
sim_stan_fits |> 
  imap(
    function(x, idx) {
      M <- sim_data[[idx]]$M
      Y <- sim_data[[idx]]$Y
      
      spread_draws(x, M[year]) %>% 
        ggplot(aes(x = year, y = M)) +
        stat_pointinterval(.width = c(0.5, 0.8, 0.95)) +
        geom_point(data = data.frame(year = 1:Y, M = M), color = "red") +
        ggtitle(paste(idx, "predicted versus true second winter mortality"))
    }
  )


