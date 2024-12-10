# Packages ----------------------------------------------------------------

pkgs <- c("here", "tidyverse", "readxl", "rstan")
#install.packages(pkgs)

library(here)
library(tidyverse); theme_set(theme_bw(base_size = 14))
library(rstan)
library(readxl)



# Load and format data for STAN input -------------------------------------


# Load stock-recruit time series by return year
run_ts <- here("1. data", "return by age time series.xlsx") |> 
  read_xlsx() |> 
  mutate(
    run = catch + escapement,
    brood_year = return_year - ttl_age
  )


# Calculate brood year returns
brood_ts <- run_ts |> 
  summarize(
    .by = c(brood_year, stock),
    return = sum(catch, escapement)
  )


# Build stock-recruit table
sr <- run_ts |> 
  mutate(
    .by = c(return_year, stock),
    spawners = sum(escapement)
  ) |> 
  pivot_wider(
    id_cols = c(return_year, stock, spawners),
    names_from = ttl_age,
    names_prefix = "N.age.",
    values_from = run,
    values_fn = sum
  ) |> 
  rowwise() |> 
  mutate(
    run = sum(c_across(contains("age"))),
    # Convert age #s to proportions
    across(contains("age"), \(x) x/run),
  ) |> 
  ungroup() |> 
  left_join(
    brood_ts,
    by = join_by(
      stock, 
      return_year == brood_year
    )
  ) |> 
  rename(
    "year" = return_year,
    "S" = spawners,
    "N" = run,
    "R" = return
  )


# Declare a function that transforms data for 1 stock into correct 
# STAN input format. Need to do this because will be required for both
# SPR and GCL.
make_stan_data <- function(stock) {
  
  fn_data <- sr[sr$stock == stock, ]
  
  A_obs <- fn_data |> 
    select(contains("age")) |> 
    as.matrix() |>
    base::`*`(1000) |> 
    round() 
    # Using 1000/year as placeholder... would need to do considerable
    # data gathering from historic files to get the time series of 
    # number of samples going back to 1977. Does this actually 
    # matter enough to be worth doing so?
    
  
  S_obs <- rnorm(
    length(fn_data$S),
    fn_data$S,
    0.05
  )
  
  H_obs <- rnorm(
    length(fn_data$S),
    (fn_data$N-fn_data$S),
    0.05
  )
  
  H_obs[H_obs<=0] <- 0.01 # Replace negative values with 0.01
  
  a_min <- min(run_ts$ttl_age) # youngest age at maturity
  a_max <- max(run_ts$ttl_age) # oldest age at maturity
  nyrs <- length(S_obs) # number of spawning years
  A <- a_max - a_min + 1 # total age classes
  nRyrs <- nyrs + A - 1 # number of recruitment years 
  
  
  stan.data <- list(
    "nyrs" = nyrs,
    "a_min" = a_min,
    "a_max" = a_max,
    "A" = A,
    "nRyrs" = nyrs + A - 1,
    "A_obs" = A_obs,
    "S_obs" = S_obs,
    "H_obs" = H_obs,
    "S_cv" = rep(0.05,length(S_obs)),
    "H_cv" = rep(0.05,length(H_obs))
  )
  
  return(stan.data)
}


# Save STAN data for GCL and SPR
GCL_stan_data <- make_stan_data("GCL")
SPR_stan_data <- make_stan_data("SPR")



# Fit STAN models for both stocks -----------------------------------------


# Embrace STAN model in a function call to iterate over stocks
fit_stan_mod <- function(stan_data) {
  stan(
    file = here(
      "2. code",
      "Stock-recruit modelling",
      "SS-SR_AR1.stan"
    ),
    model_name = "SS-SR_AR1",
    data = stan_data,
    # chains = 2,
    # iter = 500,
    # seed = 42,
    # thin = 1,
    # control = list(adapt_delta = 0.99, max_treedepth = 20)
  )
}


# Try stan model on GCL data
GCL_fit <- fit_stan_mod(GCL_stan_data)

