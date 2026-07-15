# Packages & setup --------------------------------------------------------


pkgs <- c("tidyverse", "here", "readxl", "janitor", "RColorBrewer", "mgcv", "MASS")
#install.packages(pkgs)


library(here)
library(MASS)
library(tidyverse); theme_set(theme_bw())
library(readxl)
library(janitor)
library(RColorBrewer)
library(mgcv)


# Load cleaned data -------------------------------------------------------


# These cleaned limno data are derived from:
# https://github.com/SiRE-P/Lake-Limnology
limno <- c("zpl", "ppl", "chem") |> 
  set_names() |> 
  map2(
    here(
      "1. data",
      "Limnology",
      "GCL_limno_cleaned-data_2008-2025.xlsx"
    ),
    \(x, y) read_xlsx(path = y, sheet = x)
  )


# Smolt data time series
smolts <- here(
  "3. outputs",
  "Stock-recruit data",
  "Bayesian_state-space_smolt-production_estimated_time_series.xlsx"
) |> 
  read_xlsx(sheet = "model_estimates") |> 
  filter(lake == "Great Central")


# Spawner-smolt residuals time series (freshwater survival proxy)
surv_idx <- here(
  "3. outputs",
  "Stock-recruit modelling",
  "Spawner-smolt posterior",
  "Barkley_Sockeye_spawner-smolt_residuals_summary.csv"
) |> 
  read.csv() |> 
  filter(cu == "Great Central")


# Compute annual indices of GCL productivity ------------------------------


# Sum zooplankton indices across stations within a survey
zpl_ttl <- limno$zpl |> 
  mutate(
    biomass = replace_na(biomass, 0),
    taxon = fct_lump_n(
      taxon, 
      n = 7, 
      w = biomass, 
      ties.method = "first",
      other_level = "OTHER"
    ),
    period = if_else(year > 2013, "recent", "historic"),
    taxon_period = interaction(taxon, period, sep = "_", drop = TRUE)
  ) |> 
  summarize(
    .by = c(date, year, taxon, period, taxon_period),
    across(c(biomass, density), \(x) sum(x, na.rm = TRUE))  
  ) |> 
  mutate(
    .by = period,
    day = as.integer(difftime(date, min(date), units = "days")),
    yday = as.integer(format(date, "%j")),
    survey = factor(date),
    taxon_year = interaction(taxon, year, sep = "_", drop = TRUE),
    across(c(taxon_year, taxon_period), droplevels)
  )


# Plot Zooplankton time series
zpl_ttl |> 
  ggplot(aes(date, biomass, colour = taxon)) +
  facet_wrap(
    ~period,
    scales = "free_x",
    nrow = 1
  ) +
  geom_point(position = "stack") +
  scale_color_brewer(palette = "Set1")


# Use seasonal GAM to standardize date? 
zpl_fit <- gam(
  biomass ~ taxon_period +
    s(day, by = taxon_period, bs = "cr", k = 8) +
    s(yday, by = taxon, bs = "cc", k = 8, id = 1) +
    s(taxon_year, bs = "re") +
    s(survey, bs = "re"),
  data = zpl_ttl,
  family = tw(link = "log"),
  method = "REML",
  knots = list(yday = c(0.5, 366.5))
)

gam.check(zpl_fit)


# Plot predicted versus observed values
pred_frame <- zpl_ttl |> 
  summarize(
    .by = c(period, taxon),
    min_date = min(date),
    max_date = max(date)
  ) |> 
  rowwise() |> 
  mutate(
    date = list(seq(min_date, max_date, by = "1 day")),
    .keep = "unused"
  ) |> 
  unnest(date) |> 
  mutate(
    .by = period,
    day = as.integer(difftime(date, min(date), units = "days"))
  ) |> 
  mutate(
    yday = as.integer(format(date, "%j")),
    year = as.integer(format(date, "%Y")),
    survey = unique(zpl_ttl$survey)[1],
    taxon_period = interaction(taxon, period, sep = "_", drop = TRUE),
    taxon_year = interaction(taxon, year, sep = "_", drop = TRUE)
  )


zpl_pred <- predict(
  zpl_fit, 
  pred_frame, 
  se.fit = TRUE,
  exclude = c("s(survey)", "s(taxon_year)")
) |> 
  cbind(pred_frame) |> 
  mutate(
    lci = fit - se.fit*1.96,
    uci = fit + se.fit*1.96,
    across(c(fit, lci, uci), exp)
  )


ggplot(
  data = zpl_pred,
  aes(x = day, y = fit)
) +
  facet_grid(
    taxon ~ period,
    scales = "free"
  ) +
  geom_ribbon(
    aes(ymin = lci, ymax = uci),
    alpha = 0.25
  ) +
  geom_line() +
  geom_point(
    data = zpl_ttl,
    aes(y = biomass)
  )


# Posterior simulations of annual maximum biomass
set.seed(1)

sim_max_biomass <- function(model, pred_grid, n = 2000) {
  
  X <- predict(
    model,
    newdata = pred_grid,
    type = "lpmatrix",
    exclude = "s(survey)"
  )
  
  n_draw <- n
  
  beta_draw <- MASS::mvrnorm(
    n = n_draw,
    mu = coef(model),
    Sigma = vcov(model, unconditional = TRUE)
  )
  
  eta_draw <- X %*% t(beta_draw)
  mu_draw <- exp(eta_draw)
  
  
  posterior_summary <- mu_draw |> 
    as_tibble() |> 
    cbind(pred_grid) |>
    pivot_longer(
      cols = starts_with("V"),
      names_to = "draw",
      values_to = "biomass_hat"
    ) |>
    summarize(
      total_biomass = sum(biomass_hat),
      .by = c(draw, year, date)
    ) |>
    slice_max(
      total_biomass,
      n = 1,
      with_ties = FALSE,
      by = c(draw, year)
    )
  
  annual_max <- posterior_summary |>
    summarize(
      .by = year,
      max_estimate = median(total_biomass),
      max_lwr = quantile(total_biomass, 0.025),
      max_upr = quantile(total_biomass, 0.975),
      peak_date_median = as.POSIXct(median(as.numeric(date))),
      peak_date_lwr = as.POSIXct(quantile(as.numeric(date), 0.025)),
      peak_date_upr = as.POSIXct(quantile(as.numeric(date), 0.975))
    )
  
  return(annual_max)
  
}


zpl_max_biomass_est <- sim_max_biomass(zpl_fit, pred_frame)

observed_max <- zpl_ttl |>
  summarize(
    observed_max = sum(biomass),
    .by = c(year, date)
  ) |>
  slice_max(
    observed_max,
    n = 1,
    with_ties = FALSE,
    by = year
  ) |> 
  mutate(year = as.integer(year))


left_join(zpl_max_biomass_est, observed_max, by = "year") |> 
  ggplot(aes(year, max_estimate)) +
  geom_pointrange(aes(ymin = max_lwr, ymax = max_upr)) +
  geom_point(
    aes(y = observed_max),
    colour = "red",
    size = 3
  ) +
  labs(
    y = "Annual maximum total biomass",
    caption = "Black: modeled maximum expected biomass; red: observed survey maximum"
  )



# Look at correlations with fish data -------------------------------------


# Time series of response variables (lagged to match lake year)
smolts_trim <- smolts |> 
  filter(parameter %in% c("BYO", "BYB", "w1")) |> 
  mutate(year = year + 1) |> # Convert to lake year
  pivot_wider(
    id_cols = year,
    names_from = parameter,
    values_from = `50%`
  )

surv_idx_trim <- surv_idx |> 
  filter(Rmeas == "BYO") |> 
  mutate(year = year + 1) |> # Convert to lake year
  select(year, ssr = X50.)

smolts_trim |> 
  left_join(surv_idx_trim) |> 
  pivot_longer(
    cols = c(w1, BYO, BYB, ssr),
    names_to = "response"
  ) |> 
  left_join(observed_max) |> 
  ggplot(aes(x = year, y = value)) +
  facet_wrap(
    ~ response,
    ncol = 1,
    strip.position = "left",
    scales = "free"
  ) +
  geom_point(aes(colour = observed_max)) +
  scale_colour_viridis_c(na.value = "grey") +
  theme(
    axis.title.y = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank()
  )


# Annual max zooplankton biomass
observed_max |> 
  left_join(
    smolts_trim,
    by = "year",
    relationship = "one-to-one"
  ) |> 
  left_join(
    surv_idx_trim,
    by = "year",
    relationship = "one-to-one"
  ) |> 
  pivot_longer(
    cols = c(w1, BYO, BYB, ssr),
    names_to = "response"
  ) |> 
  filter(!is.na(value)) |> 
  ggplot(aes(observed_max, value)) +
  facet_wrap(
    ~ response,
    ncol = 1,
    strip.position = "left",
    scales = "free"
  ) +
  #geom_path(colour = "grey") +
  geom_point(aes(colour = year)) +
  geom_smooth(method = "lm") +
  scale_colour_viridis_c() +
  theme(
    axis.title.y = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank()
  )


# Phosphorus max versus fish data
max_tp <- limno$chem |>
  filter(compound == "tp") |> 
  summarize(
    .by = c(date, year),
    tp = sum(measurement, na.rm = TRUE)
  ) |> 
  summarize(
    .by = year,
    max_tp = max(tp)
  )


max_tp |> 
  left_join(
    smolts_trim,
    by = "year",
    relationship = "one-to-one"
  ) |> 
  left_join(
    surv_idx_trim,
    by = "year",
    relationship = "one-to-one"
  ) |> 
  pivot_longer(
    cols = c(w1, BYO, BYB, ssr),
    names_to = "response"
  ) |> 
  filter(!is.na(value)) |> 
  ggplot(aes(max_tp, value)) +
  facet_wrap(
    ~ response,
    ncol = 1,
    strip.position = "left",
    scales = "free"
  ) +
  #geom_path(colour = "grey") +
  geom_point(aes(colour = year)) +
  geom_smooth(method = "lm") +
  scale_colour_viridis_c() +
  theme(
    axis.title.y = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank()
  )

