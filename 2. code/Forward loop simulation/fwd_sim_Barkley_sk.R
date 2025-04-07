# Packages ----------------------------------------------------------------


pkgs <- c("here", "tidyverse", "readxl", "rstan", "ggrepel")
#install.packages(pkgs)


library(here)
library(tidyverse); theme_set(theme_bw())
library(readxl)
library(ggrepel)



# Load in data and fitted models ------------------------------------------


# Spawner-recruit data time series
sr_data <- here(
  "3. outputs",
  "Stock-recruit data",
  "Barkley_Sockeye_stock-recruit_infilled.xlsx"
) |> 
  read_xlsx(sheet = "S-R data") |> 
  # Calculate harvest rates
  mutate(hr = H/N)


# Fitted state-space spawner-recruit models
AR1_frame <- list.files(
  here(
    "3. outputs",
    "Stock-recruit modelling"
  ),
  pattern = "AR1.rds",
  full.names = TRUE
) |> 
  set_names(nm = ~str_extract(.x, "(SPR|GCL|HUC).*_AR1")) |> 
  map(readRDS) |> 
  enframe(name = "spec", value = "model") 
  
  
# Historic management forecast error
fcst_err <- here(
  "1. data",
  "Somass_Sockeye_forecast_error.xlsx"
) |> 
  read_xlsx() |> 
  rename_with(tolower) |> 
  filter(forecast == "MGT") 


# Show distribution of forecast errors
ggplot(fcst_err, aes(x = err_pct)) +
  geom_density(fill = "grey") +
  scale_y_continuous(expand = c(0, 0.05))


# Bootstrap mean & SD of forecast errors
fcst_err_boot <- fcst_err |> 
  select(err_pct) |> 
  expand_grid(sample = seq(1:1000)) |> 
  nest(.by = sample) |> 
  rowwise() |> 
  mutate(boot = list(sample(data, length(data), replace = TRUE))) |> 
  select(boot) |> 
  unnest(boot) |> 
  summarize(
    mean = mean(err_pct),
    sd = sd(err_pct)
  )


# Harvest rate implementation error
hr_err <- sr_data |> 
  # Calculate annual harvest rates on Somass and Hucuktlis
  mutate(group = if_else(stock == "HUC", "Hucuktlis", "Somass")) |> 
  summarize(
    .by = c(group, year),
    across(c(N, H), sum)
  ) |> 
  mutate(hr = H/N) |> 
  pivot_wider(
    id_cols = year,
    names_from = group,
    values_from = hr,
    names_glue = "{group}_hr"
  ) |> 
  right_join(fcst_err) |> 
  # Calculate harvest rate error
  mutate(hr_err = Somass_hr-target_hr) 


# Bootstrap mean & SD of HR implementation error
hr_err_boot <- hr_err |> 
  select(hr_err) |> 
  expand_grid(sample = seq(1:1000)) |> 
  nest(.by = sample) |> 
  rowwise() |> 
  mutate(boot = list(sample(data, length(data), replace = TRUE))) |> 
  select(boot) |> 
  unnest(boot) |> 
  summarize(
    mean = mean(hr_err),
    sd = sd(hr_err)
  )


# Fit simple LM to relate Hucuktlis HR to Somass HR
Hucuktlis_hr_lm <- hr_err |> 
  distinct(year, Somass_hr, Hucuktlis_hr) %>% 
  lm(Hucuktlis_hr ~ Somass_hr, data = .)


summary(Hucuktlis_hr_lm) # Fit looks pretty strong, 
# should be good for simulation purposes


# Plot Hucuktlis versus Somass HR in recent years
(hr_corr_p <- sr_data |> 
    # Calculate annual harvest rates on Somass and Hucuktlis
    mutate(group = if_else(stock == "HUC", "Hucuktlis", "Somass")) |> 
    filter(year > 2004) |> 
    summarize(
      .by = c(group, year),
      across(c(N, H), sum)
    ) |> 
    mutate(hr = H/N) |> 
    pivot_wider(
      id_cols = year,
      names_from = group,
      values_from = hr,
      names_glue = "{group}_hr"
    ) |> 
    ggplot(aes(x = Somass_hr, y = Hucuktlis_hr)) +
    geom_abline(slope = 1, lty = 2) +
    geom_smooth(method = "lm") +
    geom_point(aes(colour = year), size = 2) +
    geom_text_repel(aes(label = year)) +
    scale_x_continuous(
      limits = c(0, 0.75),
      expand = c(0, 0),
      labels = scales::percent,
      oob = scales::oob_keep
    ) +
    scale_y_continuous(
      limits = c(0, 0.75),
      expand = c(0, 0),
      labels = scales::percent,
      oob = scales::oob_keep
    ) +
    scale_color_viridis_c() + 
    coord_fixed(1) +
    labs(
      x = "Somass harvest rate",
      y = "Hucuktlis harvest rate"
    )
)


# Save the plot
ggsave(
  plot = hr_corr_p,
  filename = here(
    "3. outputs",
    "Plots",
    "Hucuktlis_vs_Somass_HRs_2005-2024.png"
  ),
  width = 6,
  units = "in",
  dpi = "print"
)
  

# Extract data required for the forward simulation ------------------------

  

sim_params <- AR1_frame |> 
  # Start with a list of model parameters whose posteriors we need
  # to extract
  mutate(
    extract = list(
      str_subset(
        names(model), 
        "lnalpha|beta|^S|^R|^C|U"
      )
    )
  ) |> 
  unnest_longer(extract) |> 
  # Filter parameters to ensure those estimated on an annual basis
  # are extracted only for the final year in the time series
  mutate(parameter = str_remove_all(extract, "\\[\\d+\\]")) |> 
  mutate(
    .by = c(spec, parameter),
    yr = as.numeric(str_extract(extract, "\\d+")),
    max_yr = max(yr, na.rm = TRUE),
    keep = if_else(
      !is.na(yr) & yr != max_yr,
      FALSE,
      TRUE
    )
  ) |> 
  filter(keep) |> 
  select(-yr, -max_yr, -keep) |> 
  rowwise() |> 
  mutate(
    posterior = list(as.data.frame(rstan::extract(model, pars = extract))),
    stock = str_extract(spec, "GCL|SPR|HUC"),
    fert = case_when(
      str_detect(parameter, "_fert") ~ 1, 
      str_detect(parameter, "_unfert") ~ 0,
      .default = NA
    ),
    data_scope = if_else(str_detect(spec, "trim"), "trim", "full"),
    parameter = str_remove_all(parameter, "_.*")
  )
