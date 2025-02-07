# Script notes ------------------------------------------------------------


# The purpose of this script is to investigate previous estimates of annual
# smolt production made by Hyatt and try to reproduce those methods for 
# subsequent years (2008-present). There are inconsistencies in the historic
# data post-2007 that suggest an unstated change in methodology or rigour.


# Packages ----------------------------------------------------------------


pkgs <- c("here", "tidyverse", "readxl")
#install.packages(pkgs)


library(here)
library(tidyverse); theme_set(theme_bw())
library(readxl)



# Load supporting data files ----------------------------------------------


# Time series of smolt production estimates from Hyatt's work
hyatt_data <- here(
  "1. data",
  "Barkley-Sk_pre-smolt_abundances.xlsx"
) |> 
  read_xlsx() |> 
  # Trim all post-2007 data that appear unreliable
  filter(smolt_year < 2008) |> 
  # Following steps pivot the data and calculate production by brood year
  rename("huc_1" = huc_ttl) |> # assume all Hucuktlis fry smolt at age 1
  select(-contains("ttl")) |> 
  pivot_longer(
    matches("(gcl|spr|huc)_\\d+"),
    names_sep = "_",
    names_to = c("stock", "age"),
    values_to = "sockeye"
  ) |> 
  mutate(
    age = as.numeric(age),
    brood_year = smolt_year - age - 1,
    stock = toupper(stock)
  )


# Total adult returns by age for Somass CUs
Som_returns <- here(
  "1. data",
  "Barkley Sockeye time series of returns.xlsx"
) |> 
  read_xlsx(sheet = "Somass") |> 
  filter(!if_any(c(escapement, catch, ttl_age, fw_age), is.na)) |> 
  mutate(return = as.numeric(escapement) + catch) |> 
  select(stock, return_year, return, ttl_age, fw_age)


# Lazy solution for missing 2009-2010 Hucuktlis catch estimates:
# grab infilled values from formatted s-r data output
Huc_catch_infill <- here(
  "3. outputs",
  "Stock-recruit data",
  "Barkley_Sockeye_stock-recruit_infilled.xlsx"
) |> 
  read_xlsx(sheet = "S-R data") |> 
  filter(stock == "HUC") |> 
  select(year, H) |> 
  rename("catch" = H)


# Total adult returns by age for Hucuktlis
Huc_returns <- here(
  "1. data",
  "Barkley Sockeye time series of returns.xlsx"
) |> 
  read_xlsx(sheet = "Hucuktlis") |> 
  # Fill missing catch data using lookup
  left_join(
    Huc_catch_infill,
    by = "year"
  ) |> 
  select(-catch.x) |> 
  rename("catch" = catch.y) |> 
  filter(!if_any(c(escapement, catch, matches("^age_\\d+")), is.na)) |>
  mutate(
    return = escapement + catch,
    across(matches("^age_\\d+"), \(x) x*return)
  ) |> 
  pivot_longer(
    matches("^age_\\d+"),
    names_prefix = "age_",
    names_to = "age",
    values_to = "sockeye"
  ) |> 
  mutate(
    ttl_age = as.numeric(str_extract(age, "\\d")),
    fw_age = as.numeric(str_sub(age, 2, 2)),
    stock = "HUC"
  ) |> 
  select(year, stock, sockeye, ttl_age, fw_age) |> 
  rename(
    "return_year" = year,
    "return" = sockeye
  )


# Bind Somass and Hucuktlis adult return data
adult_returns <- bind_rows(
  Som_returns,
  Huc_returns
) |> 
  mutate(brood_year = return_year - ttl_age)



# Investigate brood year adult vs smolt production ----------------------


# Sum adults produced per fw age and CU by brood year
adult_by_prod <- adult_returns |> 
  summarize(
    .by = c(fw_age, stock, brood_year),
    n = sum(return)
  ) |> 
  # Calculate freshwater age composition
  mutate(
    .by = c(stock, brood_year),
    by_ttl = sum(n),
    prop = n/by_ttl,
    smolt_age = fw_age - 1,
    dataset = "adults"
  )


# Produce the same dataset for smolts
smolt_by_prod <- hyatt_data |> 
  summarize(
    .by = c(age, stock, brood_year),
    n = sum(sockeye)
  ) |> 
  # Calculate age composition
  mutate(
    .by = c(stock, brood_year),
    by_ttl = sum(n),
    prop = n/by_ttl,
    dataset = "smolts"
  ) |> 
  rename("smolt_age" = age)


# Make a dataset to show the contrast
by_prod <- bind_rows(
  adult_by_prod,
  smolt_by_prod
) |> 
  mutate(
    .by = c(stock, dataset),
    max_n = max(n, na.rm = TRUE),
    n = n/max_n
  ) |> 
  pivot_wider(
    id_cols = c(stock, brood_year, smolt_age),
    names_from = dataset,
    values_from = c(prop, n)
  )


# Merge the two datasets and plot
by_prod |> 
  ggplot(aes(x = prop_adults, y = prop_smolts)) +
  facet_grid(
    smolt_age ~ stock,
    scales = "free"
  ) +
  geom_abline(slope = 1) +
  geom_point()
# Correlation seems poor...


# Plot the data as a time series
by_prod |> 
  drop_na() |> 
  ggplot(aes(x = brood_year, y = prop_smolts)) +
  facet_grid(smolt_age ~ stock) +
  geom_segment(
    aes(
      xend = brood_year,
      yend = prop_adults
    ),
    colour = "grey"
  ) +
  geom_point(
    aes(size = prop_smolts),
    colour = "blue",
    alpha = 0.5
  ) +
  geom_point(
    aes(
      size = prop_adults,
      y = prop_adults
    ),
    colour = "red",
    alpha = 0.5
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    labels = scales::percent
  ) +
  scale_size(range = c(0.5, 3)) +
  guides(size = "none")


# Final plot just showing the discrepancies over time
by_prod |> 
  drop_na() |> 
  ggplot(aes(x = brood_year, y = prop_adults - prop_smolts)) +
  facet_grid(smolt_age ~ stock) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_segment(aes(yend = 0, xend = brood_year)) +
  geom_point() +
  scale_y_continuous(
    expand = c(0, 0.05),
    labels = scales::percent
  )


# Preliminary conclusion is that any relationship between smolt and adult age 
# compositions is loose at best...
#
# Likely the differences arise from years of favourable marine survival 
# affecting different smolt years from the same brood year. E.g. Age-2 smolt 
# survival may be much higher than age-1 smolt survival from the same brood 
# year if ocean-entry conditions are drastically different between their
# two ocean-entry years.

