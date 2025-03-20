# Packages ----------------------------------------------------------------


pkgs <- c("here", "tidyverse", "readxl")
#install.packages(pkgs)


library(here)
library(tidyverse); theme_set(theme_bw())
library(readxl)


# Load data on smolt abundances -------------------------------------------


# Smolt abundance data
smolt_abun0 <- here(
  "1. data",
  "Barkley-Sk_smolt_abundances.xlsx"
) |> 
  read_xlsx() |> 
  pivot_longer(
    !smolt_year,
    names_to = c("cu", "age"),
    names_sep = "_"
  ) |> 
  mutate(cu = toupper(cu)) |> 
  rename("year" = smolt_year)


# Load total return data for Somass stocks to infill freshwater ages
som_fw_age_props <- here(
  "1. data",
  "Barkley Sockeye time series of returns.xlsx"
) |> 
  read_xlsx(sheet = "Somass") |> 
  filter(!is.na(fw_age)) |> 
  rename("cu" = stock) |> 
  mutate(
    brood_year = return_year - ttl_age,
    across(c(catch, escapement), as.numeric)
  ) |> 
  summarize(
    .by = c(brood_year, fw_age, cu),
    return = sum(catch, escapement)
  ) |> 
  mutate(
    .by = c(cu, brood_year),
    prop = return / sum(return),
    .keep = "unused"
  )


# Add Hucuktlis fw age data
huc_by_fw_age_props <- here(
  "1. data",
  "Barkley Sockeye time series of returns.xlsx"
) |> 
  read_xlsx(sheet = "Hucuktlis") |> 
  filter(!if_all(matches("age_\\d{2}"), is.na)) |> 
  pivot_longer(
    matches("age_\\d{2}"),
    names_to = "gr_age",
    names_prefix = "age_",
    values_to = "prop"
  ) |> 
  mutate(
    fw_age = str_sub(gr_age, -1),
    ttl_age = str_sub(gr_age, 1, 1),
    across(c(fw_age, ttl_age), as.numeric),
    across(c(catch, escapement), \(x) x*prop),
    brood_year = year - ttl_age,
    by_completion = if_else(
      (brood_year < min(year) - min(ttl_age)) |
      (brood_year + max(ttl_age) > max(year)),
      "incomplete",
      "complete"
    )
  ) |> 
  filter(by_completion == "complete") |> 
  summarize(
    .by = c(brood_year, fw_age),
    return = sum(catch, escapement)
  ) |> 
  mutate(
    .by = c(brood_year),
    prop = return / sum(return),
    cu = "HED",
    .keep = "unused"
  )


# Combine the Somass & Hucuktlis data
by_fw_age_props <- bind_rows(
  som_fw_age_props,
  huc_by_fw_age_props
)


# Average recent 4-year age 2 versus 3 compositions
fw_age_4y_avg <- by_fw_age_props |> 
  filter(!is.na(prop)) |> 
  filter(brood_year > max(brood_year) - 5) |> 
  summarize(
    .by = c(cu, fw_age),
    prop_4ya = mean(prop)
  )


# Use proportions from observed returns to infill missing smolt age data
smolt_abun <- smolt_abun0 |> 
  pivot_wider(
    names_from = age,
    values_from = value,
    names_prefix = "age_"
  ) |> 
  pivot_longer(
    cols = c(age_1, age_2),
    names_prefix = "age_",
    names_to = "age",
    values_to = "smolts",
    names_transform = as.numeric
  ) |> 
  rename("ttl_smolts" = age_ttl) |> 
  mutate(
    fw_age = age + 1,
    brood_year = year - fw_age
  ) |> 
  left_join(by_fw_age_props) |> 
  left_join(fw_age_4y_avg) |> # 4 year averages
  # Infill missing age data using historical average proportions
  mutate(
    infill = if_else(is.na(smolts) & !is.na(ttl_smolts), "infill", "obs"),
    smolts = case_when(
      is.na(smolts) & !is.na(prop) ~ ttl_smolts*prop, 
      # Use recent 4-year average for incomplete recent brood years
      is.na(smolts) & is.na(prop) &!is.na(ttl_smolts) ~ ttl_smolts*prop_4ya,
      TRUE ~ smolts
    )
  ) |> 
  mutate(
    .by = c(brood_year, cu),
    infill = if_else(any(str_detect(infill, "infill")), "infill", "obs"),
    cu = factor(
      cu,
      levels = c("GCL", "SPR", "HED"),
      labels = c("Great Central", "Sproat", "Hucuktlis")
    )
  )


# Plot --------------------------------------------------------------------


# Time series per lake
(smolt_abun_p <- smolt_abun |> 
   mutate(age = factor(age, levels = c("1", "2"))) |> 
   ggplot(aes(x = year, y = smolts)) +
   facet_wrap(
     ~cu,
     strip.position = "right",
     scales = "free_y",
     ncol = 1
   ) +
   geom_col(
     aes(fill = fct_rev(age), alpha = infill),
     position = "stack",
     colour = "black"
   ) +
   scale_alpha_discrete(range = c(0.5, 1)) +
   scale_fill_viridis_d(
     option = "mako",
     end = 0.75,
     direction = -1,
     name = "Age"
   ) +
   scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
   guides(alpha = "none") +
   labs(
     x = "Survey year",
     y = "Estimated smolt abundance (millions)"
   ) +
   theme(
     panel.spacing.y = unit(1, "lines"),
     strip.background = element_rect(fill = "white"),
     legend.position = "inside",
     legend.position.inside = c(0.98, 0.28),
     legend.justification.inside = c(1, 1),
     legend.background = element_rect(colour = "black")
   )
)


# Save the plot
ggsave(
  smolt_abun_p,
  filename = here(
    "3. outputs",
    "Plots",
    "Lake_smolts_abundance_time-series.png"
  ),
  width = 6.5,
  units = "in",
  dpi = "print"
)



# Calculate average abundances by age and CU ------------------------------


# Use the non-infilled data for calculating averages
smolt_abun0 |> 
  pivot_wider(
    names_from = age,
    values_from = value
  ) |> 
  mutate(
    age_ttl = `1` + `2`,
    across(
      c(`1`, `2`), 
      \(x) x/age_ttl,
      .names = "prop_{.col}"
    )
  ) |> 
  pivot_longer(c(ttl, prop_1, prop_2)) |> 
  filter(!is.na(value)) |> 
  summarize(
    .by = c(cu, name),
    quant = list(quantile(value, c(0.25, 0.5, 0.75)))
  ) |> 
  unnest_wider(quant)
