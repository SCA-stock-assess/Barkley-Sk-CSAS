# Packages ----------------------------------------------------------------

pkgs <- c("tidyverse", "here", "readxl", "ggridges")
#install.packages(pkgs)

library(tidyverse); theme_set(theme_ridges())
library(here)
library(readxl)
library(ggridges)


# Load historical Sockeye escapement data -----------------------------


# Escapement data
escday <- read_xlsx(
  here(
    "1. data",
    "Somass_Sockeye_daily_escapement_data.xlsx"
  ),
  na = ""
) |> 
  rename(
    "adj_adults" = "Adjusted Adults.  Includes bypass since 2004 but not Biosamples.",
    "adj_jacks" = "Adjusted Jacks.  Includes bypass since 2004 but not Biosamples."
  ) |> 
  mutate(
    year = if_else(
      year < 2000,
      as.numeric(paste0(19, year)),
      year
    ),
    d_m = paste0(day, "-", month.abb[month]),
    adj_adults = case_when(
      # Preserve NAs for days not yet observed in current year
      year == max(year) & is.na(`Adjusted net Adult up count`) ~ NA_real_, 
      is.na(adj_adults) & !is.na(`Stamp Falls Adjusted Adults`) ~ `Stamp Falls Adjusted Adults`,
      is.na(adj_adults) & year < max(year) ~ 0,
      TRUE ~ adj_adults),
    adj_jacks = case_when(
      is.na(adj_jacks) ~ `Adjusted net Jack up count`,
      is.na(adj_jacks) & year < max(year) ~ 0,
      TRUE ~ adj_jacks),
    across(
      c(adj_adults, adj_jacks),
      \(x) if_else(x < 0, 0, x)
    )
  ) |> 
  select(month, day, d_m, system, year, contains("adj_")) |> 
  rowwise() |> 
  mutate(sockeye = sum(c_across(c(adj_adults, adj_jacks)), na.rm = TRUE)) 



# Ridgeline plot ---------------------------------------------------------


# Plot for current year
(ridge_p <- escday |> 
   select(year, d_m, system, sockeye) |> 
   mutate(sockeye = round(sockeye, 0)) |> 
   uncount(sockeye) |> 
   # Subsample rows down to 500000 random draws
   # Useful when toying with plot parameters
   # slice_sample(n = 5e5) |> # deactivate when saving full version of plot
   ggplot(
     aes(
       y = year, 
       x = as.Date(d_m, format = "%d-%b"),
       group = year
     )
   ) +
   facet_wrap(
     ~system,
     labeller = labeller(
       system = c(
         "GCL" = "Stamp Falls fishway",
         "SPR" = "Sproat Falls fishway"
       )
     )
   ) +
   geom_density_ridges(
     quantile_lines = TRUE,
     quantiles = 2,
     vline_colour = "red",
     #rel_min_height = 0.001
   ) +
   scale_x_date(expand = c(0, 0)) +
   scale_y_reverse(expand = c(0, 0)) +
   coord_cartesian(clip = "off") +
   labs(
     x = "Date", 
     y = "Adult return year",
     title = "Sockeye counts through Somass system fishways"
   ) +
   theme(
     strip.background = element_blank(),
     strip.text = element_text(margin = margin(t = 0.5, b = 0.3, unit = "lines"))
   )
)


# Save the plot
ggsave(
  ridge_p,
  filename = here(
    "3. outputs",
    "Plots",
    "Somass_Sockeye_escapement_ridgelines.png"
  ),
  width = 6.5,
  height = 8,
  units = "in",
  dpi = "print"
)

