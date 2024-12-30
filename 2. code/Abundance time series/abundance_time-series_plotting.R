
# Packages ----------------------------------------------------------------

pkgs <- c("tidyverse", "here", "readxl")
#install.packages(pkgs)

library(tidyverse); theme_set(theme_bw())
library(here)
library(readxl)



# Read data on full abundance time series ---------------------------------


# These data include reconstructed catch data from Henderson based on
# a relationship to Somass Sockeye harvest rate. To update these data,
# run the script under 2. code/Stock-recruit data preparation
abun_data <- here(
  "3. outputs",
  "Stock-recruit data",
  "Barkley_Sockeye_stock-recruit_infilled.xlsx"
) |> 
  read_xlsx(sheet = "S-R data") |> 
  rename(
    "Escapement" = S,
    "Catch" = H
  ) |> 
  select(year, stock, Escapement, Catch) |> 
  pivot_longer(
    c(Escapement, Catch),
    names_to = "source",
    values_to = "sockeye"
  ) |> 
  mutate(
    data_quality = case_when(
      year < 1990 ~ "moderate",
      year >= 1990 ~ "high",
      TRUE ~ "FIX"
    ),
    stock = case_when(
      stock == "GCL" ~ "Great Central",
      stock == "SPR" ~ "Sproat",
      stock == "HED" ~ "Hucuktlis"
    ) |> 
      factor(levels = c("Great Central", "Sproat", "Hucuktlis"))
  )



# Make time series plots --------------------------------------------------


# Plot separate panel for each CU
(abun_p <- ggplot(
  abun_data,
  aes(x = year, y = sockeye)
) +
  facet_grid(
    stock ~ .,
    scales = "free_y", 
    space = "free"
  ) +
  geom_col(
    aes(fill = source, alpha = data_quality)
  ) +
  scale_y_continuous(
    labels = scales::label_number(big.mark = " "), # Hair space pasted from word
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_fill_viridis_d(option = "mako", direction = -1, end = 0.8) +
  scale_alpha_discrete(range = c(1, 0.5)) +
  labs(
    y = "Number of Sockeye",
    x = "Adult return year",
    fill = "Category",
    alpha = "Data source"
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.98, 0.98),
    legend.justification.inside = c(1, 1),
    legend.box = "horizontal",
    legend.background = element_rect(colour = "black", fill = alpha("white", 0.8)),
    legend.key.size = unit(0.8, "lines"),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = NA)
  )
)


# Save plot
abun_p |> 
  ggsave(
    filename = here(
      "3. outputs",
      "Plots",
      "Barkley-Sk_abundance_time-series_CU.png"
    ),
    width = 7.5,
    height = 8,
    units = "in"
  )
