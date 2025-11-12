# Packages ----------------------------------------------------------------

pkgs <- c("tidyverse", "here", "readxl")
#install.packages(pkgs)

library(here)
library(tidyverse); theme_set(theme_bw())
library(readxl)


# Load data and build plot ------------------------------------------------


# Hucuktlis escapement time series
Huc_esc <- read.csv(
  here(
    "3. outputs",
    "Stock-recruit data",
    "Barkley_Sockeye_observed_spawners_for-SueG.csv"
  )
) |> 
  filter(
    stock == "HUC",
    year > 1980
  )


# Percentile benchmarks for two time periods 
Huc_bm <- tribble(
  ~year_start, ~year_end, ~p25, ~p50,
  -Inf, 2000, 26760, 41056,
  2007, Inf, 9576, 13948
)


# Plot
(Huc_bm_p <- ggplot(
  Huc_esc,
  aes(
    x = year,
    y = S
  )
) +
    geom_rect(
      data = Huc_bm,
      aes(
        xmin = year_start,
        xmax = year_end,
        x = year_start,
        y = p25,
        ymin = 0, 
        ymax = p25
      ),
      fill = "#F1B6DA",
      alpha = 0.65
    ) +
    geom_rect(
      data = Huc_bm,
      aes(
        xmin = year_start,
        xmax = year_end,
        x = year_start,
        y = p25,
        ymin = p25, 
        ymax = p50
      ),
      fill = "#FFFFBF",
      alpha = 0.65
    ) +
    geom_rect(
      data = Huc_bm,
      aes(
        xmin = year_start,
        xmax = year_end,
        x = year_start,
        y = p25,
        ymin = p50, 
        ymax = Inf
      ),
      fill = "#B8E086",
      alpha = 0.65
    ) +
    geom_point() +
    geom_line() +
    scale_y_continuous(
      name = "Hucuktlis escapement",
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05)),
      breaks = c(0, 3, 6, 9, 12, 15, 18)*1e4,
      labels = scales::label_number()
    ) +
    scale_x_continuous(
      name = "Return year",
      expand = expansion(mult = c(0, 0.05))
    ) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
)


# Export the plot
ggsave(
  plot = Huc_bm_p,
  filename = here(
    "3. outputs",
    "Plots",
    "Hucuktlis_Escapement_percentile-benchmarks_timeseries.png"
  ),
  width = 6,
  height = 4,
  units = "in",
  dpi = "print"
)
