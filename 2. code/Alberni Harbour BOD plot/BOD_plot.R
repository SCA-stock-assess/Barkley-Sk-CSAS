# Packages ----------------------------------------------------------------


pkgs <- c(tidyverse, here, cowplot)
#install.packages(pkgs)


library(here)
library(tidyverse) ; theme_set(theme_bw())
library(cowplot)


# Load data and recreate plot ---------------------------------------------


bod_ts <- read.csv(here("1. data", "Alberni_Harbour_BOD.csv")) |> 
  mutate(method = if_else(year <= 1969, "Estimated", "Measured"))


# Label data
labs <- data.frame(
  label = c(
    "Primary clarifier and\nASB secondary treatment\ninstalled",
    "Secondary effluent\ntreatment expanded\nto AST",
    "Aug '07 to Apr '08\ntemporary PM4\nshutdown",
    "Aug '12\naeration lagoon\neliminated"
  ),
  bod = c(41, 31, 26, 21),
  year = c(1970, 1990, 2002, 2010)
)


# Plot the full time series
(bod_p <- ggplot(
  bod_ts,
  aes(x = year, y = bod)
) +
    geom_col(
      aes(fill = method),
      colour = "grey25",
      width = 1
    ) +
    geom_segment(
      data = data.frame(
        x = c(1970.5, 1994.5, 2007.5, 2012.5),
        yend = c(40, 30, 25, 20)
      ),
      aes(
        x = x,
        xend = x,
        y = 0, 
        yend = yend
      ),
      lty = 2
    ) +
    geom_text(
      data = labs,
      aes(label = label),
      vjust = 0,
      hjust = 0,
      size = 2
    ) +
    scale_y_continuous(
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.12))
    ) +
    scale_fill_manual(values = c("grey90", "grey40")) +
    scale_x_continuous(breaks = seq(1950, 2020, by = 10)) +
    labs(
      x = NULL,
      fill = NULL,
      y = expression("Biological Oxygen Demand "~~(t%.%day^-1))
    ) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.01, 0.99),
      legend.justification.inside = c(0, 1)
    )
)


# Save plot
ggsave(
  plot = bod_p,
  filename = here(
    "3. outputs",
    "Plots",
    "Alberni_Harbour_BOD_time-series.png"
  ),
  height = 3,
  width = 7,
  units = "in",
  dpi = "print"
)
