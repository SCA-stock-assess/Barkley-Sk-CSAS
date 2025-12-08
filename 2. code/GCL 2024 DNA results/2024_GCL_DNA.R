# Packages ----------------------------------------------------------------


pkgs <- c("here", "tidyverse", "janitor")
#install.packages(pkgs)


library(here)
library(tidyverse); theme_set(theme_bw())
library(janitor)


# Load and plot DNA results -----------------------------------------------



# DNA data collected from Stamp Falls July/August 2024
strays <- here(
  "1. data",
  "Stamp 2024 Biodata With Stock ID.csv"
) |> 
  read.csv() |> 
  clean_names() |> 
  mutate(sampling_date = as.Date(sampling_date, format = "%m/%d/%Y"))


# Plot as stacked bar
(gcl_strays_p <- strays |> 
  filter(id %in% c("GCL", "Sproat")) |> 
  mutate(
    id = factor(
      id, 
      levels = c("GCL", "Sproat"), 
      labels = c("Great Central", "Sproat")
    )
  ) |> 
  count(id, sampling_date) |> 
  ggplot(aes(x = sampling_date)) +
  geom_col(
    aes(
      y = n,
      fill = fct_rev(id)
    ),
    position = "stack",
    colour = "black"
  ) +
  scale_y_continuous(
    labels = scales::label_number(),
    expand = expansion(mult = c(0, 0.05)),
    name = "Number of Sockeye sampled"
  ) +
  scale_x_date(
    date_labels = "%d %b",
    name = NULL,
    breaks = unique(strays$sampling_date)
  ) +
  scale_fill_manual(
    values = c("grey75", "grey20"),
    name = "DNA-assigned\nConservation\nUnit"
  ) +
  theme(
    #axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0,1),
    legend.justification.inside = c(0, 1),
    legend.background = element_rect(colour = "black", linewidth = 0.15)
  )
)


# Export plot
ggsave(
  plot = gcl_strays_p,
  filename = here(
    "3. outputs",
    "Plots",
    "GCL_fishway_DNA-assignments_2024.png"
  ),
  width = 7,
  height = 4,
  units = "in",
  dpi = "print"
)
