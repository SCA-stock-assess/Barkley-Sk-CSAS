# Packages ----------------------------------------------------------------

pkgs <- c("here", "tidyverse", "readxl", "geomtextpath")
#install.packages(pkgs)


library(here)
library(tidyverse)
library(readxl)
library(geomtextpath)


# Load stock-recruit data and calculate annual harvest rates --------------


# Calculate harvest rates for each stock
hr_stock_data <- here(
  "3. outputs",
  "Stock-recruit data",
  "Barkley_Sockeye_stock-recruit_infilled.xlsx"
) |> 
  read_xlsx(sheet = "S-R data") |> 
  select(year, stock, S, H) |> 
  filter(
    !if_any(S:H, is.na),
    year > 1970
  ) |> 
  mutate(hr = H/(S+H))


# Calculate total harvest rate on all stocks
hr_ttl_data <- hr_stock_data |> 
  filter(
    .by = year,
    !length(unique(stock)) < 3
  ) |> 
  summarize(
    .by = year,
    S = sum(S),
    H = sum(H)
  ) |> 
  mutate(
    hr = H/(S+H),
    stock = "ttl"
  )


# Plot harvest rate time series -------------------------------------------


# Time series plot
(hr_ts <- hr_stock_data |> 
   filter(year > 1970) |> 
   mutate(
     stock = factor(
       stock, 
       levels = c("GCL", "SPR", "HED"),
       labels = c("Great Central", "Sproat", "Hucuktlis")
     )
   ) |> 
   ggplot(aes(x = year, y = hr)) +
   geom_line(
     data = hr_ttl_data,
     aes(colour = "Barkley total"),
     linewidth = 1,
     key_glyph = "timeseries"
   ) +
   geom_line(
     aes(colour = stock),
     alpha = 0.7,
     key_glyph = "timeseries"
   ) +
   scale_y_continuous(
     labels = scales::percent,
     limits = c(0, 1),
     expand = c(0, 0)
   ) +
   scale_color_manual(values = c("mediumblue", "red3", "seagreen", "black")) +
   guides(colour = guide_legend(byrow = TRUE)) +
   labs(
     x = NULL,
     y = "Harvest rate",
     colour = "Conservation Unit"
   ) +
   theme_classic() + 
   theme(
     legend.position = "inside",
     legend.position.inside = c(0.98, 0.98),
     legend.justification.inside = c(1, 1),
     legend.background = element_rect(colour = "black", fill = alpha("white", 0.7)),
     legend.key.size = unit(0.8, "lines"),
     legend.key.spacing.y = unit(0, "lines")
   )
)


# Save the plot
ggsave(
 hr_ts,
 filename = here(
   "3. outputs",
   "Plots",
   "Barkley_Sockeye_HR_time-series.png"
 ),
 width = 6.5,
 height = 4,
 units = "in",
 dpi = "print"
)
