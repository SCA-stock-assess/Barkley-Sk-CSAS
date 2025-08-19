# Packages ----------------------------------------------------------------


pkgs <- c("here", "tidyverse", "ggpubr")
#install.packages


library(tidyverse)
library(here)
library(ggpubr)


# Create a custom theme to apply on all plots
cci_theme <- theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0, 1),
    legend.justification.inside = c(0, 1),
    legend.background = element_rect(
      colour = "black", 
      fill = alpha("white", 0.6),
      linewidth = 0.1
    ),
    legend.title = element_blank()
  )

theme_set(cci_theme)


# Load CCI data series ----------------------------------------------------


# Load all files and join into a single dataframe
cci <- list.files(
  here(
    "1. data",
    "California Current indices"
  ),
  pattern = ".csv",
  full.names = TRUE
) |> 
  map(read.csv) |> 
  map(\(x) mutate(x, date = as.Date(date))) |> 
  reduce(full_join)


# Long version of the data with proper names for ggplot
cci_long <- cci |> 
  pivot_longer(
    !date,
    names_to = "indicator"
  ) |> 
  filter(!is.na(value)) |> 
  mutate(
    long_name = case_when(
      indicator == "PDO" ~ "Pacific Decadal Oscillation",
      indicator == "NPGO" ~ "North Pacific Gyre Oscillation",
      indicator == "ONI" ~ "Ocean Niño Index (°C)",
      indicator == "CUTI" ~ "Coastal~Upwelling~Transport~Index~(45*degree*N)~(m^2%.%s^-1)",
      indicator == "BEUTI" ~ "Biologically~Effective~Upwelling~Transport~Index~(45*degree*N)~(mmol%.%s^-1%.%m^-1)",
      indicator == "MHW_intensity" ~ "Maximum intensity (std. dev.)",
      indicator == "MHW_pct_cover" ~ "Heatwave cover",
      indicator == "MHW_area" ~ "Maximum area"
    )
  )


# Replicate plots from CCIEA dashboard ------------------------------------


# Tidy versions for publication of the plots shown at:
# https://oceanview.pfeg.noaa.gov/dashboard/


# NPGO, PDO, and ONI
(large_indicators <- cci_long |> 
   filter(
     indicator %in% c("NPGO", "PDO", "ONI"),
     date > as.Date("1990-01-01")
   ) |> 
   ggplot(aes(date, value)) +
   geom_line(aes(colour = long_name)) +
   scale_y_continuous(expand = expansion(mult = c(0.05, 0.4))) +
   labs(
     x = NULL,
     y = "Anomaly score"
   )
)  


# Upwelling indices
(upwelling_indicators <- cci_long |> 
    filter(
      indicator %in% c("CUTI", "BEUTI"),
      date > as.Date("2010-01-01")
    ) |> 
    mutate(value = if_else(indicator == "BEUTI", value/13, value)) |> 
    ggplot(aes(date, value)) +
    geom_line(aes(colour = long_name)) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.4)),
      sec.axis = sec_axis(
        name = "BEUTI",
        transform = ~.x*13
      )
    ) +
    scale_colour_manual(
      values = c("turquoise3", "green3"),
      labels = function(x) {parse(text = x)}
    ) +
    labs(
      x = NULL,
      y = "CUTI"
    )
)


# Marine Heatwaves
(mhw_indicators <- cci_long |> 
    filter(
      indicator %in% c("MHW_pct_cover", "MHW_area"),
      date > as.Date("1990-01-01")
    ) |> 
    mutate(value = if_else(indicator == "MHW_area", value/10e6, value)) |> 
    ggplot(aes(date, value)) +
    geom_line(aes(colour = long_name)) +
    scale_y_continuous(
      limits = c(0, 1),
      expand = c(0, 0),
      labels = scales::percent,
      sec.axis = sec_axis(
        name = expression("Maximum area"~~(km^2)),
        transform = ~.x*10e6,
        labels = function(x) {
          ifelse(
            x == 0,
            "0",
            parse(
              text = gsub(
                "e\\+", 
                "%*%10^", 
                scales::scientific_format(digits = 1)(x)
              )
            )
          )
        }
      )
    ) +  
    labs(
      x = "Year",
      y = "Heatwave cover"
    )
)


# Assemble plots together in a grid ---------------------------------------


(aggregate_plot <- ggarrange(
  plotlist = list(large_indicators, upwelling_indicators, mhw_indicators),
  ncol = 1,
  align = "v",
  labels = "AUTO"
)  
) 


ggsave(
  plot = aggregate_plot,
  filename = here(
    "3. outputs",
    "Plots",
    "CCI_indicators_timeseries.png"
  ),
  width = 7,
  height = 7,
  units = "in",
  dpi = "print"
)
