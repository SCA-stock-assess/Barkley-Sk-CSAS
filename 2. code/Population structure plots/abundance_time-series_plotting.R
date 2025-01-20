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
  "Total return time series",
  "Barkley Sockeye total annual returns by CU_collated.xlsx"
) |> 
  read_xlsx() |> 
  mutate(
    CU = case_when(
      CU == "GCL" ~ "Great Central",
      CU == "SPR" ~ "Sproat",
      CU == "HUC" ~ "Hucuktlis",
      CU == "GCL + SPR" ~ "Somass aggregate"
    ) |> 
      factor(levels = c("Great Central", "Sproat", "Hucuktlis", "Somass aggregate"))
  )



# Make time series plot for SMU aggregate --------------------------------------


# Plot total returns broken out between catch and escapement
(abun_p <- abun_data |> 
   pivot_longer(contains("annual")) |> 
   distinct(year, name, value) |> 
   mutate(name = str_remove(name, "annual_Barkley_")) |> 
   ggplot(aes(x = year, y = value/1000)) +
   geom_vline(
     xintercept = 1969, 
     lty = 2, 
     colour = "grey50"
   ) +
   geom_col(aes(fill = name), colour = "black") +
   annotate(
     "text",
     x = 1968,
     y = 600,
     label = "Lake Fertilization Program begins",
     angle = 90,
     colour = "grey50",
     vjust = 0,
     hjust = 0
   ) +
   scale_y_continuous(
     name = "Barkley Sockeye abundance (1000s)",
     labels = scales::label_number(),
     expand = expansion(mult = c(0, 0.05))
   ) +
   scale_fill_manual(values = c("grey", "black")) +
   labs(x = "Return year") +
   theme(
     legend.position = "inside",
     legend.position.inside = c(0.02, 0.98),
     legend.justification.inside = c(0, 1),
     legend.background = element_rect(colour = "black", fill = "white"),
     legend.title = element_blank()
   )
)


# Save plot
abun_p |> 
  ggsave(
    filename = here(
      "3. outputs",
      "Plots",
      "Barkley-Sk_total_abundance_time-series.png"
    ),
    width = 7.5,
    height = 4,
    units = "in",
    dpi = "print"
  )


# Plot CU-specific time series --------------------------------------------


# Annotation data frame stating when CU-specific catch data are not provided
cu_catch <- tribble(
  ~CU, ~year_start, ~year_end, ~ymin, ~ymax, ~y, ~label,
  "Hucuktlis", 1908, 1924.5, -Inf, Inf, 125, "CU-specific catch\nunavailable",
  "Hucuktlis", 1934, 1969.5, -Inf, Inf, 125, "CU-specific catch\nunavailable",
  "Great Central", 1908, 1976.5, -Inf, Inf, 600, "CU-specific catch\nunavailable",
  "Sproat", 1908, 1976.5, -Inf, Inf, 500, "CU-specific catch\nunavailable"
) |> 
  mutate(CU = factor(CU, levels = c("Great Central", "Sproat", "Hucuktlis"))) 



# The time series plot
(cu_abun_p <- abun_data |> 
    filter(!CU == "Somass aggregate") |> 
    pivot_longer(c(catch, escapement)) |> 
    ggplot(aes(year, value/1000)) +
    facet_wrap(
      ~CU,
      strip.position = "right",
      ncol = 1,
      scales = "free_y"
    ) +
    geom_rect(
      data = cu_catch,
      aes(
        y = ymax,
        x = year_start,
        xmin = year_start,
        xmax = year_end,
        ymin = ymin,
        ymax = ymax
      ),
      fill = "grey",
      alpha = 0.5,
      colour = NA
    ) +
    geom_text(
      data = cu_catch,
      aes(
        x = (year_end - year_start)/2+year_start,
        y = y,
        label = label
      ),
      size = 2,
      colour = "grey40"
    ) +
    geom_col(
      aes(fill = name),
      colour = "black"
    ) +
    scale_y_continuous(
      name = "Sockeye abundance (1000s)",
      labels = scales::label_number(),
      expand = expansion(mult = c(0, 0.05))
    ) +
    scale_x_continuous(
      name = "Return year",
      limits = c(1908, max(abun_data$year) +2),
      expand = c(0, 0)
    ) +
    scale_fill_manual(
      values = c("grey", "black")
    ) +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.02, 0.98),
      legend.justification.inside = c(0, 1),
      legend.background = element_rect(colour = "black", fill = "white"),
      legend.title = element_blank(),
      strip.background = element_rect(colour = "black", fill = "white")
    )
)


# Save the CU-specific plot
cu_abun_p |> 
  ggsave(
    filename = here(
      "3. outputs",
      "Plots",
      "Barkley-Sk_CU_abundance_time-series.png"
    ),
    width = 7.5,
    height = 6,
    units = "in",
    dpi = "print"
  )

