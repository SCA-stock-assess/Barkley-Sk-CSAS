# Packages ----------------------------------------------------------------

pkgs <- c("here", "tidyverse", "readxl", "geomtextpath", "ggpubr")
#install.packages(pkgs)


library(here)
library(tidyverse); theme_set(theme_bw())
library(readxl)
library(geomtextpath)
library(ggpubr)


# Load time series data for each CU ---------------------------------------


# S-R data
sr_data <- here(
  "3. outputs",
  "Stock-recruit data",
  "Barkley_Sockeye_stock-recruit_infilled.xlsx"
) |> 
  read_xlsx(sheet = "S-R data") |> 
  mutate(
    stock = factor(
      stock,
      levels = c("GCL", "SPR", "HUC"),
      labels = c("Great Central", "Sproat", "Hucuktlis")
    ),
    HR = H/N
  ) |> 
  select(year, stock, H, S, HR, R)
  

# Reference points (from Working Paper)
ref_pts <- tribble(
  ~stock, ~refpt, ~q10, ~q50, ~q90, ~SEG,
  "Great Central", "Smsy", 103114, 141969, 226490, NA,
  "Great Central", "Umsy", 0.35, 0.51, 0.64, NA,
  "Great Central", "Sgen", 28735, 50729, 98582, NA,
  "Sproat", "Smsy", 79568, 102591, 155510, NA,
  "Sproat", "Umsy", 0.45, 0.61, 0.73, NA,
  "Sproat", "Sgen", 13041, 26103, 56265, NA,
  "Hucuktlis", "Smsy", 2453, 9630, 47816, 13948,
  "Hucuktlis", "Umsy", 0.07, 0.28, 0.57, NA,
  "Hucuktlis", "Sgen", 1492, 5328, 26274, 9576,
) |> 
  mutate(
    panel = if_else(
      refpt == "Umsy", 
      "Exploitation rate", 
      "Escapement (1000s)"
    ),
    across(
      c(q10, q50, q90, SEG),
      \(x) if_else(
        refpt == "Umsy",
        x,
        x/1000
      )
    ),
    stock = factor(stock, levels = c("Great Central", "Sproat", "Hucuktlis"))
  ) |> 
  pivot_wider(
    names_from = refpt,
    values_from = c(q10, q50, q90, SEG)
  )


# Data for the plots 
fourpp_data <- sr_data |> 
  pivot_longer(!c(stock, year)) |> 
  mutate(
    value = if_else(name != "HR", value/1000, value),
    long_name = case_when(
      name == "HR" ~ "Exploitation rate",
      name == "H" ~ "Harvest (1000s)",
      name == "S" ~ "Escapement (1000s)",
      name == "R" ~ "Recruitment (1000s)"
    )
  )
  


# Make the plots ----------------------------------------------------------


plots <- fourpp_data |> 
  nest(.by = c(stock, long_name), .key = "line_data") |> 
  left_join(
    ref_pts,
    by = join_by(
      stock,
      long_name == panel
    )
  ) |> 
  rowwise() |> 
  mutate(
    plot = list(      
      ggplot(
        line_data, 
        aes(
          x = year, 
          y = value
        )
      ) +
        geom_line(linewidth = 0.5) +
        # Add labelled horizontal reference line for Umsy
        geom_labelhline(
          yintercept = q50_Umsy,
          label = "*U*<sub>MSY</sub>",
          hjust = 0.1,
          rich = TRUE,
          lty = 2,
          linewidth = 0.4,
          fill = "white",
          alpha = 0.8,
          boxcolour = NA,
          label.padding = unit(0.05, "lines"),
          size = 2.5
        ) +
        # Add labelled horizontal reference line for Smsy
        geom_labelhline(
          yintercept = q50_Smsy,
          label = "*S*<sub>MSY</sub>",
          hjust = 0.1,
          rich = TRUE,
          lty = 2,
          linewidth = 0.4,
          fill = "white",
          alpha = 0.8,
          boxcolour = NA,
          label.padding = unit(0.05, "lines"),
          size = 2.5
        ) +
        #Add labelled horizontal reference line for Sgen
        geom_labelhline(
          yintercept = q50_Sgen,
          label = "*S*<sub>gen</sub>",
          hjust = 0.25,
          rich = TRUE,
          lty = 2,
          linewidth = 0.4,
          fill = "white",
          alpha = 0.8,
          boxcolour = NA,
          label.padding = unit(0.05, "lines"),
          size = 2.5
        ) +
        # Add SEGs
        geom_labelhline(
          yintercept = SEG_Smsy,
          lty = 3,
          linewidth = 0.4,
          size = 2.5,
          hjust = 0.5,
          fill = "white",
          alpha = 0.8,
          boxcolour = NA,
          label.padding = unit(0.05, "lines"),
          label = "50% SEG"
        ) +
        geom_labelhline(
          yintercept = SEG_Sgen,
          lty = 3,
          linewidth = 0.4,
          size = 2.5,
          hjust = 0.75,
          fill = "white",
          alpha = 0.8,
          boxcolour = NA,
          label.padding = unit(0.05, "lines"),
          label = "25% SEG"
        ) +
        # Add inter-quartile confidence band around Umsy 
        annotate(
          "rect",
          xmin = -Inf,
          xmax = Inf,
          ymin = q10_Umsy,
          ymax = q90_Umsy,
          fill = "black",
          alpha = 0.10
        ) +
        # Add confidence band around Smsy 
        annotate(
          "rect",
          xmin = -Inf,
          xmax = Inf,
          ymin = q10_Smsy,
          ymax = q90_Smsy,
          fill = "black",
          alpha = 0.10
        ) +
        # Add inter-quartile confidence band around Sgen
        annotate(
          "rect",
          xmin = -Inf,
          xmax = Inf,
          ymin = q10_Sgen,
          ymax = q90_Sgen,
          fill = "black",
          alpha = 0.10
        ) +
        scale_y_continuous(
          limits = if(long_name == "Exploitation rate") {
            c(0, 1)
          } else {
            c(0, NA)
          },
          expand = if(long_name == "Exploitation rate") {
            c(0, 0)
          } else {
            expansion(mult = (c(0, 0.07)))
          },
          labels = if(long_name == "Exploitation rate") {
            scales::percent
          } else {
            scales::label_number()
          }
        ) +
        labs(
          x = NULL,
          y = long_name
        ) +
        guides(colour = "none") +
        theme(
          axis.title = element_text(size = 9),
          panel.grid = element_blank(),
          axis.text.x = if(long_name %in% c("Harvest (1000s)", "Escapement (1000s)")) {
            element_blank()
          },
          axis.ticks.x = if(long_name %in% c("Harvest (1000s)", "Escapement (1000s)")) {
            element_blank()
          }
        ) 
    ),
    name = paste0(stock, "_", long_name)
  ) |> 
  pull(plot, name = name)


# Arrange plots in a 4-panel grid
(four_panel_plots <- list(
  GCL = keep_at(plots, \(x) str_detect(x, "Great Central")),
  SPR = keep_at(plots, \(x) str_detect(x, "Sproat")),
  HUC = keep_at(plots, \(x) str_detect(x, "Hucuktlis"))
) |> 
    map(
      \(x)
      ggarrange(
        plotlist = x, 
        nrow = 2, 
        ncol = 2, 
        align = "v", 
        labels = "auto",
        label.x = 0.87,
        label.y = 0.97
      )
    )
)

  
# Save the 4-panel grid plot
four_panel_plots |> 
  iwalk(
    \(x, idx) 
    ggsave(
      plot = x,
      filename = here(
        "3. outputs",
        "Plots",
        paste0("FSAR_4-panel_plot_", idx, ".png")
      ),
      height = 5,
      width = 7,
      units = "in",
      dpi = "print"
    )
  )


# Simple table summaries of the plot data ---------------------------------


# Starting point data frame with refpoints joined in
full_frame <- fourpp_data |> 
  left_join(
    ref_pts,
    by = join_by(
      stock,
      long_name == panel
    )
  ) |> 
  select(-matches("q(10|90)")) |> 
  mutate(
    .by = stock,
    nyrs = length(unique(year))
  )


# Numbers and percentages of years above/below benchmarks
above_below_ts <- full_frame |> 
  mutate(
    across(
      matches("(S|U)(gen|msy)"),
      \(x) if_else(value > x, 1, 0),
      .names = "{.col}_pass"
    )
  ) |> 
  filter(!if_all(contains("pass"), is.na)) |> 
  select(stock, year, nyrs, contains("pass")) |> 
  rename_with(\(x) str_replace_all(x, c("q50" = "SR", "_pass" = ""))) |> 
  pivot_longer(
    matches("(S|U)(gen|msy)"),
    names_sep = "_",
    names_to = c("method", "refpt"),
    values_to = "n_above"
  ) |> 
  filter(!is.na(n_above)) |> 
  pivot_wider(
    names_from = c(method, refpt),
    names_sep = "_",
    values_from = n_above
  )
  

above_below_ts |> 
  summarize(
    .by = c(stock, nyrs),
    across(matches("(S|U)(gen|msy)"), \(x) sum(x, na.rm = TRUE))
  ) |> 
  pivot_longer(
    cols = !c(stock, nyrs),
    names_to = "refpt",
    values_to = "n_above"
  ) |> 
  filter(!n_above == 0) |> # These are technically NAs
  mutate(pct_above = n_above/nyrs)
