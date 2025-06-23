# Packages ----------------------------------------------------------------

pkgs <- c("here", "tidyverse", "readxl", "geomtextpath", "ggpubr", "zoo")
#install.packages(pkgs)


library(here)
library(tidyverse); theme_set(theme_bw())
library(readxl)
library(geomtextpath)
library(ggpubr)
library(zoo)


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
  ~stock, ~refpt, ~q10, ~q50, ~q90,
  "Great Central", "Smsy", 103114, 141969, 226490,
  "Great Central", "Umsy", 0.35, 0.51, 0.64,
  "Great Central", "Sgen", 28735, 50729, 98582,
  "Sproat", "Smsy", 79568, 102591, 155510,
  "Sproat", "Umsy", 0.45, 0.61, 0.73,
  "Sproat", "Sgen", 13041, 26103, 56265,
  "Hucuktlis", "Smsy", NA, 13948, NA,
  "Hucuktlis", "Umsy", 0.07, 0.28, 0.57,
  "Hucuktlis", "Sgen", NA, 9576, NA,
) |> 
  mutate(
    panel = if_else(
      refpt == "Umsy", 
      "Exploitation rate", 
      "Escapement (1000s)"
    ),
    across(
      c(q10, q50, q90),
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
    values_from = c(q10, q50, q90)
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
    ) |> 
      factor(
        levels = c(
        "Harvest (1000s)",
        "Escapement (1000s)",
        "Exploitation rate",
        "Recruitment (1000s)"
        )
      )
  ) |> 
  arrange(stock, long_name, year) |> 
  mutate(
    .by = c(stock, long_name),
    # Calculate rolling geometric average of escapement
    gen_ma = if_else(
      name == "S", 
      exp(rollmean(x = log(value), k = 4, fill = NA, align = "right")),
      NA
    )
  )
  

# Aggregate dataset for the entire SMU
bs_agg <- fourpp_data |> 
  filter(year >= 1977) |> #Start the time series at the earliest year with all 3 CUs
  summarize(
    .by = c(year, long_name),
    value = sum(value, na.rm = TRUE)
  ) |> 
  pivot_wider(names_from = long_name) |> 
  mutate(
    `Exploitation rate` = `Harvest (1000s)` / (`Harvest (1000s)` + `Escapement (1000s)`),
    `Recruitment (1000s)` = if_else(year > max(year) - 6, NA, `Recruitment (1000s)`)
  ) |> 
  pivot_longer(
    !year,
    names_to = "long_name"
  ) |> 
  mutate(stock = "Barkley Aggregate")


# Make the plots ----------------------------------------------------------


# Lengthy code with necessary panel-specific customization to build the plots
plots <- fourpp_data |> 
  bind_rows(bs_agg) |> 
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
        # Base time series lines
        geom_line(
          aes(colour = "Annual data"),
          linewidth = 0.5
        ) +
        # Add additional line for rolling 4-year geometric mean escapement
        geom_line(
          aes(
            y = gen_ma,
            colour = "Generational\nmoving average"
          ),
          linewidth = 0.5
        ) +
        # Add labelled horizontal reference line for Umsy
        geom_labelhline(
          yintercept = q50_Umsy,
          label = "*U*<sub>MSY</sub>",
          hjust = 0.1,
          rich = TRUE,
          lty = 2,
          linewidth = 0.2,
          fill = "white",
          alpha = 0.8,
          boxcolour = NA,
          label.padding = unit(0.01, "lines"),
          size = 2.5
        ) +
        # Add labelled horizontal reference line for Smsy
        geom_labelhline(
          yintercept = q50_Smsy,
          label = if(stock != "Hucuktlis") {
            "*S*<sub>MSY</sub>"
          } else {
            "50<sup>th</sup> percentile"
          },
          hjust = 0.1,
          rich = TRUE,
          lty = 2,
          linewidth = 0.2,
          fill = "white",
          alpha = 0.8,
          boxcolour = NA,
          label.padding = unit(0.01, "lines"),
          size = 2.5
        ) +
        #Add labelled horizontal reference line for Sgen
        geom_labelhline(
          yintercept = q50_Sgen,
          label = if(stock != "Hucuktlis") {
            "*S*<sub>gen</sub>"
            } else {
              "25<sup>th</sup> percentile"
            },
          hjust = 0.45,
          rich = TRUE,
          lty = 2,
          linewidth = 0.2,
          fill = "white",
          alpha = 0.8,
          boxcolour = NA,
          label.padding = unit(0.01, "lines"),
          size = 2.5
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
        # Make panel-specific adjustments to y-axis formatting
        scale_y_continuous(
          limits = if(long_name == "Exploitation rate") {
            c(0, 1)
          } else {
            c(0, NA)
          },
          expand = if(long_name == "Exploitation rate") {
            c(0, 0)
          } else if(long_name == "Escapement (1000s)") {
            # Add more white space in escapement panel to make room for legend
            expansion(mult = (c(0, 0.15))) 
          } else {
            expansion(mult = (c(0, 0.07)))
          },
          labels = if(long_name == "Exploitation rate") {
            scales::percent
          } else {
            scales::label_number()
          }
        ) +
        # Add colour values to escapement time series lines
        scale_colour_manual(values = c("black", "grey65")) +
        labs(
          x = "Adult return year",
          y = long_name
        ) +
        # Disable legends in all panels except escapement (panel b)
        guides(
          colour = if(long_name != "Escapement (1000s)" | str_detect(stock, "Barkley")) {
            "none"
          } else {
            guide_legend()
          }
        ) +
        # Customize the theme, including removing x-axis text on top panels
        theme(
          axis.title = if(long_name %in% c("Harvest (1000s)", "Escapement (1000s)")) {
            element_blank()
          } else {
            element_text(size = 9)
          },
          panel.grid = element_blank(),
          panel.border = element_rect(
            fill = NA, 
            linewidth = unit(0.3, "lines"),
            colour = "black"
          ),
          legend.title = element_blank(),
          legend.background = element_rect(
            colour = "black", 
            fill = alpha("white", 0.8),
            linewidth = unit(0.1, "lines")
          ),
          legend.position = "inside",
          legend.position.inside = c(0.001, 0.999),
          legend.justification.inside = c(0, 1),
          # axis.ticks.x = if(long_name %in% c("Harvest (1000s)", "Escapement (1000s)")) {
          #   element_blank()
          # },
          axis.text.x = if(long_name %in% c("Harvest (1000s)", "Escapement (1000s)")) {
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
  HUC = keep_at(plots, \(x) str_detect(x, "Hucuktlis")),
  BS = keep_at(plots, \(x) str_detect(x, "Barkley"))
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

  
# Save the compiled 4-panel grid plots (one per CU/stock)
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
      height = 4.5,
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


# Total annual harvest on Barkley Sockeye
annual_catches <- full_frame |> 
  filter(name == "H") |> 
  summarize(.by = year, value = sum(value)*1000) |> 
  mutate(median = median(value))
