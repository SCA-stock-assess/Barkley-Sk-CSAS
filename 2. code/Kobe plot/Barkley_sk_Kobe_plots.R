# Packages ----------------------------------------------------------------


pkgs <- c(
  "here", "tidyverse", "readxl", "geom_textpath", "ggtext",
  "ggrepel"
  )
#install.packages(pkgs)


library(here)
library(tidyverse); theme_set(theme_bw())
library(readxl)
library(geomtextpath)
library(ggtext)
library(ggrepel)



# Import spawner & harvest time series and CU reference points ------------


# Spawner and harvest time series
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
    ER = H/N
  ) |> 
  select(year, stock, S, ER)
  

# CU specific reference points
ref_pts <- here(
  "3. outputs",
  "Stock-recruit modelling",
  "Barkley-Sk_CU_ref_pts_summary.csv"
) |> 
  read.csv() |> 
  filter(ref_pt %in% c("Smsy", "Umsy")) |> 
  pivot_wider(
    names_from = ref_pt,
    values_from = matches("q\\d+"),
    names_glue = "{ref_pt}_{.value}"
  ) |> 
  mutate(
    fert = case_when(
      stock == "Great Central" ~ "enhanced",
      stock == "Sproat" ~ "not enhanced",
      enh == 1 ~ "enhanced",
      enh == 0 ~ "not enhanced",
      stock == "Hucuktlis" & is.na(enh) ~ "mixed"
    ) |> 
      factor(levels = c("not enhanced", "mixed", "enhanced")),
    long_name = factor(
      long_name,
      levels = c(
        "Great Central Lake",
        "Sproat Lake",
        "Hucuktlis Lake (all data)",
        "Hucuktlis Lake (excl. 1993)",
        "Hucuktlis Lake (w/enhancement)",
        "Hucuktlis Lake (w/enhancement; excl. 1993)"
      )
    )
  )


# Add reference points manually from ResDoc table
ref_pts_alt <- tribble(
  ~method, ~stock, ~Umsy_q50, ~Umsy_q10, ~Umsy_q90, ~Smsy_q50, ~Smsy_q10, ~Smsy_q90,
  "S[MSY]", "Great Central", 0.58, 0.45, 0.68, 118488, 95382, 152044,
  "S[MSY]", "Sproat", 0.67, 0.54, 0.77, 89388, 73018, 113683,
  #"S-R", "Hucuktlis", 0.28, 0.07, 0.57, 9630, 2453, 47816,
  "50^th~percentile", "Hucuktlis", 0.28, 0.07, 0.57, 13948, NA_real_, NA_real_
) |> 
  mutate(stock = factor(stock, levels = c("Great Central", "Sproat", "Hucuktlis")))


# Join the Ref Pts to the abundance time series
kobe_data <- sr_data |> 
  left_join(
    ref_pts_alt, 
    by = "stock",
    relationship = "many-to-many"
  ) |> 
  mutate(
    x = S/Smsy_q50,
    y = ER/Umsy_q50,
    stock = factor(stock, levels = c("Great Central", "Sproat", "Hucuktlis"))
  ) |> 
  # Constrain x-axis range for Hucuktlis
  mutate(
    label = if_else(x > 5, round(x, 1), NA),
    x = if_else(x > 5, 5, x),
  ) |> 
  select(stock, year, x, y, method, label)



# Build Kobe plots --------------------------------------------------------


# dataframe with first and last years in the time series
end_yrs <- kobe_data |> 
  filter(
    .by = c(stock, method),
    year== min(year) | year== max(year)
  )


# Function to make kobe plot
make_kobe_p <- function(xy_data, reflines_data, endyrs_data, facet_rows) {
  
  ggplot(
    xy_data,
    aes(x = x, y = y)
    ) +
    facet_wrap(
      ~stock,
      scales = "free",
      nrow = facet_rows
    ) +
    geom_point(aes(colour = year), size = 1) +
    # Add red circles around first and last year data points
    geom_point(
      data = endyrs_data,
      colour = "red",
      shape = 21,
      size = 1,
      stroke = 1.25
    ) +
    geom_labelhline(
      data = reflines_data,
      aes(
        label = paste0("U[MSY]==", 100*round(Umsy_q50, 2), "*\'%\'"),
        yintercept = 1
      ),
      hjust = 0.95,
      colour = "grey50",
      boxcolour = NA,
      fill = "white",
      alpha = 0.5,
      label.padding = unit(0.01, "lines"),
      size = 3,
      lty = 2,
      parse = TRUE
    ) +
    geom_labelvline(
      data = reflines_data,
      aes(
        label = paste0(method, "==", round(Smsy_q50, 0)),
        xintercept = 1
      ),
      hjust = 0.95,
      colour = "grey50",
      boxcolour = NA,
      fill = "white",
      alpha = 0.5,
      label.padding = unit(0.01, "lines"),
      size = 3,
      lty = 2,
      parse = TRUE
    ) +
    # Add text labels at first and last year points
    geom_text_repel(
      data = endyrs_data,
      aes(label = paste0("'", str_sub(year, 3L, 4L))),
      size = 3
    ) +
    # Add text labels for x points > 5
    geom_text(
      aes(label = label),
      hjust = 1.15, 
      vjust = 1.25,
      size = 3
    ) +
    #add "crosshairs"
    scale_colour_viridis_c() +
    scale_y_continuous(
      limits = c(1e-6, NA),
      expand = expansion(mult = c(0, 0.06)),
      breaks = scales::pretty_breaks()
    ) +
    scale_x_continuous(
      limits = c(1e-6, NA),
      expand = expansion(mult = c(0, 0.01)),
      breaks = scales::pretty_breaks(),
      labels = scales::label_number(accuracy = 1)
    ) +
    # Fancy mathy axis labels
    labs(
      colour = "Return year",
      y = expression(frac(Exploitation~rate, U[MSY])), 
      x = expression(frac(Spawner~abundance,S[MSY])),
      parse = TRUE
    ) +
    theme_bw(base_size = 9) +
    theme(
      legend.position = "top",
      panel.spacing.x = unit(0.5, "lines"), # Move panels further apart
      panel.grid = element_blank(),
      legend.key.width = unit(dev.size()[1]/10, "in")
    )
}


# Make separate kobe plots for Hucuktlis and Somass CUs
kobe_plots <- list(
  xy_data = kobe_data,
  reflines_data = ref_pts_alt,
  endyrs_data = end_yrs
) |> 
  map(\(x) mutate(x, group = if_else(stock == "Hucuktlis", "Hucuktlis", "Somass"))) |> 
  map(\(x) split(x, x$group)) |> 
  map(\(x) enframe(x, name = "group")) |> 
  list_rbind(names_to = "name") |> 
  pivot_wider() |> 
  rowwise() |> 
  mutate(
    facet_rows = if_else(group == "Hucuktlis", 3, 1),
    kobe_p = list(make_kobe_p(xy_data, reflines_data, endyrs_data, facet_rows))
  ) |> 
  ungroup()


# View the plots
pull(kobe_plots, kobe_p, name = group)


# Save the Kobe plots
kobe_plots |> 
  mutate(height = if_else(group == "Hucuktlis", 5, 5)) |> 
  select(kobe_p, group, height) |> 
  pwalk(
    function(kobe_p, group, height) ggsave(
      plot = kobe_p,
      filename = here(
        "3. outputs",
        "Plots",
        paste0("Kobe_plots_", group,".png")
      ),
      width = 7,
      height = height,
      units = "in",
      dpi = "print"
    )
  )


# Alternatively, plot all panels together
kobe_alt <- make_kobe_p(kobe_data, ref_pts_alt, end_yrs, facet_rows = 1)


# Save the alternative version
kobe_alt |> 
  ggsave(
    filename = here(
      "3. outputs",
      "Plots",
      "Kobe_plots_all-CUs.png"
    ),
    width = 8,
    height = 3.75,
    units = "in",
    dpi = "print"
  )
