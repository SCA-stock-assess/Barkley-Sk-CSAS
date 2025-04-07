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
      stock == "Great Central" ~ "fertilized",
      stock == "Sproat" ~ "not fertilized",
      fert == 1 ~ "fertilized",
      fert == 0 ~ "not fertilized",
      stock == "Hucuktlis" & is.na(fert) ~ "mixed"
    ) |> 
      factor(levels = c("not fertilized", "mixed", "fertilized")),
    long_name = factor(
      long_name,
      levels = c(
        "Great Central Lake",
        "Sproat Lake",
        "Hucuktlis Lake (all data)",
        "Hucuktlis Lake (excl. 1993)",
        "Hucuktlis Lake (w/fertilization)",
        "Hucuktlis Lake (w/fertilization; excl. 1993)"
      )
    )
  )


# Join the Ref Pts to the abundance time series
kobe_data <- sr_data |> 
  left_join(
    ref_pts, 
    by = "stock",
    relationship = "many-to-many"
  ) |> 
  mutate(
    x = S/Smsy_q50,
    y = ER/Umsy_q50
  ) |> 
  select(stock, long_name, year, fert, x, y)



# Build Kobe plots --------------------------------------------------------


# dataframe with first and last years in the time series
end_yrs <- kobe_data |> 
  filter(
    .by = c(long_name, fert),
    year== min(year) | year== max(year)
  )


# Function to make kobe plot
make_kobe_p <- function(xy_data, reflines_data, endyrs_data, facet_rows) {
  
  ggplot(
    xy_data,
    aes(x = x, y = y)
    ) +
    facet_wrap(
      ~long_name + fert,
      scales = "free",
      nrow = facet_rows,
      labeller = labeller(.multi_line = FALSE)
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
      #fill = alpha("white", 0.75),
      size = 3,
      lty = 2,
      parse = TRUE
    ) +
    geom_labelvline(
      data = reflines_data,
      aes(
        label = paste0("S[MSY]==", round(Smsy_q50, 0)),
        xintercept = 1
      ),
      hjust = 0.95,
      colour = "grey50",
      boxcolour = NA,
      #fill = alpha("white", 0.75),
      size = 3,
      lty = 2,
      parse = TRUE
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
    # Add text labels at first and last year points
    geom_richtext(
      data = endyrs_data,
      aes(label = paste0("'", str_sub(year, 3L, 4L))),
      hjust = 0.25, 
      vjust = -0.5,
      size = 3,
      label.colour = NA,
      label.padding = unit(rep(0.01, 4), "lines"),
      fill = alpha("white", 0.75)
    ) +
    #add "crosshairs"
    scale_colour_viridis_c() +
    scale_y_continuous(
      limits = c(1e-6, NA),
      expand = expansion(mult = c(0, 0.05)),
      breaks = scales::pretty_breaks()
    ) +
    scale_x_continuous(
      limits = c(1e-6, NA),
      expand = expansion(mult = c(0, 0.05)),
      breaks = scales::pretty_breaks()
    ) +
    # Fancy mathy axis labels
    labs(
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
  reflines_data = ref_pts,
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
  mutate(height = if_else(group == "Hucuktlis", 10, 5)) |> 
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

