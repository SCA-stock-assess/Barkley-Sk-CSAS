# Packages ----------------------------------------------------------------

pkgs <- c("tidyverse", "here", "readxl", "ggridges", "janitor")
#install.packages(pkgs)


library(here)
library(readxl)
library(ggridges)
library(tidyverse); theme_set(theme_ridges())


# Load historical Sockeye escapement data -----------------------------


# Somass escapement data (historic 1943-1992)
som_hist <- read_xls(
  here(
    "1. data",
    "1943-1992 GCLESCWL.xls"
  ),
  sheet = "GCLESCWL"
) |> 
  janitor::clean_names() |> 
  mutate(
    year = as.integer(format(date, "%Y")),
    d_m = paste0(
      format(date, "%d"),
      "-",
      format(date, "%b")
    ),
    system = "GCL",
    method = "fence count"
  ) |> 
  rename("sockeye" = gc_sockeye) |> 
  select(year, d_m, sockeye, system, method)


# Somass escapement data (1975+)
som_esc <- read_xlsx(
  here(
    "1. data",
    "Somass_Sockeye_daily_escapement_data.xlsx"
  ),
  na = ""
) |> 
  rename(
    "adj_adults" = "Adjusted Adults.  Includes bypass since 2004 but not Biosamples.",
    "adj_jacks" = "Adjusted Jacks.  Includes bypass since 2004 but not Biosamples."
  ) |> 
  mutate(
    year = if_else(
      year < 2000,
      as.numeric(paste0(19, year)),
      year
    ),
    d_m = paste0(day, "-", month.abb[month]),
    adj_adults = case_when(
      # Preserve NAs for days not yet observed in current year
      year == max(year) & is.na(`Adjusted net Adult up count`) ~ NA_real_, 
      is.na(adj_adults) & !is.na(`Stamp Falls Adjusted Adults`) ~ `Stamp Falls Adjusted Adults`,
      is.na(adj_adults) & year < max(year) ~ 0,
      TRUE ~ adj_adults),
    adj_jacks = case_when(
      is.na(adj_jacks) ~ `Adjusted net Jack up count`,
      is.na(adj_jacks) & year < max(year) ~ 0,
      TRUE ~ adj_jacks),
    across(
      c(adj_adults, adj_jacks),
      \(x) if_else(x < 0, 0, x)
    )
  ) |> 
  select(d_m, system, year, contains("adj_")) |> 
  rowwise() |> 
  mutate(
    sockeye = sum(c_across(c(adj_adults, adj_jacks)), na.rm = TRUE),
    method = "video analysis",
    .keep = "unused"
  ) 


# Hucuktlis data
huc_esc <- read_xlsx(
  here(
    "1. data",
    "HucuktlisRiver_fence-counts_1997-2024.xlsx"
  )
) |> 
  pivot_longer(
    matches("\\d{4}"),
    names_to = "year",
    values_to = "sockeye"
  ) |> 
  mutate(
    d_m = paste0(
      format(Date, "%d"),
      "-",
      format(Date, "%b")
    ),
    year = as.integer(year),
    system = "HUC",
    method = "fence count",
    .keep = "unused"
  )


# Collate all escapement data
escday <- som_hist |> 
  filter(year < min(som_esc$year)) |> 
  bind_rows(
    huc_esc,
    som_esc
  ) |> 
  mutate(system = factor(system, levels = c("GCL", "SPR", "HUC")))


# Ridgeline plots --------------------------------------------------------


# Plot using calculated stats
base_p1 <- escday |> 
   filter(!is.na(sockeye)) |> 
   mutate(sockeye = round(sockeye, 0)) |> 
   uncount(sockeye) |> 
   # Subsample rows down to 500000 random draws
   # Useful when toying with plot parameters
   slice_sample(n = 5e5) |> # deactivate when saving full version of plot
   ggplot(
     aes(
       y = year, 
       x = as.Date(d_m, format = "%d-%b"),
       group = year
     )
   ) +
   geom_density_ridges(
     aes(fill = method, colour = method),
     quantile_lines = TRUE,
     quantiles = 2,
     alpha = 0.6,
     vline_colour = "red",
     rel_min_height = 1e-4
   ) 


# Ridgeline plot where values are mapped to curve height directly
base_p2 <- escday |> 
    mutate(
      .by = c(system, year),
      ttl = sum(sockeye, na.rm = TRUE),
      prop = sockeye/ttl
    ) |> 
    ggplot(
      aes(
        y = year,
        x = as.Date(d_m, format = "%d-%b"),
        height = prop,
        group = year
      )
    ) +
    geom_density_ridges(
      stat = "identity",
      aes(fill = method, colour = method),
      alpha = 0.6,
      scale = 1.5
    ) 


# Apply consistent formatting to both plots
ridge_plots <- list(
  "density_ridges" = base_p1, 
  "proportion_ridges" = base_p2
) |> 
  map(
    \(x) x +
      facet_wrap(
        ~system,
        labeller = labeller(
          system = c(
            "GCL" = "Stamp Falls fishway",
            "SPR" = "Sproat Falls fishway",
            "HUC" = "Hucuktlis River fence"
          )
        )
      ) +
      scale_x_date(
        expand = c(0, 0), 
        date_breaks = "1 month",
        date_labels = "%b"
      ) +
      scale_y_reverse(
        breaks = scales::breaks_pretty(n = 10),
        expand = c(0, 0)
      ) +
      scale_fill_manual(
        values = c("black", "grey70"),
        aesthetics = c("fill", "colour"),
        name = "Counting\nmethod"
      ) +
      coord_cartesian(clip = "off") +
      labs(
        x = "Date", 
        y = "Adult return year",
        title = "Sockeye counts through local enumeration sites"
      ) +
      theme(
        legend.position = "inside",
        legend.position.inside = c(0.98, 0.98),
        legend.justification.inside = c(1, 1),
        legend.box.background = element_rect(colour = "black", fill = alpha("white", 0.4)),
        strip.background = element_blank(),
        strip.text = element_text(margin = margin(t = 0.5, b = 0.3, unit = "lines")),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
  )


# Save the plots
ridge_plots |> 
  iwalk(
    \(x, idx) ggsave(
      x,
      filename = here(
        "3. outputs",
        "Plots",
        paste0(
          "Barkley_Sockeye_escapement_",
          idx,
          ".png"
        )
      ),
      width = 6.5,
      height = 6.5,
      units = "in",
      dpi = "print"
    )
  )
