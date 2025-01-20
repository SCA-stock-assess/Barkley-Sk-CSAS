# Packages ----------------------------------------------------------------


pkgs <- c("here", "tidyverse", "readxl", "janitor", "ggridges")
#install.packages(pkgs)


library(here)
library(tidyverse); theme_set(theme_bw())
library(readxl)
library(ggridges)


# Load smolt size data from downstream trapping 1977-2024 -----------------


# Data files come from 2 sources:
# 1) Hyatt, Rankin, and Stiff 2019-2020 data reports that collated all historic
# data from ~1977-2018 (exact years vary between lakes)
#
# 2) RST sampling program started circa 2016 and run by Hupacasath


# Begin with the Hyatt et all data
smolts0 <- c("GCL", "SPR", "HEN") |> 
  purrr::set_names() |> 
  map(
    \(x) here(
      "1. data",
      "Hyatt, Stiff, Rankin smolt data",
      paste0(x, "_Smolt_Export.xlsx")
    )
  ) |> 
  imap(
    \(x, idx) read_xlsx(
      path = x,
      sheet = paste(idx, "Smolt Data (all)")
    )
  ) |> 
  list_rbind(names_to = "cu") |> 
  janitor::clean_names() |> 
  mutate(
    cu = factor(
      cu,
      levels = c("GCL", "SPR", "HEN"),
      labels = c("Great Central", "Sproat", "Hucuktlis")
    )
  )


# now the RST data
smolts1 <- here(
  "1. data",
  "smolt size data.xlsx"
) |> 
  read_xlsx(sheet = "downstream_samples") |> 
  mutate(
    year = as.integer(format(date, "%Y")),
    location = case_when(
      cu == "GCL" ~ "Glover Pond",
      cu == "SPR" ~ "Sproat Falls",
      TRUE ~ NA
      ),
    method = case_when(
      cu == "GCL" ~ "Trap box",
      cu == "SPR" ~ "RST"
    ),
    # Use rule from Hyatt et al. papers to assign ages with reasonable certainty
    # at the extreme ends of the length distribution
    fnlage = case_when(
      len_fork_mm < 70 ~ 1,
      len_fork_mm > 100 ~ 2,
      TRUE ~ NA
    ),
    cu = factor(
      cu,
      levels = c("GCL", "SPR", "HEN"),
      labels = c("Great Central", "Sproat", "Hucuktlis")
    )
  ) |> 
  filter(
    year > max(smolts0$year),
    N == 1
  ) |> 
  select(-N, -CSR) |> 
  rename(
    "sample_date" = date,
    "fork_length" = len_fork_mm,
    "fresh_std_weight" = w_g
  )


# Combine the two datasets
smolts <- list(smolts0, smolts1) |> 
  map(\(x) mutate(x, across(everything(), as.character))) |> 
  bind_rows() |> 
  mutate(across(everything(), parse_guess))


# Build annual summary tables ---------------------------------------------


# Hyatt and all have vetted annual summary tables;
# will take those at face value and model the 2018+ data accordingly



# Plots -------------------------------------------------------------------


# Begin plotting with both age classes lumped together. 
# Hyatt et al data is apportioned by age,
# RST data will be apportioned by age when some scales received from lake
# mid-water trawl samples


# Calculate overall median values for each CU and size metric
annual_medians <- smolts |> 
  pivot_longer(c(fork_length, fresh_std_weight)) |> 
  summarize(
    .by = c(cu, name),
    median = median(value, na.rm = TRUE)
  )


# Plot of smolt sizes by year & stock
(cu_ts_agg_p <- smolts |> 
  pivot_longer(c(fork_length, fresh_std_weight)) |> 
  filter(!is.na(value)) |> 
  summarize(
    .by = c(year, cu, name),
    median = median(value),
    quant = list(quantile(value, c(0.05, 0.25, 0.75, 0.95)))
  ) |> 
  unnest_wider(quant) |> 
  rename_with(
    .cols = contains("%"),
    \(x) paste0("q_", str_remove_all(x, "%"))
  ) |> 
  ggplot(aes(x = year, y = median)) +
  facet_grid(
    name ~ cu,
    scales = "free_y",
    switch = "y",
    labeller = labeller(
      name = c(
        fork_length = "Fork length (mm)",
        fresh_std_weight = "Weight (g)"
      )
    )
  ) +
  geom_hline(
    data = annual_medians,
    aes(yintercept = median),
    lty = 2,
    colour = "grey30"
  ) +
  geom_linerange(
    aes(
      ymin = q_5,
      ymax = q_95
    ),
    colour = "grey50",
    linewidth = 0.25
  ) +
  geom_pointrange(
    aes(
      ymin = q_25,
      ymax = q_75
    ),
    size = 0.15,
    linewidth = 0.5,
    fill = "white",
    shape = 21
  ) +
  labs(
    x = "Ocean entry year",
    y = NULL
  ) +
  theme(
    strip.placement = "outside,",
    strip.background.x = element_rect(fill = "white", colour = "black"),
    strip.background.y = element_blank()
  )
)


# Save the plot
ggsave(
  cu_ts_agg_p,
  filename = here(
    "3. outputs",
    "Plots",
    "Smolt_aggregate_size_timeseries_byCU.png"
  ),
  width = 6.5,
  units = "in",
  dpi = "print"
)


# Alternative version as ridgeline plot
(cu_ridge_agg_p <- smolts |> 
  # Remove some outlier data points
  filter(
    !(fork_length > 125 & cu == "SPR"),
    fresh_std_weight < 30
  ) |> 
  pivot_longer(c(fork_length, fresh_std_weight)) |> 
  filter(!is.na(value)) |> 
  ggplot(aes(y = year, x = value, group = year)) +
  facet_grid(
    cu ~ name,
    scales = "free_x",
    switch = "x",
    labeller = labeller(
      name = c(
        fork_length = "Fork length (mm)",
        fresh_std_weight = "Weight (g)"
      )
    )
  ) +
  geom_density_ridges(
    # jittered_points = TRUE,
    # position = "raincloud"
    rel_min_height = 0.001
  ) +
  scale_y_reverse() +
  labs(
    y = "Ocean entry year",
    x = NULL,
    title = "Barkley Sockeye smolt sizes from downstream trapping"
  ) +
  theme_ridges() +
  theme(
    strip.placement = "outside,",
    strip.background.y = element_rect(fill = "white", colour = "black"),
    strip.background.x = element_blank()
  )
)


# Save the plot
ggsave(
  cu_ridge_agg_p,
  filename = here(
    "3. outputs",
    "Plots",
    "Smolt_aggregate_size_ridgeline_byCU.png"
  ),
  width = 6.5,
  height = 9,
  units = "in",
  dpi = "print"
)


# Age-specific annual medians
annual_medians_age <- smolts0 |> 
  filter(fnlage %in% c(1, 2)) |> 
  pivot_longer(c(fork_length, fresh_std_weight)) |> 
  summarize(
    .by = c(cu, name, fnlage),
    median = median(value, na.rm = TRUE)
  )


# Age-specific plots
(cu_ts_age_p <- smolts0 |> 
    filter(fnlage %in% c(1, 2)) |> 
    pivot_longer(c(fork_length, fresh_std_weight)) |> 
    filter(!is.na(value)) |> 
    summarize(
      .by = c(year, cu, fnlage, name),
      median = median(value),
      quant = list(quantile(value, c(0.05, 0.25, 0.75, 0.95)))
    ) |> 
    unnest_wider(quant) |> 
    rename_with(
      .cols = contains("%"),
      \(x) paste0("q_", str_remove_all(x, "%"))
    ) |> 
    ggplot(aes(x = year, y = median, colour = factor(fnlage))) +
    facet_grid(
      name ~ cu,
      scales = "free_y",
      switch = "y",
      labeller = labeller(
        name = c(
          fork_length = "Fork length (mm)",
          fresh_std_weight = "Weight (g)"
        )
      )
    ) +
    geom_hline(
      data = annual_medians_age,
      aes(yintercept = median, colour = factor(fnlage)),
      lty = 2,
    ) +
    geom_linerange(
      aes(
        ymin = q_5,
        ymax = q_95
      ),
      alpha = 0.8,
      linewidth = 0.25
    ) +
    geom_pointrange(
      aes(
        ymin = q_25,
        ymax = q_75
      ),
      size = 0.15,
      linewidth = 0.5,
      fill = "white",
      shape = 21
    ) +
    scale_colour_manual(values = c("black", "grey70")) +
    labs(
      x = "Ocean entry year",
      y = NULL,
      colour = "Age"
    ) +
    theme(
      strip.placement = "outside,",
      strip.background.x = element_rect(fill = "white", colour = "black"),
      strip.background.y = element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.5, 0.5),
      legend.direction = "horizontal",
      legend.background = element_rect(colour = "black")
    )
)


# Save the plot
ggsave(
  cu_ts_age_p,
  filename = here(
    "3. outputs",
    "Plots",
    "Smolt_size-at-age_timeseries_byCU.png"
  ),
  width = 6.5,
  units = "in",
  dpi = "print"
)


# Age-specific ridgelines
(cu_ridge_age_p <- smolts0 |> 
    pivot_longer(c(fork_length, fresh_std_weight)) |> 
    filter(
      !is.na(value),
      fnlage %in% c(1, 2)
    ) |> 
    mutate(fnlage = factor(fnlage)) |> 
    ggplot(
      aes(
        y = year, 
        x = value, 
        fill = fnlage, 
        colour = fnlage,
        group = interaction(year, fnlage))
    ) +
    facet_grid(
      cu ~ name,
      scales = "free_x",
      switch = "x",
      labeller = labeller(
        name = c(
          fork_length = "Fork length (mm)",
          fresh_std_weight = "Weight (g)"
        )
      )
    ) +
    geom_vline(
      data = annual_medians_age,
      aes(
        xintercept = median,
        colour = factor(fnlage)
      ),
      lty = 2
    ) +
    geom_density_ridges(
      rel_min_height = 0.001,
      alpha = 0.7
    ) +
    scale_fill_viridis_d(
      name = "Age",
      end = 0.75,
      option = "mako",
      aesthetics = c("fill", "colour")
    ) +
    scale_y_reverse() +
    labs(
      y = "Ocean entry year",
      x = NULL,
      title = "Barkley Sockeye smolt sizes from downstream trapping"
    ) +
    theme_ridges() +
    theme(
      strip.placement = "outside,",
      strip.background.y = element_rect(fill = "white", colour = "black"),
      strip.background.x = element_blank(),
      strip.text.x = element_text(margin = margin(b = 1)),
      legend.position = "inside",
      legend.position.inside = c(0.45, 0.33),
      legend.justification.inside = c(0.5, 0.5),
      legend.direction = "horizontal",
      legend.background = element_rect(fill = alpha("white", 0.7))
    )
)


# Save the plot
ggsave(
  cu_ridge_age_p,
  filename = here(
    "3. outputs",
    "Plots",
    "Smolt_size-at-age_ridgeline_byCU.png"
  ),
  width = 7,
  height = 8,
  units = "in",
  dpi = "print"
)
