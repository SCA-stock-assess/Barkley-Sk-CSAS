# Packages ----------------------------------------------------------------

pkgs <- c("tidyverse", "here", "readxl", "geomtextpath")
#install.packages(pkgs)

library(tidyverse); theme_set(theme_bw())
library(here)
library(readxl)
library(geomtextpath)



# Read data on full abundance time series ---------------------------------


# These data include reconstructed catch data from Henderson based on
# a relationship to Somass Sockeye harvest rate. To update these data,
# run the script under 2. code/Stock-recruit data preparation
abun_data0 <- here(
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


# Set colour palette for escapement and catch/stockID data
esc_pal <- "YlGnBu"
catch_pal <- "Greys"


# Create colour palette for escapement data sources
escapement_data_cols <- abun_data0 |> 
  mutate(
    escapement_data_source = factor(
      escapement_data_source,
      # Rank data sources by reliability/accuracy
      # best to worst, descending
      levels = c(
        "video counts",
        "manual transport count",
        "resistivity counter",
        "AUC survey - Clemens",
        "fence count",
        "hatchery reports",
        "fishery officer",
        "partial shoreline estimate"
      )
    )
  ) |> 
  arrange(desc(escapement_data_source)) |> 
  distinct(escapement_data_source) |> 
  drop_na() |> 
  mutate(
    .by = escapement_data_source,
    level = cur_group_id() + 0 # Add integers to level to darken the colours
  ) |> 
  mutate(escapement_colour = RColorBrewer::brewer.pal(max(level), esc_pal)[level]) |> 
  select(-level)


# Create colour palette for catch data sources
catch_data_cols <- abun_data0 |> 
  # It is really only the cannery reports that are less accurate
  # for the raw catch data
  mutate(level = if_else(catch_data_source == "cannery pack", 2, 3)) |> 
  arrange(desc(level)) |> 
  distinct(catch_data_source, level) |> 
  drop_na() |> 
  mutate(catch_colour = RColorBrewer::brewer.pal(max(level), catch_pal)[level]) |> 
  select(-level)


# Create colour palette for Stock ID sources
stockid_cols <- abun_data0 |> 
  mutate(
    stockid_source = factor(
      stockid_source,
      # Rank data sources by reliability/accuracy
      # best to worst, descending
      levels = c(
        "DNA attribution",
        "escapement",
        "MLE run timing model",
        "assumed run timing",
        "catch location"
      )
    )
  ) |> 
  arrange(desc(stockid_source)) |> 
  distinct(stockid_source) |> 
  drop_na() |> 
  mutate(
    .by = stockid_source,
    level = cur_group_id() + 1 # Add integers to level to darken the colours
  ) |> 
  mutate(stockid_colour = RColorBrewer::brewer.pal(max(level), catch_pal)[level]) |> 
  select(-level)


# Assign some colour palettes to data sources for plotting purposes
abun_data <- abun_data0 |> 
  left_join(stockid_cols) |> 
  left_join(escapement_data_cols) |> 
  left_join(catch_data_cols) |> 
  mutate(
    escapement_data_source = factor(
      escapement_data_source, 
      levels = levels(escapement_data_cols$escapement_data_source)
    ),
    catch_data_source = fct_relevel(catch_data_source, "cannery pack", after = Inf),
    stockid_source = factor(
      stockid_source, 
      levels = levels(stockid_cols$stockid_source)
    )
  )



# Make time series plot for SMU aggregate --------------------------------------


# Prepare dataframe with aggregate abundance values only
agg_abun_data <- abun_data |> 
  select(
    year,
    escapement,
    annual_Barkley_catch,
    matches("(catch|escapement)_colour"),
    matches("(catch|escapement)_data_source")
  ) |> 
  distinct(.keep_all = TRUE) |> 
  # Hacky fix to ensure aggregate Barkley catch is summed correctly within years
  add_count(year) |> 
  mutate(catch_value = annual_Barkley_catch/n) |> 
  # Rename columns to facilitate simpler pivot code
  rename("escapement_value" = escapement) |> 
  pivot_longer(
    matches("(catch|escapement)_.*"),
    names_pattern = "([[:alpha:]]*)_(.*)",
    names_to = c("name", ".value")
  ) |> 
  # Ensure catch data source information is the same within years
  group_by(year, name) |> 
  fill(colour, data_source, .direction = "updown") |> 
  ungroup() |> 
  # Ensure levels of data_source and colour appear in order in the plot
  mutate(
    data_source = factor(
      data_source,
      levels = c(
        levels(abun_data$catch_data_source), 
        levels(abun_data$escapement_data_source)
      ),
      labels = c(
        rep("catch records", 3),
        "cannery pack",
        levels(abun_data$escapement_data_source)
      )
    )
  ) |> 
  summarize(
    .by = c(year, name, colour, data_source),
    value = sum(value, na.rm = TRUE)
  ) |> 
  arrange(data_source) |> 
  mutate(colour = fct_inorder(colour)) |> 
  filter(value > 0) #remove the shoreline estimates
  

# Calculate 1970+ average annual catch, escapement, and total return
agg_abun_data |>
  #filter(year > 1969) |> 
  pivot_wider(
    id_cols = year,
    values_fn = sum,
    values_fill = 0
  ) |> 
  mutate(run = catch + escapement) |> 
  summarize(
    across(
      !year,
      .fns = list(
        "q25" = ~quantile(.x, 0.25),
        "median" = median,
        "q75" = ~quantile(.x, 0.75)
      ),
      .names = "{col}_{fn}"
    )
  ) |> 
  pivot_longer(
    everything(),
    names_sep = "_",
    names_to = c("stat", "quantile")
  ) |> 
  pivot_wider(names_from = quantile)


# Plot total returns broken out between catch and escapement
(abun_p <- agg_abun_data |> 
   ggplot(aes(x = year, y = value/1000)) +
   scale_y_continuous(
     name = "Barkley Sockeye abundance (1000s)",
     labels = scales::label_number(),
     expand = expansion(mult = c(0, 0.05))
   ) +
   labs(x = "Return year") +
   theme(
     legend.position = "inside",
     legend.position.inside = c(0.02, 0.98),
     legend.justification.inside = c(0, 1),
     legend.background = element_rect(colour = "black", fill = "white")
   )
)


# Simple version with just catch and escapement
(abun_p1 <- abun_p +
    geom_col(
      aes(fill = name),
      position = "stack",
      colour = "black",
      linewidth = 0.25
    ) +
    scale_fill_manual(values = c("grey", "black")) +
    theme(legend.title = element_blank())
)


(abun_p2 <- abun_p +
    geom_textvline(
      xintercept = 1968.5, 
      lty = 2, 
      colour = "grey50",
      label = "Great Central Lake fertilization begins",
      size = 3,
      hjust = 0.75
    ) +
    geom_col(
      aes(
        group = data_source,
        fill = colour
      ),
      position = "stack",
      colour = "black",
      linewidth = 0.25
    ) +
    scale_fill_identity(
      name = "Catch (grey) and\nescapement (colour)\ndata sources",
      guide = "legend",
      labels = levels(agg_abun_data$data_source)
    )
)


# Save plots
list(
  "simple" = abun_p1,
  "methods" = abun_p2
) |> 
  iwalk(
    \(x, idx) ggsave(
      plot = x,
      filename = here(
        "3. outputs",
        "Plots",
        paste0("Barkley-Sk_total_abundance_time-series_", idx, ".png")
      ),
      width = 7.5,
      height = 4,
      units = "in",
      dpi = "print"
    )
  )


# Plot CU-specific time series --------------------------------------------


# Annotation data frame stating when CU-specific catch data are not provided
cu_catch <- tribble(
  ~CU, ~year_start, ~year_end, ~ymin, ~ymax, ~y, ~label,
  "Hucuktlis", 1908, 1924.5, -Inf, Inf, 125, "CU-specific catch\nunavailable",
  "Hucuktlis", 1934, 1969.5, -Inf, Inf, 125, "CU-specific catch\nunavailable",
  "Great Central", 1908, 1969.5, -Inf, Inf, 600, "CU-specific catch\nunavailable",
  "Sproat", 1908, 1969.5, -Inf, Inf, 500, "CU-specific catch\nunavailable"
) |> 
  mutate(CU = factor(CU, levels = c("Great Central", "Sproat", "Hucuktlis"))) 


# Max abundances per CU
max_CU_run <- abun_data |>
  mutate(
    .by = c(CU, year),
    run = sum(catch, escapement, na.rm = TRUE)
  ) |> 
  summarize(
    .by = CU,
    max_run = max(run)
  )


# Dataframe for when fertilizer was applied to each lake
cu_fertilizer <- abun_data0 |> 
  select(year, CU, fertilized) |> 
  filter(fertilized == 1) |> 
  arrange(CU, year) |> 
  mutate(
    .by = CU,
    min_date = as.Date(paste0(year, "-01-01")),
    max_date = as.Date(paste0(year, "-12-31")),
    lag_year = lag(year, 1),
    gap = if_else(is.na(lag_year), 1, year - lag_year),
  ) |> 
  summarize(
    .by = c(CU, gap),
    min_date = min(min_date),
    max_date = max(max_date)
  ) |> 
  left_join(max_CU_run) |> 
  filter(CU != "Somass aggregate") |> 
  mutate(
    label = if_else(as.numeric(max_date - min_date) > 365, "fertilized", NA),
    across(contains("date"), \(x) as.integer(format(x, "%Y"))),
    max_date = if_else(max_date-min_date == 0, max_date + 1, max_date)
  )


# Dataframe for hatchery releases by CU
cu_releases <- abun_data0 |> 
  mutate(releases = if_else(is.na(hatchery_fry_release), hatchery_eggs_release, hatchery_fry_release)) |> 
  select(year, CU, releases) |> 
  filter(releases > 0) |> 
  left_join(max_CU_run)


# Dataframe describing fishway installation dates
fishway_dates <- tribble(
  ~CU, ~label, ~year,
  "Great Central", "Stamp Falls fishway built", 1926.5,
  "Great Central", "Major fishway upgrade", 1954.5,
  "Sproat", "Sproat Falls fishway built", 1951.5
) |> 
  mutate(CU = factor(CU, levels = c("Great Central", "Sproat", "Hucuktlis"))) 


# Prepare CU-level data for plotting
cu_abun_data <- abun_data |> 
  rename_with(.cols = c(catch, escapement), \(x) paste0(x, "_value")) |> 
  # Use Stock ID information as catch information
  mutate(
    catch_data_source = stockid_source,
    catch_colour = stockid_colour
  ) |> 
  filter(!CU == "Somass aggregate") |> 
  select(
    year, 
    CU,
    matches("(catch|escapement)_(value|data_source|colour)")
  ) |> 
  pivot_longer(
    matches("(catch|escapement)_.*"),
    names_pattern = "([[:alpha:]]*)_(.*)",
    names_to = c("name", ".value")
  ) |> 
  mutate(
    data_source = factor(
      data_source,
      levels = c(
        levels(abun_data$stockid_source), 
        levels(abun_data$escapement_data_source)
      )
    )
  ) |> 
  arrange(data_source) |> 
  mutate(colour = fct_inorder(colour)) |> 
  filter(value > 0) |> 
  droplevels()
  

# The time series plot
(cu_abun_p <- cu_abun_data |> 
    ggplot(aes(year, value/1000)) +
    facet_wrap(
      ~CU,
      strip.position = "right",
      ncol = 1,
      scales = "free_y"
    ) +
    # Add shaded area showing where catch wasn't broken out by CU
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
    # Add text to shaded years without CU-specific catch
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
    # Labeled lines showing fishway installs
    geom_textvline(
      data = fishway_dates,
      aes(
        xintercept = year,
        label = label
      ),
      lty = 2,
      size = 2,
      linewidth = 0.2,
      hjust = 0.3
    ) +
    # Add segments showing fertilization years
    geom_textsegment(
      data = filter(cu_fertilizer, !is.na(label)),
      aes(
        y = 1.02*(max_run/1000),
        yend = 1.02*(max_run/1000),
        x = min_date,
        xend = max_date,
        label = label
      ),
      remove_long = TRUE,
      linewidth = 0.2,
      size = 2
    ) +
    geom_segment(
      data = filter(cu_fertilizer, is.na(label)),
      aes(
        y = 1.02*(max_run/1000),
        yend = 1.02*(max_run/1000),
        x = min_date,
        xend = max_date,
      ),
      linewidth = 0.2
    ) +
    # Add points showing the magnitude of hatchery releases
    geom_point(
      data = cu_releases,
      aes(
        y = 0.98*(max_run/1000),
        size = releases
      ),
      shape = "|"
    ) +
    # Add the abundance data
    geom_col(
      aes(
        group = data_source,
        fill = colour
      ),
      position = "stack",
      colour = "black",
      linewidth = 0.25
    ) +
    scale_fill_identity(
      name = "Catch (grey) and\nescapement (colour)\ndata sources",
      guide = "legend",
      labels = levels(cu_abun_data$data_source)
    ) +
    scale_size_continuous(
      name = "Hatchery\nreleases",
      labels = scales::label_number(),
      breaks = seq(1e6, 9e6, length.out = 5),
      range = c(0.5,3)
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
    theme(
      legend.position = "right",
      strip.background = element_rect(colour = "black", fill = "white"),
      panel.grid.minor = element_blank()
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



# Grey scale plot with all panels shown together --------------------------


# Append aggregate abundance data to the CU-specific data
all_panels_data <- agg_abun_data |> 
  mutate(CU = "Barkley Aggregate") |> 
  bind_rows(cu_abun_data) |> 
  mutate(
    CU = factor(
      CU,
      levels = c(
        "Barkley Aggregate",
        "Great Central",
        "Sproat",
        "Hucuktlis"
      )
    ),
    name = factor(name, levels = c("catch", "escapement"))
  )


# Build grey scale plot and add annotations
(all_ts_plot <- all_panels_data |> 
  ggplot(aes(year, value/1000)) +
  facet_wrap(
    ~CU,
    strip.position = "right",
    ncol = 1,
    scales = "free_y"
  ) +
  # Add shaded area showing where catch wasn't broken out by CU
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
  # Add text to shaded years without CU-specific catch
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
  # Labeled lines showing fishway installs
  geom_textvline(
    data = fishway_dates,
    aes(
      xintercept = year,
      label = label
    ),
    lty = 2,
    size = 2,
    linewidth = 0.2,
    hjust = 0.3
  ) +
  # Add segments showing fertilization years
  geom_textsegment(
    data = filter(cu_fertilizer, !is.na(label)),
    aes(
      y = 1.02*(max_run/1000),
      yend = 1.02*(max_run/1000),
      x = min_date,
      xend = max_date,
      label = label
    ),
    remove_long = TRUE,
    linewidth = 0.2,
    size = 2
  ) +
  geom_segment(
    data = filter(cu_fertilizer, is.na(label)),
    aes(
      y = 1.02*(max_run/1000),
      yend = 1.02*(max_run/1000),
      x = min_date,
      xend = max_date,
    ),
    linewidth = 0.2
  ) +
  # Add points showing the magnitude of hatchery releases
  geom_point(
    data = cu_releases,
    aes(
      y = 0.98*(max_run/1000),
      size = releases
    ),
    shape = "|"
  ) +
  # Add the abundance data
  geom_col(
    aes(fill = name),
    position = "stack",
    colour = "black",
    linewidth = 0.25
  ) +
  scale_size_continuous(
    name = "Hatchery\nreleases",
    labels = scales::label_number(),
    breaks = seq(1e6, 9e6, length.out = 5),
    range = c(0.5,3)
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
    name = "Data type",
    values = c("grey", "black")
  ) +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank()
  )
)
  

# Save the overall grey-scale plot
all_ts_plot |> 
  ggsave(
    filename = here(
      "3. outputs",
      "Plots",
      "Barkley-Sk_all_abundance_time-series.png"
    ),
    width = 7.5,
    height = 7.5,
    units = "in",
    dpi = "print"
  )


# Make a grey-scale version with the CU-specific panels only
cu_greyscale_p <- all_ts_plot %+% 
  filter(all_panels_data,!CU == "Barkley Aggregate")
  

ggsave(
  plot = cu_greyscale_p,
  filename = here(
    "3. outputs",
    "Plots",
    "Barkley-Sk_cu_abundance_time-series_greyscale.png"
  ),
  width = 7.5,
  height = 6,
  units = "in",
  dpi = "print"
)

