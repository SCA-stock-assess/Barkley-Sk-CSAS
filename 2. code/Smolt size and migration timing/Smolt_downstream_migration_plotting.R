# Packages ----------------------------------------------------------------

pkgs <- c("here", "tidyverse", "readxl", "janitor", "cowplot")
#install.packages(pkgs)

library(here)
library(tidyverse); theme_set(theme_bw())
library(readxl)
library(cowplot)



# Load and collate RST data for the three populations ---------------------


# Start with Hucuktlis data
huc_rst <- here(
  "1. data",
  "downstream_smolt_data.xlsx"
) |> 
  read_xlsx(sheet = "Hucuktlis") |>
  rowwise() |> 
  mutate(
    count = sum(c_across(!date), na.rm = TRUE),
    d_m = format(date, "%d-%b"),
    year = format(date, "%Y"),
    stock = "Hucuktlis",
    gear = "RST",
    site = "Hucuktlis River"
  ) |> 
  ungroup()


# Great Central data
gcl_trap <-  here(
  "1. data",
  "downstream_smolt_data.xlsx"
) |> 
  read_xlsx(sheet = "GCL") |> 
  pivot_longer(
    matches("\\d{4}"),
    names_to = "year",
    values_to = "count"
  ) |> 
  mutate(
    d_m = format(Date, "%d-%b"),
    stock = "Great Central",
    gear = "trap",
    site = "Glover Pond"
  )


# Sproat data
spr_rst <- here(
  "1. data",
  "downstream_smolt_data.xlsx"
) |> 
  read_xlsx(sheet = "Sproat") |> 
  rename("count" = sockeye) |> 
  mutate(
    d_m = format(date, "%d-%b"),
    year = format(date, "%Y"),
    stock = "Sproat",
    site = "Sproat River",
    gear = "RST"
  )
  

# Collate the dataframes
ds_data <- list(huc_rst, gcl_trap, spr_rst) |> 
  map(\(x) select(x, stock, year, d_m, count, site, gear)) |> 
  list_rbind() |> 
  mutate(
    .by = c(stock, year),
    yr_ttl = sum(count, na.rm = TRUE),
    prop = count/yr_ttl,
    stock = factor(stock, levels = c("Great Central", "Sproat", "Hucuktlis")),
    julian = as.Date(paste0(d_m, "-", year), format = "%d-%b-%Y") |> 
      format("%j") |> 
      as.integer(),
    year = as.numeric(year)
  ) |> 
  filter(!is.na(count))


# Years for each stock's data in the datasets above
stk_yrs <- ds_data |> 
  distinct(stock, year) %>%
  split(.$stock) |> 
  map(\(x) pull(x, year))


# Load the data from Hyatt et al. 2019-2020 tech reports
hist_data <- here(
  "1. data",
  "Hyatt, Stiff, Rankin smolt data",
  "Smolt Sample Metadata 25.01.23.xls"
) |> 
  map2(
    c("GCL", "SPR", "HEN"),
    \(x, y) read_xls(x, sheet = y)
  ) |> 
  list_rbind() |> 
  janitor::clean_names() |> 
  mutate(
    stock = factor(
      lake,
      levels = c("GCL", "SPR", "HEN"),
      labels = levels(ds_data$stock)
    ),
    gear = tolower(gear_type),
    site = str_replace(sample_site, "Henderson", "Hucuktlis"),
    julian = as.integer(format(sample_date, "%j")),
    d_m = format(sample_date, "%d-%b"),
    count = if_else(is.na(total_catch), total_retained, total_catch)
  ) |> 
  add_count(stock, year) |> 
  filter(
    !(stock == "Great Central" & year %in% stk_yrs$`Great Central`),
    !(stock == "Sproat" & year %in% stk_yrs$Sproat),
    !(stock == "Hucuktlis" & year %in% stk_yrs$Hucuktlis),
    n > 5 # Remove years with <5 observations
  ) |> 
  mutate(
    .by = c(stock, year),
    yr_ttl = sum(count, na.rm = TRUE),
    prop = count/yr_ttl
  )


# Plot historic data quickly
hist_data |> 
  ggplot(aes(y = prop, x = as.Date(d_m, format = "%d-%b"))) +
  geom_col(aes(group = year), position = "stack") +
  facet_grid(~stock)
# Curves look plausible


# Collate both datasets together
all_data <- hist_data |> 
  select(colnames(ds_data)) |> 
  bind_rows(ds_data)



# Plot the migration timing curves ----------------------------------------


# Calculate quantiles of dates for each stock
quantiles <- all_data |> 
  add_count(stock) |> 
  summarize(
    .by = c(julian, stock, n),
    prop = sum(prop)
  ) |> 
  arrange(stock, julian) |> 
  mutate(
    .by = stock,
    cumsum = cumsum(prop),
    ttl = sum(prop),
    prop = cumsum/ttl,
    quartile = case_when(
      between(prop, 0, 0.25) ~ "q1",
      between(prop, 0.25, 0.5) ~ "q2",
      between(prop, 0.5, 0.75) ~ "q3",
      between(prop, 0.75, 1) ~ "q4",
      TRUE ~ "check"
    )
  ) |> 
  filter(cumsum > 0) |> 
  mutate(
    .by = c(stock, quartile),
    width = max(julian - min(julian))
  ) |> 
  distinct(stock, quartile, .keep_all = TRUE) |> 
  select(stock, quartile, julian, width, n) |> 
  mutate(
    start = julian,
    end = julian + width
  ) |> 
  pivot_wider(
    id_cols = c(stock, n),
    names_from = quartile,
    values_from = c(start, end)
  ) |> 
  rename(
    "min" = start_q1,
    "lwr" = start_q2,
    "median" = start_q3,
    "upr" = start_q4,
    "max" = end_q4
  ) |> 
  select(-matches("(start|end)_.*")) 


# Save summarized data for the plot 
plot_data <- all_data |> 
  mutate(
    .by = stock,
    ttl = sum(prop),
    prop = prop/ttl
  ) |> 
  complete(stock, d_m, year) |> 
  summarize(
    .by = c(d_m, stock),
    prop = sum(prop, na.rm = TRUE)
  )


# Annotations for 50% dates per CU
median_dates <- quantiles |> 
  select(stock, median) |> 
  mutate(
    d_m = format(as.Date(median, format = "%j"), "%d-%b"),
    label = paste(
      "50% date:",
      d_m
    )
  ) |> 
  # Height of the annotations = 120% of max proportion for each CU
  left_join(
    summarize(
      plot_data,
      .by = stock,
      prop = max(prop)*1.2
    )
  )


# Plot proportions as columns with density curves overlaid
(timing_p <- plot_data |> 
    ggplot(
      aes(
        x = as.Date(d_m, format = "%d-%b"),
        y = prop
      )
    ) +
    facet_wrap(
      ~stock, 
      ncol = 1,
      scales = "free_y"
    ) +
    geom_vline(
      data = median_dates,
      aes(xintercept = as.Date(d_m, format = "%d-%b")),
      colour = "red",
      lty = 2,
      linewidth = 0.6
    ) +
    geom_col(
      fill = "grey40",
      linewidth = 0.3
    ) +
    geom_density(
      data = plot_data |> 
        filter(!is.na(prop)) |> 
        mutate(prop = round(prop*100, 0)) |> 
        uncount(prop),
      aes(y = after_stat(density)),
      fill = "grey50",
      alpha = 0.5,
      bw = 5,
      colour = NA
    ) +
    geom_text(
      data = median_dates,
      aes(label = label),
      hjust = -0.05,
      vjust = 1
    ) +
    scale_y_continuous(
      labels = scales::percent,
      breaks = scales::pretty_breaks(n = 3),
      expand = expansion(mult = c(0, 0.05)),
      name = "Interannual daily average percentage of smolts counted",
    ) +
    coord_cartesian(
      xlim = c(
        as.Date("15-Mar", format = "%d-%b"),
        as.Date("20-Jun", format = "%d-%b")
      )
    ) + 
    labs(x = NULL) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
)


# Small plot showing years represented for each stock
yrs_ann <- all_data |> 
  distinct(stock, year) |> 
  ggplot(aes(y = year, x = 1)) +
  facet_wrap(
    ~stock,
    ncol = 1
  ) +
  geom_point(shape = "-", size = 5) +
  scale_y_continuous(
    position = "right",
    expand = c(0, 0)
  ) +
  # delete most plot elements to make a very simply illustration
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.line.x = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey75"),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(c(0,1,0,0), "lines")
  )


# Add year annotation panels to right side of plot
outmigration_p <- plot_grid(
  timing_p, yrs_ann, 
  align = "h",
  rel_widths = c(8,1),
  labels = c(NA, "Survey years"),
  label_size = 8,
  label_fontface = "plain",
  label_x = -0.5
)



# Save the plot
ggsave(
  outmigration_p,
  filename = here(
    "3. outputs",
    "Plots",
    "smolt_downstream_count_data.png"
  ),
  width = 6.5, 
  height = 5,
  units = "in",
  dpi = "print"
)
