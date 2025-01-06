# Packages ----------------------------------------------------------------

pkgs <- c("here", "tidyverse", "readxl")
#install.packages(pkgs)

library(here)
library(tidyverse); theme_set(theme_bw())
library(readxl)



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
    stock = "Hucuktlis"
  ) |> 
  ungroup()


# Great Central data
gcl_rst <-  here(
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
    stock = "Great Central"
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
    stock = "Sproat"
  )
  

# Collate the dataframes
rst_data <- list(huc_rst, gcl_rst, spr_rst) |> 
  map(\(x) select(x, stock, year, d_m, count)) |> 
  list_rbind() |> 
  mutate(
    .by = c(stock, year),
    yr_ttl = sum(count, na.rm = TRUE),
    prop = count/yr_ttl,
    stock = factor(stock, levels = c("Great Central", "Sproat", "Hucuktlis")),
    julian = as.Date(paste0(d_m, "-", year), format = "%d-%b-%Y") |> 
      format("%j") |> 
      as.integer()
  ) |> 
  filter(!is.na(count))


# Plot the migration timing curves ----------------------------------------


# Calculate bell curves for each stock
quantiles <- rst_data |> 
  add_count(stock) |> 
  summarize(
    .by = c(julian, stock, n),
    count = sum(count)
  ) |> 
  arrange(stock, julian) |> 
  mutate(
    .by = stock,
    cumsum = cumsum(count),
    ttl = sum(count),
    prop = cumsum/ttl,
    quartile = case_when(
      # Assume Sproat observations constitute only the back half of the curve
      stock == "Sproat" & between(prop, 0, 0.5) ~ "q3",
      stock == "Sproat" & between(prop, 0.5, 1) ~ "q4",
      between(prop, 0, 0.25) ~ "q1",
      between(prop, 0.25, 0.5) ~ "q2",
      between(prop, 0.5, 0.75) ~ "q3",
      between(prop, 0.75, 1) ~ "q4",
      TRUE ~ "check"
    )
  ) |> 
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
  # Infill Sproat q1 and q2 data assuming symmetry with last 2 quartiles
  mutate(
    end_q2 = if_else(is.na(end_q2), start_q3 - 1, end_q2),
    start_q2 = if_else(is.na(start_q2), end_q2 - (end_q3 - start_q3), start_q2),
    end_q1 = if_else(is.na(end_q1), start_q2 - 1, end_q1),
    start_q1 = if_else(is.na(start_q1), end_q1 - (end_q4 - start_q4), start_q1)
  ) |> 
  rename(
    "min" = start_q1,
    "lwr" = start_q2,
    "median" = start_q3,
    "upr" = start_q4,
    "max" = end_q4
  ) |> 
  select(-matches("(start|end)_.*")) 


# Use quantiles to estimate mean and standard deviation 
# then generate random draws from those distributions
sim_data <- quantiles |> 
  mutate(
    mean = median,
    sd = (upr - lwr)/IQR(rnorm(n))
  ) |> 
  select(stock, mean, sd) |> 
  rowwise() |> 
  mutate(
    data = list(
      rnorm(1e5, mean = mean, sd = sd) |> 
        as_tibble_col("julian")
    )
  ) |> 
  unnest(data) |> 
  mutate(
    d_m = round(julian, 0) |> 
      as.Date(format = "%j") |> 
      format("%d-%b")
  )


# Annotations for 50% date
median_dates <- quantiles |> 
  select(stock, median) |> 
  mutate(
    d_m = format(as.Date(median, format = "%j"), "%d-%b"),
    label = paste(
      "50% date:",
      d_m
    )
  ) |> 
  left_join(summarize(rst_data, .by = stock, prop = max(prop)))


# Ratio for transforming to the secondary y axis
y_trans <- 5


# Individual lines for each year with curves overlaid
(rst_p <- rst_data |> 
  complete(stock, d_m, year) |> 
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
    linewidth = 0.8
  ) +
  geom_line(
    aes(group = year),
    colour = "grey75"
  ) +
  geom_density(
    data = sim_data,
    aes(y = after_stat(density)*y_trans),
    fill = "grey50",
    alpha = 0.25,
    bw = 3,
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
    name = "Percentage of total smolts counted (lines)",
    sec.axis = sec_axis(transform = ~.*y_trans, name = "density (curves)")
  ) +
  coord_cartesian(
     xlim = c(
       as.Date("01-Mar", format = "%d-%b"),
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


# Save the plot
ggsave(
  rst_p,
  filename = here(
    "3. outputs",
    "Plots",
    "RST_smolt_count_data.png"
  ),
  width = 6.5, 
  height = 5,
  units = "in",
  dpi = "print"
)
