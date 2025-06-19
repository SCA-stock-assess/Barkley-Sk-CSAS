# Packages ----------------------------------------------------------------


pkgs <- c("tidyverse", "here", "readxl")
#install.packages(pkgs)


library(here)
library(tidyverse); theme_set(theme_bw())
library(readxl)

# Load common functions
source(here("2. code", "0. functions", "common_functions.R"))


# Load post-season run size data for Somass and Hucuktlis -----------------


# Cleaned stock-recruit data
sr_data <- here(
  "3. outputs",
  "Stock-recruit data",
  "Barkley_Sockeye_stock-recruit_infilled.xlsx"
) |> 
  read_xlsx(sheet = "S-R data")


# Run size data and harvest rates
hr_data <- sr_data |> 
  rowwise() |> 
  mutate(
    adult_prop = if_else(
      if_any(matches("N.age.\\d"), is.na),
      1,
      sum(c_across(matches("N.age.\\d")))
    ),
    run = N * adult_prop,
    catch = H * adult_prop,
    group = if_else(stock == "HUC", "Hucuktlis", "Somass")
  ) |> 
  ungroup() |> 
  summarize(
    .by = c(year, group),
    across(c(run, catch), sum)
  ) |> 
  rowwise() |> 
  mutate(
    hr = catch/run,
    target_hr = case_when(
      group == "Somass" ~ somass_mgt_rule(run),
      group == "Hucuktlis" & year > 2011 ~ hucuktlis_mgt_rule(run),
      .default = NA
    )
  )



# Plots to evaluate harvest performance -----------------------------------


# Data frame constraining years to only those when the management plan applies
hr_data_trim <- hr_data |> 
  filter(year > 2000)  |> 
  pivot_longer(contains("hr")) |> 
  mutate(
    .by = group,
    sec_axis_y = value*max(run),
    name = if_else(
      name == "hr",
      "Observed HR",
      "Target HR"
    )
  ) 


# Set up a hacky second axis to annotate the harvest rate percentages
sec_axis_txt <- hr_data_trim |> 
  mutate(year = (max(year)-min(year))*1.05 + min(year)) |> 
  summarize(
    .by = c(group, year),
    max = max(run)
  ) |> 
  expand_grid(label = seq(0, 1, by = 0.25)) |> 
  mutate(
    run = max*label,
    label = scales::percent(label)
  )


# Plot of target versus actual HR over observed run sizes
(hr_vs_abun_p <- hr_data_trim |> 
  ggplot(aes(x = year, y = run)) +
  facet_wrap(
    ~group,
    ncol = 1,
    scales = "free_y"
  ) +
  geom_col(fill = "grey", position = position_identity()) +
  geom_line(
    aes(
      y = sec_axis_y,
      colour = name
      ),
    key_glyph = "timeseries"
  ) +
  geom_text(
    data = sec_axis_txt,
    aes(label = label),
    hjust = -0.2,
    vjust = 0.25,
    size = 3
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05)),
    labels = ~./1000,
    sec.axis = dup_axis(name = "Harvest rate")
  ) +
  scale_colour_manual(values = c("red", "black")) +
  coord_cartesian(
    xlim = c(min(hr_data_trim$year), max(hr_data_trim$year)),
    clip = "off"
  ) +
  labs(
    y = "Run size (1000s)",
    x = "Return year"
  ) +
  theme(
    panel.spacing.y = unit(1, "lines"),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.05, 0.95),
    legend.justification.inside = c(0, 1),
    legend.background = element_rect(colour = "black"),
    plot.margin = unit(c(1,3,1,1), "lines"),
    axis.ticks.y = element_blank(),
    axis.text.y.right = element_blank(),
    axis.title.y.right = element_text(vjust = 13),
    strip.background = element_blank()
  )
)


# Export the plot
ggsave(
  plot = hr_vs_abun_p,
  filename = here(
    "3. outputs",
    "Plots",
    "Hucuktlis_vs_Somass_HR_w-abundance.png"
  ),
  width = 6.5,
  height = 5,
  units = "in",
  dpi = "print"
)


# Report on annual exceedances for Hucuktlis ------------------------------


# Produce a simple table summary
hr_data |> 
  ungroup() |> 
  mutate(error = hr - target_hr) |> 
  filter(
    # Include only years when current MGT plans were in effect
    case_when(
      group == "Somass" & year > 1993 ~ TRUE,
      group == "Hucuktlis" & year > 2011 ~ TRUE,
      .default = FALSE
    )
  ) |>
  add_count(group, name = "ttl_years") |> 
  filter(error > 0) |> # Focus only on exceeded years 
  summarize(
    .by = c(group, ttl_years),
    av_err = mean(error, na.rm = TRUE),
    median_err = median(error, na.rm = TRUE),
    years_exceeded = n()
  ) |> 
  mutate(pct_exceeded = years_exceeded/ttl_years)


