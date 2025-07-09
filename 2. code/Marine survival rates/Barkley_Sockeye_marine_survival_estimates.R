# Packages ----------------------------------------------------------------


pkgs <- c("tidyverse", "here")
#install.packages(pkgs)


library(here)
library(tidyverse); theme_set(theme_bw())


# Load smolt output and resulting recruitment estimates -------------------


# Load posterior values from Freshwater Life Cycle Model
posterior_df <- readRDS(
  here(
    "3. outputs",
    "Stock-recruit data",
    "Freshwater_LifeCycle_model_full-posterior.RDS"
  )
)


# Load recruitment by smolt year
r_ts <- read.csv(
  here(
    "3. outputs",
    "Stock-recruit data",
    "Barkley_Sockeye_returns_by_GR_age.csv"
  )
) |> 
  mutate(
    lake = factor(
      stock,
      levels = c("GCL", "SPR", "HUC"),
      labels = c("Great Central", "Sproat", "Hucuktlis")
    ),
    smolt_year = brood_year + fw_age,
    run = catch + escapement
  ) |> 
  summarize(
    .by = c(lake, smolt_year),
    R = sum(run)
  )


# Ocean Nino Index (ONI) time series
oni <- read.csv(
  here(
    "1. data",
    "ENSO (ONI) 250505(ONI).csv"
  )
) |> 
  mutate(
    date = as.Date(DateTime.UTC),
    start_date = date,
    end_date = lead(start_date, 1) - 1,
    oni = ONI.degrees.C
  ) |> 
  filter(
    start_date >= as.Date("1976-01-01"),
    end_date <= as.Date("2024-12-31")
  )


# Annual marine survival estimates ----------------------------------------


# Join recruit estimates with brood year smolt data and calculate survival
sas <- posterior_df |> 
  filter(parameter == "SYO") |> 
  select(year, lake, smolts = value) |> 
  nest(.by = c(lake, year), .key = "smolts") |> 
  left_join(
    r_ts, 
    by = c("year" = "smolt_year", "lake"),
    relationship = "one-to-one"
  ) |> 
  unnest(smolts) |> 
  mutate(
    survival = if_else(R/smolts > 0.999, 0.999, R/smolts), # Hacky solution
    lake = factor(lake, levels = c("Great Central", "Sproat", "Hucuktlis"))
  ) 


# Annual survival rate summary data
sas_summary <- sas |> 
  filter(!is.na(survival)) |> 
  summarize(
    .by = c(lake, year),
    surv = list(quantile(survival, probs = c(0.025, 0.1, 0.5, 0.9, 0.975)))
  ) |> 
  unnest_wider(surv) |> 
  mutate(date = as.Date(paste0(year, "-04-15")))


# Interannual summaries
sas_summary |> 
  summarize(
    .by = lake,
    median = median(`50%`),
    min = min(`50%`),
    max = max(`50%`)
  )


# Plot marine survival versus ONI -----------------------------------------


# Time series plot
(sas_ts <- sas_summary |> 
   ggplot(aes(x = date, y = `50%`)) +
   facet_wrap(
     ~lake,
     ncol = 1,
     strip.position = "right"
   ) +
   geom_rect(
     data = oni,
     aes(
       x = start_date,
       xmin = start_date -1,
       xmax = end_date +1,
       y = 0.5,
       ymin = -Inf,
       ymax = Inf,
       fill = oni
     )
   ) +
   geom_linerange(
     aes(
       ymin = `2.5%`,
       ymax = `97.5%`
     ),
     linewidth = 0.25,
     colour = "grey"
   ) +
   geom_pointrange(
     aes(
       ymin = `10%`,
       ymax = `90%`
     ),
     size = 0.25
   ) +
   scale_fill_distiller(palette = "RdYlGn") +
   scale_y_continuous(
     limits = c(0, 0.5),
     expand = c(0, 0),
     labels = scales::percent,
     oob = scales::oob_keep
   ) +
   scale_x_date(expand = c(0, 0)) +
   labs(
     x = "Brood year",
     y = "Smolt-to-adult survival",
     fill = "Ocean Ni\u00f1o\nIndex (\u00b0C)"
   ) +
   theme(
     strip.background = element_rect(fill = "white"),
     panel.grid.minor = element_blank(),
     panel.spacing.y = unit(1, "lines")
   )
)


# Plot SAS versus ONI
oni |> 
  mutate(
    month = format(date, "%b"),
    year = format(date, "%Y")
  ) |> 
  filter()


# Save the plot
ggsave(
  sas_ts,
  filename = here(
    "3. outputs",
    "Plots",
    "Smolt-to-adult_survival_with-ONI.png"
  ),
  width = 7,
  height = 5,
  units = "in",
  dpi = "print"
)


