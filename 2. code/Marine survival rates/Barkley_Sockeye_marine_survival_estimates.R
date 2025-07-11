# Packages ----------------------------------------------------------------


pkgs <- c("tidyverse", "here", "broom", "betareg")
#install.packages(pkgs)


library(here)
library(tidyverse); theme_set(theme_bw())
library(betareg)


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
    oni = ONI.degrees.C,
    month = as.character(format(date, "%b")),
    year = as.numeric(format(date, "%Y"))
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
  # Center estimate dates on 15 April for continuous date axes
  mutate(date = as.Date(paste0(year, "-04-15")))


# Inter-annual summaries by CU
sas_summary |> 
  summarize(
    .by = lake,
    median = median(`50%`),
    min = min(`50%`),
    max = max(`50%`)
  )


# Look for best-related ONI windows ---------------------------------------


# Simple first step is to look at average ONI for the first three months of 
# each year. Expectation is that ONI sets up marine food web state 3-6 
# months in advance of interactions with the incoming Barkley Sockeye
# cohort. 
oni |> 
  filter(month %in% month.abb[1:3]) |> 
  summarize(
    .by = year,
    oni = mean(oni)
  ) |> 
  left_join(
    sas_summary,
    relationship = "one-to-many"
  ) |> 
  filter(!is.na(lake)) |> 
  ggplot(aes(x = oni, y = `50%`, colour = year)) +
  facet_wrap(~lake, nrow = 1) +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`)) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "quasibinomial")) +
  scale_colour_viridis_c(option = "mako") +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, NA),
    expand = c(0, 0)
  ) +
  labs(
    x = "Average Jan-Mar ONI anomaly",
    y = "Estimated smolt-to-adult survival",
    colour = "Smolt\nyear"
  )
# Interesting that a clearer pattern isn't evident for Hucuktlis.
# Indicative of poorer quality estimates? Different marine migration timing?
# Different limiting factors?  


# Try various different aggregations of ONI data to look for best correlations
# with Marine survival rates for each lake.
oni_windows <- expand_grid(
  window_size = 3:9,
  oni_data = list(select(.data = oni, date, oni))
) |> 
  rowwise() |> 
  mutate(lead_months = list(0:(window_size-1))) |> 
  unnest_longer(lead_months) |> 
  rowwise() |> 
  mutate(oni_data = list(mutate(oni_data, oni = lead(oni, n = lead_months)))) |> 
  unnest(oni_data) |> 
  pivot_wider(
    names_from = lead_months,
    values_from = oni,
    names_prefix = "oni_lead"
  ) |> 
  nest(.by = window_size) |> 
  rowwise() |> 
  mutate(data = list(mutate(data, end_date = lead(date, n = window_size-1)))) |> 
  unnest(data) |> 
  rowwise() |> 
  mutate(mean_oni = mean(c_across(contains("oni_lead")), na.rm = TRUE)) |> 
  select(start_date = date, end_date, mean_oni, window_size)


# Join the marine survival data to the ONI windows
sas_oni_mods <- oni_windows |> 
  mutate(
    start_month = as.numeric(format(start_date, "%m")),
    end_month = as.numeric(format(end_date, "%m")),
    start_year = as.numeric(format(start_date, "%Y")),
    end_year = as.numeric(format(end_date, "%Y")),
    smolt_year = case_when(
      end_month %in% 1:10 ~ end_year,
      end_month %in% 11:12 & window_size %in% 3:5 ~ start_year,
      .default = end_year
    )
  ) |> 
  left_join(
    sas_summary, 
    by = c("smolt_year" = "year"),
    # Three rows of marine survival data per year; one for each CU
    relationship = "many-to-many"
  ) |> 
  filter(!is.na(lake)) |> # Remove years without marine survival estimates
  ungroup() |> 
  # Compartmentalize data frames by variables of interest
  nest(.by = c(start_month, end_month, window_size, lake)) |> 
  rowwise() |> 
  # Fit models for each CU
  mutate(
    model = list(betareg(`50%` ~ mean_oni, data = data)),
    glance = list(broom::glance(model)),
    window = paste0(month.abb[start_month], "-", month.abb[end_month])
  ) |> 
  unnest_wider(glance)


# Investigate which ONI windows are best for each CU
sas_oni_mods |> 
  slice_min(
    by = lake,
    order_by = AIC,
    n = 5
  )
# Interesting. Seems like very little relation exists between Hucuktlis survival
# rates and ONI.


# Plot R squared values for all models
sas_oni_mods |> 
  mutate(
    fake_start_date = as.Date(paste0("2000-", start_month, "-01")),
    fake_end_date = if_else(
      end_month < start_month,
      as.Date(paste0("2001-", end_month, "-01")),
      as.Date(paste0("2000-", end_month, "-01"))
    )
  ) |> 
  mutate(
    .by = lake,
    standardized_r_squared = pseudo.r.squared/max(pseudo.r.squared)
  ) |> 
  ggplot(aes(y = fake_start_date, x = as.factor(window_size))) +
  facet_wrap(
    ~lake,
    ncol = 1,
    strip.position = "right"
  ) +
  coord_flip() +
  annotate(
    geom = "rect",
    xmin = -Inf, 
    xmax = Inf,
    ymin = as.Date("2000-04-01"), 
    ymax = as.Date("2000-06-15"),
    fill = "grey50",
    colour = "black",
    lty = 2,
    alpha = 0.5
  ) +
  annotate(
    geom = "rect",
      xmin = -Inf, 
      xmax = Inf,
      ymin = as.Date("2001-04-01"), 
      ymax = as.Date("2001-06-15"),
    fill = "grey50",
    colour = "black",
    lty = 2,
    alpha = 0.5
  ) +
  geom_linerange(
    aes(
      ymin = fake_start_date, 
      ymax = fake_end_date,
      colour = standardized_r_squared,
      group = fct_reorder(window, start_month, .desc = TRUE)
    ),
    position = position_dodge(width = 1),
    linewidth = 1
  ) +
  scale_colour_viridis_c(option = "mako", direction = -1) +
  scale_y_date(date_labels = "%b") +
  labs(
    y = "Month",
    x = "Number of months in\nONI summary window",
    colour = expression(paste("Proportion\nof best\nmodel", R^2))
  ) +
  theme_minimal()
# Starting the windows in January seems best for GCL and SPR


# Add dataframes containing model predictions across ONI values
sas_oni_pred <- sas_oni_mods |> 
  rowwise() |> 
  mutate(
    pred_frame = list(
      tibble(
        mean_oni = seq(
          min(data$mean_oni), 
          max(data$mean_oni), 
          length.out = 100
        )
      )
    ),
    predictions = list(
      predict(
        model, 
        newdata = pred_frame,
        type = "quantile",
        at = c(0.025, 0.5, 0.975)
      ) |> 
        cbind(pred_frame)
    )
  ) |> 
  ungroup()


# Save a list of the top 5 models for each lake
best_mods <- sas_oni_pred |> 
  slice_min(
    by = lake,
    order_by = AIC,
    n = 5
  ) |> 
  select(lake, window, pseudo.r.squared, data, predictions)
# February - April model is the best for both GCL and SPR
# June - August for Hucuktlis but relationship is very weak


# Plot predictions versus survival values from top models for each lake
points_data <- best_mods |> 
  select(lake, window, pseudo.r.squared, data) |> 
  unnest(data)

lines_data <- best_mods |> 
  select(lake, window, pseudo.r.squared, predictions) |> 
  unnest(predictions)

label_data <- points_data |> 
  summarize(
    .by = c(lake, window, pseudo.r.squared),
    max_oni = max(mean_oni),
    max_surv = max(`90%`)
  ) |> 
  mutate(
    .by = lake,
    max_oni = max(max_oni),
    max_surv = max(max_surv)
  )

(best_oni_models_p <- ggplot(
  data = lines_data,
  aes(x = mean_oni, y = q_0.5)
) +
    facet_grid(
      window ~ lake,
      scales = "free"
    ) +
    geom_ribbon(
      aes(ymin = q_0.025, ymax = q_0.975),
      alpha = 0.3
    ) +
    geom_line(colour = "blue") +
    geom_linerange(
      data = points_data,
      aes(ymin = `10%`, ymax = `90%`, y = `50%`),
      alpha = 0.6
    ) +
    geom_point(
      data = points_data,
      aes(y = `50%`),
      alpha = 0.6
    ) +
    geom_text(
      data = label_data,
      aes(
        x = max_oni,
        y = max_surv,
        label = paste0("R^2==", round(pseudo.r.squared, 3))
      ),
      vjust = 1, 
      hjust = 1,
      parse = TRUE,
      size = 3
    ) +
    scale_y_continuous(
      limits = c(0, NA),
      labels = scales::percent,
      expand = expansion(mult = c(0, 0.05)),
      breaks = scales::pretty_breaks(n = 3)
    ) +
    labs(
      x = "Average ONI across time window",
      y = "Estimated marine survival"
    )
)

  
# Plot marine survival versus ONI -----------------------------------------


# Time series plot
(sas_ts_oni <- sas_summary |> 
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


# Save the time series plot
ggsave(
  sas_ts_oni,
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


# Regular time series of marine survival
(sas_ts <- sas_summary |> 
  ggplot(aes(x = date, y = `50%`)) +
  facet_wrap(
    ~lake,
    ncol = 1,
    strip.position = "right",
    scales = "free_y"
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
  scale_y_continuous(
    limits = c(0, NA),
    expand = c(0, 0),
    breaks = scales::pretty_breaks(n = 3),
    labels = scales::percent
  ) +
  labs(
    x = "Brood year",
    y = "Smolt-to-adult (marine) survival"
  ) +
  theme(
    strip.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank(),
    panel.spacing.y = unit(1, "lines")
  )
)


# Marine survival versus Feb-Apr ONI
oni_corr_p_data <- sas_oni_pred |> 
  filter(window == "Feb-Apr") |> 
  select(lake, pseudo.r.squared, data, predictions) |> 
  {
    function(x) {
      points_data <- x |> 
        select(lake, pseudo.r.squared, data) |> 
        unnest(data) |> 
        pivot_longer(
          cols = c(mean_oni, smolt_year),
          names_to = "x_var",
          values_to = "x"
        ) |> 
        mutate(
          x_var = factor(
            x_var,
            levels = c("smolt_year", "mean_oni"),
            labels = c("Ocean entry year", "Ocean Ni\u00f1o Index (\u00b0C)")
          ),
          smolt_year = if_else(
            x_var == "Ocean entry year",
            min(end_year),
            end_year
          )
        )
      
      lines_data <- x |> 
        select(lake, pseudo.r.squared, predictions) |> 
        unnest(predictions) |> 
        mutate(x_var = "Ocean Ni\u00f1o Index (\u00b0C)") |> 
        rename("x" = mean_oni)
      
      label_data <- x |> 
        select(lake, pseudo.r.squared, data) |> 
        unnest(data) |> 
        summarize(
          .by = c(lake, pseudo.r.squared),
          max_oni = max(mean_oni),
          max_surv = max(`50%`),
          x_var = "Ocean Ni\u00f1o Index (\u00b0C)"
        ) |> 
        # Bump the y value down a little for GCL
        mutate(max_surv = max(max_surv)*0.8)
      
      return(
        list(
          "points" = points_data, 
          "lines" = lines_data,
          "labels" = label_data
          )
        )
    }
  }()


(sas_corr_ts_p <- ggplot(
  data = oni_corr_p_data$points,
  aes(x = x, y = `50%`)
) +
    facet_grid(
      lake ~ x_var,
      switch = "x",
      scales = "free_x"
    ) +
    geom_ribbon(
      data = oni_corr_p_data$lines,
      aes(y = q_0.5, ymin = q_0.025, ymax = q_0.975),
      alpha = 0.2
    ) +
    geom_line(
      data = oni_corr_p_data$lines,
      aes(y = q_0.5)
    ) +
    geom_pointrange(
      aes(ymin = `10%`, ymax = `90%`, colour = smolt_year),
      size = 0.3
    ) +
    geom_text(
      data = oni_corr_p_data$labels,
      aes(
        x = max_oni,
        y = max_surv,
        label = paste0("R^2==", round(pseudo.r.squared, 3))
      ),
      vjust = 1, 
      hjust = 1,
      parse = TRUE,
      size = 3
    ) +
    scale_colour_viridis_c(option = "mako") +
    scale_y_continuous(
      name = "Estimated smolt-to-adult (marine) survival",
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05)),
      breaks = scales::pretty_breaks(n = 3),
      labels = scales::percent
    ) +
    guides(
      colour = guide_colorbar(
        title = "Ocean entry year",
        theme = theme(
          legend.direction = "horizontal",
          legend.title.position = "top",
          legend.text.position = "bottom",
          legend.key.width = unit(dev.size()[1]/4, "inches")
        )
      )
    ) +
    theme(
      strip.placement.x = "outside",
      strip.text.x = element_text(size = 11),
      strip.background.x = element_blank(),
      axis.title.x = element_blank(),
      strip.background.y = element_rect(fill = "white"),
      panel.grid.minor = element_blank(),
      legend.background = element_rect(
        colour = "grey20", 
        fill = alpha("white", 0.8),
        linewidth = 0.25
      ),
      legend.position = "inside",
      legend.position.inside = c(1, 1),
      legend.justification.inside = c(1, 1)
    )
)


ggsave(
  plot = sas_corr_ts_p,
  filename = here(
    "3. outputs",
    "Plots",
    "Marine_survival_time-series_ONI-correlation.png"
  ),
  width = 7,
  height = 7,
  units = "in",
  dpi = "print"
)
