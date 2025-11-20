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
  # Ensure only complete smolt years are used
  filter(
    between(
      smolt_year,
      min(return_year) - min(ttl_age-fw_age),
      max(return_year) - max(ttl_age-fw_age)
    )
  ) |> 
  summarize(
    .by = c(lake, smolt_year),
    R = sum(run)
  )


# Ocean Nino Index (ONI) and Pacific Decadal Oscillation (PDO) time series
cci <- list.files(
  here(
    "1. data",
    "California Current indices"
  ),
  pattern = "(?i)PDO|ONI",
  full.names = TRUE
) |> 
  map(read.csv) |> 
  reduce(full_join) |> 
  mutate(
    date = as.Date(date),
    month = as.character(format(date, "%b")),
    year = as.numeric(format(date, "%Y"))
  ) |> 
  filter(between(year, min(r_ts$smolt_year), max(r_ts$smolt_year)))
  
  
  
# Ocean temperature time series
temp <- read.csv(
  here(
    "1. data",
    "combined_temp5_1981_2024_region.csv"
  ),
  row.names = 1
) |> 
  pivot_longer(
    !Region,
    names_to = "y_m",
    values_to = "temp_5m",
    names_prefix = "X"
  ) |> 
  separate(
    y_m,
    into = c("year", "month"),
    sep = "_",
    convert = TRUE
  ) |> 
  filter(
    Region == "Somass",
    !is.na(temp_5m)
  ) |> 
  select(!Region) |> 
  mutate(month = month.abb[month])


# Merge covariate data sets
covariates <- full_join(cci, temp)


# Annual marine survival estimates ----------------------------------------


# Join recruit estimates with smolt year abundance data and calculate survival
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
    geo_mean = exp(mean(log(`50%`))), # Geometric mean of the annual medians
    min = min(`50%`),
    max = max(`50%`)
  )


# Look for best-related covariate windows ---------------------------------------


# How well are temperature and ONI correlated?
covariates |> 
  mutate(month_num = match(month, month.abb)) |> 
  arrange(date) |> 
  mutate(
    across(
      temp_5m,
      .fns = list(
        "lag0" = ~lag(.x, 0),
        "lag1" = ~lag(.x, 1),
        "lag2" = ~lag(.x, 2),
        "lag3" = ~lag(.x, 3),
        "lag4" = ~lag(.x, 4),
        "lag5" = ~lag(.x, 5),
        "lag6" = ~lag(.x, 6),
        "lag7" = ~lag(.x, 7),
        "lag8" = ~lag(.x, 8),
        "lag9" = ~lag(.x, 9)
      ),
      .names = "temp_{.fn}"
    )
  ) |> 
  pivot_longer(
    contains("temp_lag"),
    names_to = "lag",
    names_prefix = "temp_lag",
    values_to = "temp_lag"
  ) |> 
  ggplot(
    aes(
      x = ONI, 
      y = temp_lag, 
      colour = month_num,
      group = month_num
    )
  ) +
  facet_wrap(~lag) +
  geom_point() +
  geom_smooth(
    method = "lm", 
    se = FALSE
  ) +
  scale_colour_viridis_c()
# ONI seems to predict temperature during February-June with either a 0
# or 1 month lag. Considering this, it is probably fine to use month-aligned
# ONI and temperature as equivalent covariates.


# Simple first step is to look at average ONI for the first three months of 
# each year. Expectation is that ONI sets up marine food web state 3-6 
# months in advance of interactions with the incoming Barkley Sockeye
# cohort. 
covariates |> 
  filter(month %in% month.abb[1:3]) |> 
  pivot_longer(
    c(ONI, temp_5m, PDO),
    names_to = "covariate"
  ) |> 
  summarize(
    .by = c(year, covariate),
    value = mean(value)
  ) |> 
  left_join(
    sas_summary,
    relationship = "many-to-many"
  ) |> 
  filter(!is.na(lake)) |> 
  ggplot(aes(x = value, y = `50%`, colour = year)) +
  facet_grid(
    lake~covariate, 
    scales = "free",
    switch = "x"
  ) +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`)) +
  geom_point() +
  geom_smooth(
    method = "glm", 
    method.args = list(family = "quasibinomial"),
    colour = "red"
  ) +
  scale_colour_viridis_c(option = "mako") +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, NA),
    expand = c(0, 0)
  ) +
  labs(
    x = "Average Jan-Mar covariate value",
    y = "Estimated smolt-to-adult survival",
    colour = "Smolt\nyear"
  ) +
  theme(
    strip.placement = "outside",
    strip.background.x = element_blank()
  )
# Interesting that a clearer pattern isn't evident for Hucuktlis.
# Indicative of poorer quality estimates? Different marine migration timing?
# Different limiting factors? 


# Try various different aggregations of ONI data to look for best correlations
# with Marine survival rates for each lake.
covariate_windows <- expand_grid(
  window_size = 3:9,
  covariates = list(select(.data = covariates, date, ONI, temp_5m, PDO))
) |> 
  rowwise() |> 
  mutate(lead_months = list(0:(window_size-1))) |> 
  unnest_longer(lead_months) |> 
  rowwise() |> 
  mutate(
    covariates = list(
      mutate(
        covariates, 
        across(
          c(ONI, temp_5m, PDO),
          \(x) lead(x, n = lead_months)
        )
      )
    )
  ) |> 
  unnest(covariates) |> 
  pivot_longer(
    c(ONI, temp_5m, PDO),
    names_to = "covariate",
    values_to = "value"
  ) |> 
  pivot_wider(
    names_from = lead_months,
    values_from = value,
    names_prefix = "lead"
  ) |> 
  nest(.by = c(covariate, window_size)) |> 
  rowwise() |> 
  mutate(data = list(mutate(data, end_date = lead(date, n = window_size-1)))) |> 
  unnest(data) |> 
  rowwise() |> 
  mutate(cov_mean = mean(c_across(contains("lead")), na.rm = TRUE)) |> 
  select(covariate, start_date = date, end_date, cov_mean, window_size) |> 
  filter(!is.na(cov_mean), !is.nan(cov_mean))


# Join the marine survival data to the covariate windows
sas_cov_mods <- covariate_windows |> 
  mutate(
    start_month = as.numeric(format(start_date, "%m")),
    end_month = as.numeric(format(end_date, "%m")),
    start_year = as.numeric(format(start_date, "%Y")),
    end_year = as.numeric(format(end_date, "%Y")),
    # Ensure month windows align with the correct survival value
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
  nest(.by = c(covariate, start_month, end_month, window_size, lake)) |> 
  rowwise() |> 
  # Fit models for each CU
  mutate(
    model = list(betareg(`50%` ~ cov_mean, data = data)),
    glance = list(broom::glance(model)),
    window = paste0(month.abb[start_month], "-", month.abb[end_month]),
    p = broom::tidy(model) |> 
      filter(term == "cov_mean") |> 
      pull(p.value)
  ) |> 
  unnest_wider(glance)


# Investigate which covariate windows are best for each CU
sas_cov_mods |> 
  slice_max(
    by = c(lake, covariate),
    order_by = pseudo.r.squared,
    n = 5
  ) |> 
  print(n = 50)
# Interesting. Seems like very little relation exists between Hucuktlis survival
# rates and ONI or ocean temperature.


# Plot R squared values for all models
sas_cov_mods |> 
  mutate(
    fake_start_date = as.Date(paste0("2000-", start_month, "-01")),
    fake_end_date = if_else(
      end_month < start_month,
      as.Date(paste0("2001-", end_month, "-01")),
      as.Date(paste0("2000-", end_month, "-01"))
    )
  ) |> 
  mutate(
    .by = c(covariate, lake),
    standardized_r_squared = pseudo.r.squared/max(pseudo.r.squared)
  ) |> 
  ggplot(aes(y = fake_start_date, x = as.factor(window_size))) +
  facet_grid(
    lake~covariate
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
    x = "Number of months in\nsummary window",
    colour = expression(paste("Proportion\nof best\nmodel", R^2))
  ) +
  theme_minimal()
# Starting the windows in Feb-March seems best for ONI and temp;
# PDO relationships seem best in months following outmigration


# Add dataframes containing model predictions across ONI values
sas_cov_pred <- sas_cov_mods |> 
  rowwise() |> 
  mutate(
    pred_frame = list(
      tibble(
        cov_mean = seq(
          min(data$cov_mean), 
          max(data$cov_mean), 
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
best_mods <- sas_cov_pred |> 
  slice_max(
    by = c(lake, covariate),
    order_by = pseudo.r.squared,
    n = 5
  ) |> 
  select(covariate, lake, window, pseudo.r.squared, data, predictions)
# March - May model is the best for both GCL and SPR
# June - August for Hucuktlis but relationship is very weak


# Plot predictions versus survival values from top models for each lake
points_data <- best_mods |> 
  select(covariate, lake, window, pseudo.r.squared, data) |> 
  unnest(data)

lines_data <- best_mods |> 
  select(covariate, lake, window, pseudo.r.squared, predictions) |> 
  unnest(predictions)

label_data <- points_data |> 
  summarize(
    .by = c(covariate, lake, window, pseudo.r.squared),
    cov_mean = max(cov_mean),
    max_surv = max(`90%`)
  ) |> 
  mutate(
    .by = c(lake, covariate),
    cov_max = max(cov_mean),
    max_surv = max(max_surv)
  )

(best_ONI_models_p <- ggplot(
  data = lines_data,
  aes(x = cov_mean, y = q_0.5)
) +
    facet_grid(
      window ~ lake + covariate,
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
        x = cov_max,
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
      x = "Average covariate value across time window",
      y = "Estimated marine survival"
    )
)


#  Try a multivariate model structure -------------------------------------


# Arrange best time window indices together in a single dataframe per lake
sas_mr_mods <- sas_cov_mods |> 
  filter(
    case_when(
      covariate %in% c("ONI", "temp_5m") & window == "Mar-May" ~ TRUE,
      covariate == "PDO" & window == "Oct-Mar" ~ TRUE,
      .default = FALSE
    )
  ) |> 
  select(1:data) |> 
  unnest(data) |> 
  pivot_wider(
    id_cols = c(lake, smolt_year, `50%`),
    names_from = covariate,
    values_from = cov_mean
  ) |> 
  nest(.by = lake) |> 
  rowwise() |> 
  mutate(model = list(betareg(`50%` ~ PDO + temp_5m + ONI, data = data)))
  

# Individual relationships disappear--likely a result of multicollinearity
sas_mr_mods |> 
  pull(model, name = lake) |>
  map(summary)


# View correlations between indices
sas_mr_mods[1, 2] |> 
  unnest() |> 
  select(PDO, ONI, temp_5m) |> 
  cor(use = "complete.obs")
  

# check variance inflation factors
sas_mr_mods |> 
  pull(model, name = lake) |>
  map(car::vif)


  
# Plot marine survival versus covariates -----------------------------------------


# Time series plot
(sas_ts <- sas_summary |> 
   ggplot(aes(x = year, y = `50%`)) +
   facet_wrap(
     ~lake,
     ncol = 1,
     strip.position = "right"
   ) +
   geom_linerange(
     aes(
       ymin = `2.5%`,
       ymax = `97.5%`
     ),
     linewidth = 0.25,
     colour = "grey25"
   ) +
   geom_pointrange(
     aes(
       ymin = `10%`,
       ymax = `90%`
     ),
     size = 0.25
   ) +
   scale_y_continuous(
     limits = c(0, 0.5),
     expand = c(0, 0),
     labels = scales::percent,
     oob = scales::oob_keep
   ) +
   labs(
     x = "Ocean-entry year",
     y = "Smolt-to-adult survival"
   ) +
   theme(
     strip.background = element_rect(fill = "white"),
     panel.grid.minor = element_blank(),
     panel.spacing.y = unit(1, "lines")
   )
)


# Save the time series plot
ggsave(
  sas_ts,
  filename = here(
    "3. outputs",
    "Plots",
    "Smolt-to-adult_survival_time-series.png"
  ),
  width = 7,
  height = 5,
  units = "in",
  dpi = "print"
)


# Marine survival versus Mar-May ONI & temp
ONI_corr_p_data <- sas_cov_pred |> 
  filter(
    case_when(
      covariate == "PDO" & window == "Oct-Mar" ~ TRUE,
      covariate %in% c("ONI", "temp_5m") & window == "Mar-May" ~ TRUE,
      .default = FALSE
    )
  ) |> 
  select(covariate, lake, pseudo.r.squared, data, predictions, window, p) |> 
  {
    function(x) {
      points_data <- x |> 
        select(covariate, lake, pseudo.r.squared, data) |> 
        unnest(data) |> 
        pivot_wider(
          names_from = covariate,
          values_from = cov_mean
        ) |> 
        pivot_longer(
          cols = c(ONI, temp_5m, PDO),
          names_to = "x_var",
          values_to = "x"
        ) |> 
        mutate(
          x_var = factor(
            x_var,
            levels = c("PDO", "ONI", "temp_5m"),
            labels = c(
              "Pacific Decadal Oscillation", 
              "Ocean Ni\u00f1o Index (\u00b0C)", 
              "5-m temperature (\u00b0C)"
            )
          )
        )
      
      lines_data <- x |> 
        select(covariate, lake, pseudo.r.squared, predictions, p) |> 
        unnest(predictions) |> 
        mutate(
          x_var = factor(
            covariate,
            levels = c("PDO", "ONI", "temp_5m"),
            labels = c(
              "Pacific Decadal Oscillation", 
              "Ocean Ni\u00f1o Index (\u00b0C)", 
              "5-m temperature (\u00b0C)"
            )
          ),
          across(c(q_0.025, q_0.975), \(x) if_else(p < 0.05, x, NA)),
          alpha = if_else(p < 0.05, 1, 0.4)
        ) |> 
        rename("x" = cov_mean) 
      
      label_data <- x |> 
        select(!predictions) |> 
        unnest(data) |> 
        summarize(
          .by = c(covariate, lake, pseudo.r.squared, p),
          cov_max = max(cov_mean),
          max_surv = max(`50%`)
        ) |> 
        mutate(
          x_var = factor(
            covariate,
            levels = c("PDO", "ONI", "temp_5m"),
            labels = c(
              "Pacific Decadal Oscillation", 
              "Ocean Ni\u00f1o Index (\u00b0C)", 
              "5-m temperature (\u00b0C)"
            )
          ),
          alpha = if_else(p < 0.05, 1, 0.4)
        )
    
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
  data = ONI_corr_p_data$points,
  aes(x = x, y = `50%`)
) +
    facet_grid(
      lake ~ x_var,
      switch = "x",
      scales = "free"
    ) +
    geom_ribbon(
      data = ONI_corr_p_data$lines,
      aes(y = q_0.5, ymin = q_0.025, ymax = q_0.975),
      alpha = 0.2
    ) +
    geom_line(
      data = ONI_corr_p_data$lines,
      aes(y = q_0.5, alpha = alpha)
    ) +
    geom_point(aes(colour = smolt_year)) +
    # R^2 labels
    geom_text(
      data = ONI_corr_p_data$labels,
      aes(
        x = cov_max,
        y = max_surv*0.65,
        label = paste0("R^2==", round(pseudo.r.squared, 3)),
        alpha = alpha
      ),
      vjust = 1, 
      hjust = 1,
      parse = TRUE,
      size = 3
    ) +
    scale_colour_viridis_c(
      option = "mako",
      breaks = c(1985, 1995, 2005, 2015)
    ) +
    scale_y_continuous(
      name = "Estimated smolt-to-adult (marine) survival",
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05)),
      breaks = scales::pretty_breaks(n = 3),
      labels = scales::percent
    ) +
    scale_alpha_identity() +
    guides(
      colour = guide_colorbar(
        title = "Ocean-entry year",
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
  width = 8,
  height = 7,
  units = "in",
  dpi = "print"
)


# Model summaries
mods_list <- sas_cov_mods |> 
  filter(
    case_when(
      covariate == "PDO" & window == "Oct-Mar" ~ TRUE,
      covariate %in% c("ONI", "temp_5m") & window == "Mar-May" ~ TRUE,
      .default = FALSE
    )
  ) |> 
  mutate(name = paste(covariate, lake, window, sep = "_")) %>% 
  {split(.$model, .$name)} |> 
  list_flatten() 


# Likelihood ratio tests
mods_list |> 
  map(car::Anova)


mods_list |> 
  map(broom::tidy) |> 
  map(\(x) filter(x, term == "cov_mean")) |> 
  list_rbind(names_to = "name") |> 
  select(name, estimate, std.error) |> 
  mutate(
    lci = estimate - 1.96*std.error,
    uci = estimate + 1.96*std.error,
    across(c(estimate, lci, uci), \(x) -(1-exp(x)))
  )
