# Packages ----------------------------------------------------------------


pkgs <- c(
  "here", "tidyverse", "janitor", "readxl", "geom_textpath", 
  "brms", "rstan"
)
#install.packages(pkgs)


library(here)
library(tidyverse); theme_set(theme_bw())
library(geomtextpath)
library(readxl)
library(brms) # For GAM fit
library(rstan) # For state-space GAM
library(splines) # base package but must be loaded to use bs()


# Function to calculate moving average of a vector
get_mav <- function(vector, n) {

  list <- list()
  
  for(i in 1:n) {
    list[[i]] <- lag(vector, i)
  }
  
  do.call("cbind", list) |> 
    rowMeans()
  
}


# Load supporting data files ----------------------------------------------


# CU names
cu_order <- c("GCL", "SPR", "HUC")
cu_names <- c("Great Central", "Sproat", "Hucuktlis")


# Time series of smolt production estimates from Hyatt's work
hyatt_data <- here(
  "1. data",
  "Barkley-Sk_pre-smolt_abundances.xlsx"
) |> 
  read_xlsx(sheet = "smolt_production") |> 
  # Trim all post-2007 data that appear unreliable
  filter(smolt_year < 2008) |> 
  # Following steps pivot the data and calculate production by brood year
  rename("huc_1" = huc_ttl) |> # assume all Hucuktlis fry smolt at age 1
  select(-contains("ttl")) |> 
  pivot_longer(
    matches("(gcl|spr|huc)_\\d+"),
    names_sep = "_",
    names_to = c("stock", "age"),
    values_to = "sockeye"
  ) |> 
  mutate(
    age = as.numeric(age),
    brood_year = smolt_year - age - 1,
    stock = factor(toupper(stock), levels = cu_order, labels = cu_names)
  )


# Load raw ATS survey abundance estimates
ats_est <- here(
  "1. data",
  "Barkley-Sk_pre-smolt_abundances.xlsx"
) |> 
  read_xlsx(sheet = "ATS_estimates") |> 
  arrange(lake, survey_date) |> 
  group_by(lake) |> 
  # Infill missing uncertainty values with 5-y moving average CV
  mutate(
    cv = presmolt_sd/presmolt_est,
    mav = get_mav(cv, 5),
    presmolt_lwr = presmolt_est - presmolt_sd,
    presmolt_upr = presmolt_est + presmolt_sd
  ) |> 
  fill(mav, .direction = "down") |> 
  mutate(
    survey_date = as.Date(survey_date),
    cv = if_else(is.na(cv), mav, cv),
    presmolt_sd = if_else(is.na(presmolt_sd), cv*presmolt_est, presmolt_sd),
    lake = factor(lake, levels = cu_order, labels = cu_names),
    month = format(survey_date, "%b"),
    season = case_when(
      month %in% month.abb[1:3] ~ "winter",
      month %in% month.abb[4:6] ~ "spring",
      month %in% month.abb[7:9] ~ "summer",
      month %in% month.abb[10:12] ~ "fall"
    ),
    # smolt_year, as calculated below, won't be true of non-smolting 
    # age-1s (aka age_2x), but is needed to match up surveys with the 
    # corresponding smolt biosamples
    smolt_year = if_else(season == "winter", survey_year, survey_year + 1),
    smolt_year_f = factor(smolt_year),
    # Create variable that sets survey date in relation to smolt year
    # i.e. "day of smolt year"
    day_smolt_yr = as.numeric(survey_date - as.Date(paste0(smolt_year-1, "-03-20")))
  ) |> 
  ungroup() |> 
  # Remove years with no data
  filter(!is.na(presmolt_est)) |> 
  select(-mav)


# Total adult returns by age for Somass CUs
Som_returns <- here(
  "1. data",
  "Barkley Sockeye time series of returns.xlsx"
) |> 
  read_xlsx(sheet = "Somass") |> 
  filter(!if_any(c(escapement, catch, ttl_age, fw_age), is.na)) |> 
  mutate(return = as.numeric(escapement) + catch) |> 
  select(stock, return_year, return, ttl_age, fw_age)


# Lazy solution for missing 2009-2010 Hucuktlis catch estimates:
# grab infilled values from formatted s-r data output
Huc_catch_infill <- here(
  "3. outputs",
  "Stock-recruit data",
  "Barkley_Sockeye_stock-recruit_infilled.xlsx"
) |> 
  read_xlsx(sheet = "S-R data") |> 
  filter(stock == "HUC") |> 
  select(year, H) |> 
  rename("catch" = H)


# Total adult returns by age for Hucuktlis
Huc_returns <- here(
  "1. data",
  "Barkley Sockeye time series of returns.xlsx"
) |> 
  read_xlsx(sheet = "Hucuktlis") |> 
  # Fill missing catch data using lookup
  left_join(
    Huc_catch_infill,
    by = "year"
  ) |> 
  select(-catch.x) |> 
  rename("catch" = catch.y) |> 
  filter(!if_any(c(escapement, catch, matches("^age_\\d+")), is.na)) |>
  mutate(
    return = escapement + catch,
    across(matches("^age_\\d+"), \(x) x*return)
  ) |> 
  pivot_longer(
    matches("^age_\\d+"),
    names_prefix = "age_",
    names_to = "age",
    values_to = "sockeye"
  ) |> 
  mutate(
    ttl_age = as.numeric(str_extract(age, "\\d")),
    fw_age = as.numeric(str_sub(age, 2, 2)),
    stock = "HUC"
  ) |> 
  select(year, stock, sockeye, ttl_age, fw_age) |> 
  rename(
    "return_year" = year,
    "return" = sockeye
  )


# Bind Somass and Hucuktlis adult return data
adult_returns <- bind_rows(
  Som_returns,
  Huc_returns
) |> 
  mutate(
    brood_year = return_year - ttl_age,
    stock = factor(stock, levels = cu_order, labels = cu_names)
  )


# Best estimates of annual smolt size and age composition
smolt_comps <- c("GCL", "SPR", "HEN") |> 
  paste0("_Smolt_Export.xlsx") |> 
  set_names(cu_names) |> 
  map(
    \(x) read_xlsx(
      path = here("1. data", "Hyatt, Stiff, Rankin smolt data", x),
      sheet = "Weighted Size (Ages Combined)"
    )
  ) |> 
  map(janitor::clean_names) |>
  list_rbind(names_to = "lake") |> 
  select(lake, oey, matches("age\\d_(fork|pct|std_wt)$")) |> 
  rename("smolt_year" = oey)


# Annual smolt survey sample sizes
smolt_catch <- here(
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
      labels = cu_names
    ),
    gear = tolower(gear_type),
    julian = as.integer(format(sample_date, "%j")),
    d_m = format(sample_date, "%d-%b"),
    count = if_else(
      is.na(total_catch) | (total_catch < total_retained),
      total_retained, 
      total_catch
    )
  ) |> 
  add_count(stock, year) |> 
  mutate(
    .by = c(stock, year),
    yr_ttl = sum(count, na.rm = TRUE)
  )


# Smolt compositions by date and gear
smolt_ages <- c("GCL", "SPR", "HEN") |> 
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
      labels = cu_names
    )
  ) |> 
  count(cu, sample_date, fnlage, gear) |> 
  filter(fnlage %in% c(1, 2)) |> 
  pivot_wider(
    names_from = fnlage,
    names_prefix = "age",
    values_from = n,
    values_fill = 0
  ) |> 
  mutate(
    ttl = age1 + age2,
    gear = tolower(gear),
    across(c(age1, age2), \(x) x/ttl, .names = "prop_{.col}")
  )
  

# Load smolt age data received during summer 2025
new_smolt_ages <- here(
  "1. data",
  "smolt size data.xlsx"
) |> 
  read_xlsx(sheet = "post-2016_rates_of_capture") |> 
  mutate(sample_date = as.Date(sample_date))


# Investigate brood year adult vs smolt production ----------------------


# Sum adults produced per fw age and CU by brood year
adult_by_prod <- adult_returns |> 
  summarize(
    .by = c(fw_age, stock, brood_year),
    n = sum(return)
  ) |> 
  # Calculate freshwater age composition
  mutate(
    .by = c(stock, brood_year),
    by_ttl = sum(n),
    prop = n/by_ttl,
    smolt_age = fw_age - 1,
    dataset = "adults"
  )


# Produce the same dataset for smolts
smolt_by_prod <- hyatt_data |> 
  summarize(
    .by = c(age, stock, brood_year),
    n = sum(sockeye)
  ) |> 
  # Calculate age composition
  mutate(
    .by = c(stock, brood_year),
    by_ttl = sum(n),
    prop = n/by_ttl,
    dataset = "smolts"
  ) |> 
  rename("smolt_age" = age)


# Make a dataset to show the contrast
by_prod <- bind_rows(
  adult_by_prod,
  smolt_by_prod
) |> 
  mutate(
    .by = c(stock, dataset),
    max_n = max(n, na.rm = TRUE),
    n = n/max_n
  ) |> 
  pivot_wider(
    id_cols = c(stock, brood_year, smolt_age),
    names_from = dataset,
    values_from = c(prop, n)
  )


# Merge the two datasets and plot
by_prod |> 
  ggplot(aes(x = prop_adults, y = prop_smolts)) +
  facet_grid(
    smolt_age ~ stock,
    scales = "free"
  ) +
  geom_abline(slope = 1) +
  geom_point()
# Correlation seems poor...


# Plot the data as a time series
by_prod |> 
  drop_na() |> 
  ggplot(aes(x = brood_year, y = prop_smolts)) +
  facet_grid(smolt_age ~ stock) +
  geom_segment(
    aes(
      xend = brood_year,
      yend = prop_adults
    ),
    colour = "grey"
  ) +
  geom_point(
    aes(size = prop_smolts),
    colour = "blue",
    alpha = 0.5
  ) +
  geom_point(
    aes(
      size = prop_adults,
      y = prop_adults
    ),
    colour = "red",
    alpha = 0.5
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    labels = scales::percent
  ) +
  scale_size(range = c(0.5, 3)) +
  guides(size = "none")


# Final plot just showing the discrepancies over time
by_prod |> 
  drop_na() |> 
  ggplot(aes(x = brood_year, y = prop_adults - prop_smolts)) +
  facet_grid(smolt_age ~ stock) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_segment(aes(yend = 0, xend = brood_year)) +
  geom_point() +
  scale_y_continuous(
    expand = c(0, 0.05),
    labels = scales::percent
  )


# Preliminary conclusion is that any relationship between smolt and adult age 
# compositions is loose at best...
#
# Likely the differences arise from years of favourable marine survival 
# affecting different smolt years from the same brood year. E.g. Age-2 smolt 
# survival may be much higher than age-1 smolt survival from the same brood 
# year if ocean-entry conditions are drastically different between their
# two ocean-entry years.



# Plot ATS survey estimates -----------------------------------------------


# Date annotations for fertilization and hatchery releases
ann <- here(
  "3. outputs",
  "Stock-recruit data",
  "Barkley_Sockeye_stock-recruit_infilled.xlsx"
) |> 
  read_xlsx(sheet = "S-R data") |> 
  select(year, stock, fertilized, hatchery_fry_release) |> 
  mutate(hatchery_fry_release = if_else(hatchery_fry_release > 0, 1, 0)) |> 
  pivot_longer(!c(year, stock)) |> 
  filter(value == 1) |> 
  arrange(stock, name, year) |> 
  mutate(
    .by = c(stock, name),
    min_date = as.Date(paste0(year, "-01-01")),
    max_date = as.Date(paste0(year, "-12-31")),
    lag_year = lag(year, 1),
    gap = if_else(is.na(lag_year), 1, year - lag_year),
    lake = factor(stock, levels = cu_order, labels = cu_names)
  ) |> 
  summarize(
    .by = c(lake, name, gap),
    min_date = min(min_date),
    max_date = max(max_date)
  ) |> 
  left_join(
    summarize(
      .data = ats_est,
      .by = lake,
      y = max(presmolt_est + presmolt_sd, na.rm = TRUE)/1e6 + 5
    )
  ) |> 
  mutate(
    label = case_when(
      name == "fertilized" & as.numeric(max_date - min_date) > 365 ~ name,
      name == "hatchery_fry_release" & as.numeric(max_date - min_date) > 365 ~ "hatchery enhancement",
      TRUE ~ NA
    )
  ) |> 
  # pivot_longer(
  #   contains("date"),
  #   names_sep = "_",
  #   names_to = c("interval", ".value")
  # ) |> 
  select(lake, contains("date"), y, name, label) %>%
  split(.$name)



# Pointrange plot
(ats_p <- ats_est |> 
    pivot_longer(
      matches("_(est|sd)"),
      names_sep = "_",
      names_to = c("species", ".value")
    ) |> 
    mutate(
      lb = if_else(est - sd < 0, 0, est - sd),
      ub = est + sd,
      across(c(est, sd, lb, ub), \(x) x/1e6),
      species = if_else(species == "presmolt", "Sockeye fry", "stickleback")
    ) |> 
    ggplot(aes(x = survey_date, y = est)) +
    facet_wrap(
      ~lake,
      ncol = 1,
      strip.position = "right",
      scales = "free_y"
    ) +
    # Fertilization years
    geom_textsegment(
      data = ann$fertilized,
      aes(
        y = y,
        yend = y,
        x = min_date,
        xend = max_date,
        label = label
      ),
      linewidth = 0.2,
      size = 2
    ) +
    geom_segment(
      data = filter(ann$fertilized, is.na(label)),
      aes(
        y = y,
        yend = y,
        x = min_date,
        xend = max_date
      ),
      linewidth = 0.2
    ) +
    # Hatchery release years
    geom_textsegment(
      data = ann$hatchery_fry_release,
      aes(
        y = y -3,
        yend = y - 3,
        x = min_date,
        xend = max_date,
        label = label
      ),
      lty = 6,
      linewidth = 0.2,
      size = 2
    ) +
    geom_pointrange(
      aes(ymin = lb, ymax = ub, colour = species),
      size = .5,
      shape = "–",
      linewidth = 0.3
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_colour_manual(values = c("black", "red")) +
    labs(
      x = "ATS date",
      y = "Estimated pelagic fish abundance (millions)"
    ) +
    theme(
      strip.background = element_rect(fill = "white"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey85", linewidth = 0.15),
      legend.position = "inside",
      legend.position.inside = c(0.02, 0.96),
      legend.justification.inside = c(0, 1),
      legend.background = element_rect(colour = "black", fill = alpha("white", 0.6)),
      legend.title = element_blank()
    )
)


# Save the plot
ggsave(
  plot = ats_p,
  filename = here(
    "3. outputs",
    "Plots",
    "Lake_ATS_estimates_time_series.png"
  ),
  width = 6.5,
  height = 5,
  units = "in",
  dpi = "print"
)


# Standardize annual abundance estimates using GAMM ------------------------


# Model all 3 lakes together
gamm_combined <- ats_est |> 
  mutate(
    # Convert pre-smolt numbers to millions in observed data
    across(contains("presmolt_"), \(x) x/1e6),
    # Interaction term for year within lake
    lake_year = interaction(lake, smolt_year, sep = "_")
  ) |> 
  nest(.key = "obs") |> 
  rowwise() |> 
  mutate(
    model = list(
      brm(
        bf(
          as.integer(presmolt_est) ~
            s(day_smolt_yr, by = lake, bs = "fs", k = 5) +
            (1 | smolt_year_f) +
            (1 | smolt_year_f:lake),
          family = negbinomial()
        ),
        data = obs, 
        cores = 4
      )
    ),
    # Save dataframe of values to predict over
    newdata = list(
      with(
        obs,
        expand.grid(
          lake = levels(lake),
          smolt_year_f = levels(smolt_year_f),
          #smolt_year = seq.int(min(smolt_year), max(smolt_year)),
          day_smolt_yr = seq.int(1, 365)
        )
      ) |> 
        mutate(lake_year = interaction(lake, smolt_year_f, sep = "_"))
    ),
    # create dataframe of predictions
    pred = list(
      as.data.frame(
        t(
          apply(
            posterior_epred(
              model, 
              newdata = newdata,
              allow_new_levels = TRUE
            ), 
            2, 
            quantile, 
            probs = c(0.05, 0.5, 0.95)
          )
        )
      ) |> 
        bind_cols(newdata)
    )
  )


# Plot predictions from the bayesian gamm across all years
pmap(
  gamm_combined,
  \(obs, pred, ...) ggplot(
    data = obs,
    aes(day_smolt_yr, presmolt_est)
  ) +
    facet_wrap(
      ~lake,
      ncol = 1,
      strip.position = "right",
      scales = "free_y"
    ) +
    geom_line(
      data = pred,
      aes(y = `50%`, group = smolt_year_f),
      colour = "grey70"
    ) +
    geom_point() 
)


# Split data for each lake from combined gamm approach into its own row
gamm_combined_split <- gamm_combined |> 
  select(obs, pred) |> 
  as.list() |> 
  map(list_rbind) |> 
  imap(\(x, idx) nest(x, .by = lake, .key = idx)) %>% 
  {left_join(.[[1]], .[[2]])} |> 
  rowwise() |> 
  mutate(
    model = gamm_combined$model,
    # Make plot of predicted versus observed
    pred_plot = list(
      ggplot(
        data = obs,
        aes(day_smolt_yr, presmolt_est)
      ) +
        facet_wrap(
          ~smolt_year_f, 
          scales = "free_y",
          ncol = 5
        ) +
        geom_vline(
          xintercept = 365-20, # March 1st (day 365 = 20 March)
          lty = 2
        ) + 
        geom_ribbon(
          data = pred,
          aes(
            y = `50%`,
            ymin = `5%`,
            ymax = `95%`
          ),
          alpha = 0.3,
          colour = NA
        ) +
        geom_line(
          data = pred,
          aes(y = `50%`)
        ) +
        geom_text(
          aes(label = smolt_year),
          x = Inf,
          y = Inf,
          hjust = 1.5,
          vjust = 1.5,
          check_overlap = T,
          colour = "grey50",
          size = 4
        ) +
        geom_pointrange(
          aes(
            ymin = presmolt_lwr,
            ymax = presmolt_upr
          ),
          size = 0.15
        ) +
        scale_y_continuous(
          limits = c(0, NA),
          expand = expansion(mult = c(0, 0.05)),
          labels = scales::label_number(accuracy = 1),
          breaks = scales::pretty_breaks(n = 3)
        ) +
        scale_x_continuous(
          expand = c(0, 0),
          labels = scales::label_number(accuracy = 1),
          breaks = scales::pretty_breaks(n = 3)
        ) +
        labs(
          title = paste(lake, "GAMM predictions and ATS estimates by date"),
          x = "Days since 20 March",
          y = "Sockeye fry abundance (millions)"
        ) +
        theme(
          strip.text = element_blank(),
          strip.background = element_blank(),
          panel.grid.minor = element_blank()
        )
    )
  ) |> 
  ungroup()


# View annual predictions per lake
gamm_combined_split |> 
  pull(pred_plot, name = lake)


# Save the plots
gamm_combined_split |> 
  pull(pred_plot, name = lake) |> 
  iwalk(
    \(x, idx) ggsave(
      plot = x,
      filename = here(
        "3. outputs",
        "Plots",
        paste0("GAMM_fry-abun_pred_plot_", idx, ".png")
      ),
      width = 7,
      height = 8,
      units = "in",
      dpi = "print"
    )
  )


# Annual predictions per lake and stock on 1 March
ann_pred_frame <- expand.grid(
  lake = levels(ats_est$lake),
  smolt_year_f = levels(ats_est$smolt_year_f),
  day_smolt_yr = 345 # = March 1st
)

ann_gamm_pred <- posterior_epred(
  gamm_combined$model[[1]], 
  newdata = ann_pred_frame,
  allow_new_levels = TRUE
) |> 
  apply(
    2, 
    # Summarize posterior draws and convert to millions
    function(x) {
      y <- x*1e6
      
      c(
        "mean" = mean(y),
        "est" = median(y),
        "lwr" = quantile(y, 0.1),
        "upr" = quantile(y, 0.9),
        "sd" = sd(y)
      )
    }
  ) |> 
  t() |> 
  as.data.frame() |>
  bind_cols(ann_pred_frame) |> 
  mutate(
    smolt_year = as.numeric(as.character(smolt_year_f)),
    cv = sd/mean
  ) |> 
  select(lake, smolt_year, day_smolt_yr, est, lwr = `lwr.10%`, upr = `upr.90%`, cv)
    

# Export predicted pre-smolt abundance estimates
ann_gamm_pred |> 
  mutate(estimate_date = as.Date(paste0(smolt_year, "-03-01"))) |> 
  select(-day_smolt_yr) |> 
  # Ensure CI widths are explicit in the output file
  rename_with(\(x) paste0(x, "_80"), .cols = c(lwr, upr)) |> 
  # Add a flag for years where estimates are not based on an observed survey
  # (i.e. are interpolated)
  left_join(
    summarize(
      ats_est,
      .by = c(smolt_year, lake),
      interpolated = 0
    )
  ) |> 
  mutate(interpolated = if_else(is.na(interpolated), "Yes", "No")) |> 
  write.csv(
    file = here(
      "3. outputs",
      "Stock-recruit data",
      "GAMM-estimated_pre-smolt_time_series.csv"
    ),
    row.names = FALSE
  )


# Compare gamm predictions to raw values for each lake and year
(gamm_comp_p <- ats_est |> 
    filter(preferred_est == 1) |> 
    complete(lake, smolt_year) |> 
    # Screw around with column names so raw and gamm estimates can be
    # compared side-by-side in the plot
    rename_with(\(x) str_replace_all(x, "presmolt", "raw")) |> 
    rename("raw_day" = day_smolt_yr) |> 
    left_join(
      ann_gamm_pred,
      by = c("lake", "smolt_year")
    ) |> 
    rename("gamm_day" = day_smolt_yr) |> 
    rename_with(\(x) paste0("gamm_", x), .cols = c(est, lwr, upr)) |> 
    # Lengthen data so raw (ATS) and gamm estimates have their own rows
    pivot_longer(
      matches("(raw|gamm)_(est|lwr|upr|day)"),
      names_sep = "_",
      names_to = c("method", ".value")
    ) |> 
    mutate(across(c(est, lwr, upr), \(x) x/1e6)) |>  # Convert to millions
    ggplot(aes(x = smolt_year, y = est)) +
    facet_wrap(
      ~lake,
      strip.position = "right",
      scales = "free_y",
      ncol = 1
    ) +
    geom_pointrange(
      aes(
        ymin = lwr,
        ymax = upr,
        shape = method,
        colour = day 
        # Show estimate date on a continuous colour scale
        # (similar colours should have similar estimates from gamm vs. ATS)
      ),
      position = position_dodge(width = 0.5),
      linewidth = 0.2,
      size = 0.15
    ) +
    scale_colour_viridis_c(
      direction = -1,
      option = "mako",
      labels = 
    ) +
    scale_shape_manual(values = c(4, 19)) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    guides(shape = "none") +
    labs(
      x = "Smolt outmigration year",
      y = "Estimated Sockeye fry abundance (millions)",
      colour = "Days since\n20 March"
    ) +
    theme(strip.background = element_rect(fill = "white"))
)


# Save the plot
ggsave(
  plot = gamm_comp_p,
  filename = here(
    "3. outputs",
    "Plots",
    "Pre-smolt_estimates_GAMM_vs_ATS.png"
  ),
  width = 8,
  height = 5,
  units = "in",
  dpi = "print"
)



# Use GAMM predictions to estimate fry-smolt survival --------------------


# For each lake & year, determine which date corresponds to peak fry abundance
pred_pairs <- gamm_combined |> 
  unnest(pred) |> 
  select(-c(obs:newdata)) |> 
  filter(
    .by = lake_year,
    `50%` == max(`50%`)
  ) |> 
  select(lake, smolt_year_f, day_smolt_yr_max = day_smolt_yr, lake_year) |> 
  mutate(day_345 = 345) |> 
  filter(lake_year %in% paste0(ats_est$lake, "_", ats_est$smolt_year_f)) |>
  arrange(lake, smolt_year_f) |> 
  pivot_longer(
    cols = c(day_smolt_yr_max, day_345),
    names_to = "which_day",
    values_to = "day_smolt_yr"
  ) 


# Get posterior expected values (draws × rows)
# posterior_epred returns a matrix: draws x nrow(newdata)
fs_s_pred <- posterior_epred(gamm_combined[["model"]][[1]], newdata = pred_pairs) |> 
  t() |> 
  as.data.frame() |> 
  cbind(pred_pairs) |> 
  pivot_longer(
    cols = !colnames(pred_pairs),
    names_to = "draw",
    names_pattern = "V(.*)"
  ) |> 
  pivot_wider(
    id_cols = !day_smolt_yr,
    names_from = which_day,
    values_from = value
  ) |> 
  mutate(fs_s = day_345/day_smolt_yr_max)


# Multi-year summary of estimated fry-smolt survival rates
fs_s_pred |> 
  summarize(
    .by = lake,
    q50 = quantile(fs_s, 0.5),
    q05 = quantile(fs_s, 0.05),
    q95 = quantile(fs_s, 0.95)
  )


# Use smolt rates of capture to estimate annual production --------------


# Calculate smolt rate of capture by age, year, and lake
rates_of_capture <- smolt_ages |> 
  select(stock = cu, sample_date, n_age = ttl, gear, contains("prop")) |> 
  full_join(smolt_catch) |> 
  rename("n_catch" = count) |> 
  mutate(
    across(
      matches("prop_age\\d"), 
      \(x) x * n_catch, 
      .names = "catch_{str_remove(.col, 'prop_')}"
    )
  ) |> 
  select(stock, sample_date, year, matches("age\\d"), contains("n_")) |>
  # Ensure all sample dates are aggregated into a single row per lake
  # (i.e. combine data across gear types)
  summarize(
    .by = c(year, sample_date, stock),
    across(c(catch_age1, catch_age2, n_age, n_catch), \(x) sum(x, na.rm = TRUE))
  ) |> 
  bind_rows(new_smolt_ages) |> 
  summarize(
    .by = c(year, stock),
    across(c(catch_age1, catch_age2, n_age, n_catch), \(x) sum(x, na.rm = TRUE)),
    sample_days = n()
  ) |> 
  mutate(
    across(
      matches("catch_age\\d"),
      \(x) x / sample_days,
      .names = "{str_replace(.col, 'catch', 'rate')}"
    )
  ) |> 
  # Make missing values explicit
  right_join(expand_grid(stock = cu_names, year = seq.int(1977, 2023))) |> 
  arrange(stock, year) |> 
  mutate(
    .by = stock,
    `rate_age1+` = lead(rate_age2)
  ) |> 
  rowwise() |> 
  mutate(
    rate_ttl = sum(c_across(contains("rate_"))),
    across(
      matches("rate_age\\d"),
      \(x) x / rate_ttl,
      .names = "{str_replace(.col, 'rate', 'prop')}"
    )
  ) |> 
  ungroup()


# Export clean data file with smolt rates of capture and age compositions
rates_of_capture |> 
  select(year, lake = stock, sample_days, matches("(rate|prop)_age\\d")) |> 
  write.csv(
    here(
      "3. outputs",
      "Stock-recruit data",
      "Inferred_annual_smolt_age_composition.csv"
    ),
    row.names = FALSE
  )


# Export adult brood year age compositions as a potential weakly informative
# prior for infilling missing smolt age compositions
adult_by_prod |> 
  mutate(smolt_age = paste0("returns_age", smolt_age)) |> 
  pivot_wider(
    id_cols = c(stock, brood_year, by_ttl),
    names_from = smolt_age,
    values_from = prop
  ) |> 
  arrange(stock, brood_year) |> 
  mutate(
    .by = stock,
    `returns_age1+` = lead(returns_age2)
  ) |> 
  rename(
    "lake" = stock,
    "year" = brood_year,
    "ttl" = by_ttl
  ) |> 
  write.csv(
    here(
      "3. outputs",
      "Stock-recruit data",
      "Annual_adult_fw-age_composition.csv"
    ),
    row.names = FALSE
  )


# Filter ATS estimates for best annual value to use
ats_annual <- ats_est |> 
  filter(preferred_est == 1) 


# Apply smolt age estimates based on rate of capture to ATS values
ats_annual_age <- rates_of_capture |> 
  select(year, stock, contains("prop")) |> 
  left_join(
    ats_annual,
    by = join_by(
      stock == lake, 
      year == smolt_year
    )
  ) |> 
  pivot_longer(
    matches("prop"),
    names_to = "age",
    values_to = "age_prop",
    names_prefix = "prop_age"
  ) |> 
  mutate(
    across(matches("presmolt"), \(x) x * age_prop),
    brood_year = if_else(age == "2", age2_BY, age1_BY)
  )


# Calculate new brood year smolt production
by_prod2 <- ats_annual_age |> 
  mutate(age = if_else(age == "1+", "2", age)) |> 
  summarize(
    .by = c(brood_year, age, stock),
    n_smolts = sum(presmolt_est)
  ) |> 
  rename("smolt_age" = age)  |> 
  mutate(
    .by = c(brood_year, stock),
    ttl_smolts = sum(n_smolts),
    prop_smolts = n_smolts/ttl_smolts,
    smolt_age = as.numeric(smolt_age)
  ) |> 
  # Bring in the adult data
  full_join(select(by_prod, stock, brood_year, smolt_age, contains("adults")))


# Reproduce plots comparing smolt to adult comps
by_prod2 |> 
  ggplot(aes(x = prop_adults, y = prop_smolts)) +
  facet_grid(
    smolt_age ~ stock,
    scales = "free"
  ) +
  geom_abline(slope = 1) +
  geom_point()
# Correlation seems poor...
  
  
# Plot the data as a time series
by_prod2 |> 
  drop_na() |> 
  ggplot(aes(x = brood_year, y = prop_smolts)) +
  facet_grid(smolt_age ~ stock) +
  geom_segment(
    aes(
      xend = brood_year,
      yend = prop_adults
    ),
    colour = "grey"
  ) +
  geom_point(
    aes(size = prop_smolts),
    colour = "blue",
    alpha = 0.5
  ) +
  geom_point(
    aes(
      size = prop_adults,
      y = prop_adults
    ),
    colour = "red",
    alpha = 0.5
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    labels = scales::percent
  ) +
  scale_size(range = c(0.5, 3)) +
  guides(size = "none")


# Final plot just showing the discrepancies over time
by_prod2 |> 
  drop_na() |> 
  ggplot(aes(x = brood_year, y = prop_adults - prop_smolts)) +
  facet_grid(smolt_age ~ stock) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_segment(aes(yend = 0, xend = brood_year)) +
  geom_point() +
  scale_y_continuous(
    expand = c(0, 0.05),
    labels = scales::percent
  )


# Plot Bayesian model outputs alongside GAM and raw estimates -------------


# Full smolt production results from the model
bayes_smolts <- here(
  "3. outputs",
  "Stock-recruit data",
  "Bayesian_state-space_smolt-production_estimated_time_series.xlsx"
) |> 
  read_xlsx(sheet = "model_estimates") |> 
  filter(parameter == "N_lake") |> 
  mutate(
    across(matches("%"), \(x) x/1e6),
    lake = factor(lake, levels = c("Great Central", "Sproat", "Hucuktlis"))
  )


# Add Bayesian model estimates to GAM comparison plot
(gam_bayes_comp_p <- gamm_comp_p +
  geom_pointrange(
    data = bayes_smolts,
    aes(
      x = year,
      y = `50%`,
      ymin = `10%`,
      ymax = `90%`,
      colour = 345
    ),
    size = 0.15, 
    linewidth = 0.2, 
    shape = 17
  )
)


# Explort plot
ggsave(
  gam_bayes_comp_p,
  filename = here(
    "3. outputs",
    "Plots",
    "Fry_abundance_estimates_comparison.png"
  ),
  width = 8,
  height = 5,
  units = "in",
  dpi = "print"
)

