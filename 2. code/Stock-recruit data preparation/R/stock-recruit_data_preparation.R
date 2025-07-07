# Packages ----------------------------------------------------------------

pkgs <- c("here", "tidyverse", "readxl", "writexl", "broom", "geomtextpath")
#install.packages(pkgs)

library(here)
library(tidyverse); theme_set(theme_bw())
library(geomtextpath)
library(readxl)
library(writexl)
library(broom)


# Load and reformat data for Somass CUs -----------------------------------


# Load stock-recruit time series by return year
Som_run_ts <- here("1. data", "Barkley Sockeye time series of returns.xlsx") |> 
  read_xlsx(sheet = "Somass") |> 
  mutate(
    # Move categorical escapement estimates to a new column
    esc_cat = if_else(str_detect(escapement, "\\d+"), NA, escapement),
    escapement = as.double(escapement),
    run = catch + escapement,
    brood_year = return_year - ttl_age
  )


# Calculate brood year returns
Som_brood_ts <- Som_run_ts |> 
  filter(!is.na(brood_year)) |> 
  summarize(
    .by = c(brood_year, stock),
    return = sum(catch, escapement),
    n = n()
  ) |> 
  # Remove incomplete brood years
  filter(n > 5) |> 
  select(-n)


# Load data on age samples
Somass_age_samples <- here("1. data", "Barkley Sockeye time series of returns.xlsx") |> 
  read_xlsx(sheet = "Age_samples")


# Use annual proportions of GCL/SPR in the return to back-calculate 
# stock-specific age samples
Somass_age_samples_spread <- Som_run_ts |> 
  filter(return_year >= min(Somass_age_samples$year)) |> 
  summarize(
    .by = c(return_year, stock),
    run = sum(run)
  ) |> 
  mutate(
    .by = return_year,
    ttl = sum(run),
    prop = run/ttl
  ) |> 
  pivot_wider(
    id_cols = return_year,
    names_from = stock,
    values_from = prop
  ) |> 
  right_join(
    Somass_age_samples,
    by = c("return_year" = "year")
  ) |> 
  pivot_longer(
    c(GCL, SPR),
    names_to = "prob_stock",
    values_to = "prob"
  ) |> 
  mutate(
    age.samples = if_else(stock == "Somass aggregate", round(N * prob, 0), N),
    stock = if_else(stock == "Somass aggregate", prob_stock, stock)
  ) |> 
  select(return_year, stock, source, age.samples) |> 
  distinct()


# Historic average # age samples per stock to infill 1977-1979
mean_age_samples <- Somass_age_samples_spread |> 
  # Use only 1980s period when sampling was most similar
  filter(return_year < 1990) |> 
  summarize(
    .by = stock,
    age.samples = round(mean(age.samples), 0)
  ) |> 
  pull(age.samples, name = stock)


# Build stock-recruit table
Som_sr <- Som_run_ts |> 
  filter(!is.na(brood_year)) |> 
  # Add age sample number data
  left_join(select(Somass_age_samples_spread, -source)) |>
  mutate(
    .by = c(return_year, stock),
    marine_age = ttl_age - fw_age,
    spawners = sum(escapement),
    adult_spawners = sum(escapement[marine_age > 1]),
    age.samples = case_when(
      is.na(age.samples) & stock == "GCL" ~ mean_age_samples[["GCL"]],
      is.na(age.samples) & stock == "SPR" ~ mean_age_samples[["SPR"]],
      TRUE ~ age.samples
    )
  ) |> 
  pivot_wider(
    id_cols = c(return_year, stock, contains("spawners"), age.samples, fertilized),
    names_from = ttl_age,
    names_prefix = "N.age.",
    values_from = run,
    values_fn = sum
  ) |> 
  rowwise() |> 
  mutate(
    run = sum(c_across(matches("N\\.age\\.\\d"))),
    # Convert age #s to proportions
    across(matches("N\\.age\\.\\d"), \(x) x/run),
    H = run - spawners,
  ) |> 
  ungroup() |> 
  left_join(
    Som_brood_ts,
    by = join_by(
      stock, 
      return_year == brood_year
    )
  ) |> 
  rename(
    "year" = return_year,
    "S" = spawners,
    "adult_S" = adult_spawners,
    "N" = run,
    "R" = return
  )



# Load and work up Henderson/Hucuktlis data -------------------------------


# Load raw Henderson returns table
Huc_run_ts <- here("1. data", "Barkley Sockeye time series of returns.xlsx") |> 
  read_xlsx(sheet = "Hucuktlis") |> 
  rename_with(\(x) str_replace_all(x, "\\s", "_")) |> 
  mutate(
    run = catch + escapement,
    stock = "HUC"
  ) |> 
  rename("return_year" = year) |> 
  # flatten age samples across years
  mutate(across(matches("^age_\\d+"), \(x) x * age_sample_size)) |> 
  summarize(
    .by = !matches("age_(\\d+|sample)"),
    across(c(age_sample_size, matches("age_\\d+")), \(x) sum(x, na.rm = TRUE))
  ) |> 
  mutate(
    across(
      matches("^age_\\d+"), 
      \(x) if_else(
        age_sample_size > 0, 
        x / age_sample_size,
        NA_real_
      )
    ),
    adult_escapement = escapement*(age_42 + age_52 + age_53 + age_62 + age_63)
  )


# Confirm all age composition columns sum to 1
Huc_run_ts |> 
  rowwise() |> 
  mutate(sum = sum(c_across(matches("^age_\\d+")))) |> 
  count(sum)

# Ensure only 1 row per year exists
Huc_run_ts |> 
  count(return_year) |> 
  distinct(n)


# Use Somass harvest rates to hindcast Hucuktlis harvest rates ------------


# Build time series of harvest rates on both stocks
hr_ts <- list(
  "Somass" = Som_run_ts, 
  "Hucuktlis" = Huc_run_ts
  ) |> 
  imap(
    \(x, idx) summarize(
      x,
      .by = return_year,
      catch = sum(catch),
      run = sum(run),
      stock = idx
    )
  ) |> 
  list_rbind() |> 
  mutate(hr = catch/run) |> 
  filter(!is.na(hr))


# Plot Somass versus Hucuktlis harvest rate
(hr_p <- hr_ts |> 
    pivot_wider(
      id_cols = return_year,
      names_from = stock,
      values_from = hr
    ) |> 
    filter(
      !if_any(c(Somass, Hucuktlis), c(is.na, is.nan)),
      return_year > 2011 # use only years when DNA data are available
    ) |> 
    ggplot(aes(x = Somass, y = Hucuktlis)) +
    geom_abline(
      slope = 1,
      colour = "grey50",
      lty = 2
    ) +
    geom_point(
      aes(fill = return_year),
      size = 3,
      shape = 21,
      stroke = NA
    ) +
    scale_fill_viridis_c(labels = as.integer) +
    scale_x_continuous(
      labels = scales::percent,
      expand = expansion(c(0, 0.05)),
      limits = c(0, max(hr_ts$hr)),
      oob = scales::oob_keep
    ) +
    scale_y_continuous(
      labels = scales::percent,
      expand = expansion(c(0, 0.05)),
      limits = c(0, max(hr_ts$hr)),
      oob = scales::oob_keep
    ) +
    coord_fixed() +
    labs(
      title = "Somass versus Hucuktlis Sockeye harvest rates",
      colour = "Year"
    ) +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.05, 0.95),
      legend.justification.inside = c(0, 1),
      legend.background = element_rect(colour = "black")
    )
)
# Looks probable that a relationship can be established that reasonably 
# links Hucuktlis harvest rates to Somass harvest rates


# Fit various models to predict the Hucuktlis harvest rate
hr_mods <- hr_ts |> 
  filter(return_year > 2011) |> 
  mutate(esc = run - catch) |> 
  pivot_wider(
    names_from = stock,
    values_from = c(catch, run, hr, esc),
    values_fill = NA,
    names_glue = "{stock}_{.value}"
  ) |> 
  nest() |> 
  rowwise() |> 
  mutate(
    glm = list(
      glm(
        cbind(Hucuktlis_catch, Hucuktlis_esc) ~ Somass_hr,
        data = data,
        family = "binomial"
      )
    ),
    lm = list(lm(Hucuktlis_hr ~ Somass_hr, data = data)),
    loglog_lm = list(lm(log(Hucuktlis_hr) ~ log(Somass_hr), data = data)),
    loglin_lm = list(lm(log(Hucuktlis_hr) ~ Somass_hr, data = data)),
    linlog_lm = list(lm(Hucuktlis_hr ~ log(Somass_hr), data = data)),
    pred_frame = list(
      as_tibble_col(
        seq(0, max(hr_ts$hr), by = 0.01), 
        "Somass_hr"
      )
    )
  ) |> 
  pivot_longer(
    c(contains("lm")),
    names_to = "model_type",
    values_to = "model"
  ) |> 
  rowwise() |> 
  mutate(
    predict = list(
      predict(model, pred_frame, se.fit = TRUE) |> 
        as_tibble() |> 
        mutate(
          lwr = fit - 1.96*se.fit,
          upr = fit + 1.96*se.fit
        ) |> 
        cbind(pred_frame)
    ),
    glance = list(broom::glance(model))
  ) |> 
  unnest(glance)


# Pull the model predictions out as a flattened dataframe for plotting
hr_mod_preds <- hr_mods |> 
  select(model_type, predict) |> 
  unnest(predict) |> 
  mutate(
    across(
      c(fit, lwr, upr),
      \(x) case_when(
        model_type == "glm" ~ binomial()$linkinv(x),
        model_type %in% c("loglog_lm", "loglin_lm") ~ exp(x),
        TRUE ~ x
      )
    )
  )


# Plot model predictions
hr_p + 
  # Add model fit names with geom_textline?
  geom_textline(
    data = hr_mod_preds,
    aes(
      x = Somass_hr, 
      y = fit, 
      colour = model_type,
      label = model_type,
    ),
    hjust = 0.4
  ) +
  guides(colour = "none", label = "none") +
  scale_y_continuous(
    labels = scales::percent,
    expand = expansion(c(0, 0.05)),
    limits = c(0, max(hr_ts$hr)),
    oob = scales::oob_keep
  ) 
# Linear model looks best


# Stipulate which model will be used for the hindcast
Huc_hr_pred_model <- hr_mods[hr_mods$model_type=="lm",]$model[[1]]


# Function to calculate coefficient of variation from a fitted model
get_lm_cv <- function(lm_obj) {
  
  rmse <- sqrt(mean(lm_obj$residuals^2))
  
  mean_dependent_var <- mean(lm_obj$model[,1])
  
  cv = rmse/mean_dependent_var
  
  return(cv)
}


# Hindcast Henderson harvest rates using Somass harvest rates
Huc_run_ts_infill <- hr_ts |> 
  filter(stock == "Somass") |> 
  select(return_year, hr) |> 
  rename("Somass_hr" = hr) |> 
  nest() |> 
  rowwise() |> 
  mutate(
    pred = list(
      predict(
        object = Huc_hr_pred_model,
        newdata = data,
        interval = "pred",
      ) |> 
        as.data.frame()
    )
  ) |> 
  unnest(c(data, pred)) |> 
  mutate(Huc_hr_pred = if_else(fit < 0, 0, fit)) |> 
  select(return_year, contains("hr")) |> 
  right_join(Huc_run_ts) |> 
  mutate(
    # Only show predicted HRs and CVs in rows without catch data
    Huc_hr_pred_cv = if_else(
      is.na(catch), 
      get_lm_cv(Huc_hr_pred_model), 
      NA
    ),
    Huc_hr_pred = if_else(
      is.na(catch),
      Huc_hr_pred,
      NA
    ),
    # Add new catch data based on hindcast harvest rates
    catch = if_else(
      is.na(catch),
      (escapement*Huc_hr_pred)/(1-Huc_hr_pred),
      catch
    ),
    stockid_source = if_else(
      is.na(stockid_source) & !is.na(Huc_hr_pred),
      "HR model hindcast",
      stockid_source
    ),
    run = catch + escapement
  )



# Reformat and collate Somass and Hucuktlis data --------------------------


# Hucuktlis brood year returns lookup
Huc_brood_ts <- Huc_run_ts_infill |> 
  mutate(
    across(matches("^age_\\d+"), \(x) x*run, .names = "run_{.col}"),
    across(matches("^age_\\d+"), \(x) x*run, .names = "esc_{.col}")
  ) |> 
  pivot_longer(
    matches("(run|esc)_age_\\d+"),
    names_pattern = "(.*)_age_(.*)",
    names_to = c("cat", "age")
  ) |> 
  pivot_wider(
    id_cols = c(return_year, age, age_sample_size),
    names_from = cat,
    values_from = value
  ) |> 
  mutate(
    ttl_age = as.numeric(str_sub(age, 1, 1)),
    brood_year = return_year - ttl_age
  ) |> 
  filter(
    # Remove years without spawner or catch data
    !if_any(c(run, esc), is.na),
    # Don't calculate recruits when <10 age samples available
    age_sample_size >= 10 
  ) |> 
  summarize(
    .by = brood_year,
    return = sum(run),
    n = n()
  ) |> 
  # Remove brood years with incomplete returns
  filter(n > 5) |> 
  select(-n)


# Get Hucuktlis data reformatted properly
Huc_sr <- Huc_run_ts_infill |> 
  # Add recruitment data
  left_join(
    Huc_brood_ts,
    by = c("return_year" = "brood_year")
  ) |> 
  rename(
    "year" = return_year,
    "S" = escapement,
    "adult_S" = adult_escapement,
    "N" = run,
    "H" = catch,
    "R" = return,
    "age.samples" = age_sample_size
  ) |> 
  select(-Somass_hr) |> 
  rename_with(\(x) str_remove(x, "Huc_")) |> 
  # Combine age compositions by total age
  rowwise() |> 
  mutate(
    N.age.3 = sum(c_across(matches("^age_3."))),
    N.age.4 = sum(c_across(matches("^age_4."))),
    N.age.5 = sum(c_across(matches("^age_5."))),
    N.age.6 = sum(c_across(matches("^age_6."))),
    .keep = "unused"
  ) |> 
  filter(!year < 1972) # Remove years prior to 1972


# Read in CV values from excel file
esc_CV <- here(
  "1. data",
  "Barkley Sockeye time series of returns.xlsx"
  ) |> 
  read_xlsx(sheet = "Escapement_CV") |> 
  select(year, stock, S_cv = CV)

catch_CV <- here(
  "1. data",
  "Barkley Sockeye time series of returns.xlsx"
) |> 
  read_xlsx(sheet = "Catch_CV") |> 
  select(year, stock, H_cv = CV)


# Collate the two dataframes
Barkley_sk_sr <- bind_rows(Som_sr, Huc_sr) |> 
  # Add annual CV estimates for each CU's catch & escapement data
  left_join(esc_CV) |> 
  left_join(catch_CV) |> 
  mutate(H_cv = if_else(is.na(H_cv), hr_pred_cv, H_cv)) |> 
  select(
    year,
    stock, 
    S,
    adult_S,
    N,
    H,
    contains("N.age"),
    age.samples,
    fertilized,
    hatchery_fry_release,
    R,
    hr_pred,
    H_cv,
    S_cv
  ) |> 
  arrange(stock, year)


# Save metadata for solumn names
Barkley_sk_sr_metadata <- Barkley_sk_sr |> 
  colnames() |> 
  as_tibble_col("column_name") |> 
  mutate(
    definition = case_when(
      column_name == "year" ~ "Run (i.e. observation) year",
      column_name == "stock" ~ paste(
      "Conservation Unit. Either 'Great Central' (GCL),",
      "'Sproat' (SPR), or 'Hucuktlis' (formerly 'Henderson'; HUC)"
      ),
      column_name == "S" ~ "Spawners a.k.a escapement",
      column_name == "adult_S" ~ paste(
      "Spawners (S, above) with jacks (marine age = 1) removed"
      ),
      str_detect(column_name, "N.age.\\d") ~ paste(
        "Proportion of the run at the given total age (in years)"
      ),
      column_name == "N" ~ "Annual terminal run size. Equal to S +  C",
      column_name == "H" ~ "Annual terminal harvest (i.e. catch)",
      column_name == "age.samples" ~ paste(
        "Number of fish sampled to calculate the age compositions",
        "given in columns labeled 'N.age.#'"
      ),
      column_name == "fertilized" ~ paste(
        "Binary variable describing whether (1) or not (0) the CU nursary",
        "lake was fertilized in each year. Note that fertilization affects",
        "abundances of pre-smolts in year + 1 (e.g. fertilizing a lake in",
        "2010 is expected to affect the fry arising from brood year 2009)."
      ),
      column_name == "hatchery_fry_release" ~ paste(
        "Numbers of Sockeye fry released by the Hucuktlis/Henderson hatchery",
        "in each year. Pertains only to the Hucuktlis CU because outplants were",
        "discontinued in the 1930s. The hatchery ceased operations in 2007."
      ),
      column_name == "R" ~ paste(
        "Recruitment that arose from S in subsequent years. Calculated",
        "based on annual age compositions"
      ),
      column_name == "hr_pred" ~ paste(
        "Retrospective prediction of harvest rate for Hucuktlis Sockeye",
        "based on a linear model with Somass Sockeye harvest rate as a",
        "predictor. Used to estimate Hucuktlis catch in years with missing data"
      ),
      column_name == "H_cv" ~ paste(
        "Coefficient of variation on harvest data. Historical (prior to 2011)",
        "Hucuktlis Sockeye harvest rate predictions were derived from a",
        "linear model. Harvest data for Somass and Hucuktlis post-2011 are",
        "assumed to be precise."
      ),
      column_name == "S_cv" ~ paste(
        "Coefficient of variation on spawner data. Currently based on ____"
      ),
      TRUE ~ "!definition required!"
    )
  )


# Save the collated Barkley Sockeye stock-recruit data --------------------


# Save sheets as a list
list(
  "metadata" = Barkley_sk_sr_metadata,
  "S-R data" = Barkley_sk_sr
) |> 
  write_xlsx(
    path = here(
      "3. outputs", 
      "Stock-recruit data",
      "Barkley_Sockeye_stock-recruit_infilled.xlsx"
    )
  )




# Make a continuous time series of reconstructed returns ------------------


# Need Hucuktlis run by Gilbert-Rich age
HUC_rr_data <- Huc_run_ts_infill |> 
  filter(return_year >= 1972) |> 
  # Infill age compositions using historic averages
  mutate(
    across(
      matches("age_\\d{2}"),
      \(x) if_else(
        is.na(x),
        mean(x, na.rm = TRUE),
        x
      )
    )
  ) |> 
  pivot_longer(
    cols = matches("age_\\d{2}"),
    names_to = "age",
    values_to = "age_prop",
    names_prefix = "age_"
  ) |> 
  mutate(
    across(c(escapement, catch), \(x) x * age_prop),
    fw_age = as.numeric(str_sub(age, 2L, 2L)),
    ttl_age = as.numeric(str_sub(age, 1L, 1L)),
    brood_year = return_year - ttl_age
  ) |> 
  select(
    return_year, escapement, catch, fertilized, hatchery_fry_release,
    matches("(catch|escapement|age|stockid)(|_data)_source"), stock,
    fw_age, ttl_age, age_sample_size
  )


# Collate Hucuktlis and Somass time series
Som_Huc_rr_data <- Som_run_ts |>
  rename("hatchery_release" = hatchery_eggs_release) |> 
  filter(!is.na(ttl_age)) |> 
  bind_rows(rename(HUC_rr_data, "hatchery_release" = hatchery_fry_release))


# Export 
write.csv(
  Som_Huc_rr_data,
  file = here(
    "3. outputs",
    "Stock-recruit data",
    "Barkley_Sockeye_returns_by_GR_age.csv"
  ),
  row.names = FALSE
)


# Aggregate and save the Barkley Sockeye total returns by year data -------


# Load and reformat historic catch data
historic_catch <- here("1. data", "Barkley Sockeye time series of returns.xlsx") |> 
  read_xlsx(sheet = "Barkley agg") |> 
  rename("catch_data_source" = data_source) |> 
  mutate(
    total_catch = if_else(
      !is.na(Somass_catch) & !is.na(Hucuktlis_catch), 
      NA, 
      total_catch
    )
  ) |> 
  pivot_longer(
    matches(".*_catch$"),
    names_pattern = "(.*)_catch",
    names_to = "stock",
    values_to = "catch"
  ) |> 
  mutate(
    stock = case_when(
      stock == "Hucuktlis" ~ "HUC",
      stock == "Somass" ~ "GCL + SPR",
      stock == "total" ~ "TTL"
    )
  )


# Modern time series of catch and escapement for Somass stocks
Somass_catch_esc <- Som_run_ts |> 
  summarize(
    .by = c(
      return_year, stock, contains("_source"), esc_cat, 
      fertilized, hatchery_eggs_release
    ),
    across(c(catch, escapement), sum)
  )


# Collate all the data and reformat for export
Barkley_catch_esc <- Huc_run_ts |> 
  select(
    return_year, stock, contains("_source"), catch, escapement,
    fertilized, hatchery_fry_release
  ) |> 
  bind_rows(Somass_catch_esc) |> 
  rename("year" = return_year) |> 
  full_join(
    historic_catch,
    by = c("year", "stock")
  ) |> 
  mutate(
    catch = coalesce(catch.x, catch.y),
    catch_data_source = coalesce(catch_data_source.x, catch_data_source.y)
  ) |> 
  select(-matches("\\.(x|y)")) |> 
  mutate(
    .by = year,
    Barkley_catch = sum(catch, na.rm = TRUE),
    Barkley_escapement = sum(escapement, na.rm = TRUE)
  ) |> 
  filter(!stock == "TTL") |> 
  mutate(Barkley_escapement = if_else(Barkley_escapement == 0, NA, Barkley_escapement)) |> 
  rename_with(\(x) paste0("annual_", x), .cols = contains("Barkley_")) |> 
  rename(
    "qualitative_escapement" = esc_cat,
    "CU" = stock
  ) |> 
  relocate(
    escapement,
    qualitative_escapement,
    catch,
    escapement_data_source,
    catch_data_source,
    .after = CU
  )


# Save the formatted table
write_xlsx(
  list(data = Barkley_catch_esc),
  path = here(
    "3. outputs",
    "Total return time series",
    "Barkley Sockeye total annual returns by CU_collated.xlsx"
  )
)
