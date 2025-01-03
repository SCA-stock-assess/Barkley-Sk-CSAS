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
Som_run_ts <- here("1. data", "return by age time series.xlsx") |> 
  read_xlsx(sheet = "Somass") |> 
  mutate(
    run = catch + escapement,
    brood_year = return_year - ttl_age
  )


# Calculate brood year returns
Som_brood_ts <- Som_run_ts |> 
  summarize(
    .by = c(brood_year, stock),
    return = sum(catch, escapement)
  )


# Build stock-recruit table
Som_sr <- Som_run_ts |> 
  mutate(
    .by = c(return_year, stock),
    spawners = sum(escapement)
  ) |> 
  pivot_wider(
    id_cols = c(return_year, stock, spawners),
    names_from = ttl_age,
    names_prefix = "N.age.",
    values_from = run,
    values_fn = sum
  ) |> 
  rowwise() |> 
  mutate(
    run = sum(c_across(contains("age"))),
    # Convert age #s to proportions
    across(contains("age"), \(x) x/run),
    H = run - spawners,
    age.samples = 1000 # Number of age samples per stock per year. 
    # Using 1000/year as placeholder... would need to do considerable
    # data gathering from historic files to get the time series of 
    # number of samples going back to 1977. Does this actually 
    # matter enough to be worth doing so?
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
    "N" = run,
    "R" = return
  )



# Load and work up Henderson/Hucuktlis data -------------------------------


# Load raw Henderson returns table
Hed_run_ts <- here("1. data", "return by age time series.xlsx") |> 
  read_xlsx(sheet = "Henderson") |> 
  mutate(
    run = catch + escapement,
    stock = "HED"
  ) |> 
  rename("return_year" = year) |> 
  # flatten age samples across years
  mutate(across(matches("^age_\\d+"), \(x) x * age_sample_size)) |> 
  select(-sample_type) |> 
  summarize(
    .by = !contains("age"),
    across(contains("age"), \(x) sum(x, na.rm = TRUE))
  ) |> 
  mutate(
    across(
      matches("^age_\\d+"), 
      \(x) if_else(
        age_sample_size > 0, 
        x / age_sample_size,
        NA_real_
      )
    )
  )


# Confirm all age composition columns sum to 1
Hed_run_ts |> 
  rowwise() |> 
  mutate(sum = sum(c_across(matches("^age_\\d+")))) |> 
  count(sum)



# Use Somass harvest rates to hindcast Hucuktlis harvest rates ------------


# Build time series of harvest rates on both stocks
hr_ts <- list(
  "Somass" = Som_run_ts, 
  "Hucuktlis" = Hed_run_ts
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
      return_year > 2011
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
Hed_hr_pred_model <- hr_mods[hr_mods$model_type=="lm",]$model[[1]]


# Function to calculate coefficient of variation from a fitted model
get_lm_cv <- function(lm_obj) {
  
  rmse <- sqrt(mean(lm_obj$residuals^2))
  
  mean_dependent_var <- mean(lm_obj$model[,1])
  
  cv = rmse/mean_dependent_var
  
  return(cv)
}


# Hindcast Henderson harvest rates using Somass harvest rates
Hed_run_ts_infill <- hr_ts |> 
  filter(stock == "Somass") |> 
  select(return_year, hr) |> 
  rename("Somass_hr" = hr) |> 
  nest() |> 
  rowwise() |> 
  mutate(
    pred = list(
      predict(
        object = Hed_hr_pred_model,
        newdata = data,
        interval = "pred",
      ) |> 
        as.data.frame()
    )
  ) |> 
  unnest(c(data, pred)) |> 
  mutate(Hed_hr_pred = if_else(fit < 0, 0, fit)) |> 
  select(return_year, contains("hr")) |> 
  right_join(Hed_run_ts) |> 
  mutate(
    # Only show predicted HRs and CVs in rows without catch data
    Hed_hr_pred_cv = if_else(
      is.na(catch), 
      get_lm_cv(Hed_hr_pred_model), 
      NA
    ),
    Hed_hr_pred = if_else(
      is.na(catch),
      Hed_hr_pred,
      NA
    ),
    # Add new catch data based on hindcast harvest rates
    catch = if_else(
      is.na(catch),
      (escapement*Hed_hr_pred)/(1-Hed_hr_pred),
      catch
    ),
    run = catch + escapement
  )



# Reformat and collate Somass and Hucuktlis data --------------------------


# Get Hucuktlis data reformatted properly
Hed_sr <- Hed_run_ts_infill |> 
  rename(
    "year" = return_year,
    "S" = escapement,
    "N" = run,
    "H" = catch,
    "age.samples" = age_sample_size
  ) |> 
  select(-Somass_hr) |> 
  rename_with(\(x) str_remove(x, "Hed_")) |> 
  # Combine age compositions by total age
  rowwise() |> 
  mutate(
    N.age.3 = sum(c_across(matches("^age_3."))),
    N.age.4 = sum(c_across(matches("^age_4."))),
    N.age.5 = sum(c_across(matches("^age_5."))),
    N.age.6 = sum(c_across(matches("^age_6."))),
    .keep = "unused"
  )


# Collate the two dataframes
Barkley_sk_sr <- bind_rows(Som_sr, Hed_sr) |> 
  mutate(
    # State CVs for catch and escapement
    H_cv = if_else(is.na(hr_pred_cv), 0.05, hr_pred_cv),
    S_cv = case_when(
      stock != "HED" ~ 0.05,
      stock == "HED" & year < 2012 ~ 0.2,
      stock == "HED" & year >= 2012 ~ 0.1
    )
  ) |> 
  select(-contains("hr_pred"))


# Save metadata for solumn names
Barkley_sk_sr_metadata <- Barkley_sk_sr |> 
  colnames() |> 
  as_tibble_col("column_name") |> 
  mutate(
    definition = case_when(
      column_name == "year" ~ "Run (i.e. observation) year",
      column_name == "stock" ~ paste(
      "Conservation Unit. Either 'Great Central' (GCL),",
      "'Sproat' (SPR), or 'Hucuktlis' (formerly 'Henderson'; HED)"
      ),
      column_name == "S" ~ "Spawners a.k.a escapement",
      str_detect(column_name, "N.age.\\d") ~ paste(
        "Proportion of the run at the given total age (in years)"
      ),
      column_name == "N" ~ "Annual terminal run size. Equal to S +  C",
      column_name == "H" ~ "Annual terminal harvest (i.e. catch)",
      column_name == "age.samples" ~ paste(
        "Number of fish sampled to calculate the age compositions",
        "given in columns labeled 'N.age.#'"
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
        "linear model. CV for these data is calculated as RMSE of the model",
        "residuals divided by the mean of the observed Hucuktlis Sockeye harvest",
        "rates that informed the model fit (i.e. the dependent variable).",
        "Harvest data for Somass and Hucuktlis post-2011 are assumed to be precise."
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
