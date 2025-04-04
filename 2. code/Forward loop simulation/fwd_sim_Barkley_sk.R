# Packages ----------------------------------------------------------------


pkgs <- c("here", "tidyverse", "readxl", "rstan")
#install.packages(pkgs)


library(here)
library(tidyverse); theme_set(theme_bw())
library(readxl)



# Load in data and fitted models ------------------------------------------


sr_data <- here(
  "3. outputs",
  "Stock-recruit data",
  "Barkley_Sockeye_stock-recruit_infilled.xlsx"
) |> 
  read_xlsx(sheet = "S-R data")


AR1_frame <- list.files(
  here(
    "3. outputs",
    "Stock-recruit modelling"
  ),
  pattern = "AR1.rds",
  full.names = TRUE
) |> 
  set_names(nm = ~str_extract(.x, "(SPR|GCL|HUC).*_AR1")) |> 
  map(readRDS) |> 
  enframe(name = "spec", value = "model") 
  
  

# Extract data required for the forward simulation ------------------------

  
  
sim_params <- AR1_frame |> 
  # Start with a list of model parameters whose posteriors we need
  # to extract
  mutate(
    extract = list(
      str_subset(
        names(model), 
        "lnalpha|beta|^S|^R|^C|U"
      )
    )
  ) |> 
  unnest_longer(extract) |> 
  # Filter parameters to ensure those estimated on an annual basis
  # are extracted only for the final year in the time series
  mutate(parameter = str_remove_all(extract, "\\[\\d+\\]")) |> 
  mutate(
    .by = c(spec, parameter),
    yr = as.numeric(str_extract(extract, "\\d+")),
    max_yr = max(yr, na.rm = TRUE),
    keep = if_else(
      !is.na(yr) & yr != max_yr,
      FALSE,
      TRUE
    )
  ) |> 
  filter(keep) |> 
  select(-yr, -max_yr, -keep) |> 
  rowwise() |> 
  mutate(
    posterior = list(as.data.frame(rstan::extract(model, pars = extract))),
    stock = str_extract(spec, "GCL|SPR|HUC"),
    fert = case_when(
      str_detect(parameter, "_fert") ~ 1, 
      str_detect(parameter, "_unfert") ~ 0,
      .default = NA
    ),
    data_scope = if_else(str_detect(spec, "trim"), "trim", "full"),
    parameter = str_remove_all(parameter, "_.*")
  )
