# Packages ----------------------------------------------------------------


pkgs <- c("here", "readxl", "tidyverse", "broom", "gsl", "rstan")
#install.packages(pkgs)


library(here)
library(tidyverse); theme_set(theme_bw())
library(readxl)
library(broom)
library(gsl)
library(rstan)


# Load pre-smolt abundance and adult spawner data --------------


# Lake abundance estimates from Bayesian model
fry_abun <- here(
  "3. outputs",
  "Stock-recruit data",
  "Bayesian_state-space_smolt-production_estimated_time_series.xlsx"
) |> 
  read_xlsx(sheet = "model_estimates") |> 
  filter(parameter %in% c("N1", "BYO", "BYB")) |> 
  rename("cu" = lake) |> 
  mutate(brood_year = if_else(parameter == "N1", year - 2, year))


# Spawner abundance data
spwn <- here(
  "3. outputs",
  "Stock-recruit data",
  "Barkley_Sockeye_stock-recruit_infilled.xlsx"
) |> 
  read_xlsx(sheet = "S-R data") |> 
  rename("cu" = stock, "brood_year" = year) |> 
  filter(brood_year >= 1972) |> # Remove older years' data
  mutate(
    # Simplify the fertilization data for GCL and SPR
    fertilized = case_when(
      cu == "GCL" ~ 1,
      cu == "SPR" ~ 0,
      cu == "HUC" ~ fertilized
    ) |> factor(),
    # Assume 100% adults for Hucuktlis in years with missing age data
    adult_S = if_else(is.na(adult_S) & cu == "HUC", S, adult_S),
    # Use long names for CUs
    cu = case_when(
      cu == "GCL" ~ "Great Central",
      cu == "SPR" ~ "Sproat",
      cu == "HUC" ~ "Hucuktlis"
    ),
    # Align hatchery fry releases by brood year
    hatchery_fry_release = lead(hatchery_fry_release, 2)
  ) |> 
  select(brood_year, cu, S, adult_S, S_cv, R, fertilized, hatchery_fry_release)


# Join spawner data to fry data
spwn_fry <- left_join(
  fry_abun,
  select(spwn, -R)
) |>
  # Specify that "mu" and "sigma" columns refer to recruits
  rename_with(\(x) paste0(x, "_R"), .cols = c(mu, sigma)) |> 
  mutate(
    hatchery_fry_release = if_else(
      is.na(hatchery_fry_release), 
      0, 
      hatchery_fry_release
    ),
    cu = factor(cu, levels = c("Great Central", "Sproat", "Hucuktlis")),
    # Remove hatchery contributions from Hucuktlis data
    across(
      matches("%"),
      \(x) case_when(
      # When releases > number fry, assume 90% hatchery. 
      # This is following a footnote in Table 13 of the unpublished draft 
      # Henderson Stock Status PSARC report (from 2008). 
      parameter == "N1" & cu == "Hucuktlis" & hatchery_fry_release > x ~ x*0.1, 
      parameter == "N1" ~ x - hatchery_fry_release,
      parameter == "BYB" ~ x/1000, # Convert biomass to kg
      .default = x
      )
    )
  )



# Plot relationships between spawners and pre-smolts ---------------------


# Plots to examine should include `spawners` versus:
# 1) pre-smolt abundance
# 2) pre-smolt biomass
# 3) log(pre-smolt abundance per spawner)
# 4) log(pre-smolt biomass per spawner)
# 
# Also need to examine fits with & without Hucuktlis 1993 outlier removed 


sr_plot_data <- spwn_fry |> 
  # Add rows with adult recruitment as a new parameter
  bind_rows(
    mutate(
      .data = spwn, 
      parameter = "R",
      `50%` = R,
      .keep = "unused"
    )
  ) |> 
  rename(
    "Adult spawners" = adult_S,
    "Spawners" = S
  ) |> 
  mutate(
    # Give parameters informative names
    parameter = case_when(
      parameter == "N1" ~ "Age-1 fry production",
      parameter == "BYO" ~ "Smolt production",
      parameter == "BYB" ~ "Smolt biomass (kg)",
      parameter == "R" ~ "Adult returns"
    )
  ) |> 
  pivot_longer(
    matches("(?i)spawners"),
    values_to = "value"
  )


# Plot spawners versus different measures of recruitment
(sr_plots <- sr_plot_data %>%
  split(.$cu) |> 
  imap(
    \(x, idx) x |> 
      ggplot(aes(x = value, y = `50%`)) +
      facet_grid(
        parameter ~ name,
        scales = "free",
        switch = "both"
      ) +
      geom_linerange(
        aes(
          ymin = `2.5%`,
          ymax = `97.5%`
        ),
        linewidth = 0.3,
        alpha = 0.3
      ) +
      geom_errorbarh(
        aes(
          xmin = value - 1.96*value*S_cv,
          xmax = value + 1.96*value*S_cv,
        )
      ) +
      geom_pointrange(
        aes(
          ymin = `10%`,
          ymax = `90%`,
          fill = fertilized
        ),
        linewidth = 0.5,
        fatten = 3,
        shape = 21,
        stroke = 0.5
      ) +
      #scale_fill_viridis_c() +
      scale_y_continuous(
        limits = c(0, NA),
        labels = scales::label_number(),
        expand = expansion(mult = c(0, 0.05))
      ) +
      scale_x_continuous(
        limits = c(0, NA),
        labels = scales::label_number(),
        expand = expansion(mult = c(0, 0.05))
      ) +
      labs(title = idx) +
      theme(
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.placement = "outside",
        strip.background = element_blank()
      )
  )
)


# R/S version
sr_plot_data |> 
  mutate(
    across(matches("%"), \(x) x/value),
    parameter = paste0(parameter, " per spawner")
  ) %>% 
  split(.$cu) %>%
  {pmap(
    list(
      data = .,
      plot = sr_plots
    ),
    \(data, plot) plot %+% data +
      scale_y_continuous()
  )}



# Examine stock-recruit relationships ---------------------------------------


# Set up a nested dataframe containing 4 outcome variables,
# 2 predictors, and 2 different ways of filtering the data
nested_data <- expand_grid(
  response = c(
    "fry_a", "smolt_a", "smolt_b", 
    "log_frya_spwn", "log_smolta_spwn", "log_smoltb_spwn"
  ),
  predictor = c("spawners", "adult_spawners"),
  fltr = c(0, 1)
) |> 
  mutate(
    data = list(
      spwn_fry |> 
        pivot_wider(
          id_cols = c(cu, brood_year, S, adult_S, fertilized),
          names_from = parameter,
          values_from = `50%` # Use the median estimates for now
        )
    )
  ) |> 
  rowwise() |> 
  mutate(
    data = case_when(
      predictor == "spawners" ~ list(data),
      predictor == "adult_spawners" ~list(mutate(data, S = adult_S))
    ),
    data = case_when(
      response == "fry_a" ~ list(mutate(data, y = N1)),
      response == "smolt_a" ~ list(mutate(data, y = BYO)),
      response == "smolt_b" ~ list(mutate(data, y = BYB)),
      response == "log_frya_spwn" ~ list(mutate(data, y = log(N1/S))),
      response == "log_smolta_spwn" ~ list(mutate(data, y = log(BYO/S))),
      response == "log_smoltb_spwn" ~ list(mutate(data, y = log(BYB/S)))
    ),
    data = list(select(data, brood_year, cu, S, y, fertilized)),
    # Filter the Hucuktlis data to remove 1993 outlier
    data = if_else(
      fltr == 1,
      list(filter(data, !(cu == "Hucuktlis" & S > 150000))),
      list(data)
    ),
    response_long = case_when(
      response == "fry_a" ~ "Age-1 fry production",
      response == "smolt_a" ~ "Total smolt production",
      response == "smolt_b" ~ "Total smolt biomass",
      response == "log_frya_spwn" ~ "log(age-1 fry per spawner)",
      response == "log_smolta_spwn" ~ "log(smolts per spawner)",
      response == "log_smoltb_spwn" ~ "log(smolt biomass per spawner)"
    )
  ) |> 
  ungroup()
  

# Feed the different outcome variables into plots
nested_data |> 
  unnest(data) |>   
  # Filtering is only relevant for the Hucuktlis data (currently)
  filter(!(fltr == 1 & cu != "Hucuktlis")) |> 
  ggplot(
    aes(x = S, y = y, colour = predictor)
  ) +
  facet_grid(
    response_long ~ cu,
    scales = "free",
    switch = "y"
  ) +
  geom_point() +
  geom_smooth(
    aes(lty = factor(fltr)),
    method = "lm",
    alpha = 0.15
  ) +
  scale_x_continuous(
    name = "Spawners",
    labels = scales::label_number()
  ) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    strip.placement = "outside",
    strip.background.y = element_blank(),
    axis.title.y = element_blank()
  )



# Fit simple models and examine fits -------------------------------


# Use the data that were assembled above for plotting
model_fits <- nested_data |> 
  unnest(data) |> 
  # Filtering is only relevant for the Hucuktlis data (currently)
  filter(
    !(fltr == 1 & cu != "Hucuktlis"),
    !str_detect(response, "log")
  ) |> 
  nest(.by = c(cu, predictor, contains("response"), fltr)) |> 
  rowwise() |> 
  mutate(
    # Fit Ricker and Beverton-Holt models
    ricker_model = if(cu == "Hucuktlis")
      list(lm(log(y/S) ~ S + fertilized, data = data)) 
    else 
      list(lm(log(y/S) ~ S, data = data)),
    bevholt_model = if(cu == "Hucuktlis")
      list(
        glm(
          y ~ I(1/S) + fertilized, 
          family = gaussian(link = "inverse"), 
          #start = c(0, 3), # Research how to choose appropriate starting values
          data = data
        )
      )
    else
      list(
        glm(
          y ~ I(1/S), 
          family = gaussian(link = "inverse"), 
          #start = c(0, 3), # Research how to choose appropriate starting values
          data = data
        )
      ),
    # Dataframe for predictions
    newdata = case_when(
      cu == "Hucuktlis" ~ list(
        expand_grid(
          S = seq(
            0.01,
            max(data$S, na.rm = TRUE), 
            length.out = 100
          ),
          fertilized = factor(c(0, 1))
        )
      ),
      cu == "Sproat" ~ list(
        expand_grid(
          S = seq(
            0.01,
            max(data$S, na.rm = TRUE), 
            length.out = 100
          ),
          fertilized = factor(0)
        )
      ),
      cu == "Great Central" ~ list(
        expand_grid(
          S = seq(
            0.01,
            max(data$S, na.rm = TRUE), 
            length.out = 100
          ),
          fertilized = factor(1)
        )
      )
    )
  ) |> 
  mutate(
    ricker_pred = if(cu == "Hucuktlis")
      list(
      MASS::mvrnorm(
        10000, 
        mu = coef(ricker_model), 
        Sigma = as.matrix(vcov(ricker_model))
      ) |> 
        as_tibble() |> 
        rename("a0" = 1, "b" = 2, "a1" = 3) |> 
        expand_grid(newdata) |> 
        mutate(
          a1 = a0 + a1,
          across(matches("^a"), exp),
          b = -b,
          pred_y = if_else(
            fertilized == "0",
            a0*S*exp(-b*S),
            a1*S*exp(-b*S)
          )
        ) |> 
        summarize(
          .by = c(S, fertilized),
          fit = mean(pred_y),
          lwr = quantile(pred_y, 0.025),
          upr = quantile(pred_y, 0.975)
        )
      )
    else 
      list(
        MASS::mvrnorm(
          10000, 
          mu = coef(ricker_model), 
          Sigma = as.matrix(vcov(ricker_model))
        ) |> 
          as_tibble() |> 
          rename("a" = 1, "b" = 2) |> 
          expand_grid(newdata) |> 
          mutate(
            a = exp(a),
            b = -b,
            pred_y = a*S*exp(-b*S)
          ) |> 
          summarize(
            .by = c(S, fertilized),
            fit = mean(pred_y),
            lwr = quantile(pred_y, 0.025),
            upr = quantile(pred_y, 0.975)
          )
      ),
    bevholt_pred = list(
      predict(
        bevholt_model,
        newdata,
        type = "link",
        se.fit = TRUE
      ) |> 
        as_tibble() |> 
        bind_cols(newdata) |> 
        mutate(
          # Equation signs need to switch for LCI and UCI?
          lwr = fit + 1.96*se.fit,
          upr = fit - 1.96*se.fit,
          # Back-transform estimates to response scale
          across(c(fit, lwr, upr), gaussian(link = "inverse")$linkinv)
        )
    )
  ) |> 
  ungroup() |> 
  pivot_longer(
    matches("ricker|bevholt"),
    names_sep = "_",
    names_to = c("type", ".value")
  )


# Plot the various model fits
(model_plots <- model_fits %>%  
  split(.$cu) |> 
  imap(
    function(set, idx) {
      
      sub_set <- list(
        obs = unnest(set, data),
        pred = unnest(set, pred)
      )

      ggplot(
        data = sub_set$obs,
        aes(
          x = S, 
          y = y,
          colour = interaction(fltr, type)
        )
      ) +
        facet_grid(
          response_long ~ predictor,
          scales = "free",
          switch = "both"
        ) +
        geom_point(
          aes(shape = fertilized),
          colour = "black"
        ) +
        geom_line(
          data = sub_set$pred,
          aes(
            y = fit,
            lty = fertilized
          )
        ) +
        scale_shape_manual(values = c(19, 4)) +
        scale_x_continuous(
          name = NULL,
          labels = scales::label_number(),
          limits = c(0, NA),
          expand = expansion(mult = c(0, 0.05))
        ) +
        scale_y_continuous(
          name = NULL,
          limits = c(0, NA),
          expand = expansion(mult = c(0, 0.05))
        ) +
        labs(title = idx) +
        theme(
          strip.placement = "outside",
          strip.background = element_blank(),
          axis.text.x = element_text(angle = 25, hjust = 1)
        )
    }
  )
)


# Save the plotted fits
model_plots |> 
  iwalk(
    \(x, idx) ggsave(
      plot = x,
      filename = here(
        "3. outputs",
        "Plots",
        paste0(idx, "_spawner-smolt_fits.png")
      ),
      dpi = "print"
    )
  )



# Plot the recruits per spawner relationships
(model_plots <- model_fits %>%  
    split(.$cu) |> 
    imap(
      function(set, idx) {
        
        sub_set <- list(
          obs = unnest(set, data),
          pred = unnest(set, pred)
        )
        
        ggplot(
          data = sub_set$obs,
          aes(
            x = S, 
            y = log(y/S),
            colour = interaction(fltr, type)
          )
        ) +
          facet_grid(
            response_long ~ predictor,
            scales = "free",
            switch = "both",
            labeller = labeller(response_long = ~paste0("log(", .x, "/spawner)"))
          ) +
          geom_point(
            aes(shape = fertilized),
            colour = "black"
          ) +
          geom_line(
            data = sub_set$pred,
            aes(
              y = log(fit/S),
              lty = fertilized
            )
          ) +
          scale_shape_manual(values = c(19, 4)) +
          scale_x_continuous(
            name = NULL,
            limits = c(0, NA),
            labels = scales::label_number(),
            expand = expansion(mult = c(0, 0.05))
          ) +
          scale_y_continuous(
            name = NULL,
            expand = expansion(mult = c(0, 0.05))
          ) +
          labs(title = idx) +
          theme(
            strip.placement = "outside",
            strip.background = element_blank(),
            axis.text.x = element_text(angle = 25, hjust = 1)
          )
      }
    )
)

  

# Fit Beverton-Holt AR1 models to Somass data using Stan ------------------


# Pull out the Hucuktlis data, which require a different model
somass_sr <- spwn_fry |> 
  filter(cu != "Hucuktlis") |> 
  droplevels()
  

# Make input data for Stan
make_stan_data_somass <- function(cu_name) {
  
  cu_data <- somass_sr |> 
    filter(
      cu == cu_name,
      !if_any(c(mu_R, sigma_R, adult_S, S_cv), is.na),
      parameter == "BYO" # Use brood year smolt production as recruitment measure
    ) |> 
    mutate(
      S_sd = adult_S * S_cv,
      mu_S = log(adult_S^2 / sqrt(S_sd^2 + adult_S^2)),
      sigma_S = sqrt(log(1 + (S_sd^2 / adult_S^2)))
    )
  
  # Number of years in the time series
  Y <- nrow(cu_data)
  
  # Best estimates of annual spawner abundance, in log space
  S_obs <- cu_data$mu_S
  
  # Best estimates of annual smolt recruitment, in log space
  R_obs <- cu_data$mu_R
  
  # Annual standard deviation of spawner abundance, in log space
  sigma_S_obs <- cu_data$sigma_S
  
  # Annual standard deviation of smolt recruitment, in log space
  sigma_R_obs <- cu_data$sigma_R
  
  stan_data <- list(
    Y = Y,
    S_obs = S_obs,
    R_obs = R_obs,
    sigma_S_obs = sigma_S_obs,
    sigma_R_obs = sigma_R_obs
  )
  
  return(stan_data)
  
}


# Save the stan data for each CU
somass_stan_data <- levels(somass_sr$cu) |> 
  purrr::set_names() |> 
  map(make_stan_data_somass)


# Function to fit Stan model 
fit_stan_somass <- function(stan_data, cu) {
  stan(
    file = here(
      "2. code",
      "Stock-recruit modelling",
      "Stan",
      "SR_Bev-Holt_AR1.stan"
    ), 
    model_name = cu,
    data = stan_data, 
    iter = 3000, 
    chains = 4, 
    warmup = 1000
  )
}


# Fit the models for both CUs
if(
  FALSE
  #TRUE
) {
  somass_stan_fits <- somass_stan_data |> 
    imap(fit_stan_somass)}


# Save fitted models as RDS objects (toggle to TRUE to run)
if(
  #TRUE
  FALSE
) {
  somass_stan_fits |> 
    iwalk(
      \(x, idx) x |> 
        saveRDS(
          file = here(
            "3. outputs",
            "Stock-recruit modelling",
            paste("Bayesian Beverton-Holt", idx, "spawner-smolt.rds")
          )
        )
    )
}


# Load fitted models from RDS files (if not already in global environment)
somass_stan_fits <- if(exists("somass_stan_fits")) {
  somass_stan_fits
} else {
  levels(somass_sr$cu) |>
    purrr::set_names() |> 
    map(
      \(x) list.files(
        here(
          "3. outputs",
          "Stock-recruit modelling"
        ),
        pattern = paste("Bayesian Beverton-Holt", x, "spawner-smolt.rds"),
        full.names = TRUE
      )
    ) |> 
    map(readRDS)
}


# Assess convergence
worst_Rhat <- somass_stan_fits |> 
  map(\(x) summary(x)$summary) |>  
  map(as.data.frame) |> 
  list_rbind(names_to = "lake") |> 
  mutate(Rhat = round(Rhat, 3)) |> 
  arrange(desc(Rhat))

worst_Rhat %>% 
  filter(n_eff>3) %>% 
  ggplot(aes(x = n_eff, y = Rhat))+
  facet_wrap(~lake) +
  geom_point()+
  geom_hline(yintercept = 1.01, lty = 2)+
  geom_vline(xintercept = 400, lty = 2)

head(worst_Rhat, n = 20)


# Trace plots
somass_stan_fits |> 
  imap(
    \(x, idx)
    traceplot(
      x, 
      pars = rownames(filter(worst_Rhat, lake == idx))[1:20]
    ) +
      ggtitle(label = idx)
  )


# Pair plots
somass_stan_fits |> 
  map(\(x) pairs(x, pars = c("alpha", "beta")))


# Explore the posterior

# Dataframe with first year in time series for each lake
# used to properly code years in the Stan outputs
min_yrs <- somass_sr |> 
  summarize(
    .by = cu,
    min_yr = min(year)
  )


# Posterior values that are not year-specific
post_interannual <- somass_stan_fits |> 
  map(extract) |> 
  # Discard all parameters with year-specific estimates
  map(\(x) discard(.x = x, .p = \(y) is.matrix(y))) |> 
  map(\(x) map(x, \(y) as_tibble(y, rownames = "draw"))) |> 
  map(\(x) list_rbind(x, names_to = "parameter")) |> 
  list_rbind(names_to = "cu") 


# Posterior values that are year-specific
post_annual <- somass_stan_fits |> 
  map(extract) |> 
  # Keep only parameters with year-specific estimates
  map(\(x) keep(.x = x, .p = \(y) is.matrix(y))) |> 
  map(\(x) map(x, \(y) as_tibble(y, .name_repair = NULL, rownames = "draw"))) |> 
  map(\(x) list_rbind(x, names_to = "parameter")) |> 
  list_rbind(names_to = "cu") |> 
  pivot_longer(
    cols = !c(cu, parameter, draw),
    names_to = "year",
    values_to = "value",
    names_transform = \(x) str_extract(x, "\\d+")
  ) |> 
  left_join(min_yrs) |> 
  mutate(year = as.numeric(year) - 1 + min_yr) |> 
  select(-min_yr)


# All posterior values
posterior_df <- bind_rows(
  post_annual, 
  post_interannual
) 


# Plot posterior estimates for spawners and recruits
posterior_df |> 
  filter(parameter %in% c("S_true", "R_true")) |> 
  pivot_wider(
    names_from = parameter,
    values_from = value
  ) |> 
  ggplot(aes(S_true, R_true, colour = year)) +
  facet_wrap(~cu) +
  stat_summary(
    geom = "linerange",
    fun.data = mean_cl_boot,
    orientation = "x"
  ) +
  stat_summary(
    geom = "pointrange",
    fun.data = mean_cl_boot,
    orientation = "y"
  )

