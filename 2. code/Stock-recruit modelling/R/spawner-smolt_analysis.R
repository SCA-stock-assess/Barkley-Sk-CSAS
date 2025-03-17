# Packages ----------------------------------------------------------------


pkgs <- c("here", "readxl", "tidyverse", "broom")
#install.packages(pkgs)


library(here)
library(tidyverse); theme_set(theme_bw())
library(readxl)
library(broom)



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
  select(brood_year, cu, S, adult_S, R, fertilized, hatchery_fry_release)


# Join spawner data to fry data
spwn_fry <- left_join(
  fry_abun,
  select(spwn, -R)
) |>
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



# Fit simple Ricker models and examine fits -------------------------------


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
    newdata = if_else(
      cu == "Hucuktlis",
      list(
        expand_grid(
          S = seq(
            0.01,
            #min(data$S, na.rm = TRUE),
            max(data$S, na.rm = TRUE), 
            length.out = 100
          ),
          fertilized = factor(c(0, 1))
        )
      ),
      list(
        tibble(
          S = seq(
            0.01,
            #min(data$S, na.rm = TRUE),
            max(data$S, na.rm = TRUE), 
            length.out = 100
          )
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
            .by = c(S),
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
(model_plots <- model_fits |>  
  filter(cu != "Hucuktlis") %>% 
  split(.$predictor) |> 
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
          colour = type,
          fill = type
        )
      ) +
        facet_grid(
          response_long ~ cu,
          scales = "free",
          switch = "y"
        ) +
        geom_point(colour = "black") +
        geom_line(
          data = sub_set$pred,
          aes(y = fit)
        ) +
        scale_x_continuous(
          name = idx,
          labels = scales::label_number(),
          expand = expansion(mult = c(0, 0.05))
        ) +
        scale_y_continuous(
          name = NULL,
          #limits = c(0, max(sub_set$pred$upr)),
          expand = expansion(mult = c(0, 0.05))
        ) +
        theme(
          strip.placement = "outside",
          strip.background.y = element_blank(),
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
        paste0("Pre-smolts per ", idx, ".png")
      ),
      dpi = "print"
    )
  )
