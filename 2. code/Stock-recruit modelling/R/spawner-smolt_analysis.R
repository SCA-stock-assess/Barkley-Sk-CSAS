# Packages ----------------------------------------------------------------


pkgs <- c(
  "here", "readxl", "tidyverse", "broom", "gsl", "rstan", "tidybayes", 
  "bayesplot", "cowplot"
)
#install.packages(pkgs)


library(here)
library(tidyverse); theme_set(theme_bw())
library(cowplot)
library(readxl)
library(broom)
library(gsl)
library(rstan)
library(tidybayes)
library(bayesplot)


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
    # Align hatchery fry releases by brood year
    hatchery_fry_release = lead(hatchery_fry_release, 2),
    # Binary variable distinguishing enhanced years
    enhanced = case_when(
      fertilized == 1 ~ 1,
      hatchery_fry_release > 0 ~ 1,
      .default = 0
    ),
    # Assume 100% adults for Hucuktlis in years with missing age data
    adult_S = if_else(is.na(adult_S) & cu == "HUC", S, adult_S),
    # Use long names for CUs
    cu = case_when(
      cu == "GCL" ~ "Great Central",
      cu == "SPR" ~ "Sproat",
      cu == "HUC" ~ "Hucuktlis"
    )
  ) |> 
  # Ensure enhancement variable is correctly aligned to the brood year it affects
  mutate(
    .by = cu,
    enhanced = if_else(
      cu == "Great Central",
      1,
      lead(enhanced, 1L, default = 0)
    ) |> factor()
  ) |> 
  select(brood_year, cu, S, adult_S, S_cv, R, enhanced, hatchery_fry_release)


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
      #parameter == "N1" & cu == "Hucuktlis" & hatchery_fry_release > x ~ x*0.1, 
      #parameter == "N1" ~ x - hatchery_fry_release,
      parameter == "BYB" ~ x/1000, # Convert biomass to kg
      .default = x
      )
    )
  ) |> 
  arrange(cu, parameter, year)
  

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
    # Put smolts on reasonable scale,
    across(matches("%"), \(x) if_else(parameter %in% c("N1", "BYO"), x/1000, x)),
    # Give parameters informative names
    parameter = case_when(
      parameter == "N1" ~ "Age-1 fry production (1000s)",
      parameter == "BYO" ~ "Smolt production (1000s)",
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
          # 90% CIs
          xmin = value - 1.65*value*S_cv,
          xmax = value + 1.65*value*S_cv,
        )
      ) +
      geom_pointrange(
        aes(
          ymin = `10%`,
          ymax = `90%`,
          fill = enhanced
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


# plot all CUs together versus just adult spawners
(cu_sr_plots <- sr_plot_data |> 
    filter(
      name == "Adult spawners",
      parameter != "Adult returns"
    ) |> 
    mutate(group = if_else(cu == "Hucuktlis", "Hucuktlis", "Somass")) %>%
    split(.$group) |> 
    map(
      \(x) x |> 
        ggplot(aes(value, `50%`)) +
        facet_grid(
          parameter ~ cu,
          scales = "free",
          switch = "y"
        ) +
        geom_errorbarh(
          aes(
            # 90% CIs
            xmin = value - 1.65*value*S_cv,
            xmax = value + 1.65*value*S_cv,
          ),
          colour = "grey50",
          linewidth = 0.25
        ) +
        geom_pointrange(
          aes(
            ymin = `10%`,
            ymax = `90%`,
            fill = enhanced
          ),
          linewidth = 0.25,
          colour = "grey50",
          fatten = 3,
          shape = 21,
          stroke = 0.25
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
        scale_fill_manual(values = c("blue", "red")) +
        theme(
          axis.title = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.placement = "outside",
          strip.background = element_blank(),
          panel.grid = element_blank()
        )
    )
)


# Assemble plots in a grid
grid1 <- plot_grid(
  cu_sr_plots[[1]] +
    guides(fill = "none"), 
  cu_sr_plots[[2]] +
    theme(
      strip.text.y = element_blank(),
      strip.background.y = element_blank()
    ),
  rel_widths = c(1, 2)
)


# Add shared x axis label
(cu_sr_grid <- plot_grid(
  grid1,
  ggdraw() +
    draw_label(
      "Brood year adult spawners",
      size = 12
    ),
  ncol = 1,
  rel_heights = c(1, 0.075)
)
)


# Save the gridded plots 
ggsave(
  cu_sr_grid,
  filename = here(
    "3. outputs",
    "Plots",
    "Spawner-smolt_relationships_allCUs.png"
  ),
  width = 9,
  height = 7,
  units = "in",
  dpi = "print"
)


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
          id_cols = c(cu, brood_year, S, adult_S, enhanced),
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
    data = list(select(data, brood_year, cu, S, y, enhanced)),
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
      list(lm(log(y/S) ~ S + enhanced, data = data)) 
    else 
      list(lm(log(y/S) ~ S, data = data)),
    bevholt_model = if(cu == "Hucuktlis")
      list(
        glm(
          y ~ I(1/S) + enhanced, 
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
          enhanced = factor(c(0, 1))
        )
      ),
      cu == "Sproat" ~ list(
        expand_grid(
          S = seq(
            0.01,
            max(data$S, na.rm = TRUE), 
            length.out = 100
          ),
          enhanced = factor(0)
        )
      ),
      cu == "Great Central" ~ list(
        expand_grid(
          S = seq(
            0.01,
            max(data$S, na.rm = TRUE), 
            length.out = 100
          ),
          enhanced = factor(1)
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
            enhanced == "0",
            a0*S*exp(-b*S),
            a1*S*exp(-b*S)
          )
        ) |> 
        summarize(
          .by = c(S, enhanced),
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
            .by = c(S, enhanced),
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
          aes(shape = enhanced),
          colour = "black"
        ) +
        geom_line(
          data = sub_set$pred,
          aes(
            y = fit,
            lty = enhanced
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
            aes(shape = enhanced),
            colour = "black"
          ) +
          geom_line(
            data = sub_set$pred,
            aes(
              y = log(fit/S),
              lty = enhanced
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


# Simplify the recruits/spawner plots
log_rs_fits <- model_fits |> 
  filter(
    predictor == "adult_spawners",
    response != "fry_a",
    fltr == 0
  ) |> 
  (function(nested_data) {
    sub_set =  list(
      obs = unnest(nested_data, data),
      pred = unnest(nested_data, pred) |> 
        mutate(enhanced = if_else(cu != "Hucuktlis", NA, enhanced))
    )
    
    plot <- ggplot(
      data = sub_set$obs,
      aes(
        x = S, 
        y = log(y/S),
        colour = enhanced
      )
    ) +
      facet_grid(
        response_long ~ cu,
        scales = "free",
        switch = "y",
        labeller = labeller(response_long = ~paste0("log(", .x, "/spawner)"))
      ) +
      geom_point(alpha = 0.25) +
      geom_line(
        data = sub_set$pred,
        aes(
          y = log(fit/S),
          lty = type
        )
      ) +
      scale_x_continuous(
        limits = c(0, NA),
        labels = scales::label_number(),
        expand = expansion(mult = c(0, 0.05))
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      scale_colour_manual(values = c("blue", "red")) +
      labs(
        y = NULL,
        x = "Adult spawners",
        lty = "S-R fit",
        colour = "Enhancement\nstate"
      ) +
      theme(
        strip.placement = "outside",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1),
        panel.spacing.x = unit(0.75, "lines"),
        panel.grid.minor = element_blank()
      )
    
    return(plot)
  }
)()


# Save the plot for presentation in Appendix A
ggsave(
  plot = log_rs_fits,
  filename = here(
    "3. outputs",
    "Plots",
    "log_smolts-per-spanwer_exploratory_plots.png"
  ),
  width = 6.5, 
  height = 5,
  units = "in",
  dpi = "print"
)

# Fit Bayesian Beverton-Holt models to Somass data using Stan ------------------


# Prepare data from Somass CUs for modelling
somass_sr <- spwn_fry |> 
  filter(
    # Pull out the Hucuktlis data, which require a different model
    cu != "Hucuktlis",
    !if_any(c(mu_R, sigma_R, adult_S, S_cv), is.na)
  ) |> 
  # Clean up column names for recruitment values
  rename_with(
    \(x) str_remove(paste0("R_", x), "%"),
    .cols = matches("%")
  ) |> 
  # Convert biomass back to g to match mu values
  mutate(across(matches("R_\\d+"), \(x) if_else(parameter == "BYB", x*1000, x))) |> 
  droplevels()
  

# Make input data for Stan
make_stan_data_somass <- function(cu_name, R_param) {
  
  cu_data <- somass_sr |> 
    filter(
      cu == cu_name,
      parameter == R_param
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
  
  # Best estimates of annual recruitment, in log space
  R_obs <- cu_data$mu_R
  
  # Annual standard deviation of spawner abundance, in log space
  sigma_S_obs <- cu_data$sigma_S
  
  # Annual standard deviation of recruitment, in log space
  sigma_R_obs <- cu_data$sigma_R
  
  # Alpha prior: plausible maximum number of recruits per spawner
  alpha_prior <- cu_data |> 
    slice_min(
      order_by = adult_S,
      prop = 0.15
    ) |> 
    mutate(rs = R_50/adult_S) |> 
    summarize(alpha_prior = max(rs)) |> 
    pull(alpha_prior)
  
  # Beta prior: plausible maximum number of recruits
  Rmax_prior <- max(cu_data$R_50)
  
  # Alpha variability prior (on log scale)
  #sigma_alpha_prior <- 0.5
  
  # Beta variability prior (on log scale)
  #sigma_Rmax_prior <- 0.5
  
  stan_data <- list(
    Y = Y,
    S_obs = S_obs,
    R_obs = R_obs,
    sigma_S_obs = sigma_S_obs,
    sigma_R_obs = sigma_R_obs,
    #sigma_alpha_prior = sigma_alpha_prior,
    #sigma_beta_prior = sigma_beta_prior,
    alpha_prior = alpha_prior,
    Rmax_prior = Rmax_prior
  )
  
  return(stan_data)
  
}


# Save the stan data for each CU
somass_stan_data <- expand_grid(
  cu_name = levels(somass_sr$cu),
  R_param = c("BYO", "BYB")
) |> 
  # Ensure pmap captures names correctly for each model
  mutate(cu_name = set_names(cu_name, paste(cu_name, R_param, sep = "_"))) |> 
  pmap(make_stan_data_somass)


# Function to fit Stan model 
fit_stan_somass <- function(stan_data, cu) {
  stan(
    file = here(
      "2. code",
      "Stock-recruit modelling",
      "Stan",
      "SR_Bev-Holt_Somass.stan"
    ), 
    model_name = cu,
    data = stan_data, 
    #warmup = 1000,
    iter = 6000, 
    chains = 4, 
    control = list(
      max_treedepth = 15,
      adapt_delta = 0.997
    )
  )
}


# Fit the models for both CUs
if(
  FALSE
  #TRUE
) {
  somass_stan_fits <- somass_stan_data |> 
    imap(fit_stan_somass)
}


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
  names(somass_stan_data) |>
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



# Explore Somass models' posteriors ---------------------------------------


# Posterior density distributions
somass_stan_fits |> 
  map(\(x) mcmc_areas(x, pars = c("k_R", "tau_R", "sigma_R_proc", "sigma_S")))


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
  list_rbind(names_to = "cu_Rmeas") 


# Posterior values that are year-specific
post_annual <- somass_stan_fits |> 
  map(extract) |> 
  # Keep only parameters with year-specific estimates
  map(\(x) keep(.x = x, .p = \(y) is.matrix(y))) |> 
  map(\(x) map(x, \(y) as_tibble(y, .name_repair = NULL, rownames = "draw"))) |> 
  map(\(x) list_rbind(x, names_to = "parameter")) |> 
  list_rbind(names_to = "cu_Rmeas") |> 
  pivot_longer(
    cols = !c(cu_Rmeas, parameter, draw),
    names_to = "year",
    values_to = "value",
    names_transform = \(x) str_extract(x, "\\d+")
  )
  


# All posterior values
posterior_df <- bind_rows(
  post_annual, 
  post_interannual
) |> 
  separate(cu_Rmeas, sep = "_", c("cu", "Rmeas")) |> # This step takes curiously long
  left_join(min_yrs) |> 
  mutate(year = as.numeric(year) - 1 + min_yr) |> 
  select(-min_yr)


# Plot time series of AR1 process variation SD
posterior_df |> 
  filter(parameter == "v") |> 
  summarize(.by = c(cu, Rmeas, year), v_sd = sd(value)) |> 
  ggplot(aes(year, v_sd)) +
  facet_grid(Rmeas ~ cu) +
  geom_point()


# Observed values for spawners and recruits
obs_sr <- somass_sr |> 
  mutate(
    S_50 = adult_S,
    S_10 = adult_S - adult_S*S_cv*1.28,
    S_90 = adult_S + adult_S*S_cv*1.28
  ) |> 
  filter(
    parameter %in% unique(posterior_df$Rmeas),
    cu != "Hucuktlis"
  ) |> 
  select(cu, Rmeas = parameter, year, matches("(S|R)_\\d{2}")) |> 
  pivot_longer(
    matches("(S|R)_\\d{2}"),
    names_sep = "_",
    names_to = c("parameter", "quantile")
  ) |> 
  pivot_wider(names_from = quantile) |> 
  mutate(set = "obs")


# Plot posterior estimates versus observations for spawners and recruits
posterior_df |> 
  filter(parameter %in% c("S_true", "R_true")) |> 
  pivot_wider(
    names_from = parameter,
    values_from = value
  ) |> 
  rename_with(\(x) str_remove(x, "_true")) |> 
  summarize(
    .by = c(cu, Rmeas, year),
    across(
      c(S, R),
      .fns = list(
        "50" = median,
        "10" = ~quantile(.x, 0.1, na.rm = TRUE),
        "90" = ~quantile(.x, 0.9, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    )
  ) |> 
  pivot_longer(
    matches("(S|R)_\\d{2}"),
    names_sep = "_",
    names_to = c("parameter", "quantile")
  ) |> 
  pivot_wider(names_from = quantile) |> 
  mutate(set = "post") |> 
  bind_rows(obs_sr) %>% 
  split(.$Rmeas) |> 
  imap(
    \(x, idx) x |> 
      mutate(parameter = if_else(parameter == "R", idx, parameter)) |> 
      ggplot(aes(year, `50`, colour = set)) +
      facet_grid(
        parameter ~ cu, 
        switch = "y",
        scales = "free_y"
      ) +
      geom_pointrange(
        aes(
          ymin = `10`,
          ymax = `90`
        ),
        position = position_dodge(width = 0.5),
        alpha = 0.5
      ) +
      theme(
        strip.background.y = element_blank(),
        axis.title.y = element_blank(),
        strip.placement = "outside"
      )
  )
# Still some shrinkage occurring in the posterior recruitment estimates?

# Plot predicted Beverton-Holt curve versus latent states data

# Extract posterior samples of alpha and beta
somass_a_b_draws <- posterior_df |> 
  filter(parameter %in% c("alpha", "beta")) |> 
  pivot_wider(names_from = parameter)


# Create CU-specific range of spawners to predict across
S_pred <- somass_sr |> 
  summarize(
    .by = cu,
    min_S = 0,
    max_S = max(adult_S)
  ) |> 
  rowwise() |> 
  mutate(S = list(seq(min_S, max_S, length.out = 100))) |> 
  unnest(S) |> 
  select(cu, S)


# Join range of predicted spawners to alpha and beta draws
somass_pred_frame <- somass_a_b_draws |> 
  left_join(
    S_pred,
    by = "cu",
    relationship = "many-to-many"
  )

# Check that the resulting dataframe has the correct number of rows
stopifnot(
  all.equal(
    nrow(somass_pred_frame), 
    nrow(S_pred)/length(unique(S_pred$cu))*nrow(somass_a_b_draws)
    )
  )


# Latent states of spawners and recruits
latent_sr <- posterior_df |> 
  filter(
    parameter %in% c("S_true", "R_true"),
    !is.na(value)
  ) |> 
  mutate(value = value / 1000) |> 
  summarize(
    .by = c(cu, Rmeas, year, parameter),
    value = list(quantile(value, c(0.10, 0.50, 0.90)))
  ) |> 
  unnest_wider(value) |> 
  pivot_wider(
    names_from = parameter,
    values_from = matches("%"),
    names_glue = "{parameter}_{str_remove(.value, '%')}"
  ) |> 
  left_join(
    somass_sr |> 
      filter(parameter != "N1") |> 
      distinct(cu, year, enhanced)
  )


# Generate predictions for each spawner value
(somass_bevholt_fits <- somass_pred_frame |> 
    mutate(
      R_pred = ((alpha * S) / (1 + (alpha/beta)*S))/1000,
      S = S/1000
    ) |> 
    summarize(
      .by = c(cu, Rmeas, S),
      R_quant = list(quantile(R_pred, c(0.025, 0.10, 0.50, 0.90, 0.975)))
    ) |> 
    unnest_wider(R_quant) |>
    rename_with(
      \(x) str_remove(paste0("R_pred_", x), "%"),
      .cols = matches("%")
    ) |> 
    ggplot(aes(x = S, y = R_pred_50)) +
    facet_grid(
      Rmeas ~ cu, 
      scales = "free",
      switch = "y",
      labeller = labeller(
        Rmeas = c(
          "BYB" = "Smolt biomass (kg)",
          "BYO" = "Smolt abundance (1000s)"
        )
      )
    ) +
    geom_line() +
    geom_ribbon(
      aes(
        ymin = R_pred_10,
        ymax = R_pred_90
      ),
      alpha = 0.25
    ) +
    geom_linerange(
      data = latent_sr,
      aes(
        y = R_true_50,
        ymin = R_true_10,
        ymax = R_true_90,
        x = S_true_50
      ),
      colour = "grey50"
    ) +
    geom_pointrange(
      data = latent_sr,
      aes(
        y = R_true_50,
        x = S_true_50,
        xmin = S_true_10,
        xmax = S_true_90,
        fill = enhanced
      ),
      shape = 21,
      colour = "grey50"
    ) +
    scale_x_continuous(
      name = "Adult spawners (1000s)",
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05)),
      labels = scales::label_number()
    ) +
    scale_y_continuous(
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05)),
      labels = scales::label_number()
    ) +
    scale_fill_manual(values = c("blue", "red")) +
    theme(
      strip.background = element_blank(),
      axis.title.y = element_blank(),
      strip.placement = "outside",
      panel.spacing.y = unit(1, "lines"),
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
)



# Fit Bayesian Beverton-Holt models to Hucuktlis data using Stan ----------


# Prepare data from Hucuktlis for modelling
hucuktlis_sr <- spwn_fry |> 
  filter(
    cu == "Hucuktlis",
    !if_any(c(mu_R, sigma_R, adult_S, S_cv), is.na)
  ) |> 
  # Clean up column names for recruitment values
  rename_with(
    \(x) str_remove(paste0("R_", x), "%"),
    .cols = matches("%")
  ) |> 
  # Convert biomass back to g to match mu values
  mutate(across(matches("R_\\d+"), \(x) if_else(parameter == "BYB", x*1000, x))) |> 
  droplevels()


# Make input data for Stan
make_stan_data_hucuktlis <- function(R_param, ...) {
  
  cu_data <- hucuktlis_sr |> 
    filter(parameter == R_param) |> 
    mutate(
      S_sd = adult_S * S_cv, 
      mu_S = log(adult_S^2 / sqrt(S_sd^2 + adult_S^2)),
      sigma_S = sqrt(log(1 + (S_sd^2 / adult_S^2)))
    )
  
  # Number of years in the time series
  Y <- nrow(cu_data)
  
  # Binary variable indicating enhanced years
  enhanced <- as.integer(as.character(cu_data$enhanced))
  
  # Best estimates of annual spawner abundance, in log space
  S_obs <- cu_data$mu_S
  
  # Best estimates of annual recruitment, in log space
  R_obs <- cu_data$mu_R
  
  # Annual standard deviation of spawner abundance, in log space
  sigma_S_obs <- cu_data$sigma_S
  
  # Annual standard deviation of recruitment, in log space
  sigma_R_obs <- cu_data$sigma_R
  
  # Alpha prior: plausible maximum number of recruits per spawner at 
  # low spawner abundances
  alpha_prior <- cu_data |> 
    slice_min(
      order_by = adult_S,
      prop = 0.15 # Keep the lowest 15% of spawner abundances
    ) |> 
    mutate(rs = R_50/adult_S) |> 
    summarize(
      .by = enhanced,
      alpha_prior = max(rs)
    )
  
  alpha_enh_prior <- alpha_prior[alpha_prior$enhanced == 1,]$alpha_prior
  alpha_noenh_prior <- alpha_prior[alpha_prior$enhanced == 0,]$alpha_prior
  
  
  # Beta priors: plausible maximum number of recruits
  Rmax <- cu_data |> 
    filter(R_50 < max(R_50)) |> # Exclude the highest year (mitigate beta bias) 
    summarize(
      .by = enhanced,
      Rmax = max(R_50)
    )
  
  Rmax_enh_prior <- Rmax[Rmax$enhanced == 1,]$Rmax
  Rmax_noenh_prior <- Rmax[Rmax$enhanced == 0,]$Rmax
  
  
  # Alpha variability prior (on log scale)
  #sigma_alpha_prior <- 0.2
  
  # Beta variability prior (on log scale)
  #sigma_beta_prior <- 0.4
  
  stan_data <- list(
    Y = Y,
    enhanced = enhanced,
    S_obs = S_obs,
    R_obs = R_obs,
    sigma_S_obs = sigma_S_obs,
    sigma_R_obs = sigma_R_obs,
    #sigma_alpha_prior = sigma_alpha_prior,
    #sigma_beta_prior = sigma_beta_prior,
    alpha_enh_prior = alpha_enh_prior,
    alpha_noenh_prior = alpha_enh_prior,
    Rmax_enh_prior = Rmax_enh_prior,
    Rmax_noenh_prior = Rmax_noenh_prior
  )
  
  return(stan_data)
  
}


# Save the Stan data for Hucuktlis smolt abundance and biomass
hucuktlis_stan_data <- expand_grid(
  cu_name = levels(hucuktlis_sr$cu),
  R_param = c("BYO", "BYB")
) |> 
  # Ensure pmap captures names correctly for each model
  mutate(cu_name = set_names(cu_name, paste(cu_name, R_param, sep = "_"))) |> 
  pmap(make_stan_data_hucuktlis)


# Function to fit Stan model 
fit_stan_hucuktlis <- function(stan_data, cu = "hucuktlis") {
  stan(
    file = here(
      "2. code",
      "Stock-recruit modelling",
      "Stan",
      "SR_Bev-Holt_Hucuktlis.stan"
    ), 
    model_name = cu,
    data = stan_data, 
    iter = 6000, 
    chains = 4, 
    control = list(
      max_treedepth = 15,
      adapt_delta = 0.997
    )
  )
}


# Fit the models for Hucuktlis
if(
  FALSE
  #TRUE
) {
  hucuktlis_stan_fits <- hucuktlis_stan_data |> 
    imap(fit_stan_hucuktlis)
}


# Save fitted models as RDS objects (toggle to TRUE to run)
if(
  #TRUE
  FALSE
) {
  hucuktlis_stan_fits |> 
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
hucuktlis_stan_fits <- if(exists("hucuktlis_stan_fits")) {
  hucuktlis_stan_fits
} else {
  names(hucuktlis_stan_data) |>
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
worst_Rhat <- hucuktlis_stan_fits |> 
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
hucuktlis_stan_fits |> 
  imap(
    \(x, idx)
    traceplot(
      x, 
      pars = rownames(filter(worst_Rhat, lake == idx))[1:20]
    ) +
      ggtitle(label = idx)
  )


# Pair plots
hucuktlis_stan_fits |> 
  map(\(x) pairs(x, pars = c("alpha_noenh", "beta_noenh")))

# Pair plots
hucuktlis_stan_fits |> 
  map(\(x) pairs(x, pars = c("alpha_enh", "beta_enh")))


# Explore Hucuktlis models' posteriors ------------------------------------

# Posterior density distributions
hucuktlis_stan_fits |> 
  map(\(x) mcmc_areas(x, pars = c("k_R", "tau_R", "sigma_R_proc", "sigma_S")))


# Dataframe with first year in time series for each lake
# used to properly code years in the Stan outputs
min_yrs_huc <- hucuktlis_sr |> 
  summarize(
    .by = cu,
    min_yr = min(year)
  )


# Posterior values that are not year-specific
post_interannual_huc <- hucuktlis_stan_fits |> 
  map(extract) |> 
  # Discard all parameters with year-specific estimates
  map(\(x) discard(.x = x, .p = \(y) is.matrix(y))) |> 
  map(\(x) map(x, \(y) as_tibble(y, rownames = "draw"))) |> 
  map(\(x) list_rbind(x, names_to = "parameter")) |> 
  list_rbind(names_to = "cu_Rmeas") 


# Posterior values that are year-specific
post_annual_huc <- hucuktlis_stan_fits |> 
  map(extract) |> 
  # Keep only parameters with year-specific estimates
  map(\(x) keep(.x = x, .p = \(y) is.matrix(y))) |> 
  map(\(x) map(x, \(y) as_tibble(y, .name_repair = NULL, rownames = "draw"))) |> 
  map(\(x) list_rbind(x, names_to = "parameter")) |> 
  list_rbind(names_to = "cu_Rmeas") |> 
  pivot_longer(
    cols = !c(cu_Rmeas, parameter, draw),
    names_to = "year",
    values_to = "value",
    names_transform = \(x) str_extract(x, "\\d+")
  )


# All posterior values
posterior_df_huc <- bind_rows(
  post_annual_huc, 
  post_interannual_huc
) |> 
  separate(cu_Rmeas, sep = "_", c("cu", "Rmeas")) |> # This step takes curiously long
  left_join(min_yrs_huc) |> 
  mutate(year = as.numeric(year) - 1 + min_yr) |> 
  select(-min_yr) |> 
  left_join(
    hucuktlis_sr |> 
      filter(parameter != "N1") |> 
      distinct(enhanced, year)
  )


# Observed values for spawners and recruits
obs_sr_huc <- hucuktlis_sr |> 
  mutate(
    S_50 = adult_S,
    S_10 = adult_S - adult_S*S_cv*1.28,
    S_90 = adult_S + adult_S*S_cv*1.28
  ) |> 
  filter(parameter %in% unique(posterior_df_huc$Rmeas)) |> 
  select(cu, Rmeas = parameter, year, matches("(S|R)_\\d{2}")) |> 
  pivot_longer(
    matches("(S|R)_\\d{2}"),
    names_sep = "_",
    names_to = c("parameter", "quantile")
  ) |> 
  pivot_wider(names_from = quantile) |> 
  mutate(set = "obs")


# Plot posterior estimates for spawners and recruits
posterior_df_huc |> 
  filter(parameter %in% c("S_true", "R_true")) |> 
  pivot_wider(
    names_from = parameter,
    values_from = value
  ) |> 
  rename_with(\(x) str_remove(x, "_true")) |> 
  summarize(
    .by = c(cu, Rmeas, year),
    across(
      c(S, R),
      .fns = list(
        "50" = median,
        "10" = ~quantile(.x, 0.1, na.rm = TRUE),
        "90" = ~quantile(.x, 0.9, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    )
  ) |> 
  pivot_longer(
    matches("(S|R)_\\d{2}"),
    names_sep = "_",
    names_to = c("parameter", "quantile")
  ) |> 
  pivot_wider(names_from = quantile) |> 
  mutate(set = "post") |> 
  bind_rows(obs_sr_huc) %>% 
  split(.$Rmeas) |> 
  imap(
    \(x, idx) x |> 
      mutate(parameter = if_else(parameter == "R", idx, parameter)) |> 
      ggplot(aes(year, `50`, colour = set)) +
      facet_grid(
        parameter ~ cu, 
        switch = "y",
        scales = "free_y"
      ) +
      geom_pointrange(
        aes(
          ymin = `10`,
          ymax = `90`
        ),
        alpha = 0.5
      ) +
      theme(
        strip.background.y = element_blank(),
        axis.title.y = element_blank(),
        strip.placement = "outside"
      )
  )


# Plot predicted Beverton-Holt curve versus observed data
# Extract posterior samples of alpha and beta
hucuktlis_a_b_draws <- posterior_df_huc |> 
  filter(str_detect(parameter, "alpha|beta")) |> 
  pivot_wider(names_from = parameter) |> 
  select(-enhanced)


# Create CU-specific range of spawners to predict across
S_pred_huc <- hucuktlis_sr |> 
  summarize(
    .by = c(cu, enhanced),
    min_S = 0,
    max_S = max(adult_S)
  ) |> 
  rowwise() |> 
  mutate(S = list(seq(min_S, max_S, length.out = 100))) |> 
  unnest(S) |> 
  select(cu, S, enhanced)


# Join range of predicted spawners to alpha and beta draws
hucuktlis_pred_frame <- hucuktlis_a_b_draws |> 
  left_join(
    S_pred_huc,
    by = "cu",
    relationship = "many-to-many"
  )

# Check that the resulting dataframe has the correct number of rows
stopifnot(
  all.equal(
    nrow(hucuktlis_pred_frame), 
    nrow(S_pred_huc)/length(unique(S_pred_huc$cu))*nrow(hucuktlis_a_b_draws)
  )
)


# Latent states of spawners and recruits
latent_sr_huc <- posterior_df_huc |> 
  filter(
    parameter %in% c("S_true", "R_true"),
    !is.na(value)
  ) |> 
  mutate(value = value / 1000) |> 
  summarize(
    .by = c(cu, Rmeas, year, parameter),
    value = list(quantile(value, c(0.10, 0.50, 0.90)))
  ) |> 
  unnest_wider(value) |> 
  pivot_wider(
    names_from = parameter,
    values_from = matches("%"),
    names_glue = "{parameter}_{str_remove(.value, '%')}"
  ) |> 
  left_join(
    hucuktlis_sr |> 
      filter(parameter != "N1") |> 
      distinct(cu, year, enhanced)
  )


# Generate predictions for each spawner value
(hucuktlis_bevholt_fits <- hucuktlis_pred_frame |> 
  mutate(
    R_pred = if_else(
      enhanced == 1,
      ((alpha_enh * S) / (1 + (alpha_enh/beta_enh)*S))/1000,
      ((alpha_noenh * S) / (1 + (alpha_noenh/beta_noenh)*S))/1000
    ),
    S = S/1000
  ) |> 
  summarize(
    .by = c(cu, Rmeas, S, enhanced),
    R_quant = list(quantile(R_pred, c(0.025, 0.10, 0.50, 0.90, 0.975)))
  ) |> 
  unnest_wider(R_quant) |>
  rename_with(
    \(x) str_remove(paste0("R_pred_", x), "%"),
    .cols = matches("%")
  ) |> 
  ggplot(
    aes(
      x = S+1e-6, 
      y = R_pred_50,
      colour = enhanced,
      fill = enhanced
    )
  ) +
  facet_grid(
    Rmeas ~ cu, 
    scales = "free",
    switch = "y",
    labeller = labeller(
      Rmeas = c(
        "BYB" = "Smolt biomass (kg)",
        "BYO" = "Smolt abundance (1000s)"
      )
    )
  ) +
  geom_line() +
  geom_ribbon(
    aes(
      ymin = R_pred_10,
      ymax = R_pred_90
    ),
    colour = NA,
    alpha = 0.25
  ) +
  geom_linerange(
    data = latent_sr_huc,
    aes(
      y = R_true_50,
      ymin = R_true_10,
      ymax = R_true_90,
      x = S_true_50
    ),
    colour = "grey50"
  ) +
  geom_pointrange(
    data = latent_sr_huc,
    aes(
      y = R_true_50,
      x = S_true_50,
      xmin = S_true_10+1e-6,
      xmax = S_true_90,
      fill = enhanced
    ),
    shape = 21,
    colour = "grey50"
  ) +  
  scale_colour_manual(
    name = "Enhancement\nstate",
    values = c("blue", "red"),
    aesthetics = c("colour", "fill")
  ) +
  scale_x_continuous(
    name = "Adult spawners (1000s)",
    limits = c(1, NA),
    trans = "log10",
    expand = expansion(mult = c(0, 0.05)),
    labels = scales::label_number()
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05)),
    labels = scales::label_number()
  ) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    axis.title.y = element_blank(),
    strip.placement = "outside",
    panel.spacing.y = unit(1, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
)
  

# Make combined plots of B-H fits and recruitment residuals ---------------



# Start with recruitment residuals (easier)
(R_resid_p <- posterior_df |> 
  left_join(
    somass_sr |> 
      filter(parameter != "N1") |> 
      distinct(cu, year, enhanced)
  ) |> 
  bind_rows(posterior_df_huc) |> 
  filter(
    .by = cu,
    parameter == "lnresid",
    !is.na(value),
     year > min(year, na.rm = TRUE)
  ) |> 
  summarize(
    .by = c(cu, Rmeas, year, enhanced),
    resid = list(quantile(value, c(0.025, 0.10, 0.50, 0.90, 0.975)))
  ) |> 
  unnest_wider(resid) |> 
  filter(year > min(year)) |> 
  mutate(
    .by = Rmeas,
    line_y = max(`97.5%`) * 1.05,
    cu = factor(cu, levels = c("Great Central", "Sproat", "Hucuktlis"))
    ) |> 
  ggplot(aes(x = year, y = `50%`)) +
  facet_grid(
    Rmeas ~ cu, 
    scales = "free",
    switch = "y",
    labeller = labeller(
      Rmeas = c(
        "BYB" = "Smolt biomass per spawner residuals",
        "BYO" = "Smolt abundance per spawner residuals"
      )
    )
  ) +
  geom_hline(
    yintercept = 0, 
    lty = 2, 
    colour = "grey50"
  ) +
  geom_point(
    aes(
      y = line_y,
      colour = enhanced
    ),
    shape = "_",
    size = 2.5
  ) +
  geom_ribbon(
    aes(
      ymin = `2.5%`,
      ymax = `97.5%`
    ),
    alpha = 0.2
  ) +
  geom_ribbon(
    aes(
      ymin = `10%`,
      ymax = `90%`
    ),
    alpha = 0.3
  ) +
  geom_line() +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
  scale_colour_manual(values = c("blue", "red")) +
  guides(
    colour = guide_legend(
      title = "Enhancement state",
      theme = theme(
        legend.direction = "horizontal",
        legend.title.position = "top"
      )
    )
  ) +
  theme(
    strip.background = element_blank(),
    axis.title.y = element_blank(),
    strip.placement = "outside",
    panel.spacing.y = unit(1, "lines"),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "inside",
    legend.position.inside = c(0.03, 0.53),
    legend.justification.inside = c(0, 0),
    legend.background = element_rect(colour = "black", fill = alpha("white", 0.75))
  )
)
# Does it make sense to include the last year in the time series?


# Save the plot
ggsave(
  R_resid_p,
  filename = here(
    "3. outputs",
    "Plots",
    "Smolt-per-spawner_residuals_timeseries_all-CUs.png"
  ),
  width = 6.5, 
  height = 5.5,
  units = "in",
  dpi = "print"
)


# Bring B-H curve fits together in a grid plot
all_bevholt_fits <- plot_grid(
  somass_bevholt_fits +
    guides(fill = "none") +
    theme(axis.title.x = element_text(hjust = 1)),
  hucuktlis_bevholt_fits +
    theme(
      strip.text.y = element_blank(),
      axis.title.x = element_text(colour = alpha("white", 0))
      ),
  rel_widths = c(1.4, 1)
)


# Save the fitted Beverton-holt curves
ggsave(
  all_bevholt_fits,
  filename = here(
    "3. outputs",
    "Plots",
    "Smolt_production_Beverton-Holt_all-CUs.png"
  ),
  width = 11, 
  height = 7,
  units = "in",
  dpi = "print"
)


# Miscellaneous summaries and exports ----------------------


# Export posterior Beverton-Holt parameter samples for Wendell Challenger
if(
  #FALSE
  TRUE
) {
  list(posterior_df, posterior_df_huc) |> 
    map(
      \(x) filter(
        .data = x,
        case_when(
          parameter %in% c("v", "z_v", "phi", "k_R", "tau_R", "sigma_R_proc") ~ TRUE,
          str_detect(parameter, "alpha|beta") ~ TRUE,
          .default = FALSE
        )
      )
    ) |> 
    list_rbind() |> 
    write.csv(
      file = here(
        "3. outputs",
        "Stock-recruit modelling",
        "Spawner-smolt posterior",
        "Barkley_Sockeye_spawner-smolt_BevHolt-params_posterior.csv"
      ),
      row.names = FALSE
    )
}


bind_rows(
  posterior_df,
  posterior_df_huc
) |> 
  filter(str_detect(parameter, "alpha|beta")) |> 
  mutate(
    value = case_when(
      str_detect(parameter, "beta") & Rmeas == "BYO" ~ value/1e6,
      Rmeas == "BYB" ~ value/1e3,
      .default = value
    )
  ) |> 
  summarize(
    .by = c(cu, Rmeas, enhanced, parameter),
    median = median(value),
    lci = quantile(value, 0.1),
    uci = quantile(value, 0.9)
  )


# Export beta prior information for stock-recruit Ricker models
bind_rows(
  posterior_df_huc,
  posterior_df
) |> 
  filter(str_detect(parameter, "alpha|beta")) |> 
  mutate(
    enhanced = case_when(
      str_detect(parameter, "_enh") ~ 1,
      str_detect(parameter, "_noenh") ~ 0,
      .default = NA
    )
  ) |> 
  pivot_wider(names_from = parameter) |> 
  mutate(
    ricker_beta_prior = case_when(
      is.na(enhanced) ~ alpha/beta,
      enhanced == 0 ~ alpha_noenh/beta_noenh,
      enhanced == 1 ~ alpha_enh/beta_enh
    )
  ) |> 
  summarize(
    .by = c(cu, Rmeas, enhanced),
    mean_ln_b = mean(log(ricker_beta_prior)),
    sd_ln_b = sd(log(ricker_beta_prior))
  ) |> 
  saveRDS(
    file = here(
      "3. outputs",
      "Stock-recruit modelling",
      "Ricker_beta_priors_smoltBev-Holt.RDS"
    )
  )

