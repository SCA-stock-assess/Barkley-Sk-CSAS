# Packages ----------------------------------------------------------------


pkgs <- c("here", "readxl", "tidyverse", "broom")
#install.packages(pkgs)


library(here)
library(tidyverse); theme_set(theme_bw())
library(readxl)
library(broom)



# Load pre-smolt abundance and size  and adult spawner data --------------


# Lake abundance estimates
fry_abun <- here(
  "1. data",
  "Barkley-Sk_pre-smolt_abundances.xlsx"
) |> 
  read_xlsx() |> 
  pivot_longer(
    matches("gcl|spr|huc"),
    names_sep = "_",
    names_to = c("cu", ".value")
  ) |> 
  # For the purposes of this analysis, we will focus only on the total #s
  select(smolt_year, cu, ttl)


# Pre-smolt size data
fry_sizes <- here(
  "1. data",
  "smolt size data.xlsx"
) |> 
  read_xlsx(sheet = "morphometrics") |> 
  pivot_longer(
    matches("gcl|spr|huc"),
    names_sep = "_",
    names_to = c("cu", ".value")
  ) |> 
  select(-n)


# Join abundance and size data
fry <- left_join(
  fry_abun,
  fry_sizes
) |> 
  mutate(biomass = w*ttl)


# Spawner abundance data
spwn <-  here(
  "3. outputs",
  "Stock-recruit data",
  "Barkley_Sockeye_stock-recruit_infilled.xlsx"
) |> 
  read_xlsx(sheet = "S-R data") |> 
  rename("cu" = stock) |> 
  mutate(
    cu = tolower(cu),
    smolt_year = year + 2,
    # smolt_year = year for hatchery_fry_releases--lead back to align
    #hatchery_fry_release = lead(hatchery_fry_release, 2),
    hatchery_fry_release = if_else(
      is.na(hatchery_fry_release), 
      0, 
      hatchery_fry_release
    )
  ) |> 
  select(smolt_year, cu, S, adult_S, fertilized, hatchery_fry_release)


# Join spawner data to fry data
spwn_fry <- left_join(
  fry,
  spwn
) |>
  mutate(
    fertilized = factor(fertilized),
    cu = factor(
      cu,
      levels = c("gcl", "spr", "huc"),
      labels = c("Great Central", "Sproat", "Hucuktlis")
    ),
    # Remove hatchery contributions from Hucuktlis data
    ttl = if_else(
      # In 1999, assume 90% hatchery smolts. Otherwise, subtracting 10% of 
      # hatchery released fry from from the total lake ATS estimate (as is
      # done for all other years) results in a - value. This is following
      # a footnote in Tablw 13 of the unpublished draft Henderson Stock
      # Status PSARC report (from 2008). 
      cu == "Hucuktlis" & smolt_year == 1999,
      ttl*0.1, 
      ttl - (hatchery_fry_release*0.1/1e6)
      ),
    biomass = ttl * w
  )



# Plot relationships between spawners and pre-smolts ---------------------


# Plots to examine should include `spawners` versus:
# 1) pre-smolt abundance
# 2) pre-smolt biomass
# 3) log(pre-smolt abundance per spawner)
# 4) log(pre-smolt biomass per spawner)
# 
# Also need to examine fits with & without Hucuktlis 1993 outlier removed 


# Set up a nested dataframe containing those 4 outcome variables
nested_data <- tibble(
  response = c("abun", "biomass", "log_abun_spwn", "log_biom_spwn")
) |> 
  expand_grid(fltr = c(0, 1)) |> 
  mutate(data = list(spwn_fry)) |> 
  rowwise() |> 
  mutate(
    data = case_when(
      response == "abun" ~ list(mutate(data, y = ttl)),
      response == "biomass" ~ list(mutate(data, y = biomass)),
      response == "log_abun_spwn" ~ list(mutate(data, y = log(ttl/S))),
      response == "log_biom_spwn" ~ list(mutate(data, y = log(biomass/S)))
    ),
    data = list(select(data, smolt_year, cu, S, adult_S, y, fertilized)),
    # Filter the Hucuktlis data to remove 1993 outlier
    data = if_else(
      fltr == 1,
      list(filter(data, !(cu == "Hucuktlis" & S > 150000))),
      list(data)
    ),
    response_long = case_when(
      response == "abun" ~ "Smolt abundance",
      response == "biomass" ~ "Total smolt biomass",
      response == "log_abun_spwn" ~ "log(smolt abundance per spawner)",
      response == "log_biom_spwn" ~ "log(smolt biomass per spawner)"
    )
  ) |> 
  ungroup()
  

# Feed the different outcome variables into plots
nested_data |> 
  unnest(data) |>   
  # Filtering is only relevant for the Hucuktlis data (currently)
  filter(!(fltr == 1 & cu != "Hucuktlis")) %>%
  split(.$response_long) |> 
  imap(
    \(x, idx) ggplot(
      data = filter(x, fltr == 0),
      aes(x = S, y = y, colour = fertilized)
    ) +
      facet_wrap(
        ~cu,
        nrow = 1,
        scales = "free_x"
      ) +
      geom_point() +
      geom_smooth(
        data = x,
        aes(lty = factor(fltr)),
        method = "lm",
        alpha = 0.15
      ) +
      scale_x_continuous(
        name = "Spawners",
        labels = scales::label_number()
      ) +
      labs(y = idx) +
      theme(axis.text.x = element_text(angle = 25, hjust = 1))
  )



# Fit simple Ricker models and examine fits -------------------------------


# Use the data that were assembled above for plotting
model_fits <- nested_data |> 
  unnest(data) |> 
  # Filtering is only relevant for the Hucuktlis data (currently)
  filter(
    !(fltr == 1 & cu != "Hucuktlis"),
    response %in% c("abun", "biomass")
  ) |> 
  expand_grid(predictor = c("Spawners", "Adult spawners")) |> 
  nest(.by = c(cu, contains("response"), fltr, predictor)) |> 
  rowwise() |> 
  mutate(
    data = if_else(
      predictor == "Spawners",
      list(mutate(data, x = S)),
      list(mutate(data, x = adult_S))
    ),
    # Fit Ricker and Beverton-Holt models
    ricker_model = list(lm(log(y/x) ~ x, data = data)),
    bevholt_model = list(
      glm(
        y ~ I(1/x), 
        family = gaussian(link = "inverse"), 
        #start = c(0, 3), # Research how to choose appropriate starting values
        data = data
      )
    ),
    # Add predictions from Ricker and Beverton-Holt models
    newdata = list(
      tibble(
        x = seq(
          #0.01,
          min(data$x, na.rm = TRUE),
          max(data$x, na.rm = TRUE), 
          length.out = 100
        )
      )
    ),
    ricker_pred = list(
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
          pred_y = a*x*exp(-b*x)
        ) |> 
        summarize(
          .by = x,
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
  split(.$predictor) |> 
  imap(
    function(set, idx) {
      
      sub_set <- list(
        obs = unnest(set, data),
        pred = unnest(set, pred)
      )
      
      # label <- set |> 
      #   rowwise() |> 
      #   mutate(
      #     x = max(data$x, na.rm = TRUE),
      #     y = max(data$y, na.rm = TRUE)*1.1,
      #     r2 = round(adj.r.squared, 2)
      #   ) |> 
      #   ungroup()
      
      ggplot(
        data = sub_set$obs,
        aes(
          x = x, 
          y = y,
          colour = type,
          fill = type,
          lty = factor(fltr)
        )
      ) +
        facet_grid(
          response_long ~ cu,
          scales = "free",
          switch = "y"
        ) +
        geom_point() +
        # geom_ribbon(
        #   data = sub_set$pred,
        #   aes(y = fit, ymin = lwr, ymax = upr),
        #   alpha = 0.15,
        #   colour = NA
        # ) +
        geom_line(
          data = sub_set$pred,
          aes(y = fit)
        ) +
        # geom_text(
        #   data = label,
        #   aes(label = paste("R sq:", r2)),
        #   hjust = 1
        # ) +
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
