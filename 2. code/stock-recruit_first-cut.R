# Packages ----------------------------------------------------------------

pkgs <- c(
  "here", "tidyverse", "readxl", "writexl", "ggpmisc", "broom",
  "janitor" 
)
#install.packages(pkgs)

library(here)
library(tidyverse); theme_set(theme_bw(base_size = 18))
library(ggpmisc)
library(readxl)
library(writexl)
library(magrittr)


# Load and clean input data -----------------------------------------------


# Returning abundance by age time series
sk_ret0 <- read_xlsx(
  here("1. data", "return by age time series.xlsx")
) |> 
  mutate(
    brood_year = return_year - ttl_age,
    smolt_year = brood_year + fw_age,
    gr_age = paste0(ttl_age,"[", fw_age, "]"),
    return = escapement + catch
  )


# Create lookup table for spawners by brood year
by_spwns <- sk_ret0 |> 
  summarize(
    .by = c(stock, return_year),
    by_spawners = sum(escapement),
    by_adult_spawners = sum(escapement[!gr_age %in% c("3[2]", "4[3]")])
  ) |> 
  rename("brood_year" = return_year)


# Create a lookup table for sibling abundances in previous return year
sib_lu <- tribble(
  # Specify sibling ages used as predictors
  ~gr_age, ~sib_age,
  "4[2]", "3[2]",
  "5[2]", "3[2]",
  "5[2]", "4[2]",
  "5[3]", "4[3]",
  "6[3]", "4[3]",
  "6[3]", "5[3]"
) |> 
  # Add ranges of brood years and stocks observed per GR age
  left_join(
    distinct(sk_ret0, stock, brood_year, gr_age),
    relationship = "many-to-many"
  ) |> 
  # Add return data, matching age to identified sibling ages
  left_join(
    distinct(sk_ret0, stock, brood_year, gr_age, return),
    by = join_by(
      stock == stock,
      brood_year == brood_year,
      sib_age == gr_age # Pair GR age in raw data to sibling age in LU table
    )
  ) |> 
  # Sum sibling abundances per GR age group, stock, and brood year
  summarize(
    .by = c(stock, gr_age, brood_year),
    sib_return = sum(return)
  )


# Add brood year spawner abundances and sibling abundances to main data
sk_ret <- sk_ret0 |> 
  left_join(by_spwns) |> 
  left_join(sib_lu)



# Calculate expected smolt age proportions from brood year return data
by_fw_age_props <- sk_ret |> 
  summarize(
    .by = c(brood_year, fw_age, stock),
    return = sum(return)
  ) |> 
  mutate(
    .by = c(stock, brood_year),
    prop = return / sum(return),
    .keep = "unused"
  )


# Average recent 4-year age 2 versus 3 compositions
fw_age_4y_avg <- by_fw_age_props |> 
  filter(!is.na(prop)) |> 
  filter(brood_year > max(brood_year) - 5) |> 
  summarize(
    .by = c(stock, fw_age),
    prop_4ya = mean(prop)
  )


# Smolt time series data 
smolts <- read_xlsx(
  here("1. data", "lake smolts time series.xlsx"),
  sheet = "abundance"
) |> 
  pivot_longer(
    cols = !smolt_year,
    names_sep = "_",
    names_to = c("stock", "age"),
    names_transform = str_to_upper
  ) |> 
  filter(!stock == "HED") |> 
  pivot_wider(
    names_from = age,
    values_from = value,
    names_prefix = "age_"
  ) |> 
  pivot_longer(
    cols = c(age_1, age_2),
    names_prefix = "age_",
    names_to = "age",
    values_to = "smolts",
    names_transform = as.numeric
  ) |> 
  rename("ttl_smolts" = age_TTL) |> 
  mutate(
    fw_age = age + 1,
    brood_year = smolt_year - fw_age
  ) |> 
  left_join(by_fw_age_props) |> 
  left_join(fw_age_4y_avg) |> # 4 year averages
  # Infill missing age data using historical average proportions
  mutate(
    infill = if_else(is.na(smolts) & !is.na(ttl_smolts), "infill", "obs"),
    smolts = case_when(
      is.na(smolts) & !is.na(prop) ~ ttl_smolts*prop, 
      # Use recent 4-year average for incomplete recent brood years
      is.na(smolts) & is.na(prop) &!is.na(ttl_smolts) ~ ttl_smolts*prop_4ya,
      TRUE ~ smolts
    ),
    .keep = "unused"
  ) |> 
  mutate(
    .by = brood_year,
    infill = if_else(any(str_detect(infill, "infill")), "infill", "obs")
  )




# Prepare data for stock-recruit modelling --------------------------------


# Build data frame with data categories to be summarized for modelling
fw_groups <- sk_ret |> 
  summarize(
    .by = c(stock, fw_age),
    ttl_age = list(unique(ttl_age))
  ) |> 
  mutate(
    across(contains("age"), as.list),
    age_group = "fw_age"
  )
  
ttl_groups <- sk_ret |> 
  summarize(
    .by = c(stock, ttl_age),
    fw_age = list(unique(fw_age))
  ) |> 
  mutate(
    across(contains("age"), as.list),
    age_group = "ttl_age"
  )

all_groups <- sk_ret |> 
  distinct(stock, ttl_age, fw_age) |> 
  mutate(
    across(contains("age"), as.list),
    age_group = "gr_age"
  )

no_groups <- sk_ret |> 
  summarize(
    .by = stock,
    ttl_age = list(unique(ttl_age)),
    fw_age = list(unique(fw_age))
  ) |> 
  mutate(
    across(contains("age"), as.list),
    age_group = "none"
  )


data_groups <- bind_rows(fw_groups, ttl_groups, all_groups, no_groups) |> 
  # Rename columns to disambiguate during data filtering
  rename(
    "total_age" = ttl_age,
    "freshwater_age" = fw_age,
    "cu" = stock
  )


# Nested dataframe with stock-recruit data for each group
model_data <- data_groups |> 
  rowwise() |> 
  mutate(
    sr_data = list(
      sk_ret |> 
        filter(
          ttl_age %in% total_age,
          fw_age %in% freshwater_age,
          stock == cu
        ) |> 
        rename("spawners" = by_spawners) |> 
        summarize(
          .by = c(stock, brood_year, spawners),
          return = sum(return)
        ) |> 
        # Remove all rows where return or spawners is NA or <= 0
        filter(
          !if_any(c(return, spawners), is.na),
          !if_any(c(return, spawners), \(x) x <= 0)
        ) |> 
        mutate(r_s = return/spawners)
    )
  ) |> 
  ungroup()



# Model fitting -----------------------------------------------------------


# Run linear, Ricker, Beverton-Holt, Power relationships on each dataset
sr_models <- model_data |> 
  rowwise() |> 
  # Fit models
  mutate(
    linear_model = list(lm(return ~ spawners, data = sr_data)),
    ricker_model = list(lm(log(r_s) ~ spawners, data = sr_data)),
    bevholt_model = list(
      glm(
        return ~ I(1/spawners), 
        family = gaussian(link = "inverse"), 
        data = sr_data
      )
    )
  ) |> 
  # Generate predictions to compare against observed values
  mutate(
    newdata = list(
      tibble(
        spawners = seq(
          min(sr_data$spawners, na.rm = TRUE), 
          max(sr_data$spawners, na.rm = TRUE), 
          length.out = 100
        )
      )
    ),
    linear_pred = list(
      predict(
        linear_model,
        newdata, 
        interval = "confidence"
      ) |> 
        as_tibble() |> 
        bind_cols(newdata)
    ),
    # Bootstrap predictions from Ricker model
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
          pred_rec = a*spawners*exp(-b*spawners)
        ) |> 
        summarize(
          .by = spawners,
          fit = mean(pred_rec),
          lwr = quantile(pred_rec, 0.025),
          upr = quantile(pred_rec, 0.975)
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
          lwr = fit - 1.96*se.fit,
          upr = fit + 1.96*se.fit,
          # Back-transform estimates to response scale
          across(c(fit, lwr, upr), gaussian(link = "inverse")$linkinv)
        )
    )
  ) |> 
  pivot_longer(
    cols = matches("_(model|pred)$"),
    names_sep = "_",
    names_to = c("name", ".value")
  )



# Plot stock-recruit data
(sr_plots <- list(
  # find a way to make this code more elegant
  "observed" = select(sr_models, cu, contains("age"), name, sr_data) |> unnest(sr_data),
  "pred" =  select(sr_models, cu, contains("age"), name, pred) |> unnest(pred)
) |> 
    # Doesn't quite work, yet...
    pmap(
      ~ggplot(..1, aes(x = spawners, y = return)) +
        facet_grid(cu+freshwater_age+total_age~name) +
        geom_point() +
        geom_line(
          data = ..2, 
          aes(y = fit),
          colour = "blue"
        ) +
        geom_ribbon(
          data = ..2,
          aes(y = fit, ymin = lwr, ymax = upr),
          fill = "blue",
          alpha = 0.3
        ) +
        scale_x_continuous(labels = scales::comma) +
        scale_y_continuous(
          labels = scales::comma,
          breaks = seq(0, 1.5e6, by = 5e5),
          expand = expansion(mult = c(0, 0.05))
        )
    )
)


# Save stock-recruit plots
sr_plots |> 
  rlang::set_names(
    nm = ricker_data |> 
      select(population, col) |> 
      unite(col = "name") |> 
      pull()
  ) |> 
  iwalk(
    ~ggsave(
      plot = .x,
      filename = here(
        "03-output figures",
        paste0("R-PLOT_Ricker_pred", .y, ".png")
      ),
      height = 4,
      width = 7,
      units = "in"
    )
  )




# What about age-specific Ricker relationships? ---------------------------


# Look-up table of brood-year spawners by age
meta_by_spwns <- sk_ret0 |> 
  filter(!gr_age %in% c("3[2]", "4[3]")) |> 
  summarize(
    .by = c(stock, return_year, gr_age),
    by_meta_spawners = sum(escapement)
  ) |> 
  rename("brood_year" = return_year)


# Plot basic stock-recruit relationship
meta_ricker_data <- meta_by_spwns |> 
  left_join(sk_ret) |> 
  filter(!by_meta_spawners <= 0) |> 
  select(stock, gr_age, brood_year, by_meta_spawners, return) |> 
  mutate(r_s = return / by_meta_spawners) |> 
  nest(.by = c(stock, gr_age)) |> 
  rowwise() |> 
  mutate(
    model = list(lm(log(r_s) ~ by_meta_spawners, data = data)),
    r2 = summary(model)$r.squared,
    vcov_mat = list(as.matrix(vcov(model))),
    bootstrap_pred = list(
      MASS::mvrnorm(10000, mu = coef(model), Sigma = vcov_mat) |> 
        as_tibble() |> 
        rename("a" = 1, "b" = 2) |> 
        expand_grid(
          by_meta_spawners = seq(
            min(data$by_meta_spawners, na.rm = TRUE), 
            max(data$by_meta_spawners, na.rm = TRUE), 
            length.out = 100
          )
        ) |> 
        mutate(
          a = exp(a),
          b = -b,
          pred_rec = a*by_meta_spawners*exp(-b*by_meta_spawners)
        ) |> 
        summarize(
          .by = by_meta_spawners,
          mean_rec = mean(pred_rec),
          lci_rec = quantile(pred_rec, 0.025),
          uci_rec = quantile(pred_rec, 0.975)
        )
    )
  ) |> 
  ungroup() 


# Plot meta-population stock-recruit data
(meta_plots <- meta_ricker_data |> 
    select(1:3, r2, bootstrap_pred) |> 
    rowwise() |> 
    mutate(
      r2_x = max(data$by_meta_spawners, na.rm = TRUE)*0.9,
      r2_y = max(data$return, na.rm = TRUE)*0.95
    ) |> 
    ungroup() |> 
    as.list() |> 
    pmap(
      ~ggplot(..3, aes(x = by_meta_spawners, y = return)) + 
        geom_point() +
        geom_line(
          data = ..5, 
          aes(y = mean_rec),
          colour = "blue"
        ) +
        geom_ribbon(
          data = ..5,
          aes(y = mean_rec, ymin = lci_rec, ymax = uci_rec),
          fill = "blue",
          alpha = 0.3
        ) +
        annotate(
          "text", 
          label = paste0("R^2~`=`~", round(..4, 3)),
          parse = TRUE,
          x = ..6, 
          y = ..7,
          colour = "red",
          size = 8,
          hjust = 1
        ) +
        labs(
          x = parse(text = paste0("Brood~year~age~", ..2, "~spawners")),
          y = parse(text = paste0("Age~", ..2, "~return"))
        ) +
        scale_x_continuous(labels = scales::comma) +
        scale_y_continuous(labels = scales::comma) +
        ggtitle(..1)
    )
)


# Save the plots 
meta_plots |> 
  rlang::set_names(
    nm = meta_ricker_data |> 
      select(stock, gr_age) |> 
      unite(col = "name") |> 
      pull()
  ) |> 
  iwalk(
    ~ggsave(
      plot = .x,
      filename = here(
        "03-output figures",
        paste0("R-PLOT_metapop_Ricker_pred", .y, ".png")
      ),
      height = 5,
      width = 5,
      units = "in"
    )
  )




# Smolts per spawner ------------------------------------------------------


# Save stock-recruit data and Ricker models
smolt_recruit <- smolts |> 
  summarize(
    .by = c(brood_year, stock, infill),
    smolt_n = sum(smolts)*1e6
  ) |> 
  left_join(by_spwns)


# Plot smolts versus spawners
(smolt_rs_plots <- smolt_recruit |> 
    pivot_longer(
      contains("by_"),
      names_to = "spawner_units",
      names_prefix = "by_",
      values_to = "spawners"
    ) |> 
    mutate(r_s = smolt_n/spawners) |> 
    expand_grid(transform = c("log", "identity")) %>%
    split(.$transform, drop = TRUE) |> 
    imap(
      ~ .x |> 
        ggplot(aes(x = spawners, y = r_s)) +
        facet_grid(stock ~ spawner_units) +
        geom_point(aes(shape = infill)) +
        geom_smooth(method = "lm") +
        stat_poly_eq(label.x = "right", colour = "red") +
        scale_shape_manual(values = c(21, 19)) +
        scale_y_continuous(
          labels = scales::comma,
          transform = .y
        ) +
        scale_x_continuous(labels = scales::comma) +
        guides(shape = "none") +
        labs(
          y = "Smolts-per-spawner",
          x = "Brood year spawner abundance"
        ) +
        ggtitle(paste0("Y-axis transform: ", .y)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    )
)
# Looks like some negative density dependence
# Spawners generally a better predicitor than adult spawners
# Infilled data tend to be outliers...


# Save Smolt R/S plot
smolt_rs_plots$log |> 
  ggsave(
    filename = here(
      "03-output figures",
      "R-PLOT_smolts-per-spawner.png"
    ),
    height = 5,
    width = 7,
    units = "in"
  )


# Ricker model on smolts as recruits
ricker_smolt <- smolt_recruit |> 
  mutate(r_s = smolt_n/by_spawners) |> 
  nest(.by = stock) |> 
  rowwise() |> 
  mutate(
    model = list(lm(log(r_s) ~ by_spawners, data = data)),
    r2 = summary(model)$r.squared,
    vcov_mat = list(as.matrix(vcov(model))),
    bootstrap_pred = list(
      MASS::mvrnorm(10000, mu = coef(model), Sigma = vcov_mat) |> 
        as_tibble() |> 
        rename("a" = 1, "b" = 2) |> 
        expand_grid(
          by_spawners = seq(
            min(data$by_spawners, na.rm = TRUE), 
            max(data$by_spawners, na.rm = TRUE), 
            length.out = 100
          )
        ) |> 
        mutate(
          a = exp(a),
          b = -b,
          pred_smolts = a*by_spawners*exp(-b*by_spawners)
        ) |> 
        summarize(
          .by = by_spawners,
          mean_smolts = mean(pred_smolts),
          lci_smolts = quantile(pred_smolts, 0.025),
          uci_smolts = quantile(pred_smolts, 0.975)
        )
    )
  ) |> 
  ungroup() 


# Plot stock-recruit data
(ricker_smolt_plots <- ricker_smolt |> 
    select(stock, data, r2, bootstrap_pred) |> 
    as.list() |> 
    pmap(
      ~ggplot(..2, aes(x = by_spawners, y = smolt_n)) + 
        geom_point(aes(shape = infill)) +
        geom_line(
          data = ..4, 
          aes(y = mean_smolts),
          colour = "blue"
        ) +
        geom_ribbon(
          data = ..4,
          aes(y = mean_smolts, ymin = lci_smolts, ymax = uci_smolts),
          fill = "blue",
          alpha = 0.3
        ) +
        annotate(
          "text", 
          label = paste0("R^2~`=`~", round(..3, 3)),
          parse = TRUE,
          x = 4e5, 
          y = 15e6,
          colour = "red"
        ) +
        scale_shape_manual(values = c(21, 19)) +
        scale_x_continuous(labels = scales::comma) +
        scale_y_continuous(
          labels = scales::comma,
          #breaks = seq(0, 1.5e6, by = 5e5),
          expand = expansion(mult = c(0, 0.05))
        ) +
        guides(shape = "none") +
        labs(
          y = "Pre-smolt abundance",
          x = "Brood year spawner abundance"
        ) +
        ggtitle(..1)
    )
)


# Save the plots 
ricker_smolt_plots |> 
  rlang::set_names(
    nm = ricker_smolt |> 
      select(stock) |> 
      unite(col = "name") |> 
      pull()
  ) |> 
  iwalk(
    ~ggsave(
      plot = .x,
      filename = here(
        "03-output figures",
        paste0("R-PLOT_smolts_Ricker_pred", .y, ".png")
      ),
      height = 5,
      width = 7,
      units = "in"
    )
  )


# Predictors for Ricker residuals -----------------------------------------


# Calculate Ricker residuals
meta_ricker_resid <- meta_ricker_data |> 
  rowwise() |> 
  mutate(
    resid_data = list(
      select(.data = data, brood_year, by_meta_spawners, return) |> 
        mutate(
          a = exp(coef(model)[1]),
          b = -coef(model)[2],
          pred_rec = a*by_meta_spawners*exp(-b*by_meta_spawners),
          ricker_resid = return - pred_rec
        )
    )
  ) |> 
  select(stock, gr_age, resid_data) |> 
  unnest(resid_data)


# Plot predictions and residuals versus observed data 
meta_ricker_resid |> 
  ggplot(aes(by_meta_spawners, return)) +
  facet_grid(
    stock ~ gr_age,
    scales = "free",
    labeller = label_parsed
  ) +
  geom_line(
    aes(y = pred_rec),
    colour = "blue"
  ) +
  geom_segment(
    aes(
      y = pred_rec, 
      yend = pred_rec + ricker_resid
    ),
    colour = "red"
  ) +
  geom_point()
# Looks correct


# Add sibling return data to Ricker meta-population residuals
ricker_resid_pred <- meta_ricker_resid |> 
  left_join(select(sk_ret, stock, gr_age, brood_year, smolt_year, sib_return)) |> 
  mutate(
    # Looks like cube-root transformation might be appropriate
    # via: https://stats.stackexchange.com/questions/155429/how-to-transform-negative-values-to-logarithms
    ricker_resid_cuberoot = sign(ricker_resid) * (abs(ricker_resid))^(1/3),
    # Try just normal log transformation after adjusting values based on minimum 
    ricker_resid_logshift = log(ricker_resid + abs(min(ricker_resid, na.rm = TRUE)))
  ) |> 
  # Add oceanographic variables
  left_join(distinct(survival_predictors, pick(matches("sal|temp")), smolt_year))


# Compare histograms of original versus transformed Ricker residual values
ricker_resid_pred |> 
  pivot_longer(matches("ricker_resid.*")) |> 
  ggplot() +
  facet_grid(
    stock ~ name,
    switch = "y",
    scales = "free"
  ) +
  geom_histogram(aes(x = value))


# How do covariates correlate to Ricker residuals?
ricker_resid_pred |> 
  pivot_longer(
    matches("ricker_resid.*"),
    names_prefix = "ricker_resid_",
    names_to = "transformation",
    values_to = "ricker_resid"
  ) |> 
  pivot_longer(
    c(sib_return, matches("(?i)temp|sal")),
    names_to = "covariate",
    values_to = "value"
  ) |> 
  nest(.by = c(transformation, covariate)) |> 
  mutate(covariate_trans = if_else(covariate == "sib_return", "log", "identity")) |> 
  as.list() |> 
  pmap(
    ~ ..3 |> 
      ggplot(aes(x = value, y = ricker_resid)) +
      facet_grid(stock ~ gr_age, scales = "free") +
      geom_point() +
      geom_smooth(method = "lm") +
      stat_poly_eq(colour = "red") +
      scale_x_continuous(transform = ..4) +
      labs(
        x = ..2, 
        y = ..1
      )
  )


# Correlation matrix
ricker_resid_pred |> 
  mutate(across(contains("return"), ~log(.x), .names = "{.col}_log")) |> 
  pivot_longer(
    matches("ricker_resid.*|^return"),
    names_to = "outcome",
    values_to = "y"
  ) |> 
  pivot_longer(
    c(matches("(?i)temp|sal|sib_return")),
    names_to = "covariate",
    values_to = "x"
  ) |> 
  filter(!if_any(c(y, x), is.infinite)) |> 
  nest(.by = c(stock, gr_age, outcome, covariate)) |> 
  rowwise() |> 
  mutate(
    lm = list(lm(y ~ x, data = data)),
    r2 = sign(coef(lm)[2])*summary(lm)$r.squared
  ) |> 
  ggplot(aes(x = covariate, y = outcome)) +
  facet_grid(stock ~ gr_age) +
  geom_tile(aes(fill = r2)) +
  scale_fill_distiller(
    palette = "RdBu",
    limits = c(-1, 1)
  ) +
  coord_cartesian(expand = FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# At this point, I conclude that different variables appear important for the
# different ages. Sibling returns are generally strong predictors of both returns
# and Ricker residuals of returns. South Brooks temperatures appear to be 
# somewhat important for GCL 4(2)s in particular



