# Packages ----------------------------------------------------------------


pkgs <- c("here", "tidyverse", "readxl", "janitor", "ggpmisc", "magrittr")
#install.packages(pkgs)


library(here)
library(tidyverse); theme_set(theme_bw())
library(ggpmisc)
library(janitor)
library(magrittr) # For "T" pipe operators
library(readxl)



# Load historic adult biodata ---------------------------------------------


adult_bio0 <- here(
  "1. data",
  "SoxBioData1980-2022.xlsx"
) |> 
  read_xlsx(
    na = c("N/A", "NA", "#N/A", "", "0"),
    sheet = "All data",
    col_types = c(
      "numeric", "numeric", "text", "text", "numeric",
      "numeric", "text", "text", "text", "text", "text", "text",
      "text", "text", "date", "date", "text", "numeric", "numeric",
      "numeric", "numeric", "text", "text", "numeric", "numeric", 
      "text", "text"
    )
  ) |> 
  clean_names() |> 
  select(
    year, sample_type, gear, age_structure_type, gilbert_rich_age,
    sample_start_date, site, contains("length"), contains("sex")
  ) |> 
  mutate(
    # Fix a few lengths that are 10x too large
    across(contains("length"), \(x) if_else(x>1000, x/10, x)),
    # Assign fish to CU based on sampling location
    cu = case_when(
      str_detect(site, "(?i)henderson|clemens|uchuck") ~ "Hucuktlis",
      str_detect(site, "(?i)great central") ~ "Great Central",
      str_detect(site, "(?i)sproat") ~ "Sproat",
      .default = NA
    ) |> 
      factor(levels = c("Great Central", "Sproat", "Hucuktlis")),
    id = row_number()
  )



# Develop conversion for lengths to get all lengths in POH ----------


# See how many records we are working with that have known CU
adult_bio0 |> 
  filter(str_detect(sample_type, "(?i)esc")) |> 
  count(
    cu, 
    is.na(poh_length),
    is.na(fork_length),
    is.na(total_length)
  ) |> 
  View()
# In every case, it looks like total can be used to estimate POH


# Does it work to lump all stocks together?
length_comparison <- adult_bio0 |> 
  filter(
    case_when(
      # Remove all records with 1 or 0 length measurements
      is.na(poh_length) & is.na(fork_length) ~ FALSE,
      is.na(poh_length) & is.na(total_length) ~ FALSE,
      is.na(total_length) & is.na(fork_length) ~ FALSE,
      if_all(contains("length"), is.na) ~ FALSE,
      # Remove records where poh_length > fork or total length
      poh_length > fork_length ~ FALSE,
      poh_length > total_length ~ FALSE,
      .default = TRUE
    )
  ) |> 
  # POH length is most commonly recorded, use that as y
  pivot_longer(
    c(fork_length, total_length),
    names_pattern = "(.*)_.*",
    names_to = c("alt_length"),
    values_to = "x"
  ) 

# What to do about fork length? Need to ascertain how often fork length
# was collected but not total length or POH
# Could lump stocks together to create overall forklength to POH relationship?


# Plot with all CUs lumped together
(len_comp_p1 <- length_comparison |> 
    #filter(cu != "unknown") |> 
    ggplot(aes(y = poh_length, x = x)) +
    facet_wrap(
      ~alt_length,
      nrow = 1,
      strip.position = "bottom"
    ) +
    geom_point(
      aes(colour = cu),
      alpha = .1
    ) +
    geom_smooth(method = "lm") +
    stat_poly_eq(
      formula = y ~ x,
      method = "lm",
      use_label(c("eq", "R2"))
    ) +
    labs(x = NULL) +
    theme(
      strip.placement = "outside",
      strip.background = element_blank()
    )
)


# Separate out CUs
(len_comp_p2 <- length_comparison |> 
    #filter(cu != "unknown") |> 
    ggplot(aes(y = poh_length, x = x)) +
    facet_grid(
      cu ~ alt_length,
      switch =  "x"
    ) +
    geom_point(
      aes(colour = cu),
      alpha = .1
    ) +
    geom_smooth(
      aes(colour = cu),
      method = "lm"
    ) +
    stat_poly_eq(
      formula = y ~ x,
      method = "lm",
      use_label(c("eq", "R2"))
    ) +
    labs(x = NULL) +
    theme(
      strip.placement = "outside",
      strip.background = element_blank()
    )
)


# How does Hucuktlis relationship improve when removing the one outlier?
len_comp_p2 + filter(
  length_comparison, 
  !(cu == "Hucuktlis" & x - poh_length > 100)
)
# Looks much better


# Fit CU and length-specific regression models
length_models <- length_comparison |> 
  # Remove Hucuktlis outlier
  filter(!(cu == "Hucuktlis" & x - poh_length > 100)) |>
  # Put total and fork length back into their own columns
  pivot_wider(
    names_from = alt_length,
    values_from = x,
    names_glue = "{alt_length}_length"
  ) |> 
  nest(.by = cu) |> 
  rowwise() |> 
  # Fit linear models on poh_length versus total_length
  mutate(model = list(lm(poh_length ~ total_length, data = data))) |> 
  ungroup()
  

# Use models to infill missing POH lengths
adult_bio1 <- adult_bio0 |> 
  nest(.by = cu) |> 
  left_join(select(length_models, cu, model)) |> 
  rowwise() |> 
  mutate(pred_poh = list(predict(model, data))) |> 
  select(-model) |> 
  unnest(c(data, pred_poh)) |> 
  mutate(
    poh_infill = if_else(is.na(poh_length) & !is.na(pred_poh), 1, 0),
    poh_length = if_else(is.na(poh_length), pred_poh, poh_length)
  )
  

# Check consistency between predicted and observed POH
adult_bio1 |> 
  filter(poh_infill == 0) |> 
  ggplot(aes(pred_poh, poh_length)) +
  geom_point(alpha = .1) +
  geom_abline(colour = "red") +
  coord_fixed(
    ratio = 1,
    xlim = c(180, 780),
    ylim = c(180, 780)
  )
# Some outliers but it mostly looks good
# A bit surprised the outliers weren't better mitigated by fitting
# the cu-specific models

# Clean up sex assignments ------------------------------------------------


# See what the variety of sex assignments is
adult_bio1 |> 
  count(across(contains("sex"))) |> 
  print(n = 100)
# Messy!


# Use some intuitive rules to definitively assign sex
adult_bio2 <- adult_bio1 |> 
  mutate(
    sex = case_when(
      # Give "sex_resolved" assignments priority over "sex_as_sampled"
      str_detect(sex_resolved, "(?i)^(f|female)$") ~ "Female",
      str_detect(sex_resolved, "(?i)^(m|male)$") ~ "Male",
      str_detect(sex_resolved, "(?i)^j$") ~ "Male",
      str_detect(sex_resolved, "(?i)jack|jimmy") ~ "Male",
      str_detect(sex_resolved, "(?i)jill") ~ "Female",
      # Later entries in this statement are given lower priority
      str_detect(sex_as_sampled, "(?i)^(f|female)$") ~ "Female",
      str_detect(sex_as_sampled, "(?i)^(m|male)$") ~ "Male",
      str_detect(sex_as_sampled, "(?i)^j$") ~ "Male",
      str_detect(sex_as_sampled, "(?i)jack|jimmy") ~ "Male",
      str_detect(sex_as_sampled, "(?i)jill") ~ "Female",
      .default = NA_character_
    )
  ) %T>%
  # Compare new assignments against previous (i.e. ensure case_when() worked
  # as expected)
  {print(count(., across(contains("sex"))), n = 100)}
  
  

# Clean up sample_type and other cols -------------------------------------


# What are the age data like?
adult_bio2 |> 
  count(gilbert_rich_age) |> 
  print(n=100)


adult_bio3 <- adult_bio2 |> 
  mutate(
    sample_type = case_when(
      str_detect(sample_type, "(?i)^es") ~ "Escapement",
      str_detect(sample_type, "(?i)^maa|chucklesaht") ~ "Maa-nulth",
      str_detect(sample_type, "(?i)^sf|^nf|commercial - fn") ~ "Economic Opportunity",
      str_detect(sample_type, "(?i)sport") ~ "Sport",
      str_detect(sample_type, "(?i)tf|test") ~ "Test fishery",
      str_detect(sample_type, "(?i)cf|commercial") ~ "Commercial",
      .default = sample_type
    ),
    gilbert_rich_age = if_else(
      gilbert_rich_age %in% paste0(rep(2:7, each = 4), 1:4),
      gilbert_rich_age,
      NA
    ),
    total_age = as.numeric(str_extract(gilbert_rich_age, "^\\d{1}")),
    brood_year = year - total_age
  ) %T>%
  {print(count(., sample_type))}



# Plots -------------------------------------------------------------------


# Use escapement data only for now
plot_data <- adult_bio3 |> 
  filter(sample_type == "Escapement")


# Calculate interannual mean POH lengths
poh_grand_means <- plot_data |> 
  filter(!if_any(c(cu, sex, poh_length), is.na)) |> 
  summarize(
    .by = c(cu, sex),
    poh_length_sum = sum(poh_length),
    n = n(),
    q99 = quantile(poh_length, 0.997)
  ) |> 
  mutate(mean_poh = poh_length_sum/n)


# Plot time series of POH length by CU and sex
(poh_ts <- plot_data |> 
    summarize(
      .by = c(cu, year, sex),
      across(
        poh_length,
        .fns = list(
          "q5" = ~quantile(.x, 0.05, na.rm = TRUE),
          "q25" = ~quantile(.x, 0.25, na.rm = TRUE),
          "q50" = ~quantile(.x, 0.5, na.rm = TRUE),
          "q75" = ~quantile(.x, 0.75, na.rm = TRUE),
          "q95" = ~quantile(.x, 0.95, na.rm = TRUE)
        ),
        .names = "{.fn}_{.col}"
      ),
      n = n()
    ) |> 
    filter(!if_any(c(cu, sex), is.na)) |> 
    ggplot(aes(x = year, y = q50_poh_length)) +
    facet_grid(
      cu ~ sex,
      scales = "free_y"
    ) +
    geom_hline(
      data = poh_grand_means,
      aes(yintercept = mean_poh),
      linetype = 2
    ) +
    geom_linerange(
      aes(
        ymin = q5_poh_length,
        ymax = q95_poh_length
      ),
      colour = "grey50",
      linewidth = 0.25
    ) +
    geom_pointrange(
      aes(
        size = n,
        ymin = q25_poh_length,
        ymax = q75_poh_length
      ),
      linewidth = 0.5,
      fill = "white"
    ) +
    # geom_text(
    #   data = poh_grand_means,
    #   aes(
    #     label = paste0("Interannual mean POH = ", round(mean_poh, 0), " mm"),
    #     y = q99,
    #     x = 2003
    #   ),
    #   vjust = 1
    # ) +
    scale_size(range = c(0.03,0.3)) +
    labs(
      x = "Return year",
      y = "Post-orbital to hypural length (mm)",
      size = "Sample\nsize"
    ) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "white")
    )
)


# Save the plot
poh_ts |> 
  ggsave(
    filename = here(
      "3. outputs",
      "Plots",
      "Barkley_Sockeye_Escapement_POH-lengths_1980-2022.png"
    ),
    height = 4,
    width = 7,
    units = "in",
    dpi = "print"
  )


# Interannual mean sex ratios
sex_grand_means <- plot_data |> 
  filter(
    !if_any(c(cu, sex), is.na),
    between(brood_year, min(year) - 3, max(year) - 6)
  ) |> 
  summarize(
    .by = c(cu, sex),
    n = n()
  ) |> 
  pivot_wider(
    names_from = sex,
    values_from = n
  ) |> 
  mutate(mean_f = Female/(Female+Male))


# Plot annual sex ratios by brood year
(sex_comp_ts <- plot_data |> 
    filter(
      !if_any(c(cu, sex), is.na),
      between(brood_year, min(year) - 3, max(year) - 6)
    ) |> 
    summarize(
      .by = c(brood_year, cu, sex),
      n = n()
    ) |> 
    mutate(
      .by = c(brood_year, cu),
      N = sum(n)
    ) |> 
    filter(N > 49) |> 
    ggplot(aes(x = brood_year)) +
    facet_wrap(
      ~cu,
      ncol = 1,
      strip.position = "right"
    ) +
    geom_col(
      aes(
        fill = sex, 
        y = n
      ),
      colour = "black",
      width = 1,
      position = "fill"
    ) +
    geom_label(
      data = sex_grand_means,
      aes(
        label = paste0("Long-term mean = ", round(mean_f*100, 0), "% female"),
        y = 0.9,
        x = 1995
      ),
      vjust = 1
    ) +
    scale_fill_manual(
      name = "Sex",
      values = c("grey90", "grey20")
    ) +
    scale_y_continuous(
      name = "Percentage in escapement samples",
      expand = c(0, 0),
      labels = scales::percent
    ) +
    scale_x_continuous(
      name = "Brood year",
      expand = c(0, 0)
    ) +
    theme(
      panel.grid = element_blank(),
      panel.spacing.y = unit(1, "line"),
      strip.background = element_rect(fill = "white")
    )
)

# Save the plot
sex_comp_ts |> 
  ggsave(
    filename = here(
      "3. outputs",
      "Plots",
      "Barkley_Sockeye_Escapement_sex-composition_1977-2016.png"
    ),
    height = 4,
    width = 7,
    units = "in",
    dpi = "print"
  )


# Age composition by sex
(age_sex_ts <- plot_data |> 
    filter(
      !if_any(c(cu, gilbert_rich_age, sex), is.na),
      between(brood_year, min(year) - 3, max(year) - 6),
      gilbert_rich_age %in% c(32, 42, 43, 52, 53, 62, 63)
    ) |> 
    mutate(
      gilbert_rich_age = paste0(
        str_sub(gilbert_rich_age, 1, 1),
        "[",
        str_sub(gilbert_rich_age, 2, 2),
        "]"
      )
    ) |> 
    summarize(
      .by = c(brood_year, cu, sex, gilbert_rich_age),
      n = n()
    ) |> 
    mutate(
      .by = c(brood_year, cu, sex),
      N = sum(n)
    ) |> 
    filter(N > 25) |> 
    ggplot(aes(x = brood_year)) +
    facet_grid(cu ~ sex) +
    geom_col(
      aes(
        fill = fct_rev(gilbert_rich_age), 
        y = n
      ),
      width = 1,
      position = "fill"
    ) +
    scale_fill_brewer(
      palette = "Paired",
      labels = scales::parse_format(),
      aesthetics = c("colour", "fill"),
      name = "Gilbert-\nRich age"
    ) +
    scale_y_continuous(
      name = "Percentage in escapement samples",
      expand = c(0, 0),
      labels = scales::percent
    ) +
    scale_x_continuous(
      name = "Brood year",
      expand = c(0, 0)
    ) +
    theme(
      panel.grid = element_blank(),
      panel.spacing.y = unit(1, "line"),
      strip.background = element_rect(fill = "white")
    )
)


# Save the plot
age_sex_ts |> 
  ggsave(
    filename = here(
      "3. outputs",
      "Plots",
      "Barkley_Sockeye_Escapement_sex-age-composition_1977-2016.png"
    ),
    height = 4,
    width = 7,
    units = "in",
    dpi = "print"
  )



# Table for ResDoc that has interannual median age compositions -----------


adult_bio3 |> 
  filter(
    !is.na(cu),
    !if_any(c(poh_length, gilbert_rich_age), is.na),
    gilbert_rich_age %in% c(32, 42, 43, 52, 53, 54, 62, 63, 64, 73, 74)
  ) |> 
  summarize(
    .by = c(gilbert_rich_age, cu),
    value = paste0(round(median(poh_length), 0), " (", n(), ")")
  ) |> 
  pivot_wider(names_from = cu) |> 
  arrange(gilbert_rich_age) |> 
  write.table(
    "clipboard",
    sep = ";",
    row.names = FALSE,
    quote = FALSE
  )

