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
  "SoxBioData1980-2021.xlsx"
) |> 
  read_xlsx(sheet = "All data") |> 
  clean_names() |> 
  mutate(
    # Fix a few lengths that are 10x too large
    across(contains("length"), \(x) if_else(x>1000, x/10, x)),
    # Assign fish to CU based on sampling location
    cu = case_when(
      str_detect(site, "(?i)henderson") ~ "Hucuktlis",
      str_detect(site, "(?i)great central") ~ "Great Central",
      str_detect(site, "(?i)sproat") ~ "Sproat",
      .default = "unknown"
    )
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
len_comp_p2 %+% filter(
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
      str_detect(sex_resolved, "(?i)^f$") ~ "female",
      str_detect(sex_resolved, "(?i)^m$") ~ "male",
      str_detect(sex_resolved, "(?i)^j$") ~ "male",
      str_detect(sex_resolved, "(?i)jack|jimmy") ~ "male",
      str_detect(sex_resolved, "(?i)jill") ~ "female",
      # Later entries in this statement are given lower priority
      str_detect(sex_as_sampled, "(?i)^f$") ~ "female",
      str_detect(sex_as_sampled, "(?i)^m$") ~ "male",
      str_detect(sex_as_sampled, "(?i)^j$") ~ "male",
      str_detect(sex_as_sampled, "(?i)jack|jimmy") ~ "male",
      str_detect(sex_as_sampled, "(?i)jill") ~ "female",
      .default = NA_character_
    )
  ) %T>%
  # Compare new assignments against previous (i.e. ensure case_when() worked
  # as expected)
  {print(count(., across(contains("sex"))), n = 100)}
  
  
  
  
  