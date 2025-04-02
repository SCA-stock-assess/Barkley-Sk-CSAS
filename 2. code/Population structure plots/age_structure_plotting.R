# Packages ----------------------------------------------------------------

pkgs <- c("here", "tidyverse", "readxl", "ggstream")
#install.packages(pkgs)

library(here)
library(tidyverse); theme_set(theme_bw())
library(readxl)
library(ggstream)


# Load stock-recruit data -------------------------------------------------


# Need raw data, which contains fw and marine ages
som_data <- here(
  "1. data",
  "Barkley Sockeye time series of returns.xlsx"
) |> 
  read_xlsx(sheet = "Somass") |> 
  rename("year" = return_year) |> 
  mutate(
    age = paste0(
      ttl_age,
      "[",
      fw_age,
      "]"
    ),
    sockeye = as.numeric(escapement) + catch
  ) |> 
  mutate(
    .by = c(year, stock),
    ttl = sum(sockeye),
    prop = sockeye/sum(sockeye)
  ) |> 
  select(year, stock, prop, age)
  

huc_data <- here(
  "1. data",
  "Barkley Sockeye time series of returns.xlsx"
) |> 
  read_xlsx(sheet = "Hucuktlis") |> 
  mutate(across(matches("age_\\d+"), \(x) x * age_sample_size)) |> 
  summarize(
    .by = c(year, escapement, catch),
    across(matches("age_(sample_size|\\d+)"), sum)
  ) |> 
  mutate(across(matches("age_\\d+"), \(x) x / age_sample_size)) |> 
  pivot_longer(
    matches("age_\\d+"),
    names_prefix = "age_",
    names_to = "age",
    values_to = "prop"
  ) |> 
  select(year, age, prop) |> 
  filter(
    .by = year, 
    !any(is.na(prop))
  ) |> 
  mutate(    
    stock = "HUC",
    age = paste0(
      str_sub(age, 1, 1),
      "[",
      str_sub(age, 2, 2),
      "]"
    )
  ) 


run_ages <- bind_rows(huc_data, som_data) |> 
  mutate(stock = factor(stock, levels = c("GCL", "SPR", "HUC"))) |> 
  filter(!str_detect(age, "NA"))


# Plot age composition time series for each stock -------------------------


# Formatted base plot
(ages_base_p <- run_ages |> 
   complete(stock, age, year, fill = list(prop = 0)) |> 
   ggplot(
     aes(
       year, 
       prop,
       fill = fct_rev(age)
     )
   ) +
   facet_wrap(
     ~stock,
     ncol = 1,
     labeller = labeller(
       stock = c(
         "GCL" = "Great Central",
         "SPR" = "Sproat",
         "HUC" =  "Hucuktlis"
       )
     ),
     strip.position = "right"
   ) +
   scale_fill_brewer(
     palette = "Paired",
     labels = scales::parse_format(),
     aesthetics = c("colour", "fill"),
     name = "Gilbert-\nRich age"
   ) +
   scale_y_continuous(labels = scales::percent) +
   coord_cartesian(expand = FALSE) +
   labs(
     x = "Adult return year", 
     y = "Age composition"
   ) +
   theme(
     panel.grid = element_blank(),
     strip.background = element_rect(colour = "black", fill = "white"),
     panel.spacing.y = unit(1, "lines")
   )
)


# Version using stream chart
 # (ages_stream_p <- ages_base_p +
 #     geom_stream(
 #       aes(y = prop*100),
 #       type = "proportional",
 #       bw = 0.4
 #     )
 # )
# Can't get Hucuktlis data to fill correctly--ggstream glitch?

# Version using column chart
(ages_col_p <- ages_base_p +
    geom_col(
      position = "fill",
      width = 1
    ) 
)


# Version using stacked area
(ages_area_p <- ages_base_p +
  geom_area(
    position = "fill"
  )
)


# Save the column and area charts
list(
  "column" = ages_col_p,
  "area" = ages_area_p
) |> 
  iwalk(
    \(x, idx) ggsave(
      plot = x,
      filename = here(
        "3. outputs",
        "Plots",
        paste0(
          "Stock_age-composition_",
          idx,
          ".png"
        )
      ),
      dpi = "print",
      width = 6.5,
      height = 5,
      units = "in"
    )
  )


# Plot age compositions by brood year -------------------------------------



# Need raw data, which contains fw and marine ages
som_data_brood <- here(
  "1. data",
  "Barkley Sockeye time series of returns.xlsx"
) |> 
  read_xlsx(sheet = "Somass") |> 
  mutate(
    age = paste0(
      ttl_age,
      "[",
      fw_age,
      "]"
    ),
    sockeye = as.numeric(escapement) + catch,
    brood_year = return_year - ttl_age
  ) |> 
  rename("year" = brood_year) |> 
  mutate(
    .by = c(year, stock),
    ttl = sum(sockeye),
    prop = sockeye/sum(sockeye)
  ) |> 
  select(year, stock, prop, age)


huc_data_brood <- here(
  "1. data",
  "Barkley Sockeye time series of returns.xlsx"
) |> 
  read_xlsx(sheet = "Hucuktlis") |> 
  mutate(across(matches("age_\\d+"), \(x) x * age_sample_size)) |> 
  summarize(
    .by = c(year, escapement, catch),
    across(matches("age_(sample_size|\\d+)"), sum)
  ) |> 
  mutate(across(matches("age_\\d+"), \(x) x / age_sample_size)) |> 
  pivot_longer(
    matches("age_\\d+"),
    names_prefix = "age_",
    names_to = "age",
    values_to = "prop"
  ) |> 
  filter(
    .by = year, 
    !any(is.na(prop))
  ) |> 
  mutate(
    catch = if_else(is.na(catch), 100, catch), # Hacky fix
    stock = "HUC",
    age = paste0(
      str_sub(age, 1, 1),
      "[",
      str_sub(age, 2, 2),
      "]"
    ),
    ttl_age = as.numeric(str_sub(age, 1, 1)),
    brood_year = year - ttl_age,
    sockeye = (catch + escapement) * prop
  ) |> 
  summarize(
    .by = c(brood_year, age, stock),
    sockeye = sum(sockeye)
  ) |> 
  mutate(
    .by = brood_year,
    ttl = sum(sockeye),
    prop = sockeye/ttl
  ) |> 
  select(year = brood_year, age, stock, prop)


run_ages_brood <- bind_rows(huc_data_brood, som_data_brood) |> 
  mutate(stock = factor(stock, levels = c("GCL", "SPR", "HUC"))) |> 
  filter(!str_detect(age, "NA")) |>
  # Remove incomplete brood years
  add_count(stock, year) |> 
  filter(
    .by = stock,
    !n < max(n)
  )


# Update plot with brood year data
(ages_base_p2 <- ages_base_p %+% run_ages_brood +
  labs(x = "Brood year")
) 


# Add stacked columns
(ages_brood_p <- ages_base_p2 +
  geom_col(
    position = "fill",
    width = 1
  ) 
)


# Save the stacked column chart
ggsave(
  plot = ages_brood_p,
  filename = here(
    "3. outputs",
    "Plots",
    "Stock_brood-year_age-compositions.png"
    ),
  dpi = "print",
  width = 6.5,
  height = 5,
  units = "in"
)
