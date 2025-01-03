
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
  "return by age time series.xlsx"
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
    sockeye = escapement + catch
  ) |> 
  select(year, stock, sockeye, age)
  

huc_data <- here(
  "1. data",
  "return by age time series.xlsx"
) |> 
  read_xlsx(sheet = "Henderson") |> 
  mutate(across(matches("age_\\d+"), \(x) x * age_sample_size)) |> 
  summarize(
    .by = c(year, escapement, catch),
    across(contains("age_"), sum)
  ) |> 
  mutate(
    across(matches("age_\\d+"), \(x) x / age_sample_size),
    run = escapement + catch,
    across(matches("age_\\d+"), \(x) x * run)
  ) |> 
  pivot_longer(
    matches("age_\\d+"),
    names_prefix = "age_",
    names_to = "age",
    values_to = "sockeye"
  ) |> 
  select(year, age, sockeye) |> 
  filter(
    .by = year, 
    !if_any(sockeye, is.na)
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
  mutate(stock = factor(stock, levels = c("GCL", "SPR", "HUC")))


# Plot age composition time series for each stock -------------------------


# Formatted base plot
(ages_base_p <- run_ages |> 
   filter(sockeye > 1) |> 
   ggplot(
     aes(
       year, 
       sockeye,
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
(ages_stream_p <- ages_base_p +
    geom_stream(
      type = "proportional",
      bw = 0.4
    )
)


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


# Save the column and stream charts
list(
  "column" = ages_col_p,
  "stream" = ages_stream_p
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
