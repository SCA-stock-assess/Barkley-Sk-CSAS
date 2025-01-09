# Packages ----------------------------------------------------------------

pkgs <- c("tidyverse", "MBA","reshape2","magrittr","RColorBrewer", "readxl", "here")
#install.packages(pkgs)

library(here)
library(readxl)
library(tidyverse); theme_set(theme_bw(base_size=15))
library(fuzzyjoin)
library(gganimate)
library(RColorBrewer)
library(magrittr)
library(MBA)
library(reshape2)

# Colour palette for the temp-oxy index
idx_pal <- c("#252525",brewer.pal(n = 11, name = "RdYlGn")[-c(6,11)])

# Analysis year
curr_yr <- 2021

# Load and manipulate data ---------------------------------------------------


# Order of sites
site_order <- c("5km Center", "Polly's Pt.", "Hohm Isle", "Estuary")


# Load the raw data
hs0 <- here(
  "1. data",
  "Catalyst_HarbourSurvey_2021_all.xlsx"
) |> 
  read_xlsx() |> 
  select(-pH) |> 
  mutate(
    date = as.Date(date),
    julian = date |> format("%j") |> as.integer(),
    site = case_when(
      site == "River/Outfall" & str_detect(depth, "(?i)river") ~ "River",
      site == "River/Outfall" ~ "Estuary",
      site == "Hohm Island 2" ~ "Hohm Isle",
      site == "Polly's Point 2" ~ "Polly's Pt.",
      TRUE ~ site
    ),
    # Add categorical index based on Howard Stiff's recommendations
    idx_cat = case_when(
      temp < 12 & do_mgl > 4 ~ 5,
      between(temp, 12,16) & do_mgl > 4 ~ 4,
      between(temp, 16,18) & do_mgl > 4 ~ 3,
      temp < 18 & between(do_mgl,3,4) ~ 2,
      between(temp, 18,24) | between(do_mgl,2,3) ~ 1,
      temp > 24 | do_mgl < 2 ~ 0 # Catastrophic levels
    )
  ) |> 
  filter(site != "River")  |>  
  mutate(
    across(!c(date, site), as.double),
    site = factor(site, levels = site_order)
  ) 


# Add data to expanded dataframe with all depths for each site to infill missing observations
# (assumes uniform water column below bottom measurement)
hs <- hs0 |> 
  complete(nesting(date, julian), nesting(site, depth)) |> 
  mutate(
    dist_km = case_when(
      site == "5km Center" ~ 5,
      site == "Hohm Isle" ~ 2.24,
      site == "Polly's Pt." ~ 3.45,
      site == "Estuary" ~ 0),
    # Shift the measurements at river and 5km mark for plotting
    plotting_shift = case_when(
      dist_km == 5 ~ 4.95, 
      dist_km == 0 ~ 0.05, 
      TRUE ~ dist_km
    ),
    # Add lat/long data for each site
    lat = case_when(
      site == "5km Center" ~ 49.201260,
      site == "Polly's Pt." ~ 49.215729,
      site == "Hohm Isle" ~ 49.227984,
      site == "Estuary" ~ 49.247812
    ),
    lon = case_when(
      site == "5km Center" ~ -124.820737,
      site == "Polly's Pt." ~ -124.821066,
      site == "Hohm Isle" ~ -124.823543,
      site == "Estuary" ~ -124.816844
    )
  ) |> 
  arrange(site, julian, depth) #|> 
  # Infill missing data for deeper depths with closest measurement
  # fill(
  #   c(salinity, temp, contains("do"), idx_cat), 
  #   .direction = "down"
  # )
  

# Extract max depths
max_depths <- hs0 |> 
  group_by(site) |> 
  summarize(max_depth = max(depth))



# Define plotting function and interpolate time series --------------------



# set date range for plots.
max_jul <- max(hs$julian) 
min_jul <- min(hs$julian)


# Create function that interpolates values between surveys
intp_fn <- function(station, var) {
  
  no <- max_jul - min_jul
  
  xyz <- hs |> 
    filter(site == station, !is.na({{var}})) |> 
    select(julian, depth, {{var}})
  
  df <- xyz |> 
    mba.surf(
      no.X = no, 
      no.Y = no, 
      extend = T,
      # Parameters to control the stretch along the x axis,
      # Higher n:m ratio means more x stretch, less y stretch
      n = 1, 
      m = (max(xyz$julian)-min(xyz$julian))/max(xyz$depth)
    ) 
  
  dimnames(df$xyz.est$z) <- list(df$xyz.est$x, df$xyz.est$y)
  
  melt(
    df$xyz.est$z, 
    varnames = c('julian', 'depth')
    ) |> 
    mutate(
      #julian = round(julian, 0),
      date = as.POSIXct.numeric(julian*24*3600, origin = paste0(curr_yr-1, "-12-31"))
    ) |>  
    summarize(
      .by = c(date, julian, depth),
      value = mean(value)
    )
}


# Raster plot for time series data (based on temp data, tweaks needed for DO and salinity)
ts_p_fn <- function(df, hs_var, min_julian, max_julian) {
  
  point_data <- hs0 |> 
    filter(between(julian, min_julian, max_julian)) |> # Truncate raw data
    drop_na(salinity:do_sat)
  
  plot <- df |> 
    filter(between(julian, min_julian, max_julian)) |> # Truncate interpolated data to desired time period
    ggplot(aes(x = date, y = depth)) +
    facet_grid(site~., scales = "free_y", space = "free_y") +
    geom_raster(aes(fill = value), interpolate = TRUE, alpha = 0.75) +
    geom_point(
      data = point_data, 
      aes(x = as.POSIXct(date), y = depth, colour = idx_cat),
      size = 1, 
      #alpha = 0.15,
      #shape = 3
    ) +
    # Custom colour scale based on the index and Howard Stiff's suggestions
    scale_fill_gradientn(
      colours = idx_pal,
      name = "Temp-oxy\nindex",
      limits = c(0,5),
      aesthetics = c("colour", "fill")
    ) +
    geom_contour(
      aes(z = value), 
      binwidth = 1, 
      colour = "black", 
      alpha = 0.2
    ) +
    labs(y = "Depth (m)", x = NULL) +
    scale_y_reverse(expand = c(0,0), labels = as.integer) +
    scale_x_datetime(
      expand = c(0,0), 
      breaks = "2 weeks", 
      date_labels = "%d %b"
    ) +
    theme(
      strip.background = element_rect(fill = "white"),
      legend.position = "inside",
      legend.position.inside = c(0.98,0.60),
      legend.justification.inside = c("right", "bottom"),
      legend.direction = "horizontal",
      legend.background = element_rect(colour = "black",fill = alpha("white",0.75)),
      panel.spacing.y = unit(1, "lines"),
      axis.title.y = element_text(vjust = 2)
    )
  
  return(plot)
}


# Apply function to all sites in the data and combine into single dfs
sites <- purrr::set_names(site_order) # List of site names


# Index
idx <- sites |> 
  map(\(x) intp_fn(x, idx_cat)) |> 
  list_rbind(names_to = "site") |> 
  mutate(
    site = factor(site, levels = site_order),
    # Constrain values to [5,0]
    # Equation below from: https://stats.stackexchange.com/a/281164
    value = ((5-0)*((value - min(value))/(max(value)-min(value)))) 
  )


# Plot from the continuous index 
(idx_p <- ts_p_fn(idx, "hold_idx", min_jul, max_jul))


# Save the plot
ggsave(
  filename = here(
    "3. outputs",
    "Plots",
    "Temp-oxy-index_time-series_2021.png"
  ),
  plot = idx_p,
  height = 8,
  width = 6.5,
  units = "in",
  dpi = "print"
)



  

# Plot holding conditions index ------------------------------------------------


# Interpolate time series of the continuous index
idx <- map_df(sites, ~intp_fn(.x, hold_idx, stretch = (1/1)), .id = "site") |> 
  mutate(site = factor(site, levels = site_order))









# Plot from the categorical index
idx_c <- map_df(sites, ~intp_fn(.x, idx_cat, stretch = (1/1)), .id = "site") |> 



idxcat_p <- ts_p_fn(idx_c, "idx_cat", min_jul, max_jul) + 
  

idxcat_p[["layers"]][[3]][["stat_params"]][["binwidth"]] <- 1 # Change contour spacing

idxcat_p


# Save time series plot
ggsave(
  filename = paste0(
    here("Harbour Survey", "plots"),
    "/R-PLOT_time series temp-oxy index",
    Sys.Date(),
    ".png"
  ),
  plot = idxcat_p,
  height = 7,
  width = 9,
  units = "in"
)


# Clean up the workspace
rm(idx, idx_c, idx_p, idxcat_p)
gc()



# Cross-section of the Inlet from most recent survey ----------------------


# Create function that interpolates values spatially between sites
intp_dist_fn <- function(var, day_of_year) {
  df <- hs |> 
    filter(julian == day_of_year) |> 
    #Add some rows to extend interpolated values down to 60 m. 
    #Assumes linear water column below bottom deepest measurement at 5km site
    add_row(depth = c(55,60), site = rep("5km Center", 2)) |> 
    #add_row(depth = c(25,30), site = rep("Polly's Pt.", 2)) |> # Hacky fix for 7 July 
    #filter(!(site %in% site_order[-1] & is.na(.data[[var]]))) |> #.data[[]] for when variable names are given in quotes
    arrange(site, depth) |> 
    group_by(site) |> 
    fill(everything(), .direction = "down") |> 
    ungroup() |> 
    # Now keep only the columns needed for interpolation
    select(dist_km, depth, .data[[var]]) |> 
    #add_row(depth = c(30,40), idx_cat = c(2, 0.5), dist_km = c(3.45, 3.45)) |> # Hacky fix for 16 June
    filter(!is.na(.data[[var]])) |>  # Added step on 21 June 2023. Might help when one site missing data?
    mba.surf(
      no.X = 61, no.Y = 61, 
      extend = T,
      n = 6, m = 1 # Parameters to control the stretch along the x axis
    ) 
  
  dimnames(df$xyz.est$z) <- list(df$xyz.est$x, df$xyz.est$y)
  
  melt(df$xyz.est$z, varnames = c('dist_km', 'depth')) |> 
    mutate(julian = day_of_year)
}



# Manually input bathymetric data from web application: https://data.chs-shc.ca/map
bathy <- data.frame(
  dist_km = c(-0.1,-0.1,0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5,4.75,4.8,4.85, 5.01),
  depth = c(61,10,12, 15, 18, 20, 24, 30, 33, 33, 35, 35,46,55,61)
)

# Plotting function
xs_p_fn <- function(df, date_tag, ...) {
  df %>% 
    # Merge in columns delineating where actual observations occurred
    difference_left_join(
      hs0 |> 
        drop_na(salinity:do_sat) |> 
        select(julian, depth, dist_km, plotting_shift) |> 
        mutate(depth_sample = depth),
      by = c("julian", "depth", "dist_km"),
      max_dist = 0.2
    ) |>
    select(-contains(".y")) |> 
    rename_with(~ gsub(".x", "", .x)) |> 
    ggplot(aes(x = dist_km, y = depth), ...) +
    geom_raster(aes(fill = value), interpolate = TRUE, alpha = 0.75) +
    geom_point(
      aes(x = plotting_shift, y = depth_sample),
      size = 0.2, 
      alpha = 0.05,
      shape = 8
    ) +
    #geom_contour(aes(z = value), binwidth = 1, colour = "black", alpha = 0.2) +
    geom_polygon(data = bathy, fill = "grey80", colour = "black") +
    labs(
      y = "Depth (m)", 
      x = "Distance from river mouth (km)",
      title = paste("Survey date:", date_tag)
    ) +
    coord_cartesian(
      xlim = c(0, max(filter(hs, !is.na(salinity), julian == max(julian))$dist_km)) + 0.01, 
      ylim = c(60, 0), 
      expand = F
    ) +
    scale_y_reverse(labels = as.integer) +
    theme(
      strip.background = element_rect(fill = "white"),
      legend.position = c(0.02,0.02),
      legend.direction = "horizontal",
      legend.justification = c("left", "bottom"),
      legend.background = element_rect(colour = "black"),
      panel.spacing.y = unit(1, "lines"),
      axis.title.y = element_text(vjust = 2)
    )
}

# Take maximum date where all sites have data:
complete_survey <- hs |> 
  group_by(julian, site) |>
  filter(!all(is.na(idx_cat))) |> 
  group_by(julian) |> 
  filter(!length(unique(site)) < 4) |> 
  ungroup() %>% 
  select(julian) |> 
  max()


# Plot for temp-oxy index
(xs_plot_idx <- xs_p_fn(
  intp_dist_fn("idx_cat", max_jul) %>% 
    # Squish the interpolated index values into [0,5]
    mutate(value = (5-0)*((value - min(value))/(max(value)-min(value)))+0),
  date_tag = as.Date(
    max_jul, 
    format = "%j", 
    origin = as.Date(paste0(curr_yr-1, "-12-31"))
  ) |> 
    format("%d-%b")
) +
    scale_fill_gradientn(
      colours = idx_pal,
      name = "Temp-oxy\nindex", 
      limits = c(0,5)
    ) 
)


# Plot for Dissolved Oxygen
(xs_plot_do <- xs_p_fn(
  intp_dist_fn("do_mgl", 208),
  date_tag = as.Date(
    208, 
    format = "%j", 
    origin = as.Date(paste0(curr_yr-1, "-12-31"))
  ) |> 
    format("%d-%b")
) +
    scale_fill_gradientn(
      colours = do_pal$colour, 
      name = "DO (mg/L)",
      labels = as.integer, 
      breaks = c(0, 4, 8, 12),
      values = scales::rescale(do_pal$value, to = c(0, 1)),
    )
)


# Save plots
list(xs_plot_do, xs_plot_idx) |> 
  set_names(c("DO", "Index")) |> 
  iwalk(
    ~ggsave(
      plot = .x, 
      filename = paste0(
        here("Harbour Survey", "plots"),
        "/R-PLOT_Inlet cross section",
        .y,
        Sys.Date(),
        ".png"
      ),
      height = 3,
      width = 8,
      units = "in"
    )
  )


# Animated time series of Inlet cross-section -----------------------------



# Apply function to all dates in the data and combine into single df
anim_data <- hs |> 
  filter(
    between(
      date, 
      as.Date(paste0(curr_yr, "-03-31")),
      as.Date(paste0(curr_yr, "-10-31"))
    )
  ) |> 
  distinct(julian) |> 
  pull() |> 
  as.list() |> 
  purrr::set_names() |> 
  map_dfr(
    ~ intp_dist_fn("idx_cat", .x) |> 
      mutate(value = (5-0) * ((value - min(value)) / (max(value) - min(value))) + 0),
    .id = "julian"
  ) |> 
  mutate(
    julian = as.numeric(julian),
    date = as.Date(
      julian, 
      format = "%j", 
      origin = as.Date(paste0(curr_yr-1, "-12-31"))
    )
  ) |> 
  filter(is.finite(value)) # Drop NaN values


# Save the plot as an animation with frames for each date
anim <- xs_p_fn(
  anim_data, 
  aes(group = date, frame = date)
) +
  scale_fill_gradientn(
    colours = idx_pal,
    name = "Temp-oxy\nindex", 
    limits = c(0,5)
  ) +
  transition_time(date) +
  labs(title = "{frame_time}")


# Customize the animation
animate(
  anim, 
  width = 800, height = 400,
  fps = 20, 
  duration = 45,
  renderer = gifski_renderer()
)


# Save
anim_save(
  here(
    "Harbour Survey", 
    "plots", 
    "inlet_cross-section_animation.gif"
    )
  )

