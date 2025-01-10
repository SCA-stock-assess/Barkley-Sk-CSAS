# Packages ----------------------------------------------------------------

pkgs <- c(
  "tidyverse", "MBA","reshape2", "cowplot", "magick",
  "magrittr","RColorBrewer", "readxl", "here"
)
#install.packages(pkgs)

library(here)
library(readxl)
library(tidyverse); theme_set(theme_bw(base_size=12))
library(cowplot)
library(magick)
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



# Interpolate values across dates and depths  --------------------


# set date range for plots.
max_jul <- max(hs$julian) 
min_jul <- 80


# Create function that interpolates values between surveys
intp_fn <- function(station, var) {
  
  resolution <- max_jul - min_jul
  
  xyz <- hs |> 
    filter(site == station, !is.na({{var}})) |> 
    select(julian, depth, {{var}})
  
  df <- xyz |> 
    mba.surf(
      no.X = resolution, 
      no.Y = resolution, 
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
      date = as.POSIXct.numeric(
        julian*24*3600, 
        origin = paste0(curr_yr-1, "-12-31")
      )
    )
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




  

# Plot holding conditions index ------------------------------------------------


# Data frame with observed points from the raw data
point_data <- hs0 |> 
  filter(between(julian, min_jul, max_jul)) |> # Truncate raw data
  drop_na(salinity:do_sat)


# Plot with interpolated values overlaid with raw points
(idx_ts_p <- idx |> 
  filter(between(julian, min_jul, max_jul)) |> # Truncate interpolated data to desired time period
  ggplot(aes(x = date, y = depth)) +
  facet_grid(site~., scales = "free_y", space = "free_y") +
  geom_raster(aes(fill = value), interpolate = TRUE, alpha = 0.75) +
  geom_point(
    data = point_data, 
    aes(x = as.POSIXct(date), y = depth, colour = idx_cat),
    size = 0.3, 
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
    breaks = "1 month", 
    date_labels = "%b"
  ) +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.position = "inside",
    legend.position.inside = c(0.98,0.62),
    legend.justification.inside = c("right", "bottom"),
    legend.direction = "horizontal",
    legend.background = element_rect(colour = "black",fill = alpha("white",0.75)),
    panel.spacing.y = unit(1, "lines"),
    axis.title.y = element_text(vjust = 2)
  )
)


# Save the plot
ggsave(
  filename = here(
    "3. outputs",
    "Plots",
    "Temp-oxy-index_time-series_2021.png"
  ),
  plot = idx_ts_p,
  height = 7,
  width = 6.5,
  units = "in",
  dpi = "print"
)


# Load in-river migration data for 2021
esc_data <- here(
  "1. data",
  "Somass_Sockeye_daily_escapement_data.xlsx"
) |> 
  read_xlsx() |> 
  filter(year == 2021) |> 
  mutate(ttl = `Adjusted net Adult up count` + `Adjusted net Jack up count`) |> 
  select(year, month, day, system, ttl) |> 
  mutate(date = as.POSIXct(paste(year, month, day, sep = "-"))) |> 
  summarize(
    .by = date,
    ttl = sum(ttl, na.rm = TRUE)
  ) |> 
  filter(
    date <= as.POSIXct.numeric(
      max_jul*24*3600, 
      origin = paste0(curr_yr-1, "-12-31")
    )
  )

comp_site <- "Estuary"
ratio <- max(idx[idx$site == comp_site,]$depth)/max(esc_data$ttl)


# Plot just Estuary for comparison of in-river migration curves
heat_dome_plot <- idx |> 
  filter(
    between(julian, min_jul, max_jul),
    site == comp_site
  ) |> 
  mutate(depth = -depth + max(depth)) |> 
  ggplot(aes(x = date, y = depth)) +
  geom_raster(aes(fill = value), interpolate = TRUE, alpha = 0.75) +
  geom_contour(
    aes(z = value), 
    binwidth = 1, 
    colour = "black", 
    alpha = 0.2
  ) +
  annotate(
    geom = "rect",
    xmin = as.POSIXct("2021-06-25"),
    xmax = as.POSIXct("2021-07-07"),
    ymin = -Inf,
    ymax = Inf,
    lty = 2,
    fill = "grey",
    alpha = 0.3,
    colour = "grey60"
  ) +
  annotate(
    "text",
    x = as.POSIXct("2021-07-02"),
    y = 8,
    label = "Heat dome",
    angle = 270
  ) +
  geom_area(
    data = esc_data,
    aes(
      x = date - 4,
      y = ttl*ratio
    ),
    colour = "black",
    fill = "black",
    alpha = 0.3,
    linewidth = 0.8
  ) +
  # Custom colour scale based on the index and Howard Stiff's suggestions
  scale_fill_gradientn(
    colours = idx_pal,
    name = "Temp-oxy\nindex",
    limits = c(0,5),
    aesthetics = c("colour", "fill")
  ) +
  labs(
    title = "Somass Estuary temperature-oxygen profile in 2021",
    y = "Depth (m)",
    x = NULL
    ) +
  scale_y_continuous(
    expand = c(0,0), 
    labels = function(x) -(x - max(idx[idx$site == comp_site,]$depth)),
    sec.axis = sec_axis(
      transform = ~./ratio,
      name = "Daily Somass Sockeye Escapement",
      labels = scales::label_number()
    )
  ) +
  scale_x_datetime(
    expand = c(0,0), 
    breaks = "1 month", 
    date_labels = "%b"
  ) +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.position = "inside",
    legend.position.inside = c(0.98,0.98),
    legend.justification.inside = c("right", "top"),
    legend.direction = "horizontal",
    legend.background = element_rect(colour = "black",fill = alpha("white",0.75)),
    panel.spacing.y = unit(1, "lines"),
    axis.title.y = element_text(vjust = 2)
  )


# Add heat dome map image
image <- image_read( 
  here(
    "1. data",
    "WNA_Heat_Wave_Temp_Anomaly.jpg"
  )
)


# Add image
heat_dome_plot2 <- heat_dome_plot +
  annotate(
    "rect",
    xmin = as.POSIXct("2021-03-28"),
    xmax = as.POSIXct("2021-05-25"),
    ymin = 6, 
    ymax = 11,
    colour = "black",
    linewidth = 2
  ) +
  annotation_raster(
    image,
    xmin = as.POSIXct("2021-03-28.5"),
    xmax = as.POSIXct("2021-05-25.5"),
    ymin = 6, 
    ymax = 11
  ) +
  annotate(
    "text",
    label = "Image source: NASA Earth Observatory",
    x = as.POSIXct("2021-03-29"),
    hjust = 0,
    y = 6.15,
    size = 2,
    colour = "grey30"
  )



# Cross-section of the Inlet during heat dome ----------------------


# Day to interpolate cross section data for
intp_day <- 181


# Interpolates values spatially between sites
xs_intp <- hs |> 
  filter(julian == intp_day) |> 
  # Remove the 5km site, which was not surveyed due to the extreme heat
  filter(
    .by = site,
    !all(is.na(idx_cat))
  ) |> 
  arrange(site, depth) |> 
  group_by(site) |> 
  # Fill down observed values to deeper waters (i.e. assume uniform conditions)
  # Might not be accurate but it's the best option we have
  fill(everything(), .direction = "down") |> 
  ungroup() |> 
  select(dist_km, depth, idx_cat) |> 
  mba.surf(
    no.X = 61, no.Y = 61, 
    extend = T,
    n = 6, m = 1 
  ) 

dimnames(xs_intp$xyz.est$z) <- list(xs_intp$xyz.est$x, xs_intp$xyz.est$y)

xs_xyz <- melt(xs_intp$xyz.est$z, varnames = c('dist_km', 'depth')) |> 
  mutate(julian = intp_day)


# Manually input bathymetric data from web application: https://data.chs-shc.ca/map
bathy <- data.frame(
  dist_km = c(-0.1,-0.1,0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5,4.75,4.8,4.85, 5.01),
  depth = c(61,10,12, 15, 18, 20, 24, 30, 33, 33, 35, 35,46,55,61)
)


# Plot cross section data
(heat_dome_xs_p <- xs_xyz |> 
    mutate(value = (5-0)*((value - min(value))/(max(value)-min(value)))+0) |> 
    ggplot(aes(x = dist_km, y = depth)) +
    geom_raster(aes(fill = value), interpolate = TRUE, alpha = 0.75) +
    geom_point(
      data = filter(hs, julian == intp_day, !is.na(idx_cat)),
      aes(x = plotting_shift),
      size = 0.5,
      alpha = 0.6,
      shape = 8
    ) +
    geom_contour(
      aes(z = value), 
      binwidth = 1, 
      colour = "black", 
      alpha = 0.2
    ) +
    geom_polygon(data = bathy, fill = "grey80", colour = "black") +
    labs(
      y = "Depth (m)", 
      x = "Distance from river mouth (km)",
      title = paste(
        "Inlet cross-section", 
        as.POSIXct.numeric(
          intp_day*24*3600, 
          origin = paste0(curr_yr-1, "-12-31")
        ) |> 
          format("%d %b %Y")
      )
    ) +
    coord_cartesian(
      xlim = c(0, max(xs_xyz$dist_km)) + 0.01, 
      ylim = c(max(xs_xyz$depth), 0), 
      expand = F
    ) +
    scale_fill_gradientn(
      colours = idx_pal,
      name = "Temp-\noxy\nindex", 
      limits = c(0,5)
    ) +
    scale_y_reverse(labels = as.integer) +
    guides(colour = "none") +
    theme(axis.title.y = element_text(vjust = 2))
)





# Build combined time series and cross section plot from heat dome --------


# Lay plots in grid
(layout1 <- plot_grid(
  heat_dome_plot2 + theme(legend.position = "none"),
  heat_dome_xs_p,
  labels = c("A", "B"),
  rel_heights = c(1, 0.7),
  ncol = 1
)
)


# Save the heat dome plot
ggsave(
  plot = layout1,
  filename = here(
    "3. outputs",
    "Plots",
    "Heat_Dome_temp-oxy_plot.png"
  ),
  width = 8,
  height = 7,
  units = "in",
  dpi = "print"
)
