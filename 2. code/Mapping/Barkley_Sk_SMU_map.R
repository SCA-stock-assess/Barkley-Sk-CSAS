# Packages ----------------------------------------------------------------

pkgs <- c(
  "tidyverse", "ggrepel", "here", 
  "sf", "ggspatial", "ggOceanMaps", "cowplot", "geomtextpath"
  )
#install.packages(pkgs)
#remotes::install_github("MikkoVihtakari/ggOceanMaps") # Alt version


# Load packages
library(here)
library(tidyverse); theme_set(theme_bw(base_size = 14))
library(ggrepel)
library(sf)
library(ggspatial)
library(ggOceanMaps)
library(geomtextpath)
library(cowplot)


# Load spatial features ---------------------------------------------------


# Coordinates for the zoomed-out map
coords_big <- c(xmin = -150, xmax = -120, ymin = 46, ymax = 63)

# Bounding box for zoomed-out coastline map
bb_big <- st_bbox(
  coords_big,
  crs = "NAD83"
)


# Coordinates for the Barkley Sound map
coords_small <- c(xmin = -126.25, xmax = -124.7, ymin = 48.7, ymax = 49.45)

# Bounding box for Barkley Sound (plus a little extra room)
bb_small <- st_bbox(
  coords_small,
  crs = "NAD83"
)


# High resolution coastline data from Freshwater Atlas
vi_coastline <- read_sf(
  here(
    "1. data",
    "FWA_COASTLINES_SP", 
    "FWCSTLNSSP_line.shp" # This is the too-large file (157 MB)
  )
) |> 
  st_transform(crs = "NAD83") |> 
  # Constrain data to include only Van Isle watersheds
  filter(
    WTRSHDGRPC %in% c(
      "NEVI", "NIMP", "COMX", "TAHS", "COWN", "VICT", 
      "SANJ", "PARK", "ALBN", "CLAY", "GOLD", "BRKS",
      "TSIT", "HOLB", "CAMB", "SALM"
    )
  ) |> 
  # Unite all the line ends to form a polygon for the whole of Vancouver Island
  st_union() |> 
  st_polygonize()


# Shapefile with lake polygons from Freshwater Atlas
lakes <- read_sf(
  here(
    "1. data",
    "FWA_LAKES_POLY",
    "FWLKSPL_polygon.shp"
  )
) |> 
  rename_with(str_to_lower) |> 
  filter(str_detect(gnsnm1, "(?i)hucuktlis|great\\scentral|sproat")) |> 
  st_transform(crs = "NAD83") |> 
  mutate(CU = str_remove_all(gnsnm1, "\\sLake"))


# Vancouver island stream linestrings from BC Freshwater Atlas.
# Filter to include only Sproat/Stamp/Somass
stream_lines <- read_sf(
  here(
    "1. data",
    "FWA_STREAM_NETWORKS_SP",
    "FWSTRMNTWR_line.shp"
  )
) |> 
  rename_with(str_to_lower) |> 
  filter(str_detect(gnis_name, "(?i)sproat|stamp|somass|hucuktlis|clemens")) |> 
  mutate(
    CU = case_when(
      str_detect(gnis_name, "(?i)sproat") ~ "Sproat",
      str_detect(gnis_name, "(?i)stamp") ~ "Great Central",
      str_detect(gnis_name, "(?i)hucuktlis|clemens") ~ "Hucuktlis",
      TRUE ~ "Multi"
    )
  )
  


# Build plot with small inset --------------------------------------------


# Annotation data for "Gulf of Alaska" text
goa_ann <- data.frame(
  lon = -147:-137,
  lat = 55
) |> 
  as.matrix() |> 
  st_linestring() |> 
  list() |> 
  as_tibble_col("geometry") |> 
  mutate(text = "Gulf of\nAlaska") |> 
  st_as_sf(crs = 4326)


# Zoomed out map showing the bounding box context
(big_map <- basemap(
  limits = coords_big,
  land.col = "grey70",
  land.border.col = NA,
  grid.col = "grey95",
  rotate = TRUE
  ) +
   geom_sf(
     data = st_as_sfc(bb_small),
     colour = "red",
     fill = NA,
     linewidth = 0.5
   ) +
    geom_textsf(
      data = goa_ann,
      aes(label = text),
      linecolour = NA,
      size = 3
    ) +
   theme(
     axis.ticks.x = element_blank(),
     axis.ticks.y = element_blank(),
     axis.text = element_blank(),
     axis.title = element_blank(),
     panel.background = element_rect(fill = "lightblue1"),
     panel.ontop = FALSE,
     panel.border = element_rect(
       fill = NA, 
       colour = "black", 
       linewidth = 1
     )
   )
)


# Small map showing Barkley Sound
(small_map <- ggplot(vi_coastline) +
    geom_sf(
      fill = "grey80",
      colour = "grey60"
    ) +
    geom_sf(
      data = stream_lines,
      aes(colour = CU)
    ) +
    geom_sf(
      data = lakes,
      aes(fill = CU),
      colour = "grey60"
    ) +
    geom_text_repel(
      data = lakes,
      aes(
        label = CU,
        colour = CU,
        geometry = geometry
      ),
      stat = "sf_coordinates",
      min.segment.length = Inf,
      hjust = "right",
      direction = "y",
      nudge_x = -0.2,
      nudge_y = 0.05,
      bg.color = "white",
      bg.r = 0.05
    ) +
    scale_colour_manual(
      values = c("deeppink3", "seagreen3", "black", "chocolate3"),
      aesthetics = c("colour", "fill")
    ) +
    scale_x_continuous(breaks = c(-125, -126)) +
    scale_y_continuous(breaks = c(48.8, 49, 49.2, 49.4)) +
    annotation_north_arrow(
      location = "tr",
      style = ggspatial::north_arrow_nautical(),
      height = unit(3, "lines"),
      width = unit(3, "lines")
    ) +
    annotation_scale(location = "br") +
    coord_sf(
      expand = FALSE,
      xlim = coords_small[1:2],
      ylim = coords_small[3:4]
    ) +
    guides(fill = "none", colour = "none") +
    theme(
      axis.title = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "lightblue2"),
      panel.ontop = FALSE
    )
)


# Combined map adding the zoomed out map to the small map
(smu_map <- ggdraw(small_map) +
  draw_plot(
    big_map,
    x = 0.15,
    y = -0.25,
    width = 0.25
  )
)


# Save the map
ggsave(
  plot = smu_map,
  filename = here(
    "3. outputs",
    "Barkley_Sockeye_SMU_map.png"
  ),
  dpi = "print",
  scale = 1.8
)


  

# Build plot with large inset ---------------------------------------------


(pnw <- basemap(
  limits = c(-165, -125, 40, 62),
  rotate = TRUE,
  land.border.col = NA
) +
  geom_sf(
    data = st_as_sfc(bb_small),
    colour = "red",
    linewidth = 1,
    fill = NA
  ) +
  theme_minimal()
)


# Lakes to show on new map
lakes2 <- # Shapefile with lake polygons from Freshwater Atlas
  lakes <- read_sf(
    here(
      "1. data",
      "FWA_LAKES_POLY",
      "FWLKSPL_polygon.shp"
    )
  ) |> 
  rename_with(str_to_lower) |> 
  filter(
    str_detect(
      gnsnm1, 
      "(?i)hucuktlis|great\\scentral|sproat|dickson|ash|elsie|oshinow|nahmint"
      )
    ) |> 
  st_transform(crs = "NAD83")


# Stream lines for new map
stream_lines2 <- read_sf(
  here(
    "1. data",
    "FWA_STREAM_NETWORKS_SP",
    "FWSTRMNTWR_line.shp"
  )
) |> 
  rename_with(str_to_lower) |> 
  filter(
    str_detect(
      gnis_name, 
      "(?i)sproat|stamp|somass|hucuktlis|clemens|^ash\\s|nahmint"
      )
    )



# Barkley Sound Map
(bs <- ggplot(vi_coastline) +
    geom_sf(
      fill = NA,
      colour = "grey60"
    ) +
    geom_textsf(
      data = stream_lines2,
      aes(label = gnis_name),
      vjust = 1.5,
      hjust = 0.9,
      text_smoothing = 65,
      remove_long = FALSE,
      size = 3
    ) +
     geom_textsf(
       data = lakes2,
       aes(label = gnsnm1),
       #fill = NA,
       size = 3,
       vjust = 2
     ) +
    # scale_colour_manual(
    #   values = c("deeppink3", "seagreen3", "black", "chocolate3"),
    #   aesthetics = c("colour", "fill")
    # ) +
    annotation_north_arrow(
      location = "tr",
      style = ggspatial::north_arrow_nautical(),
      height = unit(3, "lines"),
      width = unit(3, "lines")
    ) +
    annotation_scale(location = "bl") +
    coord_sf(
      expand = FALSE,
      xlim = c(-125.7, -124.75),
      ylim = c(48.8, 49.48)
    ) +
    guides(fill = "none", colour = "none") +
    theme_void()
)
