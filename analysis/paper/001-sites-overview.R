############################################################
# Purpose: Filter and plot fossil insect sites in Europe
#          sampled from natural deposits and categorize
#          them by region.
#
# Methods:
#   1. Filter sites from natural deposits and time range
#   2. Convert to sf objects
#   3. Plot sites on a map of Europe
############################################################

# ---- Load packages ----
# pacman::p_load() loads the packages and installs them if missing
pacman::p_load(data.table, tidyverse, sf, ggrepel, spatstat, rnaturalearth, rnaturalearthdata)

# ---- 1. Import site and species data ----
bugs <-
  fread(here::here("analysis", "data", "raw_data", "bugs_europe_extraction_samples_20250612.csv"),
        na.strings = c("", "NA", "NULL"),
        encoding = "UTF-8")

# ---- 2. Prepare spatial data ----
# Switch off spherical geometry
sf_use_s2(use_s2 = FALSE)

# Define study area bounding box
bbox <- st_bbox(
  c(xmin= -25, ymin= 30, xmax= 40, ymax= 75), crs = 4326) %>%
  st_as_sfc()

# ---- 3. Filter and categorize sites ----
# Filter for stratigraphic sequences (natural deposits), exclude non-presence samples, age limits, and Greenland
sites <- bugs %>%
  filter(context == "Stratigraphic sequence",
         sample != "BugsPresence",
         between(age_older, -500, 16000),
         between(age_younger, -500, 16000),
         (age_older - age_younger) <= 2000,
         country != "Greenland") %>%
  distinct(site, age_older, age_younger, latitude, longitude) %>%
  # Assign sites to regions based on latitude/longitude
  mutate(region = case_when(
    between(latitude, 53.5, 74.8) & between(longitude, 2.8, 41) ~ "Scandinavia/Baltic",
    between(latitude, 49.8, 62.6) & between(longitude, -12.6, 1.8) ~ "British/Irish Isles",
    between(latitude, 63, 66.3) & between(longitude, -25, -12) ~ "Iceland",
    latitude < 45.5 ~ "Meridional",
    TRUE ~ "Continental"
  ))

# ---- 4. Convert sites to sf objects ----
# Split sites into British/Irish Isles vs other regions
sites.sf.list <-
  split(sites, sites$region == "British/Irish Isles")

# Convert to sf objects
sites.gb.sf <-
  st_as_sf(sites.sf.list[["TRUE"]], coords = c("longitude", "latitude"), crs = 4326)
sites.sf <-
  st_as_sf(sites.sf.list[["FALSE"]], coords = c("longitude", "latitude"), crs = 4326)

# ---- 5. Fetch country polygons and crop ----
countries.sf <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(crs = 4326) %>% # Ensure same CRS as sites
  st_crop(bbox)

# Switch spherical geometry back on for plotting
sf_use_s2(TRUE)

# ---- 6. Prepare data for plotting ----
# Combine all sites into one sf object and assign point sizes
# Smaller size for dense British/Irish Isles points
# Add a size column based on region
sites.all.sf <- bind_rows(sites.sf, sites.gb.sf) %>%
  mutate(point_size = ifelse(region == "British/Irish Isles", 2, 3))

# ---- 7. Plot the sites ----
fig <- ggplot() +
  # Draw boundaries of non-GB sites
  geom_sf(data = st_boundary(sites.sf), color = "black") +
  # Plot countries
  geom_sf(data = countries.sf, fill = "grey70", linewidth = 0.8, col = "black", alpha = 0.2) +
  # Plot study area bounding box
  geom_sf(data = st_as_sf(bbox), linewidth = 0.8, col = "grey", fill = NA) +
  # Plot all sites with fill by region and variable point size
  geom_sf(
    data = sites.all.sf,
    aes(fill = region, size = point_size),
    shape = 21,
    color = "black",
    alpha = 0.9,
    show.legend = c(fill = TRUE)
  ) +
  # Assign colors to regions and enlarge legend points
  scale_fill_manual(
    values = c(
      "Scandinavia/Baltic" = "#E69F00",
      "British/Irish Isles" = "#CC79A7",
      "Iceland" = "#0073C2FF",
      "Continental" = "#009E73",
      "Meridional" = "#F0E442"
    ),
    name = "",
    guide = guide_legend(
      override.aes = list(size = 4)
    )
  ) +
  scale_size_identity() + # Use actual point sizes for the plot
  coord_sf(datum = st_crs(4326), label_graticule = "NW", clip = "on") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 14, face = "bold", colour = "black"),
    axis.title.y = element_text(size = 14, face = "bold", colour = "black"),
    legend.text = element_text(size = 8, colour = "black"),
    legend.position = "bottom",
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
  )

# ---- 8. Save figure ----
ggsave(
  "001-sites-overview.jpg",
  fig,
  device = "jpg",
  path = here("analysis", "figures"),
  width=2400,
  height=2400,
  units = "px",
  dpi = 300
)

message("✅ Sites have been plotted and figure saved: '001-sites-overview.jpg'")
