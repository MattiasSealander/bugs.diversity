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
pacman::p_load(cowplot, data.table, here, tidyverse, sf, ggrepel, spatstat, rnaturalearth, rnaturalearthdata)


# ---- 1. Import site and species data ----
bugs <-
  fread(here::here("analysis", "data", "raw_data", "bugs_europe_extraction_samples_20250612.csv"),
        na.strings = c("", "NA", "NULL"),
        encoding = "UTF-8")


# ---- 2. Prepare spatial data ----
# Switch off spherical geometry
sf_use_s2(use_s2 = FALSE)

# Define Europe bounding box
bbox <- st_bbox(
  c(xmin= -25, ymin= 30, xmax= 40, ymax= 75), crs = 4326) %>%
  st_as_sfc()

# Define study area bounding box (United Kingdom)
gbbox <-
  st_bbox(
    c(xmin= -11.5, ymin= 49.8, xmax= 1.8, ymax= 61), crs = 4326) %>%
  st_as_sfc()


# ---- 3. Filter and categorize sites ----
# Filter to Stratigraphic sequences (natural deposits) within age limits, exclude temporary samples (BugsPresence) and Greenland
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
    between(latitude, 49.8, 62.6) & between(longitude, -12.6, 1.8) ~ "British/Irish Isles",
    TRUE ~ ""
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


# ---- 5a. Fetch Europe country polygons and crop ----
countries.sf <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(crs = 4326) %>% # Ensure same CRS as sites
  st_crop(bbox)

# ---- 5b. Fetch country polygons for the study area and crop ----
countries.gb.sf <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(name %in% c("Ireland", "United Kingdom")) %>%
  st_transform(crs = 4326) %>% # Ensure same CRS as sites
  st_crop(gbbox)

# Switch spherical geometry back on for plotting
sf_use_s2(TRUE)


# ---- 6. Prepare data for plotting ----
# Combine all sites into one sf object and assign point sizes
# Smaller size for dense British/Irish Isles points
# Add a size column based on region
sites.all.sf <- bind_rows(sites.sf, sites.gb.sf) %>%
  mutate(point_size = ifelse(region == "British/Irish Isles", 3, 5))


# ---- 7. Plot the figure ----

# 7a. Plot the European sites
plot.1 <- ggplot() +
  # Draw boundaries of non-GB sites
  geom_sf(data = st_boundary(sites.sf), color = "black") +
  # Plot countries
  geom_sf(data = countries.sf, fill = "grey70", linewidth = 0.8, col = "black", alpha = 0.2) +
  # Plot study area bounding box
  geom_sf(data = st_as_sf(bbox), linewidth = 0.8, col = "grey", fill = NA) +
  # Plot study area bounding box
  geom_sf(data = st_as_sf(gbbox), linewidth = 0.6, col = "black", fill = "pink2", alpha = .4) +
  # Plot all sites with fill by region and variable point size
  geom_sf(
    data = sites.all.sf,
    aes(size = point_size, fill = age_younger),
    shape = 21,
    color = "white",
    alpha = 0.9,
    show.legend = c(fill = TRUE)
  ) +
  scale_size_identity() + # Use actual point sizes for the plot
  scale_fill_viridis(
    name = "Time (Years BP)",
    direction = -1
  ) +
  coord_sf(datum = st_crs(4326), label_graticule = "NW", clip = "on") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 14, face = "bold", colour = "black"),
    axis.title.y = element_text(size = 14, face = "bold", colour = "black"),
    legend.text = element_text(size = 8, colour = "black"),
    legend.position = "none",
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
  )

# 7b. Plot the study area sites (inset map)
plot.2 <- ggplot() +
  # Plot countries
  geom_sf(
    data = countries.gb.sf,
    fill = "grey70",
    linewidth = 0.8,
    col = "black",
    alpha = 0.2
  ) +
  # Plot all sites with fill by region and variable point size
  geom_sf(
    data = sites.gb.sf,
    aes(size = 6, fill = age_younger),
    shape = 21,
    color = "white",
    alpha = 0.8,
    lwd = 3
  ) +
  scale_size_identity() + # Use actual point sizes for the plot
  scale_fill_viridis(
    name = "Time (Years BP)",
    direction = -1
  ) +
  coord_sf(datum = NA, clip = "on", expand = FALSE) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
    panel.grid.major = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
    legend.position = "none",
  )

# ---- 7c. Generate legend plot ----
legend_plot <- ggplot() +
  geom_sf(
    data = sites.gb.sf,
    aes(fill = age_younger),
    shape = 21,
    color = "white",
    alpha = 0.8,
    lwd = 3
  ) +
  scale_fill_viridis(
    name = "Time (Years BP)",
    direction = -1
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.width = unit(0.10, 'npc')
  )

# Extract the legend
legend <- get_legend(legend_plot)

# 7d. Draw the study area as an inset map on the European distribution map
fig <-
  ggdraw(plot.1) +
  draw_plot({plot.2},
    x = 0.6,
    y = 0.4,
    width = 0.46,
    height = 0.46)

# Add legend at the bottom
fig_with_legend <- plot_grid(
  fig,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.05)
)


# ---- 8. Save figure ----
ggsave(
  "001-sites-overview-britain.jpg",
  fig_with_legend,
  device = "jpg",
  path = here("analysis", "figures"),
  width=40,
  height=40,
  units = "cm",
  dpi = 300
)

message("✅ Sites have been plotted and figure saved: '001-sites-overview-britain.jpg'")
