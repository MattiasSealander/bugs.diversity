# ==============================================================
# Purpose: Filter and plot fossil insect sites in Europe
#          sampled from natural deposits and categorize
#          them by region.
#
# Methods:
#   1. Load fossil insect occurrence data.
#   2. Filter sites by region, context, and age constraints.
#   3. Import country vectors
#   4. Convert to sf objects
#   5. Plot sites on a map of Europe with UK/Ireland as inset map
#
# Key Notes:
#   - ggdraw is used to add UK/Ireland map as inset
#   - Harmonised sample age ranges are visualised using viridis
#
# Output:
#   - Plot: "001-sites-overview-britain.jpg"
# ==============================================================

# ---- 0. Load packages ----
pacman::p_load(cowplot, data.table, here, tidyverse, sf, ggspatial, spatstat, rnaturalearth, rnaturalearthdata, viridis
               )

# ==============================================================
# 1. Import species occurrence data
# ==============================================================
bugs <- fread(here("analysis", "data", "raw_data", "bugs_europe_extraction_samples_20250612.csv")) %>%
  mutate(sample_id = paste(sample, sample_group, site, sep = "@")) %>%
  as_tibble()

# ==============================================================
# 2. Prepare spatial data
# --------------------------------------------------------------
# Define bounding boxes for Europe and the British Isles.
# ==============================================================
# Switch off spherical geometry for spatial operations.
sf_use_s2(use_s2 = FALSE)

# Europe bounding box
bbox <- st_bbox(c(xmin = -25, ymin = 30, xmax = 40, ymax = 75), crs = 4326) %>%
  st_as_sfc()

# British Isles bounding box
gbbox <- st_bbox(c(xmin = -11.5, ymin = 49.8, xmax = 1.8, ymax = 61), crs = 4326) %>%
  st_as_sfc()

# ==============================================================
# 3. Filter and categorize sites
# --------------------------------------------------------------
# Filter to natural deposits (Stratigraphic sequences) according to selection criteria
# Convert to sf objects and categorize by region
# ==============================================================
# Filter data
sites <- bugs %>%
  mutate(age_range = age_older - age_younger) %>%
  filter(
    !family_name %in% c("ACANTHOSOMATIDAE","ANTHOCORIDAE", "APHALARIDAE", "APHIDOIDEA", "AUCHENORRHYNCHA", "BIBIONIDAE",
                        "BRACHYCENTRIDAE", "CALOPTERYGIDAE", "CERCOPIDAE", "CHIRONOMIDAE", "CICADELLIDAE", "CICADOMORPHA",
                        "CIXIIDAE", "CORIXIDAE", "CYCLORRHAPHA", "CYDNIDAE", "DELPHACIDAE", "DERMAPTERA", "DIPTERA", "FORFICULIDAE",
                        "FORMICIDAE", "FULGOROMORPHA", "GERRIDAE", "GLOSSOSOMATIDAE", "GOERIDAE", "HEBRIDAE", "HEMIPTERA",
                        "HETEROPTERA", "HOMOPTERA", "HYDROPSYCHIDAE", "HYDROPTILIDAE", "HYMENOPTERA", "LEPIDOPTERA", "LEPIDOSTOMATIDAE",
                        "LEPTOCERIDAE", "LIMNEPHILIDAE", "LYGAEIDAE", "MEMBRACIDAE", "MICROPHYSIDAE", "MICROSPORIDAE", "MIRIDAE",
                        "MOLANNIDAE", "NABIDAE", "NEMATOCERA", "NEMOURIDAE", "ODONATA", "PARASITICA", "PENTATOMIDAE", "PHRYGANEIDAE",
                        "POLYCENTROPIDAE", "PSYCHOMYIIDAE", "PSYLLIDAE", "RAPHIDIIDAE", "SALDIDAE", "SCUTELLERIDAE", "SERICOSTOMATIDAE",
                        "SIALIDAE", "TRICHOPTERA", "THYREOCORIDAE", "TINGIDAE", "TIPULIDAE", "TRIOPSIDAE", "TRIOZIDAE"),
    context == "Stratigraphic sequence",
    sample != "BugsPresence",
    age_range <= 2000,
    country != "Greenland"
  ) %>%
  mutate(mid_age = (age_older + age_younger) / 2) %>%
  filter(between(mid_age, -500, 16000)) %>%
  distinct(site, age_older, age_younger, latitude, longitude) %>%
  mutate(region = case_when(
    between(latitude, 49.8, 62.6) & between(longitude, -12.6, 1.8) ~ "British/Irish Isles",
    TRUE ~ ""
  ))

# Split sites by region and convert to sf objects
sites.sf.list <- split(sites, sites$region == "British/Irish Isles")
sites.gb.sf <- st_as_sf(sites.sf.list[["TRUE"]], coords = c("longitude", "latitude"), crs = 4326)
sites.sf <- st_as_sf(sites.sf.list[["FALSE"]], coords = c("longitude", "latitude"), crs = 4326)

# ==============================================================
# 4. Prepare country polygons
# --------------------------------------------------------------
# Fetch country polygons for Europe and British Isles and crop to bounding boxes.
# ==============================================================
# Europe
countries.sf <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(crs = 4326) %>%
  st_crop(bbox)
# British Isles
countries.gb.sf <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(name %in% c("Ireland", "United Kingdom")) %>%
  st_transform(crs = 4326) %>%
  st_crop(gbbox)

# Switch spherical geometry back on for plotting
sf_use_s2(TRUE)

# ==============================================================
# 5. Prepare data for plotting
# --------------------------------------------------------------
# Combine and plot data, using an inset map for the UK
# ==============================================================
# Combine sites into one sf object
sites.all.sf <- bind_rows(sites.sf, sites.gb.sf) %>%
  mutate(point_size = ifelse(region == "British/Irish Isles", 3, 5))

# 5a. Plot European sites
plot.1 <- ggplot() +
  geom_sf(data = countries.sf, fill = "grey70", linewidth = 0.8, col = "black", alpha = 0.2) +
  geom_sf(data = st_as_sf(bbox), linewidth = 0.8, col = "grey", fill = NA) +
  geom_sf(data = st_as_sf(gbbox), linewidth = 0.6, col = "black", fill = "pink2", alpha = .4) +
  geom_sf(data = sites.all.sf, aes(size = point_size, fill = age_younger),
          shape = 21, color = "black", alpha = 0.9, show.legend = c(fill = TRUE)) +
  scale_size_identity() +
  scale_fill_viridis(name = "Time (Years BP)", direction = -1) +
  coord_sf(datum = st_crs(4326), label_graticule = "NW", clip = "on") +
  theme_minimal() +
  theme(legend.position = "none")

# 5b. Plot British Isles (inset map)
plot.2 <- ggplot() +
  geom_sf(data = countries.gb.sf, fill = "grey70", linewidth = 0.8, col = "black", alpha = 0.2) +
  geom_sf(data = sites.gb.sf, aes(size = 6, fill = age_younger),
          shape = 21, color = "white", alpha = 0.8, lwd = 3) +
  scale_size_identity() +
  scale_fill_viridis(name = "Time (Years BP)", direction = -1) +
  annotation_north_arrow(location = "tl", which_north = "true",
                                    height = unit(1.5, "cm"),
                                    width = unit(1.5, "cm")) +
  coord_sf(datum = NA, clip = "on", expand = FALSE) +
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        plot.background = element_blank())

# 5c. Create legend
legend_plot <- ggplot() +
  geom_sf(data = sites.gb.sf, aes(fill = age_younger),
          shape = 21, color = "white", alpha = 0.8, lwd = 3) +
  scale_fill_viridis(name = "Time (Years BP)", direction = -1, breaks=c(2000, 6000, 10000, 14000)) +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.key.width= unit(1.5, 'cm'))
legend <- get_legend(legend_plot)

# 5d. Combine main map and inset
fig <- ggdraw(plot.1) +
  draw_plot(plot.2, x = 0.6, y = 0.4, width = 0.46, height = 0.46)

# 5e. Add legend
fig_with_legend <- plot_grid(fig, legend, ncol = 1, rel_heights = c(1, 0.05))

# ==============================================================
# 6. Save figure
# --------------------------------------------------------------
# Save the combined figure to the figures directory.
# ==============================================================
ggsave(
  "001-sites-overview.jpg",
  fig_with_legend,
  device = "jpg",
  path = here("analysis", "figures"),
  width = 5900,
  height = 5900,
  units = "px",
  dpi = 500
)

message("✅ Sites have been plotted and figure saved: '001-sites-overview.jpg'")
