pacman::p_load(readr, tidyverse, sf, ggrepel, ggplot2, spatstat, rnaturalearth, rnaturalearthdata)
# Import descriptive metadata
bugs.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "bugs_europe_extraction_samples_20250612.csv"), sep = ";", dec = ".", header = TRUE, na = c("", "NA", "NULL"), encoding = "UTF-8")

# Spherical geometry (s2) switched off
sf_use_s2(use_s2 = FALSE)

# Bounding box for study area
bbox <- st_bbox(
  c(xmin= -25, ymin= 30, xmax= 40, ymax= 75), crs = 4326) %>%
  st_as_sfc()

# Filter sites sampled in natural contexts between modern day and 16 000 BP, with max age range of 2 000 years
sites <-
  bugs.csv %>%
  mutate(age_range=age_older - age_younger) %>%
  filter(context == "Stratigraphic sequence", sample != "BugsPresence", between(age_older, -500, 16000) & between(age_younger, -500, 16000),
         age_range <= 2000, !country == "Greenland") %>%
  select(site, age_older, age_younger, latitude, longitude) %>%
  distinct() %>%
  mutate(place = case_when(
    between(latitude, 53.5, 74.8) & between(longitude, 2.8, 41) ~ "Scandinavia/Baltic",
    between(latitude, 49.8, 62.6) & between(longitude, -12.6, 1.8) ~ "British/Irish Isles",
    between(latitude, 63, 66.3) & between(longitude, -25, -12) ~ "Iceland",
    latitude < 45.5 ~ "Meridional")) %>%
  mutate(region = ifelse(is.na(place),
                         "Continental",
                         place)) %>%
  select(-place)

# Store all British/Irish sites
sites.gb.sf <-
  sites %>%
  filter(region == "British/Irish Isles") %>%
  st_as_sf(coords = c("longitude", "latitude"), crs=4326) #as geographic (sf) data with WGS 84 projection

# Store remaining sites, with region
sites.sf <-
  sites %>%
  filter(!region == "British/Irish Isles") %>%
  st_as_sf(coords = c("longitude", "latitude"), crs=4326) #as geographic (sf) data with WGS 84 projection

# Fetch polygons for countries from Natural Earth R data, crop to study area bounding box
countries.sf <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(crs = 4326) %>%
  st_crop(bbox)

# Spherical geometry (s2) switched on
sf_use_s2(use_s2 = TRUE)

# Setup figure
fig <-
  ggplot() +
  geom_sf(data = st_boundary(sites.sf)) +
  geom_sf(data = countries.sf, fill = "grey70", linewidth=.8, col= "black", alpha = .2) +
  geom_sf(data = bbox, linewidth=.8, col="grey", fill = NA) +
  geom_sf(data = sites.sf, aes(fill = region), shape = 21, size = 3, color = "black", alpha = .9, show.legend = 'point') +
  geom_sf(data = sites.gb.sf, aes(fill = "British/Irish Isles"), shape = 21, size = 2, color = "black", alpha = .8, show.legend = FALSE) +
  scale_fill_manual(values = c("Scandinavia/Baltic" = "#E69F00",  "British/Irish Isles" = "#CC79A7", "Iceland" = "#0073C2FF", "Continental" = "#009E73", "Meridional" = "#F0E442"),
                    name = "") +
  coord_sf(datum=st_crs(4326),
           label_graticule = "NW",
           clip = "on") +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 14, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 14, face = "bold", colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        legend.position = "bottom",
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0)) # Left margin

# Save figure
ggsave("001-sites-overview.jpg",
       fig,
       device = "jpg",
       here::here("analysis", "figures"),
       width=2400,
       height=2400,
       units = "px",
       dpi = 300)
