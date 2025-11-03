############################################################
# Plot British/Irish Isles fossil insect sites per 500-yr bin
# Samples are assigned to every temporal bin they intersect
############################################################

pacman::p_load(tidyverse, data.table, IRanges, sf, rnaturalearth, rnaturalearthdata, cowplot, here)

# ---- 1. Load fossil insect data ----
bugs <- fread(
  here::here("analysis" "data" "raw_data" "bugs_europe_extraction_samples_20250612.csv"),
  na.strings = c("", "NA", "NULL"),
  encoding = "UTF-8"
)

# ---- 2. Filter to British/Irish Isles and natural context ----
samples <- bugs %>%
  filter(
    context == "Stratigraphic sequence",
    sample != "BugsPresence",
    between(age_older, -500, 16000),
    between(age_younger, -500, 16000),
    (age_older - age_younger) <= 2000,
    between(latitude, 49.8, 62.6),
    between(longitude, -12.6, 1.8)
  ) %>%
  mutate(sample_id = paste(sample, sample_group, site, sep = "|")) %>%
  distinct(sample, sample_group, site, latitude, longitude, age_older, age_younger, sample_id)

# ---- 3. Define 500-year temporal bins ----
bins <- tibble(
  start = c(-500, seq(1, 15501, by = 500)),
  end   = seq(0, 16000, by = 500)
)

# ---- 4. Assign samples to all overlapping bins ----
sample_ranges <- IRanges(start = samples$age_younger, end = samples$age_older)
bin_ranges    <- IRanges(start = bins$start, end = bins$end)

overlaps <- findOverlaps(sample_ranges, bin_ranges, type = "any")

sites_bin_map <- tibble(
  site      = samples$site[queryHits(overlaps)],
  latitude  = samples$latitude[queryHits(overlaps)],
  longitude = samples$longitude[queryHits(overlaps)],
  bin_end   = bins$end[subjectHits(overlaps)]
) %>%
  distinct()


# ---- 5. Convert to sf object ----
sites_sf <- st_as_sf(
  sites_bin_map,
  coords = c("longitude", "latitude"),
  crs = 4326
)

# ---- 6. Fetch UK/Ireland polygons ----
uk_sf <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(name %in% c("United Kingdom", "Ireland")) %>%
  st_transform(crs = 4326)

# ---- 7. Plot sites facetted by bin ----
fig <-
  ggplot() +
  geom_sf(data = uk_sf, fill = "grey80", color = "black") +
  geom_sf(data = sites_sf, shape = 21, color = "white", fill = "#1f78b4", size = 2, alpha = 0.8) +
  facet_wrap(~bin_end, nrow = 4) +
  coord_sf(xlim = c(-12.6, 1.8), ylim = c(49.8, 62.6), expand = FALSE) +
  labs(title = "Fossil insect sites in the British/Irish Isles (per 500-year bin)") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

# ==============================================================
# 8. Save figure
# ==============================================================
ggsave(
  filename = "supplementary-S1.jpg",
  plot = fig,
  path = here("analysis", "figures"),
  units = "cm",
  width = 30,
  height = 30,
  dpi = 300
)

message("✅ Data has been summarized and figure saved: 'supplementary-S1.jpg'")
