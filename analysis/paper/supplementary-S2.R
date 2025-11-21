# ---- 0. Load required packages ----
pacman::p_load(
  grid, gridExtra, data.table, tidyverse, IRanges, here, magick
)

# ==============================================================
# 1. Import species occurrence data
# ==============================================================
bugs <- fread(
  here::here("analysis", "data", "raw_data", "bugs_europe_extraction_samples_20250612.csv"),
  na.strings = c("", "NA", "NULL"),
  encoding = "UTF-8"
)

# ==============================================================
# 2. Filter data
# ==============================================================
# Unique samples which has a defined age range within criteria
samples <- bugs %>%
  select(country, latitude, longitude, sample, site, sample_group, age_older, age_younger, context) %>%
  mutate(age_range = age_older - age_younger,
         region = case_when(
           between(latitude, 49.8, 62.6) & between(longitude, -12.6, 1.8) ~ "British/Irish Isles",
           TRUE ~ "")
  ) %>%
  filter(
    region == "British/Irish Isles",
    context == "Stratigraphic sequence",
    sample != "BugsPresence",
    age_range <= 2000,
    country != "Greenland"
  ) %>%
  mutate(mid_age = (age_older + age_younger) / 2) %>%
  filter(between(mid_age, -500, 16000)) %>%
  mutate(sample_id = paste(sample, sample_group, site, sep = "@")) %>%
  distinct(sample_id, sample, sample_group, site, age_younger, age_older)

# ==============================================================
# 3. Define and assign samples to 500-year bins
# ==============================================================
# Define time bins
bins <- tibble(
  start = c(-500, seq(1, 15501, by = 500)),
  end   = seq(0, 16000, by = 500)
)

# Assign samples to time bins
sample_ranges <- IRanges(start = samples$age_younger, end = samples$age_older)
bin_ranges <- IRanges(start = bins$start, end = bins$end)

overlaps <- findOverlaps(sample_ranges, bin_ranges, type = "any")


# ==============================================================
# 4. Build results table
# ==============================================================
# Matrix with sample-bin assignments
sample_bin_map <- tibble(
  sample_id = samples$sample_id[queryHits(overlaps)],
  bin_id    = subjectHits(overlaps)
)

# ==============================================================
# 5. Summarise data
# ==============================================================
# Count number of unique bins per sample
sample_bin_counts <- sample_bin_map %>%
  group_by(sample_id) %>%
  summarise(n_bins = n_distinct(bin_id), .groups = "drop")

# Summarise frequency of samples per number of bins
bin_distribution_table <- sample_bin_counts %>%
  count(n_bins, name = "Samples") %>%
  dplyr::rename("Number of Bins" = n_bins) %>%
  arrange(1)

# ==============================================================
# 6. Create table
# ==============================================================
# Black-and-white table theme
bw_theme <- ttheme_default(
  core = list(
    fg_params = list(cex = 1.2, fontfamily = "sans"),   # normal text
    bg_params = list(fill = "white", col = "black")     # white background, black borders
  ),
  colhead = list(
    fg_params = list(cex = 1.3, fontface = "bold", fontfamily = "sans"), # bold headers
    bg_params = list(fill = "white", col = "black")     # same style for header
  ),
  padding = unit(c(4, 4), "mm") # some padding for readability
)

# Create table graphic
table_plot <- tableGrob(
  bin_distribution_table,
  rows = NULL,
  theme = bw_theme
)

# Render table to PNG
png(here("analysis", "figures", "supplementary-S2.png"),
    units = "cm",
    width = 20,
    height = 20,
    res = 300)
# Display table and save
grid.newpage()
grid.draw(table_plot)
dev.off()

# Read image in order to crop white space
img <- image_read(here("analysis", "figures", "supplementary-S2.png"))

# Trim white space using magick package
img_trimmed <- image_trim(img)

# Add a small border back
img_final <- image_border(img_trimmed, color = "white", geometry = "5x5")

# Save cropped image
image_write(img_final, path = here("analysis", "figures", "supplementary-S2.png"), format = "png")

message("✅ Data has been summarized and table saved: 'supplementary-S2jpg'")
