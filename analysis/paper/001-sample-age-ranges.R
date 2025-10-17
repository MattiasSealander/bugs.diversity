############################################################
# Purpose: Summarize sample age ranges for fossil insect data from
#          the United Kingdom and Ireland during the Late Glacial - Holocene period.
#
# Methods:
#   1. Filter fossil insect data to period and region of study
#   2. Order samples in ascending order of age
#   3. Plot data
############################################################

# ---- 0. Load packages ----
pacman::p_load(tidyverse, IRanges, cowplot, here, data.table)


# ==============================================================
# 1. Import species occurrence data
# --------------------------------------------------------------
# The raw dataset includes sample-level fossil insect occurrences across Europe from SEAD.
# ==============================================================
bugs <- fread(
  here::here("analysis/data/raw_data/bugs_europe_extraction_samples_20250612.csv"),
  na.strings = c("", "NA", "NULL"),
  encoding = "UTF-8"
)


# ==============================================================
# 2. Filter and prepare data
# --------------------------------------------------------------
# Filter for samples from natural contexts (Stratigraphic sequences) and valid age ranges, from the United Kingdom and Ireland.
# Samples are identified by concatenation of site, sample_group, and sample.
# Add row numbers for figure axis labels.
# ==============================================================
data <- bugs %>%
  mutate(
    region = case_when(
      between(latitude, 49.8, 62.6) & between(longitude, -12.6, 1.8) ~ "British/Irish Isles")) %>%
  filter(
    region == "British/Irish Isles",
    context == "Stratigraphic sequence",
    sample != "BugsPresence",
    between(age_older, -500, 16000),
    between(age_younger, -500, 16000),
    (age_older - age_younger) <= 2000,
    country != "Greenland"
  ) %>%
  distinct(sample, sample_group, site, age_older, age_younger) %>%
  mutate(sample_uid = interaction(site, sample_group, sample, sep = "_")) %>%
  arrange((age_younger)) %>%  # order oldest to youngest
  mutate(sample_row = row_number())  # row number for y-axis


# ==============================================================
# 3. Plot figure
# ==============================================================
fig <-
  ggplot(data, aes(y = sample_row)) +
  geom_segment(aes(x = age_older, xend = age_younger, y = sample_row, yend = sample_row, color = age_older),
               size = 3) +
  scale_x_reverse(limits = c(16000, -500), breaks = scales::pretty_breaks(n = 5)) +
  scale_color_gradient(low = "lightgreen", high = "darkgreen") +
  labs(
    x = "Time (Years BP)",
    y = "Sample"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none"
  )


# ==============================================================
# 4. Save figure
# ==============================================================
ggsave(
  "001-sample-age-ranges.jpg",
  fig,
  device = "jpg",
  path = here("analysis", "figures"),
  units = "cm",
  width = 10,
  height = 20,
  dpi = 300
)

message("✅ Data has been summarized and figure saved: '001-sample-age-ranges.jpg'")
