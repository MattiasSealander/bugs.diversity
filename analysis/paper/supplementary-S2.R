library(tidyverse)
library(IRanges)
library(data.table)
library(here)

# ---- 1. Load and filter samples ----
samples <- fread(
  here::here("analysis/data/raw_data/bugs_europe_extraction_samples_20250612.csv"),
  na.strings = c("", "NA", "NULL"),
  encoding = "UTF-8"
) %>%
  filter(
    context == "Stratigraphic sequence",
    sample != "BugsPresence",
    between(age_older, -500, 16000),
    between(age_younger, -500, 16000),
    (age_older - age_younger) <= 2000,
    between(latitude, 49.8, 62.6),
    between(longitude, -12.6, 1.8)
  ) %>%
  mutate(sample_id = paste(sample, sample_group, site, sep = "@")) %>%
  distinct(sample_id, sample, sample_group, site, age_younger, age_older)

# ---- 2. Define 500-year bins ----
bins <- tibble(
  start = c(-500, seq(1, 15501, by = 500)),
  end   = seq(0, 16000, by = 500)
)

# ---- 3. Assign samples to overlapping bins ----
sample_ranges <- IRanges(start = samples$age_younger, end = samples$age_older)
bin_ranges <- IRanges(start = bins$start, end = bins$end)

overlaps <- findOverlaps(sample_ranges, bin_ranges, type = "any")

# ---- 4. Build a table of sample -> bin assignments ----
sample_bin_map <- tibble(
  sample_id = samples$sample_id[queryHits(overlaps)],
  bin_id    = subjectHits(overlaps)
)

# ---- 5. Count number of unique bins per sample ----
sample_bin_counts <- sample_bin_map %>%
  group_by(sample_id) %>%
  summarise(n_bins = n_distinct(bin_id), .groups = "drop")

# ---- 6. Summarise frequency of samples per number of bins ----
bin_distribution <- sample_bin_counts %>%
  count(n_bins, name = "n_samples")

# ---- 7. Plot with text labels ----
fig <-
  ggplot(bin_distribution, aes(x = factor(n_bins), y = n_samples)) +
  geom_bar(stat = "identity", color = "black", width = 0.7, fill = "green4") +
  geom_text(
    aes(label = n_samples),
    vjust = -0.4,
    size = 4,
    fontface = "bold"
  ) +
  labs(
    x = "Number of 500-year bins intersected by sample age range",
    y = "Number of samples",
    title = "Distribution of samples by number of overlapping 500-year bins"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )


# ==============================================================
# 8. Save figure
# ==============================================================
ggsave(
  filename = "supplementary-S2.jpg",
  plot = fig,
  path = here("analysis", "figures"),
  units = "cm",
  width = 20,
  height = 20,
  dpi = 300
)

message("✅ Data has been summarized and figure saved: 'supplementary-S2jpg'")
