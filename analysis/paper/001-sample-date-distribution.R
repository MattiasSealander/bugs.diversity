############################################################
# Purpose: Summarize fossil insect data for Europe during
#          the Late Glacial - Holocene period.
#
# Methods:
#   1. Filter fossil insect data to period and region of study
#   2. Summarize sites, samples and total abundance of insects
#   3. Plot summarized data over time
############################################################

# ---- Load packages ----
pacman::p_load(tidyverse, IRanges, cowplot, here, data.table)

# ---- 1. Import site and species data ----
bugs <- fread(
  here::here("analysis/data/raw_data/bugs_europe_extraction_samples_20250612.csv"),
  na.strings = c("", "NA", "NULL"),
  encoding = "UTF-8"
)

# ---- 2. Filter and prepare data ----
# Filter for stratigraphic sequences, valid age ranges, and natural contexts
# Deduplicate by sample, sample_group, and site to ensure unique entries
# Compute start and end ages for each sample
nat <- bugs %>%
  filter(
    context == "Stratigraphic sequence",
    sample != "BugsPresence",
    between(age_older, -500, 16000),
    between(age_younger, -500, 16000),
    (age_older - age_younger) <= 2000,
    country != "Greenland"
  ) %>%
  distinct(sample, sample_group, site, age_older, age_younger) %>%
  mutate(
    start = age_younger,
    end   = age_older
  ) %>%
  select(sample, sample_group, site, start, end)

# ---- 3. Prepare 500-year age bins ----
age_bins <- tibble(
  start = c(-500, seq(1, 15501, by = 500)),
  end   = seq(0, 16000, by = 500)
)

# ---- 4. Find overlaps using IRanges ----
nat_ranges <- IRanges(start = nat$start, end = nat$end)
bin_ranges <- IRanges(start = age_bins$start, end = age_bins$end)

# Find overlaps between sample age ranges and bins
overlaps <- findOverlaps(nat_ranges, bin_ranges, type = "any")

# ---- 5. Build overlap hits tibble ----
# Each row corresponds to a sample overlapping a specific age bin
overlap_hits <- tibble(
  sample  = nat$sample[queryHits(overlaps)],
  sample_group = nat$sample_group[queryHits(overlaps)],
  site    = nat$site[queryHits(overlaps)],
  bin_start = age_bins$start[subjectHits(overlaps)],
  bin_end   = age_bins$end[subjectHits(overlaps)]
)

# ---- 6. Summarise hits ----
# 6a. Number of unique samples per bin
hits.samples <- overlap_hits %>%
  group_by(bin_end) %>%
  summarise(samples = n_distinct(interaction(sample, sample_group, site)), .groups = "drop")

# 6b. Number of unique sites per bin
# Deduplicate to ensure each site is counted only once per bin
hits.sites <- overlap_hits %>%
  distinct(bin_end, site) %>%
  group_by(bin_end) %>%
  summarise(sites = n(), .groups = "drop")

# 6c. Total abundance per bin
# Join back to original abundance data and sum per bin
hits.abund <- overlap_hits %>%
  left_join(bugs, by = c("sample", "sample_group", "site")) %>%
  group_by(bin_end) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop")

# ---- 7. Plotting function ----
# Generalized function to create bar plots for samples, sites, or abundance
# 'show_x' controls whether to display the x-axis label (used to avoid repetition)
plot_bar <- function(data, y_var, y_label, show_x = TRUE) {
  ggplot(data, aes(x = bin_end, y = !!sym(y_var))) +
    geom_bar(stat = "identity", fill = "green4", color = "black", alpha = 0.8, just = 0) +
    scale_x_reverse(limits = c(16000, -500), breaks = scales::pretty_breaks(n = 10)) +
    labs(x = if(show_x) "Age bin (BP)" else "", y = y_label) +
    theme_minimal() +
    theme(
      axis.title.y = element_text(size = 10, face = "bold", colour = "black"),
      axis.text.y = element_text(size = 10, colour = "black"),
      axis.title.x = element_text(size = 12, face = "bold", colour = "black")
    )
}

# ---- 8. Generate plots ----
fig.sites <- plot_bar(hits.sites, "sites", "Sites", show_x = FALSE)
fig.samples <- plot_bar(hits.samples, "samples", "Samples", show_x = FALSE)
fig.abundance <- plot_bar(hits.abund, "abundance", "Abundance", show_x = TRUE)

# ---- 9. Combine plots ----
fig <- cowplot::plot_grid(
  fig.sites, fig.samples, fig.abundance,
  align = "v",
  labels = c('A', 'B', 'C'),
  ncol = 1,
  label_size = 12
)

# ---- 10. Save figure ----
ggsave(
  "001-sample-date-distribution.jpg",
  fig,
  device = "jpg",
  path = here("analysis", "figures"),
  width = 8,
  height = 8,
  dpi = 300
)

message("✅ Data has been summarized and figure saved: '001-sample-date-distribution.jpg'")
