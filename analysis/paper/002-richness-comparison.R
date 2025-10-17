############################################################
# Purpose:
#   Compare three richness estimates by 500-year time bins:
#   - Raw richness (unadjusted)
#   - Rarefied richness (sample-size standardized, Cmin = 10)
#   - SRS standardized richness (Scaling with Ranked Subsampling, Cmin = 10)
#
# Description:
#   The script processes subfossil beetle data from Europe,
#   binning samples into 500-year intervals between 16,000 and -500 BP.
#   Each bin’s mean taxon richness and standard error are calculated.
#
#   Temporal bins are visualized following the identified periods
#   of biotic transition defined by Pilotto et al. (2022).
#   The resulting figure compares the three richness measures side-by-side.
############################################################

# ---- Load packages ----
pacman::p_load(
  data.table, tidyverse, IRanges, vegan, SRS, ggh4x, ggsci
)

# ---- 1. Define parameters ----
Cmin <- 10            # Rarefaction / SRS threshold
SRS_seed <- 1988      # reproducible seed
age_min <- -500
age_max <- 16000
bin_width <- 500

# ---- 2. Import data ----
bugs <- fread(
  here::here("analysis/data/raw_data/bugs_europe_extraction_samples_20250612.csv"), na.strings = c("", "NA", "NULL"), encoding = "UTF-8"
)

# ---- 3. Filter relevant records ----
# Keep only stratigraphic, non-Greenland samples with valid ages and narrow ranges
bugs_eu <- bugs %>%
  mutate(age_range = age_older - age_younger,
         region = case_when(
           between(latitude, 49.8, 62.6) & between(longitude, -12.6, 1.8) ~ "British/Irish Isles")) %>%
  filter(
    region == "British/Irish Isles",
    context == "Stratigraphic sequence",
    sample != "BugsPresence",
    between(age_older, age_min, age_max),
    between(age_younger, age_min, age_max),
    age_range <= 2000,
    country != "Greenland"
  ) %>%
  distinct()

# ---- 4. Create unique sample identifiers ----
# Combine 'sample', 'sample_group', and 'site' to ensure uniqueness
sample_meta <- bugs_eu %>%
  distinct(sample, sample_group, site, .keep_all = TRUE) %>%
  transmute(
    sample_id = paste(sample, sample_group, site, sep = "@"),
    sample, sample_group, site,
    start = age_younger, end = age_older
  )

# ---- 5. Define 500-year time bins ----
age_bins <- tibble(
  start = c(age_min, seq(1, age_max - bin_width + 1, by = bin_width)),
  end   = seq(0, age_max, by = bin_width)
)

# ---- 6. Assign samples to time bins (IRanges overlap) ----
ir_samples <- IRanges(start = sample_meta$start, end = sample_meta$end)
ir_bins    <- IRanges(start = age_bins$start, end = age_bins$end)
ov <- findOverlaps(ir_samples, ir_bins, type = "any")

hits <- tibble(
  sample_id = sample_meta$sample_id[queryHits(ov)],
  bin_start = age_bins$start[subjectHits(ov)],
  bin_end   = age_bins$end[subjectHits(ov)]
) %>%
  separate(sample_id, into = c("sample", "sample_group", "site"), sep = "@", remove = FALSE) %>%
  inner_join(bugs, by = c("sample", "sample_group", "site"), relationship = "many-to-many") %>%
  select(sample, sample_group, site, sample_id, end = bin_end, taxon, abundance) %>%
  distinct()

# ---- 7. Pivot data to a sample × taxon matrix ----
raw_wide <- hits %>%
  select(sample_id, end, taxon, abundance) %>%
  distinct() %>%
  pivot_wider(
    id_cols = c(sample_id, end),
    names_from = taxon,
    values_from = abundance,
    values_fill = 0
  ) %>%
  arrange(end)

species_cols <- setdiff(colnames(raw_wide), c("sample_id", "end"))

# Ensure abundance columns are numeric
raw_wide <- raw_wide %>%
  mutate(across(all_of(species_cols), ~ as.numeric(.)))

# ---- 8. Compute raw richness ----
raw.rich <- raw_wide %>%
  rowwise() %>%
  mutate(richness = specnumber(c_across(all_of(species_cols)))) %>%
  ungroup() %>%
  group_by(end) %>%
  summarise(
    mean_rich = round(mean(richness, na.rm = TRUE), 2),
    err = sd(richness, na.rm = TRUE) / sqrt(n()),
    type = "Raw",
    .groups = "drop"
  )

# ---- 9. Compute rarefied richness ----
# Exclude samples with fewer than Cmin individuals
rare_wide <- raw_wide %>%
  rowwise() %>%
  mutate(total_abund = sum(c_across(all_of(species_cols)))) %>%
  ungroup() %>%
  filter(total_abund >= Cmin)

if (nrow(rare_wide) > 0) {
  rare_rich <- tibble(
    end = rare_wide$end,
    richness = rarefy(select(rare_wide, all_of(species_cols)) %>% as.matrix(), sample = Cmin)
  ) %>%
    group_by(end) %>%
    summarise(
      mean_rich = round(mean(richness), 2),
      err = sd(richness) / sqrt(n()),
      type = "Rarefaction",
      .groups = "drop"
    )
} else {
  rare_rich <- tibble(end = numeric(), mean_rich = numeric(), err = numeric(), type = character())
}

# ---- 10. Compute SRS standardized richness ----
# Prepare taxa × samples matrix
srs_prep <- hits %>%
  mutate(sample_label = paste(sample, sample_group, site, end, sep = "@")) %>%
  distinct(taxon, sample_label, abundance) %>%
  pivot_wider(
    names_from = sample_label,
    values_from = abundance,
    values_fill = 0
  ) %>%
  column_to_rownames("taxon")

# Remove empty samples
srs_prep <- srs_prep[, colSums(srs_prep) > 0, drop = FALSE]

# Keep as data.frame (required for SRS internal behavior)
srs_prep <- as.data.frame(srs_prep)

# Remove samples with total < Cmin
sample_sums <- colSums(srs_prep)
if (any(sample_sums < Cmin)) {
  warning(sum(sample_sums < Cmin), " samples have fewer individuals than Cmin — excluded.")
  keep <- sample_sums >= Cmin
  srs_prep <- srs_prep[, keep, drop = FALSE]
}

# Run SRS algorithm
set.seed(SRS_seed)
srs_out <- SRS(srs_prep, Cmin = Cmin, set_seed = TRUE, seed = SRS_seed)

# Transpose to samples × taxa
srs_mat <- as.data.frame(t(srs_out))

# Compute richness per bin
srs.rich <- srs_mat %>%
  rownames_to_column("sample_id") %>%
  separate(sample_id, into = c("sample", "sample_group", "site", "end"), sep = "@", fill = "right") %>%
  mutate(end = as.numeric(end)) %>%
  rowwise() %>%
  mutate(richness = specnumber(c_across(-c(sample, sample_group, site, end)))) %>%
  ungroup() %>%
  group_by(end) %>%
  summarise(
    mean_rich = round(mean(richness, na.rm = TRUE), 2),
    err = sd(richness, na.rm = TRUE) / sqrt(n()),
    type = "SRS",
    .groups = "drop"
  )

# ---- 11. Combine all richness estimates ----
rich.mat <- bind_rows(raw.rich, rare_rich, srs.rich)
rich.mat$type <- factor(rich.mat$type, levels = c("Raw", "Rarefaction", "SRS"))

# ---- 12. Define time-period background rectangles (Pilotto et al. 2022) ----
rects.n <- data.frame(
  xstart = c(12000,9500,8000,6500,4500,4000,3500,-500),
  xend   = c(16000,12000,9500,8000,6500,4500,4000,3500),
  col    = c("TM 1","TM 2","TM 3","TM 4","TM 5","TM 6","TM 7","TM 8")
)
rects.n$col <- factor(rects.n$col, levels = paste("TM", 1:8))

# ---- 13. Generate figure ----
fig <- ggplot() +
  geom_rect(
    data = rects.n,
    aes(ymin = xstart, ymax = xend, xmin = -Inf, xmax = Inf, fill = col),
    alpha = 0.5
  ) +
  geom_ribbon(
    data = rich.mat,
    aes(y = end, xmin = mean_rich - err, xmax = mean_rich + err),
    alpha = 0.5
  ) +
  geom_path(
    data = rich.mat,
    aes(y = end, x = mean_rich),
    linewidth = 1, linetype = "dashed"
  ) +
  geom_point(
    data = rich.mat,
    aes(y = end, x = mean_rich),
    shape = 21, size = 3, fill = "green4", col = "black"
  ) +
  scale_y_reverse(limits = c(16000, -500), breaks = scales::pretty_breaks(n = 10)) +
  scale_fill_jco(name = "Time Periods") +
  facet_grid2(~ type, scales = "free", axes = "margins", independent = "x") +
  labs(y = "Age (BP)", x = "Mean richness", color = NULL) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.x = element_text(size = 14, face = "bold", colour = "black"),
    plot.title = element_text(size = 16, face = "bold", colour = "black", hjust = 0.5),
    strip.text.x = element_text(size = 14, face = "bold", colour = "black"),
    strip.text.y = element_text(size = 12, face = "bold", colour = "black")
  )

# ---- 15. Save figure ----
ggsave("002-richness-comparison-britain.jpg",
       fig,
       device = "jpg",
       here::here("analysis", "figures"),
       width=20,
       height=30,
       units = "cm",
       dpi = 300)

message("✅ Species richness comparison completed and figure saved: 002-richness-comparison-britain.jpg")
