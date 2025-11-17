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

# ---- 0. Load required packages ----
pacman::p_load(
  cowplot, data.table, entropart, ggh4x, ggplot2,
  ggsci, tidyverse, IRanges, here, SRS
)

# Define parameters
Cmin <- 10            # Rarefaction / SRS threshold
SRS_seed <- 1988      # reproducible seed

# ==============================================================
# 1. Import species occurrence data
# ==============================================================

bugs <- fread(
  here::here("analysis", "data", "raw_data", "bugs_europe_extraction_samples_20250612.csv"),
  na.strings = c("", "NA", "NULL"),
  encoding = "UTF-8"
)

# ==============================================================
# 2. Define temporal periods (TM1–TM8)
# ==============================================================

rects <- data.frame(
  ystart = c(12000, 9500, 8000, 6500, 4500, 4000, 3500, -500),
  yend   = c(16000, 12000, 9500, 8000, 6500, 4500, 4000, 3500),
  col    = factor(
    c("TM 1", "TM 2", "TM 3", "TM 4", "TM 5", "TM 6", "TM 7", "TM 8"),
    levels = c("TM 1", "TM 2", "TM 3", "TM 4", "TM 5", "TM 6", "TM 7", "TM 8")
  )
)

# ==============================================================
# 3. Filter and prepare temporal range data
# ==============================================================
time_mat <- bugs %>%
  select(country, latitude, longitude, sample, site, sample_group, age_older, age_younger, context) %>%
  mutate(age_range = age_older - age_younger,
         region = case_when(
           between(latitude, 49.8, 62.6) & between(longitude, -12.6, 1.8) ~ "British/Irish Isles",
           TRUE ~ "")
         )%>%
  filter(
    region == "British/Irish Isles",
    context == "Stratigraphic sequence",
    sample != "BugsPresence",
    age_range <= 2000,
    country != "Greenland"
  ) %>%
  mutate(mid_age = (age_older + age_younger) / 2) %>%
  filter(between(mid_age, -500, 16000)) %>%
  distinct() %>%
  select(-age_range, -context, -mid_age) %>%
  mutate(sample = paste(sample, sample_group, site, sep = "@")) %>%
  column_to_rownames("sample") %>%
  select(age_older, age_younger) %>%
  dplyr::rename(start = age_younger, end = age_older)

# ==============================================================
# 4. Define 500-year temporal bins
# ==============================================================

range <- data.frame(
  start = c(-500, seq(1, 15501, by = 500)),
  end = seq(0, 16000, by = 500)
)

# ==============================================================
# 5. Identify overlapping samples with time bins
# ==============================================================

intersection <- findOverlaps(
  query = do.call(IRanges, time_mat),
  subject = do.call(IRanges, range),
  type = "any"
)

# ==============================================================
# 6. Merge sample metadata with bin assignments
# ==============================================================

hits <- data.frame(time_mat[queryHits(intersection),], range[subjectHits(intersection),]) %>%
  as_tibble(rownames = "sample") %>%
  separate(sample, into = c("sample", "sample_group", "site"), sep = "@") %>%
  inner_join(bugs, by = c("sample", "sample_group")) %>%
  select(-start, -end)

# ==============================================================
# 7. Construct raw species abundance matrix
# ==============================================================
raw_mat <- hits %>%
  select(site.x, latitude, longitude, sample_group, sample, end = end.1, taxon, abundance) %>%
  distinct() %>%
  mutate(sample_id = paste(sample, sample_group, site.x, end, sep = "@")) %>%
  mutate(
    place = case_when(
      between(latitude, 53.5, 74.8) & between(longitude, 2.8, 41) ~ "Scandinavia/Baltic",
      between(latitude, 49.8, 62.6) & between(longitude, -12.6, 1.8) ~ "British/Irish Isles",
      between(latitude, 63, 66.3) & between(longitude, -25, -12) ~ "Iceland",
      latitude < 45 ~ "Meridional"
    ),
    region = ifelse(is.na(place), "Continental", place)
  ) %>%
  select(-place, -latitude, -longitude)

# ==============================================================
# 8. Pivot data to a sample × taxon matrix
# ==============================================================
raw_mat <- hits %>%
  mutate(sample_id = paste(sample, sample_group, site.x, end.1, sep = "@")) %>%
  select(sample_id, end = end.1, taxon, abundance) %>%
  distinct() %>%
  pivot_wider(
    id_cols = c(sample_id, end),
    names_from = taxon,
    values_from = abundance,
    values_fill = 0
  ) %>%
  arrange(end)

species_cols <- setdiff(colnames(raw_mat), c("sample_id", "end"))

# Ensure abundance columns are numeric
raw_mat <- raw_mat %>%
  mutate(across(all_of(species_cols), ~ as.numeric(.)))

# ==============================================================
# 9. Compute raw richness
# ==============================================================
raw.rich <- raw_mat %>%
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

# ==============================================================
# 10. Compute rarefied richness
# ==============================================================
# Exclude samples with fewer than Cmin individuals
raw_mat <- raw_mat %>%
  rowwise() %>%
  mutate(total_abund = sum(c_across(all_of(species_cols)))) %>%
  ungroup() %>%
  filter(total_abund >= Cmin)

if (nrow(raw_mat) > 0) {
  rare_rich <- tibble(
    end = raw_mat$end,
    richness = rarefy(select(raw_mat, all_of(species_cols)) %>% as.matrix(), sample = Cmin)
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

# ==============================================================
# 11. Compute SRS standardized richness
# ==============================================================
# Prepare taxa × samples matrix
srs_prep <- hits %>%
  mutate(sample_id = paste(sample, sample_group, site.x, end.1, sep = "@")) %>%
  distinct(taxon, sample_id, abundance) %>%
  pivot_wider(
    names_from = sample_id,
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

# ==============================================================
# 12. Combine all richness estimates
# ==============================================================
rich.mat <- bind_rows(raw.rich, rare_rich, srs.rich)
rich.mat$type <- factor(rich.mat$type, levels = c("Raw", "Rarefaction", "SRS"))

# ==============================================================
# 13. Generate figure
# ==============================================================
fig <- ggplot() +
  geom_rect(
    data = rects,
    aes(ymin = ystart, ymax = yend, xmin = -Inf, xmax = Inf, fill = col),
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
  facet_wrap(~type, scales = "free_x") +
  labs(y = "Time (Years BP)", x = "Mean richness", color = NULL) +
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

# ==============================================================
# 14. Save figure
# ==============================================================
ggsave("002-richness-comparison-britain.jpg",
       fig,
       device = "jpg",
       here::here("analysis", "figures"),
       width=20,
       height=30,
       units = "cm",
       dpi = 300)

message("✅ Species richness comparison completed and figure saved: 002-richness-comparison-britain.jpg")
