# ==============================================================
# Script: Alpha Diversity Analysis (SRS Standardized)
#
# Purpose:
#   Compute alpha diversity metrics (Hill numbers: Richness, Shannon, Simpson)
#   for Holocene fossil insect data from the British/Irish Isles,
#   using SRS (Scaling with Ranked Subsampling) standardized abundances.
#
# Workflow:
#   1. Load fossil insect occurrence data.
#   2. Filter samples by region, context, and age constraints.
#   3. Assign samples to 500-year temporal bins using IRanges.
#   4. Standardize abundances using SRS:
#        - Deduplicate taxon × sample combinations.
#        - Remove taxa with zero abundance and samples below Cmin threshold.
#        - Apply SRS to equalize sampling effort across samples.
#   5. Build list of community matrices per time bin (taxon × sample).
#   6. Compute alpha diversity metrics for each bin using entropart:
#        - Richness (q = 0)
#        - Shannon (q = 1)
#        - Simpson (q = 2)
#      Correction method: "UnveilJ" (accounts for unseen species).
#   7. Plot diversity metrics across time with Holocene time slices.
#
# Key Notes:
#   - Uses IRanges for robust interval overlap detection.
#   - SRS ensures fair comparison by standardizing sample sizes.
#   - Bins with <2 taxa or <2 samples are skipped (metrics not meaningful).
#
# Output:
#   - Plot: "002-alpha-diversity-britain-srs.jpg"
# ==============================================================

# ---- 0. Load packages ----
pacman::p_load(
  data.table, tidyverse, here, IRanges, SRS, entropart, ggsci
)

# ==============================================================
# 1. Import species occurrence data
# ==============================================================
bugs <- fread(here("analysis", "data", "raw_data", "bugs_europe_extraction_samples_20250612.csv")) %>%
  mutate(sample_id = paste(sample, sample_group, site, sep = "@")) %>%
  as_tibble()

# ==============================================================
# 2. Define temporal periods (TM1–TM8) after Pilotto et al. (2022)
# ==============================================================
rects <- tibble(
  ystart = c(12000, 9500, 8000, 6500, 4500, 4000, 3500, -500),
  yend   = c(16000, 12000, 9500, 8000, 6500, 4500, 4000, 3500),
  col    = factor(c("TM 1", "TM 2", "TM 3", "TM 4", "TM 5", "TM 6", "TM 7", "TM 8"),
                  levels = c("TM 1", "TM 2", "TM 3", "TM 4", "TM 5", "TM 6", "TM 7", "TM 8"))
)

# ==============================================================
# 3. Filter and prepare temporal range data
# ==============================================================
time_mat <- bugs %>%
  select(country, latitude, longitude, sample_id, sample, age_older, age_younger, context) %>%
  mutate(age_range = age_older - age_younger,
         mid_age = (age_older + age_younger) / 2,
         region = case_when(
           between(latitude, 49.8, 62.6) & between(longitude, -12.6, 1.8) ~ "British/Irish Isles",
           TRUE ~ "")
  ) %>%
  filter(region == "British/Irish Isles",
         context == "Stratigraphic sequence",
         sample != "BugsPresence",
         age_range <= 2000,
         between(mid_age, -500, 16000)) %>%
  distinct() %>%
  select(sample_id, age_younger, age_older) %>%
  dplyr::rename(start = age_younger, end = age_older)

# Fix reversed ranges
if (any(time_mat$start > time_mat$end)) {
  time_mat <- time_mat %>%
    mutate(start = pmin(start, end), end = pmax(start, end))
}

# ==============================================================
# 4. Define 500-year temporal bins
# ==============================================================
range <- tibble(
  start = c(-500, seq(1, 15501, by = 500)),
  end   = seq(0, 16000, by = 500)
)

# ==============================================================
# 5. Identify overlapping samples with time bins
# ==============================================================
query <- IRanges(start = time_mat$start, end = time_mat$end)
subject <- IRanges(start = range$start, end = range$end)

intersection <- findOverlaps(query, subject, type = "any")
if (length(intersection) == 0) stop("No overlaps found. Check age ranges and bin definitions.")

# ==============================================================
# 6. Merge sample metadata with bin assignments
# ==============================================================
hits <- tibble(
  sample_id = time_mat$sample_id[queryHits(intersection)],
  bin_start = range$start[subjectHits(intersection)],
  bin_end   = range$end[subjectHits(intersection)]
) %>%
  inner_join(bugs, by = "sample_id", relationship = "many-to-many")

# ==============================================================
# 7. Prepare raw matrix for SRS
# ==============================================================
raw_mat <- hits %>%
  select(sample_id, bin_end, taxon, abundance) %>%
  distinct()

# Deduplicate taxon × sample_id
srs_agg <- raw_mat %>%
  group_by(taxon, sample_id) %>%
  summarise(abundance = sum(abundance), .groups = "drop")

# Pivot to taxa × samples
srs_prep <- srs_agg %>%
  pivot_wider(names_from = sample_id, values_from = abundance, values_fill = 0) %>%
  column_to_rownames("taxon")

# Drop taxa with all zeros
Cmin=10
initial_taxa <- nrow(srs_prep)
srs_prep <- srs_prep[rowSums(srs_prep) > 0, , drop = FALSE]
dropped_taxa <- initial_taxa - nrow(srs_prep)
message("✅ Taxa retained: ", nrow(srs_prep), " (Dropped: ", dropped_taxa, ")")

# Drop samples below threshold
initial_samples <- ncol(srs_prep)
srs_prep <- srs_prep[, colSums(srs_prep) >= Cmin, drop = FALSE]
dropped_samples <- initial_samples - ncol(srs_prep)
message("✅ Samples retained: ", ncol(srs_prep), " (Dropped: ", dropped_samples, ")")

# Perform SRS standardization with reproducible seed
SRS_seed <- 1988
set.seed(SRS_seed)
srs_out <- SRS(as.data.frame(srs_prep), Cmin = Cmin, set_seed = TRUE, seed = SRS_seed)
rownames(srs_out) <- rownames(srs_prep)

# Convert to long format and merge bin info
srs_long <- as_tibble(srs_out, rownames = "taxon") %>%
  pivot_longer(-taxon, names_to = "sample_id", values_to = "abundance") %>%
  filter(abundance > 0) %>%
  inner_join(raw_mat %>% select(sample_id, bin_end), by = "sample_id", relationship = "many-to-many") %>%
  distinct()

# ==============================================================
# 8. Build list of matrices per bin (same as raw script)
# ==============================================================
build_Listdf <- function(df) {
  split(df, df$bin_end) %>%
    map(~ .x %>%
          select(taxon, sample_id, abundance) %>%
          pivot_wider(names_from = sample_id, values_from = abundance, values_fill = 0) %>%
          column_to_rownames("taxon") %>%
          as.matrix())
}

Listdf_SRS <- build_Listdf(srs_long)

# ==============================================================
# 9. Compute alpha diversity per bin
# ==============================================================
alpha_results <- map_dfr(names(Listdf_SRS), function(time_name) {
  mat <- Listdf_SRS[[time_name]]

  # Initial counts
  initial_taxa <- nrow(mat)
  initial_samples <- ncol(mat)

  # Drop empty taxa/samples
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  mat <- mat[, colSums(mat) > 0, drop = FALSE]

  dropped_taxa <- initial_taxa - nrow(mat)
  dropped_samples <- initial_samples - ncol(mat)

  message("Bin ", time_name, ": Taxa retained = ", nrow(mat),
          " (Dropped: ", dropped_taxa, "), Samples retained = ", ncol(mat),
          " (Dropped: ", dropped_samples, ")")

  # Skip bins with too few taxa/samples
  if (nrow(mat) < 2 || ncol(mat) < 2) {
    message("⚠️ Bin ", time_name, " skipped (too few taxa or samples).")
    return(NULL)
  }

  # Compute diversity using abundance weighting
  Weights <- colSums(mat) / sum(colSums(mat))
  MC <- MetaCommunity(mat, Weights)

  tibble(
    Time = as.numeric(time_name),
    Richness = AlphaDiversity(MC, q = 0, Correction = "UnveilJ")$Total,
    Shannon  = AlphaDiversity(MC, q = 1, Correction = "UnveilJ")$Total,
    Simpson  = AlphaDiversity(MC, q = 2, Correction = "UnveilJ")$Total
  )
})

# ==============================================================
# 10. Plot alpha diversity
# ==============================================================
alpha_results <- alpha_results %>%
  pivot_longer(cols = c(Richness, Shannon, Simpson), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric, levels = c("Richness", "Shannon", "Simpson")))

plot_min <- min(alpha_results$Time, na.rm = TRUE)
plot_max <- max(alpha_results$Time, na.rm = TRUE)
rects_filtered <- rects %>% filter(ystart <= plot_max, yend >= plot_min)
time_breaks <- pretty(alpha_results$Time, n = 10)

p <- ggplot(alpha_results, aes(x = Value, y = Time)) +
  geom_rect(data = rects_filtered,
            aes(ymin = ystart, ymax = yend, xmin = -Inf, xmax = Inf, fill = col),
            inherit.aes = FALSE, alpha = 0.5) +
  geom_path(aes(group = Metric), color = "black", linetype = "dashed") +
  geom_point(size = 3, color = "black") +
  facet_wrap(~Metric, scales = "free_x") +
  scale_fill_jco(guide = guide_legend(reverse = TRUE)) +
  scale_y_reverse(breaks = time_breaks, labels = time_breaks) +
  labs(title = "Alpha Diversity", subtitle = "SRS standardized abundances",
       x = "Hill number", y = "Time (Years BP)", fill = "Time Periods") +
  coord_cartesian(ylim = c(plot_max, plot_min)) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(size = 12, colour = "black", angle = 45, vjust = .5),
    axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
    axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
    strip.background = element_rect(fill = "#f0f0f0", color = NA),
    strip.text.x = element_text(size = 14, face = "bold", colour = "black"),
    plot.title = element_text(size = 16, face = "bold", colour = "black"),
    legend.title = element_text(size = 16, face = "bold", colour = "black"),
    legend.text = element_text(size = 12, colour = "black"),
    legend.position = "right"
  )

# ==============================================================
# 11. Save figure
# ==============================================================
ggsave(filename = "003-alpha-diversity-britain-srs.jpg",
       plot = p, path = here("analysis", "figures"),
       width = 3300, height = 4200, units = "px", dpi = 300)

message("✅ Alpha diversity plot (SRS) generated successfully.")
