# ==============================================================
# Script: Richness and Coverage Analysis (Raw Data)
#
# Purpose:
#   Compute richness and coverage metrics (Observed Richness,
#   Coverage, Estimated Richness) for Holocene fossil insect data
#   from British/Irish Isles using entropart.
#
# Workflow:
#   1. Load fossil insect occurrence data.
#   2. Filter samples by region, context, and age constraints.
#   3. Assign samples to 500-year temporal bins using IRanges.
#   4. Prepare raw abundance matrix.
#   5. Build list of matrices per bin.
#   6. Compute richness and coverage metrics for each bin.
#   7. Plot results.
#
# Key Notes:
#   - Uses IRanges for robust interval overlap detection.
#   - Diagnostic messages report taxa and samples dropped during prep
#     and bins skipped due to insufficient data.
#
# Output:
#   - Plot: "002-richness-estimation-britain.jpg"
# ==============================================================

# ---- 0. Load packages ----
pacman::p_load(
  data.table, tidyverse, here, IRanges, entropart, ggsci, ggplot2
)

# ==============================================================
# 1. Import species occurrence data
# ==============================================================
bugs <- fread(here("analysis", "data", "raw_data", "bugs_europe_extraction_samples_20250612.csv"),
              na.strings = c("", "NA", "NULL"), encoding = "UTF-8") %>%
  mutate(sample_id = paste(sample, sample_group, site, sep = "@")) %>%
  as_tibble()

# ==============================================================
# 2. Define temporal periods (TM1–TM8)
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
# 7. Prepare raw abundance matrix
# ==============================================================
raw_mat <- hits %>%
  select(sample_id, bin_end, taxon, abundance) %>%
  distinct() %>%
  mutate(bin_end = as.numeric(bin_end))

if (nrow(raw_mat) == 0) stop("raw_mat is empty after filtering.")

# ==============================================================
# 8. Build lists of matrices per bin
# ==============================================================
build_Listdf <- function(df) {
  split(df, df$bin_end) %>%
    map(~ .x %>%
          select(taxon, sample_id, abundance) %>%
          pivot_wider(names_from = sample_id, values_from = abundance, values_fill = 0) %>%
          column_to_rownames("taxon") %>%
          as.matrix())
}

Listdf_raw <- build_Listdf(raw_mat)

# ==============================================================
# 9. Compute richness and coverage metrics
# ==============================================================
message("Running richness and coverage analysis...")

results <- map_dfr(names(Listdf_raw), function(time_name) {
  mat <- Listdf_raw[[time_name]]

  # Initial counts
  initial_taxa <- nrow(mat)
  initial_samples <- ncol(mat)

  # Drop empty taxa/samples
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  mat <- mat[, colSums(mat) > 0, drop = FALSE]

  dropped_taxa <- initial_taxa - nrow(mat)
  dropped_samples <- initial_samples - ncol(mat)

  message("[Raw] Bin ", time_name, ": Taxa retained = ", nrow(mat),
          " (Dropped: ", dropped_taxa, "), Samples retained = ", ncol(mat),
          " (Dropped: ", dropped_samples, ")")

  # Skip bins with too few taxa/samples
  if (nrow(mat) < 2 || ncol(mat) < 2) {
    message("⚠️ [Raw] Bin ", time_name, " skipped (too few taxa or samples).")
    return(NULL)
  }

  # Compute metrics
  obs_rich <- sum(rowSums(mat) > 0)
  cov <- tryCatch(Coverage(mat), error = function(e) NA)

  Weights <- colSums(mat)
  Weights <- Weights / sum(Weights)
  MC <- tryCatch(MetaCommunity(mat, Weights), error = function(e) NULL)

  richness <- NA
  if (!is.null(MC)) {
    richness <- tryCatch(Richness(MC$Ns, Correction = "Chao1"), error = function(e) NA)
  }

  tibble(
    Time = as.numeric(time_name),
    Metric = c("Observed Richness", "Coverage", "Estimated Richness"),
    Value = c(obs_rich, cov, richness)
  )
})

if (nrow(results) == 0) stop("No valid bins for richness estimation.")

# ==============================================================
# 10. Plot richness and coverage metrics
# ==============================================================
results <- results %>%
  mutate(Metric = factor(Metric, levels = c("Observed Richness", "Coverage", "Estimated Richness")))

plot_min <- min(results$Time)
plot_max <- max(results$Time)
rects_filtered <- rects %>% filter(ystart <= plot_max, yend >= plot_min)
time_breaks <- pretty(results$Time, n = 10)

p <- ggplot(results, aes(x = Value, y = Time)) +
  geom_rect(data = rects_filtered,
            aes(ymin = ystart, ymax = yend, xmin = -Inf, xmax = Inf, fill = col),
            inherit.aes = FALSE, alpha = 0.5) +
  geom_path(aes(group = Metric), color = "black", linetype = "dashed") +
  geom_point(size = 3, color = "black") +
  facet_wrap(~Metric, scales = "free_x") +
  scale_fill_jco(guide = guide_legend(reverse = TRUE)) +
  scale_y_reverse(breaks = time_breaks, labels = time_breaks) +
  labs(title = "Richness and Coverage",
       x = "Value", y = "Time (Years BP)", fill = "Time Periods") +
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
ggsave(filename = "002-richness-estimation-britain.jpg",
       plot = p, path = here("analysis", "figures"),
       width = 3300, height = 4200, units = "px", dpi = 300)

message("✅ Richness and coverage plot generated successfully.")
