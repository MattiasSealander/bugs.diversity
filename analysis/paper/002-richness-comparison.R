# ==============================================================
# Script: Richness Comparison (Raw, Rarefaction, SRS)
#
# Purpose:
#   Compare three richness estimates by 500-year time bins:
#     - Raw richness (unadjusted)
#     - Rarefied richness (sample-size standardized, Cmin = 10)
#     - SRS standardized richness (Scaling with Ranked Subsampling, Cmin = 10)
#
# Workflow:
#   1. Load fossil insect occurrence data.
#   2. Filter samples by region, context, and age constraints.
#   3. Assign samples to 500-year temporal bins using IRanges.
#   4. Prepare raw abundance matrix.
#   5. Compute raw richness per bin.
#   6. Compute rarefied richness per bin (exclude samples < Cmin).
#   7. Standardize abundances using SRS and compute richness per bin.
#   8. Combine all richness estimates and plot.
#
# Key Notes:
#   - Diagnostic messages report taxa and samples dropped during SRS prep
#     and samples excluded during rarefaction filtering.
#   - Bins with too few samples are skipped.
#
# Output:
#   - Plot: "002-richness-comparison-britain.jpg"
# ==============================================================

# ---- 0. Load packages ----
pacman::p_load(
  data.table, tidyverse, cowplot, here, IRanges, SRS, vegan, ggsci
)

# Parameters
Cmin <- 10
SRS_seed <- 1988

# ==============================================================
# 1. Import species occurrence data
# ==============================================================
bugs <- fread(here("analysis", "data", "raw_data", "bugs_europe_extraction_samples_20250612.csv")) %>%
  mutate(sample_id = paste(sample, sample_group, site, sep = "@")) %>%
  as_tibble()

# ==============================================================
# 2. Define periods
# ==============================================================
rects <- tibble(
  ystart = c(11700, 7000, 5000, 0),
  yend   = c(16000, 11700, 7000, 5000),
  col    = factor(c("Lateglacial", "Early Holocene", "Mid-Holocene", "Late Holocene"),
                  levels = c("Lateglacial", "Early Holocene", "Mid-Holocene", "Late Holocene"))
)

# ==============================================================
# 3. Filter and prepare temporal range data
# ==============================================================
time_mat <- bugs %>%
  select(country, latitude, longitude, sample_id, sample, age_older, age_younger, context, family_name) %>%
  mutate(age_range = age_older - age_younger,
         mid_age = (age_older + age_younger) / 2,
         region = case_when(
           between(latitude, 49.8, 62.6) & between(longitude, -12.6, 1.8) ~ "British/Irish Isles",
           TRUE ~ "")
  ) %>%
  filter(!family_name %in% c("ACANTHOSOMATIDAE","ANTHOCORIDAE", "APHALARIDAE", "APHIDOIDEA", "AUCHENORRHYNCHA", "BIBIONIDAE",
                             "BRACHYCENTRIDAE", "CALOPTERYGIDAE", "CERCOPIDAE", "CHIRONOMIDAE", "CICADELLIDAE", "CICADOMORPHA",
                             "CIXIIDAE", "CORIXIDAE", "CYCLORRHAPHA", "CYDNIDAE", "DELPHACIDAE", "DERMAPTERA", "DIPTERA", "FORFICULIDAE",
                             "FORMICIDAE", "FULGOROMORPHA", "GERRIDAE", "GLOSSOSOMATIDAE", "GOERIDAE", "HEBRIDAE", "HEMIPTERA",
                             "HETEROPTERA", "HOMOPTERA", "HYDROPSYCHIDAE", "HYDROPTILIDAE", "HYMENOPTERA", "LEPIDOPTERA", "LEPIDOSTOMATIDAE",
                             "LEPTOCERIDAE", "LIMNEPHILIDAE", "LYGAEIDAE", "MEMBRACIDAE", "MICROPHYSIDAE", "MICROSPORIDAE", "MIRIDAE",
                             "MOLANNIDAE", "NABIDAE", "NEMATOCERA", "NEMOURIDAE", "ODONATA", "PARASITICA", "PENTATOMIDAE", "PHRYGANEIDAE",
                             "POLYCENTROPIDAE", "PSYCHOMYIIDAE", "PSYLLIDAE", "RAPHIDIIDAE", "SALDIDAE", "SCUTELLERIDAE", "SERICOSTOMATIDAE",
                             "SIALIDAE", "TRICHOPTERA", "THYREOCORIDAE", "TINGIDAE", "TIPULIDAE", "TRIOPSIDAE", "TRIOZIDAE"),
         region == "British/Irish Isles",
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
  inner_join(bugs, by = "sample_id")

# ==============================================================
# 7. Prepare raw abundance matrix
# ==============================================================
raw_mat <- hits %>%
  select(sample_id, bin_end, taxon, abundance) %>%
  distinct() %>%
  mutate(bin_end = as.numeric(bin_end))

if (nrow(raw_mat) == 0) stop("raw_mat is empty after filtering.")

# ==============================================================
# 8. Compute raw richness per bin
# ==============================================================
raw_rich <- raw_mat %>%
  group_by(sample_id, bin_end) %>%
  summarise(richness = n_distinct(taxon), .groups = "drop") %>%
  group_by(bin_end) %>%
  summarise(mean_rich = round(mean(richness), 2),
            err = sd(richness) / sqrt(n()),
            type = "Raw", .groups = "drop")

# ==============================================================
# 9. Compute rarefied richness
# ==============================================================
# Pivot to sample × taxon matrix
rare_prep <- raw_mat %>%
  pivot_wider(names_from = taxon, values_from = abundance, values_fill = 0)

# Filter samples below Cmin
rare_prep <- rare_prep %>%
  mutate(total_abund = rowSums(across(-c(sample_id, bin_end)))) %>%
  filter(total_abund >= Cmin)

message("✅ Rarefaction: Samples retained = ", nrow(rare_prep),
        " (Dropped: ", nrow(raw_mat) - nrow(rare_prep), ")")

rare_rich <- if (nrow(rare_prep) > 0) {
  tibble(
    bin_end = rare_prep$bin_end,
    richness = rarefy(select(rare_prep, -c(sample_id, bin_end)), sample = Cmin)
  ) %>%
    group_by(bin_end) %>%
    summarise(mean_rich = round(mean(richness), 2),
              err = sd(richness) / sqrt(n()),
              type = "Rarefaction", .groups = "drop")
} else {
  tibble(bin_end = numeric(), mean_rich = numeric(), err = numeric(), type = character())
}

# ==============================================================
# 10. Compute SRS standardized richness
# ==============================================================
# Prepare taxa × samples matrix
srs_prep <- raw_mat %>%
  distinct(taxon, sample_id, abundance) %>%
  pivot_wider(names_from = sample_id, values_from = abundance, values_fill = 0) %>%
  column_to_rownames("taxon")

# Drop taxa and samples with diagnostics
initial_taxa <- nrow(srs_prep)
initial_samples <- ncol(srs_prep)
srs_prep <- srs_prep[rowSums(srs_prep) > 0, , drop = FALSE]
message("✅ SRS: Taxa retained = ", nrow(srs_prep), " (Dropped: ", initial_taxa - nrow(srs_prep), ")")
srs_prep <- srs_prep[, colSums(srs_prep) >= Cmin, drop = FALSE]
message("✅ SRS: Samples retained = ", ncol(srs_prep), " (Dropped: ", initial_samples - ncol(srs_prep), ")")

# Perform SRS
set.seed(SRS_seed)
srs_out <- SRS(as.data.frame(srs_prep), Cmin = Cmin, set_seed = TRUE, seed = SRS_seed)
rownames(srs_out) <- rownames(srs_prep)

# Compute richness per bin
srs_long <- as_tibble(srs_out, rownames = "taxon") %>%
  pivot_longer(-taxon, names_to = "sample_id", values_to = "abundance") %>%
  filter(abundance > 0) %>%
  inner_join(raw_mat %>% select(sample_id, bin_end), by = "sample_id", relationship = "many-to-many") %>%
  group_by(sample_id, bin_end) %>%
  summarise(richness = n_distinct(taxon), .groups = "drop") %>%
  group_by(bin_end) %>%
  summarise(mean_rich = round(mean(richness), 2),
            err = sd(richness) / sqrt(n()),
            type = "SRS", .groups = "drop")

# ==============================================================
# 11. Combine all richness estimates
# ==============================================================
rich.mat <- bind_rows(raw_rich, rare_rich, srs_long)
rich.mat$type <- factor(rich.mat$type, levels = c("Raw", "Rarefaction", "SRS"))

# ==============================================================
# 12. Plot richness comparison
# ==============================================================
plot_min <- min(rich.mat$bin_end)
plot_max <- max(rich.mat$bin_end)
rects_filtered <- rects %>% filter(ystart <= plot_max, yend >= plot_min)
time_breaks <- pretty(rich.mat$bin_end, n = 10)

fig <- ggplot(rich.mat, aes(x = mean_rich, y = bin_end)) +
  geom_rect(data = rects_filtered,
            aes(ymin = ystart, ymax = yend, xmin = -Inf, xmax = Inf, fill = col),
            inherit.aes = FALSE, alpha = 0.5) +
  geom_ribbon(aes(y = bin_end, xmin = mean_rich - err, xmax = mean_rich + err), alpha = 0.5) +
  geom_path(linewidth = 1, linetype = "dashed") +
  geom_point(shape = 21, size = 4, fill = "green4", color = "black") +
  facet_wrap(~type, scales = "free_x") +
  scale_fill_jco(guide = guide_legend(reverse = TRUE)) +
  scale_y_reverse(breaks = time_breaks, labels = time_breaks) +
  coord_cartesian(ylim = c(plot_max, plot_min)) +
  labs(title = "Species Richness",
       x = "Mean richness", y = "Time (Years BP)", fill = "Time Periods") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(size = 12, colour = "black"),
    axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
    axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
    strip.background = element_rect(fill = "#f0f0f0", color = NA),
    strip.text.x = element_text(size = 14, face = "bold", colour = "black"),
    plot.title = element_text(size = 16, face = "bold", colour = "black"),
    legend.position = "right"
  )

ggsave("002-richness-comparison.jpg",
       fig, path = here("analysis", "figures"),
       width = 3300, height = 4200, units = "px", dpi = 300)

message("✅ Richness comparison completed and figure saved: 002-richness-comparison.jpg")
