
# ==============================================================
# Script: Gamma Diversity (SRS only) + δ18O (GISP2) — Side-by-side
#
# Purpose:
#   Compute gamma diversity (Shannon, Simpson) using SRS-standardized
#   abundances for Holocene fossil insect data (British/Irish Isles),
#   and compare visually with GISP2 δ18O on aligned time axes.
#
# Output:
#   - Combined plot: "005-gamma-srs-vs-d18o-coleoptera.jpg"
# ==============================================================

# ---- 0. Load packages ----
pacman::p_load(
  data.table, ggsci, tidyverse, here, IRanges, SRS, entropart, patchwork, viridis
)

# Ensure figure directory exists
fig_dir <- here("analysis", "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# ---- 1. Import species occurrence data ----
bugs <- fread(here("analysis", "data", "raw_data", "bugs_europe_extraction_samples_20250612.csv")) %>%
  mutate(sample_id = paste(sample, sample_group, site, sep = "@")) %>%
  as_tibble()

# ---- 2. Filter and prepare temporal range data ----
time_mat <- bugs %>%
  select(country, latitude, longitude, sample_id, sample, age_older, age_younger, context, family_name) %>%
  mutate(
    age_range = age_older - age_younger,
    mid_age   = (age_older + age_younger) / 2,
    region = case_when(
      between(latitude, 49.8, 62.6) & between(longitude, -12.6, 1.8) ~ "British/Irish Isles",
      TRUE ~ ""
    )
  ) %>%
  filter(
    !family_name %in% c("ACANTHOSOMATIDAE","ANTHOCORIDAE", "APHALARIDAE", "APHIDOIDEA", "AUCHENORRHYNCHA", "BIBIONIDAE",
                        "BRACHYCENTRIDAE", "CALOPTERYGIDAE", "CERCOPIDAE", "CHIRONOMIDAE", "CICADELLIDAE", "CICADOMORPHA",
                        "CIXIIDAE", "CORIXIDAE", "CYCLORRHAPHA", "CYDNIDAE", "DELPHACIDAE", "DERMAPTERA", "DIPTERA", "FORFICULIDAE",
                        "FORMICIDAE", "FULGOROMORPHA", "GERRIDAE", "GLOSSOSOMATIDAE", "GOERIDAE", "HEBRIDAE", "HEMIPTERA",
                        "HETEROPTERA", "HOMOPTERA", "HYDROPSYCHIDAE", "HYDROPTILIDIDAE", "HYDROPTILIDAE", "HYMENOPTERA", "LEPIDOPTERA",
                        "LEPIDOSTOMATIDAE", "LEPTOCERIDAE", "LIMNEPHILIDAE", "LYGAEIDAE", "MEMBRACIDAE", "MICROPHYSIDAE", "MICROSPORIDAE",
                        "MIRIDAE", "MOLANNIDAE", "NABIDAE", "NEMATOCERA", "NEMOURIDAE", "ODONATA", "PARASITICA", "PENTATOMIDAE",
                        "PHRYGANEIDAE", "POLYCENTROPIDAE", "PSYCHOMYIIDAE", "PSYLLIDAE", "RAPHIDIIDAE", "SALDIDAE", "SCUTELLERIDAE",
                        "SERICOSTOMATIDAE", "SIALIDAE", "TRICHOPTERA", "THYREOCORIDAE", "TINGIDAE", "TIPULIDAE", "TRIOPSIDAE", "TRIOZIDAE"),
    region == "British/Irish Isles",
    context == "Stratigraphic sequence",
    sample != "BugsPresence",
    age_range <= 2000,
    between(mid_age, -500, 16000)
  ) %>%
  distinct() %>%
  select(sample_id, age_younger, age_older) %>%
  dplyr::rename(start = age_younger, end = age_older)

# Fix reversed ranges if any
if (any(time_mat$start > time_mat$end)) {
  time_mat <- time_mat %>%
    mutate(start = pmin(start, end), end = pmax(start, end))
}

# ==============================================================
# 2. Define periods
# ==============================================================
rects <- tibble(
  ystart = c(11700, 7000, 5000, 0),
  yend   = c(16000, 11700, 7000, 5000),
  col    = factor(c("Lateglacial", "Early Holocene", "Mid-Holocene", "Late Holocene"),
                  levels = c("Lateglacial", "Early Holocene", "Mid-Holocene", "Late Holocene"))
)

# ---- 3. Define 500-year temporal bins ----
range <- tibble(
  start = c(-500, seq(1, 15501, by = 500)),
  end   = seq(0, 16000, by = 500)
)

# ---- 4. Identify overlapping samples with time bins ----
query   <- IRanges(start = time_mat$start, end = time_mat$end)
subject <- IRanges(start = range$start, end = range$end)

intersection <- findOverlaps(query, subject, type = "any")
if (length(intersection) == 0) stop("No overlaps found. Check age ranges and bin definitions.")

# ---- 5. Merge sample metadata with bin assignments ----
hits <- tibble(
  sample_id = time_mat$sample_id[queryHits(intersection)],
  bin_start = range$start[subjectHits(intersection)],
  bin_end   = range$end[subjectHits(intersection)]
) %>%
  inner_join(bugs, by = "sample_id", relationship = "many-to-many")

# ---- 6. Prepare raw abundance matrix (pre-SRS) ----
raw_mat <- hits %>%
  select(sample_id, bin_end, taxon, abundance) %>%
  distinct() %>%
  mutate(bin_end = as.numeric(bin_end))

if (nrow(raw_mat) == 0) stop("raw_mat is empty after filtering.")

# ---- 7. SRS standardization ----
Cmin    <- 10
SRS_seed <- 1988
message("Running SRS standardization...")

# Aggregate abundances at taxon x sample
srs_agg <- raw_mat %>%
  group_by(taxon, sample_id) %>%
  summarise(abundance = sum(abundance), .groups = "drop")

# Pivot to taxa x samples
srs_prep <- srs_agg %>%
  pivot_wider(names_from = sample_id, values_from = abundance, values_fill = 0) %>%
  column_to_rownames("taxon")

# Drop empty taxa and samples (with diagnostics)
initial_taxa    <- nrow(srs_prep)
initial_samples <- ncol(srs_prep)

srs_prep <- srs_prep[rowSums(srs_prep) > 0, , drop = FALSE]
message("✅ Taxa retained after zero-filter: ", nrow(srs_prep), " (Dropped: ", initial_taxa - nrow(srs_prep), ")")

srs_prep <- srs_prep[, colSums(srs_prep) >= Cmin, drop = FALSE]
message("✅ Samples retained after Cmin filter: ", ncol(srs_prep), " (Dropped: ", initial_samples - ncol(srs_prep), ")")

# Perform SRS
set.seed(SRS_seed)
srs_out <- SRS(as.data.frame(srs_prep), Cmin = Cmin, set_seed = TRUE, seed = SRS_seed)
rownames(srs_out) <- rownames(srs_prep)

# Long format + attach bin info
srs_long <- as_tibble(srs_out, rownames = "taxon") %>%
  pivot_longer(-taxon, names_to = "sample_id", values_to = "abundance") %>%
  filter(abundance > 0) %>%
  inner_join(raw_mat %>% select(sample_id, bin_end), by = "sample_id", relationship = "many-to-many") %>%
  distinct()

# ---- 8. Build list of matrices per bin (SRS only) ----
build_Listdf <- function(df) {
  split(df, df$bin_end) %>%
    purrr::map(~ .x %>%
                 select(taxon, sample_id, abundance) %>%
                 pivot_wider(names_from = sample_id, values_from = abundance, values_fill = 0) %>%
                 column_to_rownames("taxon") %>%
                 as.matrix())
}

Listdf_SRS <- build_Listdf(srs_long)

# ---- 9. Compute gamma diversity (SRS only; Shannon & Simpson) ----
compute_gamma <- function(data_list, dataset_type) {
  purrr::map_dfr(names(data_list), function(time_name) {
    mat <- data_list[[time_name]]

    # Initial counts
    initial_taxa    <- nrow(mat)
    initial_samples <- ncol(mat)

    # Drop empty taxa/samples
    mat <- mat[rowSums(mat) > 0, , drop = FALSE]
    mat <- mat[, colSums(mat) > 0, drop = FALSE]

    dropped_taxa    <- initial_taxa - nrow(mat)
    dropped_samples <- initial_samples - ncol(mat)

    message("[", dataset_type, "] Bin ", time_name, ": Taxa retained = ", nrow(mat),
            " (Dropped: ", dropped_taxa, "), Samples retained = ", ncol(mat),
            " (Dropped: ", dropped_samples, ")")

    # Skip bins with too few taxa/samples
    if (nrow(mat) < 2 || ncol(mat) < 2) {
      message("⚠️ [", dataset_type, "] Bin ", time_name, " skipped (too few taxa or samples).")
      return(NULL)
    }

    # Compute diversity (q = 1: Shannon; q = 2: Simpson)
    MC <- MetaCommunity(mat)
    tibble(
      Time   = as.numeric(time_name),
      Metric = c("Richness", "Shannon", "Simpson"),
      Value  = c(GammaDiversity(MC, q = 0, Correction = "UnveiliC"),
                 GammaDiversity(MC, q = 1, Correction = "UnveiliC"),
                 GammaDiversity(MC, q = 2, Correction = "UnveiliC")),
      Type   = dataset_type
    )
  })
}

gamma_srs <- compute_gamma(Listdf_SRS, "SRS")

# ---- 10. Read δ18O (GISP2) and harmonize time span ----
d18o <- fread(here("analysis", "data", "raw_data", "d18o.csv"), na.strings = "NaN") %>%
  as_tibble() %>%
  filter(age_BP_SEAD <= 16000)

# ---- 11. Build plots: SRS diversity (no rects) + δ18O; side-by-side ----
time_min <- 0
time_max <- 16000
time_breaks <- seq(time_min, time_max, by = 1000)

# Diversity data for plotting (only Shannon/Simpson, restrict to 0–16 ka BP)
gamma_srs_plot <- gamma_srs %>%
  dplyr::filter(Metric %in% c("Shannon", "Simpson"),
                Time >= time_min, Time <= time_max)

# SRS gamma diversity plot (no shaded rects)
plot_min <- min(gamma_srs_plot$Time)
plot_max <- max(gamma_srs_plot$Time)
rects_filtered <- rects %>% filter(ystart <= plot_max, yend >= plot_min)

p_div <- ggplot(gamma_srs_plot, aes(x = Value, y = Time)) +
  geom_rect(data = rects_filtered,
            aes(ymin = ystart, ymax = yend, xmin = -Inf, xmax = Inf, fill = col),
            inherit.aes = FALSE, alpha = 0.5) +
  geom_hline(yintercept = seq(time_min, time_max, by = 500),  colour = "grey", linetype = "dashed") +
  geom_hline(yintercept = seq(time_min, time_max, by = 1000), colour = "grey") +
  geom_path(aes(group = Metric), linewidth = 1, color = "black", linetype = "dashed") +
  geom_point(size = 3, color = "black") +
  scale_fill_jco(guide = guide_legend(reverse = TRUE, ncol = 1), name = "Time Periods") +
  facet_wrap(~Metric, scales = "free_x", ncol = 3) +
  scale_y_reverse(breaks = time_breaks, labels = time_breaks) +
  coord_cartesian(ylim = c(time_max, time_min)) +
  labs(
    title = "Gamma Diversity (SRS)",
    subtitle = "Shannon (q = 1) and Simpson (q = 2)",
    x = "Hill number",
    y = "Years BP"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y   = element_text(size = 12, colour = "black"),
    axis.text.x   = element_text(size = 12, colour = "black", angle = 45, vjust = 0.5),
    axis.title.y  = element_text(size = 12, face = "bold", colour = "black"),
    axis.title.x  = element_text(size = 12, face = "bold", colour = "black"),
    strip.background = element_rect(fill = "#f0f0f0", color = NA),
    strip.text.x  = element_text(size = 14, face = "bold", colour = "black"),
    plot.title    = element_text(size = 16, face = "bold", colour = "black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "none"
  )


d18o_line <- d18o %>%
  dplyr::filter(!is.na(d18O_GISP2_per_mille), !is.na(age_BP_SEAD)) %>%
  dplyr::arrange(age_BP_SEAD)

p_d18 <- ggplot(d18o_line, aes(x = age_BP_SEAD, y = d18O_GISP2_per_mille)) +
  geom_line(linewidth = .5, color = "black", alpha = 1, aes(group = 1)) +
  geom_vline(xintercept = seq(time_min, time_max, by = 500),  colour = "grey", linetype = "dashed") +
  geom_vline(xintercept = seq(time_min, time_max, by = 1000), colour = "grey") +
  scale_x_reverse(breaks = time_breaks, labels = time_breaks) +
  coord_cartesian(xlim = c(time_max, time_min)) +
  coord_flip() +
  labs(title = "GISP2 δ\u00b9\u2078O (‰)", y = "δ\u00b9\u2078O (‰)", x = "") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y   = element_blank(),
    axis.text.x   = element_text(size = 12, colour = "black", angle = 45, vjust = 0.5),
    axis.title.y  = element_text(size = 12, face = "bold", colour = "black"),
    axis.title.x  = element_text(size = 12, face = "bold", colour = "black"),
    plot.title    = element_text(size = 16, face = "bold", colour = "black"),
    legend.position = "none"
  )


# Combine side-by-side on one row
combined <- p_div + p_d18 + plot_layout(ncol = 2, widths = c(1, .5))

# Save combined figure
ggsave(
  filename = "007-gamma-srs-vs-d18o-coleoptera.jpg",
  plot     = combined,
  path     = fig_dir,
  width = 30,
  height = 40,
  units = "cm",
  dpi = 300
)

message("✅ Combined SRS diversity vs δ18O saved as 005-gamma-srs-vs-d18o-coleoptera.jpg")



# ---- 12. Correlation: SRS Diversity (Shannon/Simpson) vs δ18O ----

# Harmonize analysis window (same as plots)
time_min <- 0
time_max <- 16000
edges    <- seq(time_min, time_max, by = 500)     # bin edges: 0, 500, ..., 16000
ends   <- edges[-1]  # 500, 1000, ..., 16000

# (A) Bin δ18O to 500-year bins (mean per bin_end)
d18o_bins <- d18o %>%
  dplyr::filter(age_BP_SEAD >= time_min, age_BP_SEAD <= time_max) %>%
  mutate(
    # bin index (1 for (0,500], 2 for (500,1000], ...)
    bin_idx = findInterval(age_BP_SEAD, vec = edges, rightmost.closed = TRUE),
    # map index to true bin_end values (500..16000)
    bin_end = ends[bin_idx]
  ) %>%
  group_by(bin_end) %>%
  summarise(
    d18O_mean = mean(d18O_GISP2_per_mille, na.rm = TRUE),
    d18O_sd   = sd(d18O_GISP2_per_mille,  na.rm = TRUE),
    n_obs     = sum(!is.na(d18O_GISP2_per_mille)),
    .groups   = "drop"
  )

# (B) Prepare SRS diversity (only Shannon & Simpson; restrict to 0–16 ka)
gamma_srs_plot <- gamma_srs %>%
  dplyr::filter(Metric %in% c("Richness", "Shannon", "Simpson"),
                Time >= time_min, Time <= time_max)

# (C) Join diversity with δ18O by bin_end (Time)
gamma_join <- gamma_srs_plot %>%
  inner_join(d18o_bins, by = c("Time" = "bin_end"))

# Quick diagnostic
message("Joined rows (per metric):")
print(gamma_join %>% count(Metric))

# (D) Correlation tests per metric
cor_results <- gamma_join %>%
  group_by(Metric) %>%
  summarise(
    n                 = n(),
    pearson_r         = cor(Value, d18O_mean, use = "complete.obs", method = "pearson"),
    pearson_p         = cor.test(Value, d18O_mean, method = "pearson")$p.value,
    spearman_rho      = cor(Value, d18O_mean, use = "complete.obs", method = "spearman"),
    spearman_p        = cor.test(Value, d18O_mean, method = "spearman")$p.value,
    .groups           = "drop"
  )

message("✅ Correlation results (SRS diversity vs δ18O):")
print(cor_results)

# (E) Save results
readr::write_csv(cor_results, file = file.path(fig_dir, "005-correlation-srs-d18o.csv"))
message("✅ Saved correlation table to 005-correlation-srs-d18o.csv")

# (F) Scatter plot with linear fit (one panel per metric)
p_corr <- ggplot(gamma_join, aes(x = d18O_mean, y = Value)) +
  geom_smooth(method = "lm", se = TRUE, color = "steelblue") +
  geom_point(aes(fill=Time), shape = 21, color = "black", size = 3, alpha = 0.8) +
  scale_fill_viridis(name = "Time (Years BP)", direction = -1, breaks=c(2000, 6000, 10000, 14000)) +
  facet_wrap(~Metric, scales = "free_y") +
  labs(
    title = "SRS Gamma Diversity v. GISP2 δ¹⁸O",
    x = "δ¹⁸O (‰)",
    y = "Hill number"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y   = element_text(size = 12, colour = "black"),
    axis.text.x   = element_text(size = 12, colour = "black"),
    strip.text    = element_text(size = 14, face = "bold"),
    plot.title    = element_text(size = 16, face = "bold"),
    legend.position = "bottom",
    legend.key.width= unit(1.5, 'cm')
  )

ggsave(
  "007-correlation-gamma-d18o.jpg",
  p_corr,
  device = "jpg",
  path = here("analysis", "figures"),
  width = 40,
  height = 40,
  units = "cm",
  dpi = 300
)

