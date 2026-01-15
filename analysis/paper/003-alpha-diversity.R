# ==============================================================
# Script: Alpha Diversity Analysis
#
# Purpose:
#   Compute alpha diversity metrics (Hill numbers: Richness, Shannon, Simpson)
#   for Lateglacial-Holocene subfossil insect data from the UK/Ireland,
#   using raw abundances and SRS (Scaling with Ranked Subsampling) standardized data
#
# Workflow:
#   1. Load fossil insect occurrence data.
#   2. Filter samples by region, context, and age constraints.
#   3. Assign samples to 500-year temporal bins using IRanges.
#   4. Standardize abundances using SRS:
#        - Deduplicate taxon × sample combinations.
#        - Remove taxa with zero abundance and samples below Cmin threshold.
#        - Apply SRS to equalize sampling effort across samples.
#   5. Build lists of community matrices per time bin (taxon × sample).
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
#   - Plot (SRS): "supplementary-S3-alpha-diversity-srs.jpg"
#   - Plot (RAW): "supplementary-S3-alpha-diversity-raw.jpg"
# ==============================================================

# ---- 0. Load packages ----
pacman::p_load(
  data.table, tidyverse, here, IRanges, SRS, entropart, ggsci, purrr, patchwork
)

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
                        "HETEROPTERA", "HOMOPTERA", "HYDROPSYCHIDAE", "HYDROPTILIDAE", "HYMENOPTERA", "LEPIDOPTERA", "LEPIDOSTOMATIDAE",
                        "LEPTOCERIDAE", "LIMNEPHILIDAE", "LYGAEIDAE", "MEMBRACIDAE", "MICROPHYSIDAE", "MICROSPORIDAE", "MIRIDAE",
                        "MOLANNIDAE", "NABIDAE", "NEMATOCERA", "NEMOURIDAE", "ODONATA", "PARASITICA", "PENTATOMIDAE", "PHRYGANEIDAE",
                        "POLYCENTROPIDAE", "PSYCHOMYIIDAE", "PSYLLIDAE", "RAPHIDIIDAE", "SALDIDAE", "SCUTELLERIDAE", "SERICOSTOMATIDAE",
                        "SIALIDAE", "TRICHOPTERA", "THYREOCORIDAE", "TINGIDAE", "TIPULIDAE", "TRIOPSIDAE", "TRIOZIDAE"),
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
# 4. Define 500-year temporal bins
# ==============================================================
range <- tibble(
  start = c(-500, seq(1, 15501, by = 500)),
  end   = seq(0, 16000, by = 500)
)

# ==============================================================
# 5. Identify overlapping samples with time bins
# ==============================================================
query   <- IRanges(start = time_mat$start, end = time_mat$end)
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
# 7. Prepare raw matrix for SRS and raw alpha
# ==============================================================
raw_mat <- hits %>%
  select(sample_id, bin_end, taxon, abundance) %>%
  distinct()

# Deduplicate taxon × sample_id for both SRS and RAW
srs_agg <- raw_mat %>%
  group_by(taxon, sample_id) %>%
  summarise(abundance = sum(abundance), .groups = "drop")

# ---------- SRS standardization ----------
# Pivot to taxa × samples
srs_prep <- srs_agg %>%
  pivot_wider(names_from = sample_id, values_from = abundance, values_fill = 0) %>%
  column_to_rownames("taxon")

# Drop taxa with all zeros
Cmin <- 10
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

# Convert SRS result to long format and attach bin info
srs_long <- as_tibble(srs_out, rownames = "taxon") %>%
  pivot_longer(-taxon, names_to = "sample_id", values_to = "abundance") %>%
  filter(abundance > 0) %>%
  inner_join(raw_mat %>% select(sample_id, bin_end), by = "sample_id", relationship = "many-to-many") %>%
  distinct()

# ==============================================================
# 8. Build list of matrices
# ==============================================================
build_Listdf <- function(df) {
  split(df, df$bin_end) %>%
    purrr::map(~ .x %>%
                 select(taxon, sample_id, abundance) %>%
                 pivot_wider(names_from = sample_id, values_from = abundance, values_fill = 0) %>%
                 column_to_rownames("taxon") %>%
                 as.matrix())
}

Listdf_SRS <- build_Listdf(srs_long)
Listdf_RAW <- build_Listdf(raw_mat)

# ==============================================================
# 9. Compute alpha diversity
# ==============================================================
compute_alpha <- function(list_of_mats) {
  purrr::map_dfr(names(list_of_mats), function(time_name) {
    mat <- list_of_mats[[time_name]]

    # Initial counts
    initial_taxa    <- nrow(mat)
    initial_samples <- ncol(mat)

    # Drop empty taxa/samples
    mat <- mat[rowSums(mat) > 0, , drop = FALSE]
    mat <- mat[, colSums(mat) > 0, drop = FALSE]

    dropped_taxa    <- initial_taxa - nrow(mat)
    dropped_samples <- initial_samples - ncol(mat)

    message("Bin ", time_name, ": Taxa retained = ", nrow(mat),
            " (Dropped: ", dropped_taxa, "), Samples retained = ", ncol(mat),
            " (Dropped: ", dropped_samples, ")")

    # Skip bins with too few taxa/samples
    if (nrow(mat) < 2 || ncol(mat) < 2) {
      message("⚠️ Bin ", time_name, " skipped (too few taxa or samples).")
      return(NULL)
    }

    # Compute alpha diversity with abundance weighting (entropart)
    Weights <- colSums(mat) / sum(colSums(mat))
    MC <- MetaCommunity(mat, Weights)

    tibble(
      Time     = as.numeric(time_name),
      Richness = AlphaDiversity(MC, q = 0, Correction = "UnveilJ")$Total,
      Shannon  = AlphaDiversity(MC, q = 1, Correction = "UnveilJ")$Total,
      Simpson  = AlphaDiversity(MC, q = 2, Correction = "UnveilJ")$Total
    )
  })
}

alpha_srs <- compute_alpha(Listdf_SRS) %>%
  mutate(Type = "SRS")

alpha_raw <- compute_alpha(Listdf_RAW) %>%
  mutate(Type = "RAW")

# ==============================================================
# 10. Prepare transform and dual-axis plotting functions
# ==============================================================

# 10.1 Long-format binder
prep_alpha_long <- function(alpha_raw, alpha_srs) {
  raw_l <- alpha_raw %>%
    pivot_longer(cols = c(Richness, Shannon, Simpson),
                 names_to = "Metric", values_to = "Value") %>%
    mutate(Type = "RAW")

  srs_l <- alpha_srs %>%
    pivot_longer(cols = c(Richness, Shannon, Simpson),
                 names_to = "Metric", values_to = "Value") %>%
    mutate(Type = "SRS")

  bind_rows(raw_l, srs_l) %>%
    mutate(Metric = factor(Metric, levels = c("Richness","Shannon","Simpson")))
}

# 10.2 Per-metric linear transform function (SRS -> RAW (x_raw = a * x_srs + b) )
compute_transforms <- function(df_long) {
  df_long %>%
    group_by(Metric) %>%
    summarise(
      min_raw = min(Value[Type == "RAW"], na.rm = TRUE),
      max_raw = max(Value[Type == "RAW"], na.rm = TRUE),
      min_srs = min(Value[Type == "SRS"], na.rm = TRUE),
      max_srs = max(Value[Type == "SRS"], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      # handle edge cases (constant ranges) safely
      a = ifelse(max_srs > min_srs,
                 (max_raw - min_raw) / (max_srs - min_srs),
                 1),
      b = min_raw - a * min_srs
    )
}

# ==============================================================
# 11. Plotting helper
# ==============================================================
plot_metric_dual <- function(df_long, metric_name, rects) {
  # filter to the metric
  d <- df_long %>% filter(Metric == metric_name)

  # compute transform for this metric
  tr <- compute_transforms(d) %>% filter(Metric == metric_name)
  a <- tr$a[1]; b <- tr$b[1]

  # transform SRS values into RAW x-space for plotting
  d_plot <- d %>%
    mutate(Value_plot = ifelse(Type == "SRS", a * Value + b, Value))

  # axis range on RAW scale
  x_min <- min(d_plot$Value_plot, na.rm = TRUE)
  x_max <- max(d_plot$Value_plot, na.rm = TRUE)
  time_min <- min(d_plot$Time, na.rm = TRUE)
  time_max <- max(d_plot$Time, na.rm = TRUE)
  rects_filtered <- rects %>% filter(ystart <= time_max, yend >= time_min)
  time_breaks <- pretty(d_plot$Time, n = 10)

  # Build plot
  ggplot(d_plot, aes(x = Value_plot, y = Time, shape = Type, color = Type, alpha = Type)) +
    geom_rect(
      data = rects_filtered,
      aes(ymin = ystart, ymax = yend, xmin = -Inf, xmax = Inf, fill = col),
      inherit.aes = FALSE, alpha = 0.35
    ) +
    geom_path(aes(group = Type, linetype = Type), linewidth = 0.5) +
    geom_point(size = 2) +
    scale_fill_jco(name = "Time Periods", guide = guide_legend(reverse = TRUE)) +
    scale_color_manual(values = c(RAW = "black", SRS = "black")) +
    scale_shape_manual(values = c(RAW = 16, SRS = 17)) +
    scale_linetype_manual(values = c(RAW = "solid", SRS = "dashed")) +
    scale_alpha_manual(name = "Type", values = c(1, .5)) +
    # Primary RAW axis (bottom), secondary SRS axis (top) using inverse transform
    scale_x_continuous(
      name = "Hill number",
      sec.axis = sec_axis(~ (.-b)/a, name = "Hill number (SRS)")
    ) +
    scale_y_reverse(breaks = time_breaks, labels = time_breaks, name = "Time (Years BP)") +
    coord_cartesian(ylim = c(time_max, time_min)) +
    labs(title = paste("Alpha Diversity:", metric_name)) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y  = element_text(size = 12, colour = "black"),
      axis.text.x  = element_text(size = 11, colour = "black"),
      axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
      axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
      plot.title   = element_text(size = 16, face = "bold", colour = "black"),
      legend.title = element_text(size = 14, face = "bold", colour = "black"),
      legend.text  = element_text(size = 12, colour = "black"),
      legend.position = "right",
      axis.ticks.x.top = element_line(),   # show ticks on top axis
      axis.text.x.top  = element_text(size = 11, colour = "black")
    )
}

# ==============================================================
# 11. Function to generate and save dual-axis alpha diversity plot
# ==============================================================

plot_alpha_dual <- function(alpha_raw, alpha_srs, rects, filename) {
  df_long <- prep_alpha_long(alpha_raw, alpha_srs)

  #p_rich <- plot_metric_dual(df_long, "Richness", rects) +
  #  labs(title = "Richness")
  p_shan <- plot_metric_dual(df_long, "Shannon", rects) +
    labs(title = "Shannon")
  p_simp <- plot_metric_dual(df_long, "Simpson", rects) +
    labs(title = "Simpson") +
    theme(
      axis.title.y = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank()
    )

  # Arrange horizontally, collect legend once, and harmonize plot margins
  #p_all <- (p_rich | p_shan | p_simp) +
  p_all <- (p_shan | p_simp) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right",
          plot.title = element_text(size = 16, face = "bold"))

  # Add a global title/subtitle
  p_all <- p_all + plot_annotation(
    title = "Alpha Diversity"
  )

  # Save as a wide, landscape image
  ggsave(
    filename = filename,
    plot = p_all,
    path = here("analysis", "figures"),
    width = 3900, height = 4900, units = "px", dpi = 500
  )
  message("✅ Dual-axis alpha diversity (horizontal) plot saved: ", filename)
}

# ==============================================================
# 11. Generate and save figure
# ==============================================================
plot_alpha_dual(
  alpha_raw  = alpha_raw,
  alpha_srs  = alpha_srs,
  rects      = rects,
  filename   = "003-alpha-diversity.jpg"
)

message("✅ Alpha diversity generated and figure saved: '003-alpha-diversity.jpg'")
