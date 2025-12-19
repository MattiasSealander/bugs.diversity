# ==============================================================
# Script: Gamma Diversity Analysis (Raw vs SRS)
#
# Purpose:
#   Compute gamma diversity metrics (Hill numbers: Richness, Shannon, Simpson)
#   for Holocene fossil insect data from British/Irish Isles,
#   comparing raw abundances and SRS-standardized abundances.
#
# Workflow:
#   1. Load fossil insect occurrence data.
#   2. Filter samples by region, context, and age constraints.
#   3. Assign samples to 500-year temporal bins using IRanges.
#   4. Prepare raw abundance matrix.
#   5. Standardize abundances using SRS.
#   6. Build list of matrices per bin for raw and SRS.
#   7. Compute gamma diversity metrics for each bin using entropart.
#   8. Plot results comparing raw vs SRS.
#
# Key Notes:
#   - Uses IRanges for robust interval overlap detection.
#   - Diagnostic messages report taxa and samples dropped during SRS prep
#     and during bin-level filtering.
#   - Bins with <2 taxa or <2 samples are skipped.
#
# Output:
#   - Plot: "004-gamma-diversity-raw-srs.jpg"
# ==============================================================

# ---- 0. Load packages ----
pacman::p_load(
  data.table, tidyverse, here, BiocManager, IRanges, SRS, entropart, ggsci, patchwork
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
# 8. SRS standardization
# ==============================================================
Cmin <- 10
SRS_seed <- 1988
message("Running SRS standardization...")

# Aggregate abundances
srs_agg <- raw_mat %>%
  group_by(taxon, sample_id) %>%
  summarise(abundance = sum(abundance), .groups = "drop")

# Pivot to taxa × samples
srs_prep <- srs_agg %>%
  pivot_wider(names_from = sample_id, values_from = abundance, values_fill = 0) %>%
  column_to_rownames("taxon")

# Drop taxa and samples with diagnostics
initial_taxa <- nrow(srs_prep)
initial_samples <- ncol(srs_prep)
srs_prep <- srs_prep[rowSums(srs_prep) > 0, , drop = FALSE]
message("✅ Taxa retained after zero-filter: ", nrow(srs_prep), " (Dropped: ", initial_taxa - nrow(srs_prep), ")")
srs_prep <- srs_prep[, colSums(srs_prep) >= Cmin, drop = FALSE]
message("✅ Samples retained after Cmin filter: ", ncol(srs_prep), " (Dropped: ", initial_samples - ncol(srs_prep), ")")

# Perform SRS
set.seed(SRS_seed)
srs_out <- SRS(as.data.frame(srs_prep), Cmin = Cmin, set_seed = TRUE, seed = SRS_seed)
rownames(srs_out) <- rownames(srs_prep)

# Convert to long format and merge bin info
srs_long <- as_tibble(srs_out, rownames = "taxon") %>%
  pivot_longer(-taxon, names_to = "sample_id", values_to = "abundance") %>%
  filter(abundance > 0) %>%
  inner_join(raw_mat %>% select(sample_id, bin_end), by = "sample_id", relationship = "many-to-many") %>%
  distinct()

# ---- 9. Build lists of matrices per bin ----
build_Listdf <- function(df) {
  split(df, df$bin_end) %>%
    map(~ .x %>%
          select(taxon, sample_id, abundance) %>%
          pivot_wider(names_from = sample_id, values_from = abundance, values_fill = 0) %>%
          column_to_rownames("taxon") %>%
          as.matrix())
}

Listdf_raw <- build_Listdf(raw_mat)
Listdf_SRS <- build_Listdf(srs_long)

# ---- 10. Compute gamma diversity for Raw and SRS ----


compute_gamma <- function(data_list, dataset_type) {
  purrr::map_dfr(names(data_list), function(time_name) {
    mat <- data_list[[time_name]]

    # Initial counts and filtering
    initial_taxa    <- nrow(mat)
    initial_samples <- ncol(mat)
    mat <- mat[rowSums(mat) > 0, , drop = FALSE]
    mat <- mat[, colSums(mat) > 0, drop = FALSE]
    dropped_taxa    <- initial_taxa    - nrow(mat)
    dropped_samples <- initial_samples - ncol(mat)

    message("[", dataset_type, "] Bin ", time_name, ": Taxa retained = ", nrow(mat),
            " (Dropped: ", dropped_taxa, "), Samples retained = ", ncol(mat),
            " (Dropped: ", dropped_samples, ")")

    if (nrow(mat) < 2 || ncol(mat) < 2) {
      message("⚠️ [", dataset_type, "] Bin ", time_name, " skipped (too few taxa or samples).")
      return(NULL)
    }

    # Gamma diversity with correction (same as your script)
    MC <- entropart::MetaCommunity(mat)

    tibble::tibble(
      Time  = as.numeric(time_name),
      Metric = c("Shannon", "Simpson"),
      Value  = c(
        entropart::GammaDiversity(MC, q = 1, Correction = "UnveiliC"),
        entropart::GammaDiversity(MC, q = 2, Correction = "UnveiliC")
      ),
      Type = dataset_type
    )
  })
}


gamma_raw <- compute_gamma(Listdf_raw, "Raw")
gamma_srs <- compute_gamma(Listdf_SRS, "SRS")


# Long-format binder (only Shannon & Simpson)
prep_gamma_long <- function(gamma_raw, gamma_srs) {
  raw_l <- gamma_raw %>%
    dplyr::mutate(Metric = factor(Metric, levels = c("Shannon","Simpson"))) %>%
    dplyr::mutate(Type = "RAW")

  srs_l <- gamma_srs %>%
    dplyr::mutate(Metric = factor(Metric, levels = c("Shannon","Simpson"))) %>%
    dplyr::mutate(Type = "SRS")

  dplyr::bind_rows(raw_l, srs_l)
}

# Per-metric linear transform mapping SRS -> RAW  (x_raw = a * x_srs + b)
compute_transforms <- function(df_long) {
  df_long %>%
    dplyr::group_by(Metric) %>%
    dplyr::summarise(
      min_raw = min(Value[Type == "RAW"], na.rm = TRUE),
      max_raw = max(Value[Type == "RAW"], na.rm = TRUE),
      min_srs = min(Value[Type == "SRS"], na.rm = TRUE),
      max_srs = max(Value[Type == "SRS"], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      a = ifelse(max_srs > min_srs, (max_raw - min_raw)/(max_srs - min_srs), 1),
      b = min_raw - a * min_srs
    )
}


# ---- 11. Plot function ----
plot_gamma <- function(gamma_data, dataset_type, filename) {
  gamma_data <- gamma_data %>%
    mutate(Metric = factor(Metric, levels = c("Richness", "Shannon", "Simpson")))

  plot_min <- min(gamma_data$Time)
  plot_max <- max(gamma_data$Time)
  rects_filtered <- rects %>% filter(ystart <= plot_max, yend >= plot_min)
  time_breaks <- pretty(gamma_data$Time, n = 10)

  p <- ggplot(gamma_data, aes(x = Value, y = Time)) +
    geom_rect(data = rects_filtered,
              aes(ymin = ystart, ymax = yend, xmin = -Inf, xmax = Inf, fill = col),
              inherit.aes = FALSE, alpha = 0.5) +
    geom_path(aes(group = Metric), linewidth = 1, color = "black", linetype = "dashed") +
    geom_point(size = 3, color = "black") +
    facet_wrap(~Metric, scales = "free_x") +
    scale_fill_jco(guide = guide_legend(reverse = TRUE)) +
    scale_y_reverse(breaks = time_breaks, labels = time_breaks) +
    labs(title = "Gamma Diversity", subtitle = paste(dataset_type, "values"),
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

  ggsave(filename = filename,
         plot = p, path = here("analysis", "figures"),
         width = 3300, height = 4200, units = "px", dpi = 300)

  message("✅ Gamma diversity plot for ", dataset_type, " saved as ", filename)
}

# ---- 12. Generate and save plots ----
plot_gamma(gamma_raw, "Raw (unstandardized) abundances", "003-gamma-diversity-raw.jpg")
plot_gamma(gamma_srs, "SRS standardized abundances", "003-gamma-diversity-srs.jpg")



plot_metric_dual_gamma <- function(df_long, metric_name, rects) {
  d <- df_long %>% dplyr::filter(Metric == metric_name)

  tr <- compute_transforms(d) %>% dplyr::filter(Metric == metric_name)
  a <- tr$a[1]; b <- tr$b[1]

  d_plot <- d %>%
    dplyr::mutate(Value_plot = ifelse(Type == "SRS", a * Value + b, Value))

  # Axis/time helpers (use your rects/time span)
  time_min <- min(d_plot$Time, na.rm = TRUE)
  time_max <- max(d_plot$Time, na.rm = TRUE)
  rects_filtered <- rects %>% dplyr::filter(ystart <= time_max, yend >= time_min)
  time_breaks <- pretty(d_plot$Time, n = 10)

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


plot_gamma_dual_horizontal <- function(gamma_raw, gamma_srs, rects, filename) {
  df_long <- prep_gamma_long(gamma_raw, gamma_srs)

  p_shan <- plot_metric_dual_gamma(df_long, "Shannon", rects) + ggplot2::labs(title = "Shannon")
  p_simp <- plot_metric_dual_gamma(df_long, "Simpson", rects) + ggplot2::labs(title = "Simpson") +
    ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )

  p_all <- (p_shan | p_simp) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "right",
                   plot.title = ggplot2::element_text(size = 16, face = "bold"))

  p_all <- p_all + patchwork::plot_annotation(
    title = "Gamma Diversity",
  )

  ggplot2::ggsave(
    filename = filename,
    plot = p_all,
    path = here::here("analysis", "figures"),
    width = 3900, height = 4900, units = "px", dpi = 500
  )
  message("✅ Dual-axis gamma diversity (horizontal) plot saved: ", filename)
}



# 12. Generate and save the combined figure (Shannon | Simpson)
gamma_raw <- compute_gamma(Listdf_raw, "RAW")
gamma_srs <- compute_gamma(Listdf_SRS, "SRS")

plot_gamma_dual_horizontal(
  gamma_raw = gamma_raw,
  gamma_srs = gamma_srs,
  rects     = rects,
  filename  = "004-gamma-diversity-RAW-vs-SRS-dual-axis-horizontal.jpg"
)








