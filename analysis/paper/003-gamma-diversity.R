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
  data.table, tidyverse, here, BiocManager, IRanges, SRS, entropart, ggsci, purrr, patchwork
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

# ==============================================================
# 9. Compute gamma diversity metrics for Raw and SRS
# ==============================================================

# Function for computing gamma diversity
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

# Compute gamma diversity (Raw and SRS)
gamma_raw <- compute_gamma(Listdf_raw, "Raw")
gamma_srs <- compute_gamma(Listdf_SRS, "SRS")

# ==============================================================
# 10. Prepare transform and dual-axis plotting functions
# ==============================================================

# 10.1 Long-format binder
prep_gamma_long <- function(gamma_raw, gamma_srs) {
  df_long <- dplyr::bind_rows(
    Raw = gamma_raw,
    SRS = gamma_srs,
    .id = "Type"
  ) %>%
    dplyr::mutate(
      Type   = factor(Type, levels = c("Raw", "SRS")),
      Metric = factor(Metric, levels = c("Shannon","Simpson"))
    )
  df_long
}

# 10.2 Per-metric linear transform function (SRS -> RAW (x_raw = a * x_srs + b) )
compute_transforms <- function(df_long) {
  sums <- df_long %>%
    dplyr::group_by(Metric, Type) %>%
    dplyr::summarise(
      min = min(Value, na.rm = TRUE),
      max = max(Value, na.rm = TRUE),
      .groups = "drop_last"
    ) %>%
    tidyr::pivot_wider(
      names_from  = Type,
      values_from = c(min, max),
      names_sep   = "_"
    ) %>%
    dplyr::ungroup()

  # Compute a, b for each Metric
  sums %>%
    dplyr::mutate(
      denom = pmax(max_SRS - min_SRS, 0),
      a = dplyr::if_else(denom > 0, (max_Raw - min_Raw) / denom, 1),
      b = min_Raw - a * min_SRS
    ) %>%
    dplyr::select(Metric, a, b)
}


# 10.3 Attach transform parameters and pre-compute Value_plot
attach_dual_transform <- function(df_long, tr_df) {
  df_long %>%
    dplyr::left_join(tr_df, by = "Metric") %>%
    dplyr::mutate(
      Value_plot = dplyr::if_else(Type == "SRS", a * Value + b, Value)
    )
}

# 10.4 Plot function (assumes Value_plot, a, b are already present)
plot_metric_dual_gamma <- function(df_metric, rects, a, b,
                                   base_title = "Alpha Diversity") {
  # Axis/time helpers (using df_metric span)
  time_min <- 0
  time_max <- 16000
  time_breaks <- seq(0, 16000, by = 1000)

  rects_filtered <- rects %>%
    dplyr::filter(ystart <= time_max, yend >= time_min)

  ggplot2::ggplot(
    df_metric,
    ggplot2::aes(x = Value_plot, y = Time, shape = Type, linetype = Type)
  ) +
    ggplot2::geom_rect(
      data = rects_filtered,
      ggplot2::aes(ymin = ystart, ymax = yend, xmin = -Inf, xmax = Inf, fill = col),
      inherit.aes = FALSE, alpha = 0.35
    ) +
    ggplot2::geom_path(linewidth = 0.5, color = "black") +
    ggplot2::geom_point(size = 2, color = "black") +
    ggsci::scale_fill_jco(name = "Time Periods", guide = ggplot2::guide_legend(reverse = TRUE)) +
    ggplot2::scale_shape_manual(values = c(Raw = 16, SRS = 17)) +
    ggplot2::scale_linetype_manual(values = c(Raw = "solid", SRS = "dashed")) +
    # Primary Raw axis (bottom), secondary SRS axis (top) using inverse transform
    ggplot2::scale_x_continuous(
      name = "Hill number",
      sec.axis = ggplot2::sec_axis(~ (.-b)/a, name = "Hill number (SRS)")
    ) +
    ggplot2::scale_y_reverse(breaks = time_breaks, labels = time_breaks, name = "Time (Years BP)") +
    ggplot2::coord_cartesian(ylim = c(time_max, time_min)) +
    ggplot2::labs(title = base_title) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_text(size = 12, colour = "black"),
      axis.text.x  = ggplot2::element_text(size = 11, colour = "black"),
      axis.title.y = ggplot2::element_text(size = 12, face = "bold", colour = "black"),
      axis.title.x = ggplot2::element_text(size = 12, face = "bold", colour = "black"),
      plot.title   = ggplot2::element_text(size = 16, face = "bold", colour = "black"),
      legend.title = ggplot2::element_text(size = 14, face = "bold", colour = "black"),
      legend.text  = ggplot2::element_text(size = 12, colour = "black"),
      legend.position = "right",
      axis.ticks.x.top = ggplot2::element_line(),   # show ticks on top axis
      axis.text.x.top  = ggplot2::element_text(size = 11, colour = "black")
    )
}

# 10.5 Compose the two panels (Shannon | Simpson) — return separate plots
plot_gamma_dual_horizontal <- function(gamma_raw, gamma_srs, rects,
                                       strip_y_right = TRUE,
                                       titles = c(Shannon = "Shannon", Simpson = "Simpson")) {
  df_long <- prep_gamma_long(gamma_raw, gamma_srs)
  tr      <- compute_transforms(df_long)
  df_pre  <- attach_dual_transform(df_long, tr)

  # Helper to build each panel
  build_panel <- function(metric, title_text, strip_y = FALSE) {
    d_m <- df_pre %>% dplyr::filter(Metric == metric)
    # a,b are constant per metric after compute_transforms()
    a <- d_m$a[1]; b <- d_m$b[1]

    p <- plot_metric_dual_gamma(
      df_metric  = d_m,
      rects      = rects,
      a          = a, b = b,
      base_title = title_text
    )

    if (strip_y) {
      p <- p + ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        axis.text.y  = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank()
      )
    }
    p
  }

  p_shan <- build_panel("Shannon", titles[["Shannon"]], strip_y = FALSE)
  p_simp <- build_panel("Simpson", titles[["Simpson"]], strip_y = isTRUE(strip_y_right))

  # Return both plots so you can customize them later
  # (a named list is convenient for downstream use)
  list(
    p_shan = p_shan,
    p_simp = p_simp
  )
}

# ==============================================================
# 11. Generate plot
# ==============================================================

# Generate dual-axis gamma diversity plots
plots <- plot_gamma_dual_horizontal(
  gamma_raw = gamma_raw,
  gamma_srs = gamma_srs,
  rects     = rects
)

# Extract individual plots to be called in manuscript figure
p_shan <- plots$p_shan
p_simp <- plots$p_simp

# Compose dual-figure layout with shared legend and title
p_all <- (p_shan | p_simp) +
  patchwork::plot_layout(guides = "collect") +
  plot_annotation(title = "Gamma Diversity") &
  ggplot2::theme(
    legend.position = "right",
    plot.title = ggplot2::element_text(size = 16, face = "bold"),
  )

# ==============================================================
# 12. Save plots
# ==============================================================

# Save individual plots as RDS for manuscript figure assembly
saveRDS(p_shan, file = here::here("analysis", "data", "derived_data", "gamma_shan.rds"))
saveRDS(p_simp, file = here::here("analysis", "data", "derived_data", "gamma_simp.rds"))

# Save dual-plot
ggsave(
  filename = "003-gamma-diversity.jpg",
  plot = p_all,
  path = here::here("analysis", "figures"),
  width = 3900, height = 4900, units = "px", dpi = 500
)

message("✅ Dual-axis gamma diversity (horizontal) plot saved: 003-gamma-diversity.jpg")


