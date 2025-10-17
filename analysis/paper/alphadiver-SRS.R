# ==============================================================
# Description:
# This script performs alpha diversity estimation on fossil insect data
# extracted from the BugsCEP database. It applies the SRS (Scaling with Ranked
# Subsampling) method to standardize abundance matrices,
# then computes alpha diversity metrics (Observed richness, Coverage,
# Shannon, and Inverse Simpson) for each 500-year time slice.
#
# Workflow overview:
#   1. Load and clean fossil insect occurrence data.
#   2. Define 500-year temporal bins (16 ka – present).
#   3. Assign each sample to one or more bins based on its dated range.
#   4. Scale sample abundance data using SRS.
#   5. Construct community (taxon × sample) matrices per time bin.
#   6. Compute alpha diversity metrics using entropart.
#   7. Produce time-series plots for each metric.
#
# Outputs:
#   - Alpha diversity summary table (data.frame)
#   - Multi-panel ggplot of alpha metrics over time
# ==============================================================


# ---- 0. Load required packages ----
pacman::p_load(
  cowplot, data.table, entropart, ggh4x, ggplot2,
  ggsci, tidyverse, IRanges, here, SRS
)


# ---- 1. Import data ----
# Import fossil insect data extracted from BugsCEP (Europe)
bugs <- fread(
  here::here("analysis/data/raw_data/bugs_europe_extraction_samples_20250612.csv"),
  na.strings = c("", "NA", "NULL"),
  encoding = "UTF-8"
)


# ---- 2. Define temporal periods (after Pilotto et al. 2022) ----
# These intervals define key Holocene time slices for plotting background rectangles.
rects <- data.frame(
  ystart = c(12000, 9500, 8000, 6500, 4500, 4000, 3500, -500),
  yend   = c(16000, 12000, 9500, 8000, 6500, 4500, 4000, 3500),
  col    = factor(
    c("TM 1", "TM 2", "TM 3", "TM 4", "TM 5", "TM 6", "TM 7", "TM 8"),
    levels = c("TM 1", "TM 2", "TM 3", "TM 4", "TM 5", "TM 6", "TM 7", "TM 8")
  )
)


# ---- 3. Prepare time intervals ----
# Filter to focus on samples from natural contexts (Stratigraphic sequence) and Late Glacial - Holocene period
time_mat <- bugs %>%
  select(country, sample, site, sample_group, age_older, age_younger, context) %>%
  mutate(age_range = age_older - age_younger) %>%
  filter(
    context == "Stratigraphic sequence",
    sample != "BugsPresence",
    between(age_older, -500, 16000),
    between(age_younger, -500, 16000),
    age_range <= 2000,
    country != "Greenland"
  ) %>%
  distinct() %>%
  select(-age_range, -context) %>%
  mutate(sample = paste(sample, sample_group, site, sep = "@")) %>%
  column_to_rownames("sample") %>%
  select(-country, -site, -sample_group) %>%
  dplyr::rename(start = age_younger, end = age_older)


# ---- 4. Define 500-year temporal bins ----
range <- data.frame(
  start = c(-500, seq(1, 15501, by = 500)),
  end   = seq(0, 16000, by = 500)
)


# ---- 5. Identify sample–bin overlaps ----
# Samples overlapping each time bin are assigned accordingly.
intersection <- findOverlaps(
  query = do.call(IRanges, time_mat),
  subject = do.call(IRanges, range),
  type = "any"
)


# ---- 6. Merge sample and bin info ----
# Combine sample information with corresponding time bins.
hits <- data.frame(time_mat[queryHits(intersection),], range[subjectHits(intersection),]) %>%
  as_tibble(rownames = "sample") %>%
  separate(sample, into = c("sample", "sample_group", "site"), sep = "@") %>%
  inner_join(bugs, by = c("sample", "sample_group"), relationship = "many-to-many") %>%
  select(-start, -end)


# ---- 7. Filter to region and construct sample IDs ----
region_name_focus <- "British/Irish Isles"

raw_region <- hits %>%
  select(site.x, latitude, longitude, sample_group, sample, end.1, taxon, abundance) %>%
  distinct() %>%
  mutate(sample_id = paste(sample, sample_group, site.x, sep = "@")) %>%  # exclude end.1
  filter(
    between(latitude, 49.8, 62.6),
    between(longitude, -12.6, 1.8)
  ) %>%
  mutate(region = region_name_focus) %>%
  select(-latitude, -longitude)

if (nrow(raw_region) == 0) stop("No records found for region: ", region_name_focus)

# ---- 8. Deduplicate taxon × sample rows safely (do NOT sum) ----
# For each taxon x sample_id:
#  - if duplicates have identical abundance -> keep one row
#  - if duplicates disagree on abundance    -> raise an error for inspection

# 8a. Count duplicates and collect unique abundance values per group
agg_check <- raw_region %>%
  group_by(taxon, sample_id) %>%
  summarise(
    n_rows = n(),
    unique_abundances = list(sort(unique(abundance))),
    .groups = "drop"
  )

# 8b. Find problematic groups where there are >1 distinct abundance values
problems <- agg_check %>%
  filter(sapply(unique_abundances, length) > 1)

if (nrow(problems) > 0) {
  # Provide a concise diagnostic to inspect the first few problematic groups
  message("ERROR: Found taxon x sample_id groups with conflicting abundance values.\n",
          "Please inspect these groups in the raw data. Showing up to 10 examples:\n")
  print(head(problems, 10))
  stop("Abundance conflicts detected for taxon x sample_id groups. Resolve before continuing.")
}

# 8c. If there are duplicates but they all agree, we deduplicate safely;
#    otherwise we have already stopped above.
srs_agg <- raw_region %>%
  distinct(taxon, sample_id, abundance) %>%   # collapse identical rows to a single row
  select(taxon, sample_id, abundance)

# ---- 9. Prepare abundance matrix (taxa × samples) ----
srs_prep <- srs_agg %>%
  pivot_wider(names_from = sample_id, values_from = abundance, values_fill = 0) %>%
  column_to_rownames("taxon")

# Drop taxa with all zeros
srs_prep <- srs_prep[rowSums(srs_prep) > 0, , drop = FALSE]

# ---- 10. Perform SRS standardization ----
sample_totals <- colSums(srs_prep)
Cmin <- 10  # minimum total abundance threshold (change if needed)
keep_samples <- names(sample_totals)[sample_totals >= Cmin]
srs_prep <- srs_prep[, keep_samples, drop = FALSE]

SRS_seed <- 1988
set.seed(SRS_seed)
srs_out <- SRS(
  data = as.data.frame(srs_prep),
  Cmin = Cmin,
  set_seed = TRUE,
  seed = SRS_seed
)

rownames(srs_out) <- rownames(srs_prep)

# ---- 11. Add metadata back and build per–time-slice abundance tables ----
# Map SRS-scaled samples to their time slices
sample_time_lookup <- raw_region %>%
  distinct(sample_id, end.1)

# Convert to long format for merging
srs_long <- as.data.table(srs_out, keep.rownames = "taxon") %>%
  melt(id.vars = "taxon", variable.name = "sample_id", value.name = "abundance") %>%
  left_join(sample_time_lookup, by = "sample_id") %>%
  mutate(region = region_name_focus)

# ---- Build community matrices(flat: time slice -> taxon × sample) ----
# Each matrix: rows = taxa, columns = samples.
Listdf <- list()
et_order <- sort(unique(srs_long$end.1))

for (et in et_order) {
  dt_slice <- srs_long[end.1 == et]
  if (nrow(dt_slice) == 0) next

  tab <- dcast(dt_slice, taxon ~ sample_id, value.var = "abundance", fill = 0)
  setDF(tab)
  rownames(tab) <- tab$taxon
  tab$taxon <- NULL

  tab <- tab[rowSums(tab) > 0, , drop = FALSE]
  tab <- tab[, colSums(tab) > 0, drop = FALSE]

  if (nrow(tab) > 0 && ncol(tab) > 0)
    Listdf[[as.character(et)]] <- tab
}

message("✅ Listdf constructed for ", region_name_focus,
        " with ", length(Listdf), " valid time slices.")


# ==============================================================
# Function: prepare_for_MetaCommunity
# ==============================================================
# Prepares a taxon × sample abundance matrix for use with entropart package.
# Ensures numeric-only data, removes empty taxa/samples,
# and constructs a MetaCommunity object safely.
# ==============================================================

prepare_for_MetaCommunity <- function(df, verbose = FALSE) {
  if ("taxon" %in% names(df)) df <- tibble::column_to_rownames(df, "taxon")
  df_num <- suppressWarnings(as.data.frame(lapply(df, as.numeric)))
  df_num[is.na(df_num)] <- 0

  df_num <- df_num[rowSums(df_num) > 0, , drop = FALSE]
  df_num <- df_num[, colSums(df_num) > 0, drop = FALSE]

  if (nrow(df_num) < 2 || ncol(df_num) < 2) return(NULL)

  Abundances <- as.matrix(df_num)
  Weights <- colSums(Abundances)
  if (sum(Weights) == 0) return(NULL)
  Weights <- Weights / sum(Weights)

  tryCatch(
    entropart::MetaCommunity(Abundances, Weights),
    error = function(e) {
      if (verbose) message("⚠️ MetaCommunity failed: ", e$message)
      NULL
    }
  )
}


# ==============================================================
# Function: generate_alphaDiversity_SRS
# ==============================================================
# Description:
#   Computes alpha diversity metrics (Observed Richness, Coverage,
#   Shannon, Inverse Simpson) for SRS-scaled abundance data.
#
# Arguments:
#   Listdf: list of taxon × sample matrices per time slice
#   rects:  data frame of temporal intervals for background shading
#   output_dir: directory to save figure
#
# Returns:
#   List containing:
#     - $plot: ggplot object of alpha diversity metrics
#     - $results: tidy data.frame of metric values per time slice
# ==============================================================

generate_alphaDiversity_SRS <- function(Listdf, rects,
                                        output_dir = here::here("analysis", "figures")) {

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  results <- data.frame(Time = character(), Metric = character(), Value = numeric())

  for (time_name in names(Listdf)) {
    comm_df <- Listdf[[time_name]]
    if (is.null(comm_df) || nrow(comm_df) == 0) next

    # ---- Ensure numeric input ----
    comm_num <- as.data.frame(comm_df)
    comm_num <- comm_num[sapply(comm_num, is.numeric)]
    comm_num[is.na(comm_num)] <- 0

    # ---- Drop empty taxa/samples ----
    comm_num <- comm_num[rowSums(comm_num) > 0, , drop = FALSE]
    comm_num <- comm_num[, colSums(comm_num) > 0, drop = FALSE]

    # 🟡 Filtering rule:
    # Skip bins that have only one taxon or one sample,
    # since diversity indices are not meaningful in those cases.
    if (nrow(comm_num) < 2 || ncol(comm_num) < 2) next

    # ---- Compute diversity metrics ----
    obs_rich <- sum(rowSums(comm_num) > 0)
    cov <- tryCatch(entropart::Coverage(as.matrix(comm_num)), error = function(e) NA)

    # ---- Compute abundance-weighted meta-community ----
    # 🟢 Note:
    #   Weights are based on each sample total abundance
    Weights <- colSums(comm_num)
    Weights <- Weights / sum(Weights)
    MC <- tryCatch(entropart::MetaCommunity(as.matrix(comm_num), Weights), error = function(e) NULL)

    shannon_div <- simpson_div <- NA
    if (!is.null(MC)) {
      shannon_div <- tryCatch(entropart::AlphaDiversity(MC, q = 1, Correction = "Best")$Total,
                              error = function(e) NA)
      simpson_div <- tryCatch(entropart::AlphaDiversity(MC, q = 2, Correction = "Best")$Total,
                              error = function(e) NA)
    }

    # ---- Append results ----
    add_metric <- function(name, value) {
      if (!is.na(value)) {
        results <<- rbind(results,
                          data.frame(Time = time_name, Metric = name, Value = value))
      }
    }
    add_metric("Observed Richness", obs_rich)
    add_metric("Coverage", cov)
    add_metric("Shannon", shannon_div)
    add_metric("Inverse Simpson", simpson_div)
  }

  if (nrow(results) == 0) {
    message("No valid time slices found for alpha diversity.")
    return(NULL)
  }

  # ---- Prepare for plotting ----
  results$TimeNumeric <- suppressWarnings(as.numeric(results$Time))
  results$Metric <- factor(results$Metric,
                           levels = c("Coverage", "Observed Richness", "Shannon", "Inverse Simpson"))

  plot_min <- min(results$TimeNumeric, na.rm = TRUE) - 500
  plot_max <- max(results$TimeNumeric, na.rm = TRUE) + 500
  rects_filtered <- rects[rects$ystart <= plot_max & rects$yend >= plot_min, ]
  time_breaks <- pretty(results$TimeNumeric, n = 10)

  # ---- Plot ----
  p <- ggplot() +
    geom_rect(data = rects_filtered,
              aes(ymin = ystart, ymax = yend, xmin = -Inf, xmax = Inf, fill = col),
              alpha = 0.5) +
    geom_point(data = results,
               aes(x = Value, y = TimeNumeric), color = "black", size = 3) +
    geom_path(data = results,
              aes(x = Value, y = TimeNumeric, group = Metric),
              color = "black", linetype = "dashed") +
    facet_wrap(~Metric, scales = "free_x") +
    scale_fill_jco(guide = guide_legend(reverse = TRUE)) +
    scale_y_reverse(breaks = time_breaks, labels = time_breaks) +
    labs(
      title = "Alpha diversity metrics (British/Irish Isles)",
      subtitle = "Scaled with SRS",
      x = "Value",
      y = "Time (Years BP)",
      fill = "Time Periods"
    ) +
    coord_cartesian(ylim = c(plot_max, plot_min)) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y = element_text(size = 12, colour = "black"),
      axis.text.x = element_text(size = 12, colour = "black", angle = 45, vjust = .5),
      axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
      strip.background = element_rect(fill = "#f0f0f0", color = NA),
      strip.text.x = element_text(size = 14, face = "bold", colour = "black"),
      legend.title = element_text(size = 16, face = "bold", colour = "black"),
      legend.text = element_text(size = 12, colour = "black"),
      legend.position = "right"
    )

  # ---- Save plot ----
  ggsave(
    filename = "003-alpha-diversity-Britain-SRS.jpg",
    plot = p,
    path = output_dir,
    width = 3300, height = 4200, units = "px", dpi = 300
  )

  message("✅ Alpha diversity plot generated successfully.")
  return(list(plot = p, results = results))
}


# ---- 14. Run analysis ----
res_alpha <- generate_alphaDiversity_SRS(
  Listdf,
  rects,
  output_dir = here::here("analysis", "figures")
)

# Inspect results
# head(res_alpha$results)
# res_alpha$plot
