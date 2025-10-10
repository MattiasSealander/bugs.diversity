# ==============================================================
# Description:
#   Computes gamma diversity metrics (Richness, Shannon, Inverse Simpson,
#   Coverage) for Holocene fossil insect data using Scaling with Ranked Subsampling
#   (SRS). The script performs temporal binning (500-year bins),
#   applies SRS scaling across samples, and summarizes gamma diversity
#   trajectories by region and time.
#
# Outputs:
#   - Listdf2: nested list structure (region → time bin → community matrix)
#   - Gamma diversity plots and result tables (saved in analysis/figures/)
# ==============================================================

# ---- 0. Load required packages ----
pacman::p_load(
  cowplot, data.table, entropart, ggh4x, ggplot2,
  ggsci, tidyverse, IRanges, here, SRS
)

# ==============================================================
# 1. Import species occurrence data
# --------------------------------------------------------------
# The raw dataset includes sample-level fossil insect occurrences
# across Europe. Only samples within stratigraphic contexts are used.
# ==============================================================

bugs <- fread(
  here::here("analysis", "data", "raw_data", "bugs_europe_extraction_samples_20250612.csv"),
  na.strings = c("", "NA", "NULL"),
  encoding = "UTF-8"
)

# ==============================================================
# 2. Define temporal periods (TM1–TM8)
# --------------------------------------------------------------
# These are visually used to shade plot backgrounds, representing period
# of biotic transition in coleoptera during Holocene, after Pilotto et al. (2022)
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
# --------------------------------------------------------------
# - Excludes Greenland, and filters to samples from natural contexts (stratigraphic)
# - Restricts to 16 ka – present and ≤2000 yr time ranges
# - Produces sample-level start and end ages for IRanges binning
# ==============================================================

time.mat <- bugs %>%
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

# ==============================================================
# 4. Define 500-year temporal bins
# --------------------------------------------------------------
# Bins are used for assigning samples to time intervals.
# ==============================================================

range <- data.frame(
  start = c(-500, seq(1, 15501, by = 500)),
  end = seq(0, 16000, by = 500)
)

# ==============================================================
# 5. Identify overlapping samples with time bins
# --------------------------------------------------------------
# Uses IRanges to assign each sample to all 500-year
# bins that overlap its age range.
# ==============================================================

intersection <- findOverlaps(
  query = do.call(IRanges, time.mat),
  subject = do.call(IRanges, range),
  type = "any"
)

# ==============================================================
# 6. Merge sample metadata with bin assignments
# --------------------------------------------------------------
# Produces a table linking each sample to its temporal bin.
# ==============================================================

hits <- data.frame(time.mat[queryHits(intersection),], range[subjectHits(intersection),]) %>%
  as_tibble(rownames = "sample") %>%
  separate(sample, into = c("sample", "sample_group", "site"), sep = "@") %>%
  inner_join(bugs, by = c("sample", "sample_group")) %>%
  select(-start, -end)

# ==============================================================
# 7. Construct raw species abundance matrix
# --------------------------------------------------------------
# Combines site and sample-level metadata and assigns regions
# based on latitude/longitude thresholds.
# ==============================================================

raw.mat <- hits %>%
  select(site.x, latitude, longitude, sample_group, sample, end.1, taxon, abundance) %>%
  distinct() %>%
  mutate(sample = paste(sample, sample_group, site.x, end.1, sep = "@")) %>%
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
# 7b. Standardized Rarefaction by Subsampling (SRS)
# --------------------------------------------------------------
# Purpose:
#   Standardizes sample abundances to a common total count (Cmin)
#   while preserving relative frequencies as closely as possible.
#
# Notes:
#   - Cmin = 10 ensures comparability across low-abundance samples
#   - Reproducible via a fixed random seed
# ==============================================================

Cmin <- 10
SRS_seed <- 1988
message("Running SRS standardization...")

# 1. Aggregate by taxon and unique sample ID
srs_agg <- raw.mat %>%
  group_by(taxon, sample) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop")

# 2. Pivot to wide matrix (taxa x samples)
srs_prep <- srs_agg %>%
  pivot_wider(names_from = sample, values_from = abundance, values_fill = 0) %>%
  column_to_rownames("taxon")

# 3. Remove empty taxa/samples
srs_prep <- srs_prep[rowSums(srs_prep) > 0, , drop = FALSE]
srs_prep <- srs_prep[, colSums(srs_prep) >= Cmin, drop = FALSE]
if (nrow(srs_prep) == 0 || ncol(srs_prep) == 0)
  stop("No data left after filtering for SRS.")

# 4. Preserve taxon names for reassignment
taxa_names <- rownames(srs_prep)

# 5. Run SRS normalization
set.seed(SRS_seed)
srs_out <- SRS(srs_prep, Cmin = Cmin, set_seed = TRUE, seed = SRS_seed)
rownames(srs_out) <- taxa_names

# 6. Convert to long format for downstream use
dt <- as.data.table(srs_out, keep.rownames = "taxon")
srs_long <- melt(dt, id.vars = "taxon",
                 variable.name = "sample", value.name = "abundance")


# ==============================================================
# 7c. Merge regional and temporal information
# --------------------------------------------------------------
# Restores sample metadata and reattaches region assignments
# from raw.mat, after decomposing the concatenated sample IDs.
# ==============================================================

setDT(srs_long)
srs_long[, c("sample", "sample_group", "site.x", "end.1") :=
           tstrsplit(sample, "@", fixed = TRUE)]

# Extract region info from raw.mat
region_info <- as.data.table(raw.mat)[,
                                      c("sample", "sample_group", "site.x", "end.1") :=
                                        tstrsplit(sample, "@", fixed = TRUE)
][, end.1 := as.numeric(end.1)][
  , .(sample, sample_group, site.x, end.1, region)
]

# Remove duplicates per sample to prevent Cartesian expansion
region_info <- unique(region_info, by = c("sample", "sample_group", "site.x", "end.1"))

# Ensure data types match
srs_long[, end.1 := as.numeric(end.1)]

# Merge region assignments
srs_long <- merge(
  srs_long,
  region_info,
  by = c("sample", "sample_group", "site.x", "end.1"),
  all.x = TRUE,
  sort = FALSE
)


# ==============================================================
# 8. Prepare SRS data for diversity list construction
# --------------------------------------------------------------
# - Drops spatial coordinates (no longer needed)
# - Ensures unique taxon × sample × time combinations
# - Creates stable sample identifiers
# ==============================================================

srs_long[, sample_id := paste(sample, sample_group, site.x, sep = "@")]
srs_long <- unique(srs_long, by = c("region", "end.1", "taxon", "sample_id"))
srs_long[, end.1 := as.numeric(end.1)]

# ==============================================================
# 9. Build Listdf2 structure
# --------------------------------------------------------------
# Output format:
#   Listdf2[[region]][[time_bin]] → taxon x sample abundance matrix
# ==============================================================

Listdf2 <- srs_long[, {
  et_order <- sort(unique(end.1))
  slices <- vector("list", length = length(et_order))
  names(slices) <- as.character(et_order)

  for (et in et_order) {
    dt_slice <- .SD[end.1 == et]
    if (nrow(dt_slice) == 0) next
    tab <- dcast(dt_slice, taxon ~ sample_id, value.var = "abundance", fill = 0)
    setDF(tab)
    rownames(tab) <- tab$taxon
    tab$taxon <- NULL
    slices[[as.character(et)]] <- tab
  }
  .(region_list = list(slices))
}, by = region]

Listdf2 <- setNames(Listdf2$region_list, Listdf2$region)
message("✅ Listdf2 successfully constructed. Regions: ",
        paste(names(Listdf2), collapse = ", "))

# ==============================================================
# 10. Gamma Diversity Calculation & Plotting Function
# --------------------------------------------------------------
# Calculates coverage and Hill numbers (q = 0, 1, 2)
# for each region and time interval, then produces
# publication-ready diversity plots.
# ==============================================================

generate_gammaDiversity_plots <- function(Listdf, rects,
                                          output_dir = here::here("analysis", "figures")) {
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  requireNamespace("ggsci")
  requireNamespace("entropart")

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  plot_list <- list()
  results_list <- list()

  for (region_name in names(Listdf)) {
    region_data <- Listdf[[region_name]]
    region_results <- data.frame(Time = character(), Metric = character(), Value = numeric())

    for (time_name in names(region_data)) {
      sublist <- region_data[[time_name]]
      if (is.null(sublist) || nrow(sublist) == 0) next

      comm_mat <- tryCatch(as.matrix(sublist), error = function(e) NULL)
      if (is.null(comm_mat) || ncol(comm_mat) == 0) next

      meta <- tryCatch(entropart::MetaCommunity(Abundances = comm_mat),
                       error = function(e) NULL)
      if (is.null(meta)) next

      cov_time <- tryCatch(entropart::Coverage(rowSums(comm_mat)), error = function(e) NA_real_)
      richness <- tryCatch(entropart::GammaDiversity(meta, q = 0), error = function(e) NA_real_)
      shannon  <- tryCatch(entropart::GammaDiversity(meta, q = 1), error = function(e) NA_real_)
      simpson  <- tryCatch(entropart::GammaDiversity(meta, q = 2), error = function(e) NA_real_)

      region_results <- bind_rows(
        region_results,
        data.frame(Time = time_name, Metric = "Coverage", Value = cov_time),
        data.frame(Time = time_name, Metric = "Richness", Value = richness),
        data.frame(Time = time_name, Metric = "Shannon", Value = shannon),
        data.frame(Time = time_name, Metric = "Inverse Simpson", Value = simpson)
      )
    }

    if (nrow(region_results) == 0) next

    # Prepare for plotting
    region_results$TimeNumeric <- suppressWarnings(as.numeric(region_results$Time))
    region_results <- region_results[!is.na(region_results$TimeNumeric), ]
    region_results$Metric <- factor(region_results$Metric,
                                    levels = c("Coverage", "Richness", "Shannon", "Inverse Simpson"))

    # Restrict rectangles to visible range
    plot_min <- min(region_results$TimeNumeric)
    plot_max <- max(region_results$TimeNumeric)
    rects_filtered <- rects[rects$ystart <= plot_max & rects$yend >= plot_min, , drop = FALSE]
    if (nrow(rects_filtered) > 0) {
      rects_filtered$ystart <- pmin(rects_filtered$ystart, plot_max)
      rects_filtered$yend <- pmax(rects_filtered$yend, plot_min)
    }

    time_breaks <- pretty(region_results$TimeNumeric, n = 10)

    # Plot diversity metrics
    p <- ggplot() +
      geom_rect(data = rects_filtered,
                aes(ymin = ystart, ymax = yend, xmin = -Inf, xmax = Inf, fill = col),
                alpha = 0.5) +
      geom_point(data = region_results,
                 aes(x = Value, y = TimeNumeric),
                 color = "black", size = 3) +
      geom_path(data = region_results,
                aes(x = Value, y = TimeNumeric, group = Metric),
                color = "black", linetype = "dashed") +
      facet_wrap(~Metric, scales = "free_x") +
      scale_fill_jco(guide = guide_legend(reverse = TRUE)) +
      scale_y_reverse(breaks = time_breaks, labels = time_breaks) +
      labs(title = paste("Gamma diversity metrics for", region_name),
           x = "Diversity Value", y = "Time (Years BP)", fill = "Time Periods") +
      coord_cartesian(ylim = c(plot_max, plot_min)) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, vjust = .5),
        strip.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 10)
      )

    safe_name <- gsub("[^A-Za-z0-9_]", "_", region_name)
    ggsave(filename = paste0("004-gamma-diversity-SRS2-", safe_name, ".jpg"),
           plot = p, path = output_dir, width = 3300, height = 4200, units = "px", dpi = 300)

    plot_list[[region_name]] <- p
    results_list[[region_name]] <- region_results
  }

  message("✅ Gamma diversity plots generated and saved successfully.")
  return(list(plots = plot_list, results = results_list))
}

# ==============================================================
# 11. Run gamma diversity analysis
# --------------------------------------------------------------
# Executes gamma diversity calculations and plotting for all regions.
# ==============================================================

generate_gammaDiversity_plots(
  Listdf2,
  rects,
  output_dir = here::here("analysis", "figures")
)
