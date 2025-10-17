# ==============================================================
# Description:
#   This script computes alpha diversity metrics (Observed Richness,
#   Shannon, Inverse Simpson, and Coverage) for Holocene fossil insect data
#   from the British/Irish Isles region.
#
# Workflow overview:
#   1. Load and clean fossil insect occurrence data.
#   2. Define 500-year temporal bins (16 ka – present).
#   3. Assign each sample to one or more bins based on its dated range.
#   4. Construct community (taxon × sample) matrices per time bin.
#   5. Compute alpha diversity metrics using entropart.
#   6. Produce time-series plots for each metric.
#
# Outputs:
#   - `Listdf2`: list of taxon × sample matrices per 500-year time bin.
#   - Figure: "003-alpha-diversity-Britain.jpg"
#     (saved in ./analysis/figures)
# ==============================================================


# ---- 0. Load required packages ----
pacman::p_load(
  cowplot, data.table, entropart, ggh4x, ggplot2,
  ggsci, tidyverse, IRanges, here
)


# ---- 1. Import species occurrence data ----
# Load fossil insect dataset (CSV format).
# The file includes metadata such as site, sample age, and taxon abundance.
bugs <- fread(
  here::here("analysis/data/raw_data/bugs_europe_extraction_samples_20250612.csv"),
  na.strings = c("", "NA", "NULL"),
  encoding = "UTF-8"
)


# ---- 2. Define temporal periods (after Pilotto et al. 2022) ----
# Used to color-code time intervals on plots.
rects <- data.frame(
  ystart = c(12000, 9500, 8000, 6500, 4500, 4000, 3500, -500),
  yend   = c(16000, 12000, 9500, 8000, 6500, 4500, 4000, 3500),
  col    = factor(
    c("TM 1", "TM 2", "TM 3", "TM 4", "TM 5", "TM 6", "TM 7", "TM 8"),
    levels = c("TM 1", "TM 2", "TM 3", "TM 4", "TM 5", "TM 6", "TM 7", "TM 8")
  )
)


# ---- 3. Prepare sample metadata for temporal binning ----
# Filter samples to retain only samples from natural contexts (Stratigraphic sequence) within the Late Glacial - Holocene window.
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
  rename(start = age_younger, end = age_older)


# ---- 4. Define 500-year temporal bins ----
range <- data.frame(
  start = c(-500, seq(1, 15501, by = 500)),
  end   = seq(0, 16000, by = 500)
)


# ---- 5. Identify sample–bin overlaps ----
# Samples overlapping each time bin are assigned accordingly.
intersection <- findOverlaps(
  query = do.call(IRanges, time.mat),
  subject = do.call(IRanges, range),
  type = "any"
)


# ---- 6. Merge sample–bin relationships ----
hits <- data.frame(
  time.mat[queryHits(intersection),],
  range[subjectHits(intersection),]
) %>%
  as_tibble(rownames = "sample") %>%
  separate(sample, into = c("sample", "sample_group", "site"), sep = "@") %>%
  inner_join(bugs, by = c("sample", "sample_group")) %>%
  select(-start, -end)


# ---- 7. Construct raw species abundance data and assign regions ----
# Define regions based on lat/lon; focus only on the British/Irish Isles.
raw.mat <- hits %>%
  select(site.x, latitude, longitude, sample_group, sample, end.1, taxon, abundance) %>%
  distinct() %>%
  mutate(sample = paste(sample, sample_group, site.x, sep = "@")) %>%
  mutate(
    place = case_when(
      between(latitude, 53.5, 74.8) & between(longitude, 2.8, 41) ~ "Scandinavia/Baltic",
      between(latitude, 49.8, 62.6) & between(longitude, -12.6, 1.8) ~ "British/Irish Isles",
      between(latitude, 63, 66.3) & between(longitude, -25, -12) ~ "Iceland",
      latitude < 45 ~ "Meridional"
    ),
    region = ifelse(is.na(place), "Continental", place)
  ) %>%
  filter(region == "British/Irish Isles") %>%  # Focus analysis on this region
  select(-place, -latitude, -longitude)


# ---- 8. Build community matrices (one per time bin) ----
# Each matrix: rows = taxa, columns = samples.
Listdf2 <- raw.mat %>%
  group_by(end.1) %>%
  group_split() %>%
  setNames(map_chr(., ~ unique(.$end.1))) %>%
  map(~ {
    df <- .x %>%
      select(taxon, sample, abundance) %>%
      distinct() %>%  # Keep unique taxon–sample pairs, avoid summing duplicates
      pivot_wider(names_from = sample, values_from = abundance, values_fill = 0)
    if ("taxon" %in% names(df)) column_to_rownames(df, "taxon") else df
  })


# ---- 9. Helper functions ----

# Check if all numeric values are zero or NA.
is_all_zero_numeric_df <- function(df) {
  num_data <- df[sapply(df, is.numeric)]
  if (ncol(num_data) == 0) return(TRUE)
  num_vals <- unlist(num_data)
  num_vals <- num_vals[!is.na(num_vals)]
  if (length(num_vals) == 0) return(TRUE)
  all(num_vals == 0)
}

# Remove samples (columns) that contain no taxa.
clean_df_remove_empty_communities <- function(df) {
  numeric_cols <- sapply(df, is.numeric)
  numeric_data <- df[, numeric_cols, drop = FALSE]
  col_sums <- colSums(numeric_data, na.rm = TRUE)
  keep_cols <- col_sums > 0
  if (!any(keep_cols)) return(NULL)

  df_clean <- cbind(
    df[, !numeric_cols, drop = FALSE],
    numeric_data[, keep_cols, drop = FALSE]
  )
  return(df_clean)
}


# ---- 10. Clean and filter community matrices ----
Listdf2 <- Listdf2 %>%
  discard(is_all_zero_numeric_df) %>%
  map(clean_df_remove_empty_communities) %>%
  discard(is.null)


# ---- 11. Compute and visualize alpha diversity metrics ----
# Computes Observed Richness, Coverage, Shannon, and Inverse Simpson per bin.
# Returns a list with a results table and ggplot object.

generate_alphaDiversity_plots_BI <- function(Listdf2, rects,
                                             output_dir = here::here("analysis", "figures")) {

  region_results <- data.frame()

  for (time_name in names(Listdf2)) {
    sublist <- Listdf2[[time_name]]
    if (is.null(sublist) || nrow(sublist) == 0) next

    comm_data <- as.matrix(sublist)
    if (ncol(comm_data) == 0) next

    # Compute alpha diversity metrics
    obs_rich <- tryCatch(sum(rowSums(comm_data) > 0), error = function(e) NA)
    cov       <- tryCatch(entropart::Coverage(comm_data), error = function(e) NA)
    MC        <- tryCatch(entropart::MetaCommunity(Abundances = comm_data), error = function(e) NULL)

    shannon_div <- if (!is.null(MC)) tryCatch(
      as.numeric(entropart::AlphaDiversity(MC, q = 1, Correction = "Best")$Total),
      error = function(e) NA
    ) else NA

    simpson_div <- if (!is.null(MC)) tryCatch(
      as.numeric(entropart::AlphaDiversity(MC, q = 2, Correction = "Best")$Total),
      error = function(e) NA
    ) else NA

    # Append to results
    add_metric <- function(name, value) {
      if (!is.na(value))
        region_results <<- bind_rows(
          region_results,
          data.frame(Time = time_name, Metric = name, Value = value)
        )
    }

    add_metric("Observed Richness", obs_rich)
    add_metric("Coverage", cov)
    add_metric("Shannon", shannon_div)
    add_metric("Inverse Simpson", simpson_div)
  }

  # ---- Plotting ----
  region_results$TimeNumeric <- suppressWarnings(as.numeric(region_results$Time))
  region_results$Metric <- factor(
    region_results$Metric,
    levels = c("Coverage", "Observed Richness", "Shannon", "Inverse Simpson")
  )

  plot_min <- min(region_results$TimeNumeric, na.rm = TRUE) - 500
  plot_max <- max(region_results$TimeNumeric, na.rm = TRUE) + 500
  rects_filtered <- rects[rects$ystart <= plot_max & rects$yend >= plot_min, ]
  time_breaks <- pretty(region_results$TimeNumeric, n = 10)

  p <- ggplot() +
    geom_rect(
      data = rects_filtered,
      aes(ymin = ystart, ymax = yend, xmin = -Inf, xmax = Inf, fill = col),
      alpha = 0.5
    ) +
    geom_point(data = region_results, aes(x = Value, y = TimeNumeric), color = "black", size = 3) +
    geom_path(data = region_results, aes(x = Value, y = TimeNumeric, group = Metric),
              color = "black", linetype = "dashed") +
    facet_wrap(~Metric, scales = "free_x") +
    scale_fill_jco(guide = guide_legend(reverse = TRUE)) +
    scale_y_reverse(breaks = time_breaks, labels = time_breaks) +
    labs(
      title = "Alpha Diversity Metrics – British/Irish Isles",
      subtitle = "Raw abundance data (500-year bins)",
      x = "Diversity Value", y = "Time (Years BP)", fill = "Time Periods"
    ) +
    coord_cartesian(ylim = c(plot_max, plot_min)) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, vjust = .5),
      strip.text = element_text(size = 14, face = "bold"),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 10)
    )

  ggsave(
    filename = "003-alpha-diversity-British_Irish_Isles.jpg",
    plot = p, path = output_dir,
    width = 3300, height = 4200, units = "px", dpi = 300
  )

  message("✅ Alpha diversity plot generated successfully for British/Irish Isles.")
  list(plot = p, results = region_results)
}


# ---- 12. Run analysis ----
generate_alphaDiversity_plots_BI(
  Listdf2,
  rects,
  output_dir = here::here("analysis", "figures")
)
