# ==============================================================
# Title: Holocene Fossil Insect Gamma Diversity Analysis
# Author: [Your Name]
# Date: [YYYY-MM-DD]
#
# Description:
#   This script processes fossil insect occurrence data, bins it into
#   500-year intervals, constructs regional community matrices, and
#   computes gamma-diversity metrics (Richness, Shannon, Inverse Simpson,
#   and Sample Coverage) across time.
#
#   Outputs:
#     - Regional gamma-diversity time-series plots
#     - Summary data for each diversity metric
#
#   Temporal periods follow Pilotto et al. (2022) Holocene definitions.
# ==============================================================


# ---- 0. Load required packages ----
pacman::p_load(
  cowplot, data.table, entropart, ggh4x, ggplot2,
  ggsci, tidyverse, IRanges, here
)


# ---- 1. Import species occurrence data ----
#' @description Load fossil insect dataset containing site, sample,
#'              age, and taxonomic abundance data.
bugs <- fread(
  here::here("analysis", "data", "raw_data", "bugs_europe_extraction_samples_20250612.csv"),
  na.strings = c("", "NA", "NULL"),
  encoding = "UTF-8"
)


# ---- 2. Define temporal periods (Pilotto et al. 2022) ----
#' @description Define Holocene time intervals for plot shading.
rects <- data.frame(
  ystart = c(12000, 9500, 8000, 6500, 4500, 4000, 3500, -500),
  yend   = c(16000, 12000, 9500, 8000, 6500, 4500, 4000, 3500),
  col    = factor(
    c("TM 1", "TM 2", "TM 3", "TM 4", "TM 5", "TM 6", "TM 7", "TM 8"),
    levels = c("TM 1", "TM 2", "TM 3", "TM 4", "TM 5", "TM 6", "TM 7", "TM 8")
  )
)


# ---- 3. Prepare data for temporal binning ----
#' @description Filter valid stratigraphic records, constrain to 16 ka–present,
#'              and remove samples exceeding a 2000-year age range.
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


# ---- 4. Define 500-year age bins ----
#' @description Create equally spaced 500-year bins for temporal overlap detection.
range <- data.frame(
  start = c(-500, seq(1, 15501, by = 500)),
  end = seq(0, 16000, by = 500)
)


# ---- 5. Identify overlapping samples with each time bin ----
#' @description Determine which samples fall within or overlap each 500-year bin.
intersection <- findOverlaps(
  query = do.call(IRanges, time.mat),
  subject = do.call(IRanges, range),
  type = "any"
)


# ---- 6. Merge sample and bin information ----
#' @description Join age-bin intersections with original species data.
hits <- data.frame(time.mat[queryHits(intersection),], range[subjectHits(intersection),]) %>%
  as_tibble(rownames = "sample") %>%
  separate(sample, into = c("sample", "sample_group", "site"), sep = "@") %>%
  inner_join(bugs, by = c("sample", "sample_group")) %>%
  select(-start, -end)


# ---- 7. Construct raw species abundance matrix ----
#' @description Build sample × species abundance matrix and assign
#'              samples to biogeographic regions.
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
  select(-place, -latitude, -longitude) %>%
  group_by(end.1, region, sample) %>%
  pivot_wider(
    id_cols = c(taxon, region, end.1),
    names_from = sample,
    values_from = abundance,
    values_fill = 0
  ) %>%
  ungroup() %>%
  distinct() %>%
  select(taxon, end.1, region, everything())


# ---- 8. Split data into nested regional lists ----
#' @description Create nested list by region → time bin.
Listdf <- split(raw.mat[, -3], raw.mat$region)


# ---- 9. Helper: Check if data frame is all zeros ----
#' @title is_all_zero_numeric_df
#' @description Check if all numeric values in a data frame are 0 or NA.
is_all_zero_numeric_df <- function(df) {
  num_data <- df[sapply(df, is.numeric)]
  if (ncol(num_data) == 0) return(TRUE)
  num_vals <- unlist(num_data)
  num_vals <- num_vals[!is.na(num_vals)]
  if (length(num_vals) == 0) return(TRUE)
  all(num_vals == 0)
}


# ---- 10. Helper: Remove empty communities ----
#' @title clean_df_remove_empty_communities
#' @description Removes columns (communities) where total abundance = 0.
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


# ---- 11. Clean, split, and prepare data for MetaCommunity analysis ----
#' @description Clean and structure data for region × time metacommunity analysis.
Listdf2 <- lapply(Listdf, function(df) {
  splits <- split(df[, -2], df$end.1)
  splits_clean <- Filter(function(subdf) !is_all_zero_numeric_df(subdf), splits)
  splits_clean <- lapply(splits_clean, function(subdf) {
    cleaned <- clean_df_remove_empty_communities(subdf)
    if (is.null(cleaned)) return(NULL)
    cleaned
  })
  splits_clean <- Filter(Negate(is.null), splits_clean)
  lapply(splits_clean, function(subdf) {
    if ("taxon" %in% names(subdf)) column_to_rownames(subdf, "taxon") else subdf
  })
})


# ---- 12. Helper: Prepare data for MetaCommunity ----
#' @title prepare_for_MetaCommunity
#' @description Clean and validate abundance matrices for entropart::MetaCommunity.
#' @param df Data frame with species (rows) and communities (columns)
#' @param verbose Logical; print diagnostic messages
#' @return MetaCommunity object or NULL
prepare_for_MetaCommunity <- function(df, verbose = TRUE) {
  if ("taxon" %in% names(df)) df <- tibble::column_to_rownames(df, "taxon")
  df_num <- as.data.frame(lapply(df, function(x) as.numeric(as.character(x))))
  df_num <- df_num[, colSums(is.na(df_num)) == 0, drop = FALSE]
  df_num <- df_num[rowSums(df_num) > 0, , drop = FALSE]
  df_num <- df_num[, colSums(df_num) > 0, drop = FALSE]

  if (nrow(df_num) < 2 || ncol(df_num) < 2) {
    if (verbose) message("Insufficient data for MetaCommunity.")
    return(NULL)
  }

  Abundances <- as.matrix(df_num)
  Weights <- colSums(Abundances)
  if (any(Weights == 0)) return(NULL)
  Weights <- Weights / sum(Weights)

  tryCatch(
    entropart::MetaCommunity(Abundances = Abundances, Weights = Weights),
    error = function(e) {
      if (verbose) message("MetaCommunity() failed: ", e$message)
      NULL
    }
  )
}


# ---- 13. Main Function: Compute and Plot Gamma Diversity ----
#' @title generate_gammaCoverage_plots
#' @description
#' Computes gamma diversity metrics (Richness, Shannon, Inverse Simpson,
#' and Sample Coverage) across time for each biogeographic region.
#' Produces time-series plots shaded by Holocene time intervals.
#'
#' @param Listdf Nested list (region → time bin) of community matrices
#' @param rects Data frame of period rectangles for plot background
#' @param output_dir Directory path for saving figures
#'
#' @return List with two elements:
#'   \item{plots}{Named list of ggplot objects by region}
#'   \item{results}{Named list of data frames of diversity metrics}
#'
#' @examples
#' generate_gammaCoverage_plots(Listdf2, rects)
#'
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
    region_results <- data.frame(Time = character(), Metric = character(), Value = numeric(),
                                 stringsAsFactors = FALSE)

    for (time_name in names(region_data)) {
      sublist <- region_data[[time_name]]
      if (is.null(sublist) || nrow(sublist) == 0) next

      comm_mat <- tryCatch(as.matrix(sublist), error = function(e) NULL)
      if (is.null(comm_mat) || ncol(comm_mat) == 0) next

      # Compute MetaCommunity and diversity metrics
      meta <- tryCatch({
        entropart::MetaCommunity(Abundances = comm_mat)
      }, error = function(e) {
        warning(sprintf("MetaCommunity failed for %s %s: %s", region_name, time_name, e$message))
        NULL
      })
      if (is.null(meta)) next

      cov_time <- tryCatch(entropart::Coverage(rowSums(comm_mat)), error = function(e) NA_real_)
      richness <- tryCatch(entropart::GammaDiversity(meta, q = 0), error = function(e) NA_real_)
      shannon <- tryCatch(entropart::GammaDiversity(meta, q = 1), error = function(e) NA_real_)
      simpson <- tryCatch(entropart::GammaDiversity(meta, q = 2), error = function(e) NA_real_)

      region_results <- bind_rows(
        region_results,
        data.frame(Time = time_name, Metric = "Coverage", Value = cov_time),
        data.frame(Time = time_name, Metric = "Richness", Value = richness),
        data.frame(Time = time_name, Metric = "Shannon", Value = shannon),
        data.frame(Time = time_name, Metric = "Inverse Simpson", Value = simpson)
      )
    }

    if (nrow(region_results) == 0) next

    # ---- Plot preparation ----
    region_results$TimeNumeric <- suppressWarnings(as.numeric(region_results$Time))
    region_results <- region_results[!is.na(region_results$TimeNumeric), ]
    region_results$Metric <- factor(region_results$Metric,
                                    levels = c("Coverage", "Richness", "Shannon", "Inverse Simpson"))

    plot_min <- min(region_results$TimeNumeric)
    plot_max <- max(region_results$TimeNumeric)
    rects_filtered <- rects[rects$ystart <= plot_max & rects$yend >= plot_min, , drop = FALSE]
    if (nrow(rects_filtered) > 0) {
      rects_filtered$ystart <- pmin(rects_filtered$ystart, plot_max)
      rects_filtered$yend <- pmax(rects_filtered$yend, plot_min)
    }

    # Time axis breaks
    time_breaks <- pretty(region_results$TimeNumeric, n = 10)

    # ---- Plot ----
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
      scale_y_reverse(breaks = time_breaks, labels = time_breaks,
                      expand = expansion(mult = c(0.01, 0.01))) +
      labs(title = paste("Gamma diversity metrics for", region_name),
           x = "Diversity Value", y = "Time (Years BP)", fill = "Time Periods") +
      coord_cartesian(ylim = c(plot_max, plot_min)) +
      theme_minimal() +
      theme(
        axis.text = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(angle = 45, vjust = .5),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
        strip.background = element_rect(fill = "#f0f0f0", color = NA),
        strip.text = element_text(size = 14, face = "bold", colour = "black"),
        legend.title = element_text(size = 16, face = "bold", colour = "black"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.position = "right"
      )

    # ---- Save Plot ----
    safe_name <- gsub("[^A-Za-z0-9_]", "_", region_name)
    file_name <- paste0("004--gamma-diversity-", safe_name, ".jpg")
    ggsave(file_name, plot = p, path = output_dir,
           width = 3300, height = 4200, units = "px", dpi = 300)

    plot_list[[region_name]] <- p
    results_list[[region_name]] <- region_results
  }

  message("✅ Gamma diversity plots generated and saved successfully.")
  return(list(plots = plot_list, results = results_list))
}


# ---- 14. Run analysis ----
generate_gammaDiversity_plots(
  Listdf2,
  rects,
  output_dir = here::here("analysis", "figures")
)
