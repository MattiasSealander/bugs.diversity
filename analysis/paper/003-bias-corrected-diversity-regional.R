pacman::p_load(cowplot, entropart, ggh4x, ggplot2, ggsci, tidyverse, IRanges)

# Import site and species data
bugs.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "bugs_europe_extraction_samples_20250612.csv"), sep = ";", dec = ".", header = TRUE, na = c("", "NA", "NULL"), encoding = "UTF-8")

# Prepare data for age bin intersection analysis
time.mat <-
  bugs.csv %>%
  select (country, sample, site, sample_group, age_older, age_younger, context) %>%
  mutate(age_range=age_older - age_younger) %>%
  filter(context == "Stratigraphic sequence", sample != "BugsPresence", between(age_older, -500, 16000) & between(age_younger, -500, 16000), age_range <= 2000,
         country != "Greenland") %>%
  distinct() %>%
  select(-age_range, -context) %>%
  mutate(sample = paste(sample, sample_group, site, sep = "@")) %>%
  column_to_rownames("sample") %>%
  select(-country, -site, -sample_group) %>%
  dplyr::rename(start = age_younger, end = age_older)

# Select age bin(s)
range <- data.frame(
  start = c(-500, seq(1, 15501, by = 500)),
  end = seq(0, 16000, by = 500))

# Find what samples overlap in time with the specified age bin(s)
intersection <-
  findOverlaps(query = do.call(IRanges, time.mat), subject = do.call(IRanges, range), type = "any")

# Retrieve the hits from the intersection and name of samples
hits <-
  data.frame(time.mat[queryHits(intersection),], range[subjectHits(intersection),]) %>%
  # Retrieve concatenated row name
  as_tibble(rownames = "sample") %>%
  # Separate concatenated sample name into individual columns
  separate(col = sample,
           sep = "@",
           into = c("sample", "sample_group", "site")) %>%
  # Join the age bins to fossil insect data on sample + sample group name
  inner_join(bugs.csv, by = c("sample", "sample_group")) %>%
  select(-start, -end)

# Species matrix of raw values
# Filter data and pivot insect species with abundance values, grouped by (site + date + sample group + sample)
raw.mat <-
  hits %>%
  select (site.x, latitude, longitude, sample_group, sample, end.1, taxon, abundance) %>%
  distinct() %>%
  mutate(sample = paste(sample, sample_group, site.x, sep = "@")) %>%
  mutate(place = case_when(
    between(latitude, 53.5, 74.8) & between(longitude, 2.8, 41) ~ "Scandinavia/Baltic",
    between(latitude, 49.8, 62.6) & between(longitude, -12.6, 1.8) ~ "British/Irish Isles",
    between(latitude, 63, 66.3) & between(longitude, -25, -12) ~ "Iceland",
    latitude < 45 ~ "Meridional")) %>%
  mutate(region = ifelse(is.na(place),
                         "Continental",
                         place)) %>%
  select(-place, -latitude, -longitude) %>%
  group_by(end.1, region, sample) %>%
  pivot_wider(id_cols = c(taxon, region, end.1),
              names_from = sample,
              values_from = abundance,
              values_fill = 0) %>%
  ungroup() %>%
  distinct() %>%
  select(taxon, end.1, region, everything())

# List with sublists by region
Listdf <- split(raw.mat[,-3], raw.mat$region)

# ---- CHECK IF DATA FRAME IS ALL ZEROS ----
# This function checks whether all numeric values in a dataframe are zero or NA.
# Returns TRUE if all numeric entries are zero or missing.
is_all_zero_numeric_df <- function(df) {
  num_data <- df[sapply(df, is.numeric)]  # Keep only numeric columns

  if (ncol(num_data) == 0) return(TRUE)   # No numeric columns → treat as all zero
  num_vals <- unlist(num_data)            # Flatten to vector
  num_vals <- num_vals[!is.na(num_vals)]  # Remove NA values

  if (length(num_vals) == 0) return(TRUE) # All values were NA
  all(num_vals == 0)                      # Are all remaining values zero?
}

# ---- REMOVE EMPTY COMMUNITIES FROM DATA FRAME ----
# This function removes columns (communities) that sum to zero across all taxa.
# It preserves non-numeric columns (like "taxon") and only filters numeric ones.
# If all communities are empty, it returns NULL to signal removal.
clean_df_remove_empty_communities <- function(df) {
  numeric_cols <- sapply(df, is.numeric)                       # Identify numeric columns
  numeric_data <- df[, numeric_cols, drop = FALSE]             # Subset numeric data

  col_sums <- colSums(numeric_data, na.rm = TRUE)              # Sum of each community
  keep_cols <- col_sums > 0                                    # Logical vector for filtering

  if (!any(keep_cols)) return(NULL)                            # If no valid communities remain

  df_clean <- cbind(
    df[, !numeric_cols, drop = FALSE],                         # Keep non-numeric (e.g., taxon)
    numeric_data[, keep_cols, drop = FALSE]                    # Keep only non-zero communities
  )

  return(df_clean)
}

# ---- SPLIT, FILTER, AND CLEAN DATA BY TIME INTERVAL ----
# This pipeline:
# 1. Splits each top-level data frame in Listdf by the `end.1` time column.
# 2. Filters out sub-dataframes that contain only zeros.
# 3. Removes communities (columns) that have zero abundance across all taxa.
# 4. Removes sublists where all data was filtered out.
# 5. Sets "taxon" column as rownames, if present.

Listdf2 <- lapply(Listdf, function(df) {

  # ---- STEP 1: SPLIT BY TIME INTERVAL ----
  splits <- split(df[,-2], df$end.1)  # Remove column 2 (end.1) to avoid duplication

  # ---- STEP 2: FILTER OUT ZERO-ONLY DATAFRAMES ----
  splits_clean <- Filter(function(subdf) {
    !is_all_zero_numeric_df(subdf)
  }, splits)

  # ---- STEP 3: REMOVE EMPTY COMMUNITIES (ZERO-COLUMNS) ----
  splits_clean <- lapply(splits_clean, function(subdf) {
    cleaned <- clean_df_remove_empty_communities(subdf)

    # ---- STEP 4: DISCARD IF CLEANING REMOVED EVERYTHING ----
    if (is.null(cleaned)) return(NULL)
    cleaned
  })

  # ---- STEP 5: REMOVE NULL ELEMENTS (EMPTY AFTER CLEANING) ----
  splits_clean <- Filter(Negate(is.null), splits_clean)

  # ---- STEP 6: SET TAXON AS ROW NAMES (IF AVAILABLE) ----
  lapply(splits_clean, function(subdf) {
    if ("taxon" %in% names(subdf)) {
      column_to_rownames(subdf, "taxon")
    } else {
      subdf
    }
  })
})


#### ---- PREPARE DATA FOR MetaCommunity() ---- ####
# This function converts a data frame of taxon abundances into a MetaCommunity object
# compatible with the `entropart` package. It handles:
# - Taxon name assignment
# - Numeric coercion and cleaning
# - Filtering of empty rows/columns
# - Generation of weights
# - Error handling for minimal input or invalid structure

prepare_for_MetaCommunity <- function(df, verbose = TRUE) {

  # ---- HANDLE TAXON COLUMN ----
  # If there is a "taxon" column, convert it to rownames
  if ("taxon" %in% names(df)) {
    df <- tibble::column_to_rownames(df, "taxon")
  }

  # ---- CONVERT ALL VALUES TO NUMERIC ----
  # Apply numeric coercion to all columns
  df_num <- as.data.frame(lapply(df, function(x) as.numeric(as.character(x))))

  # ---- REMOVE COLUMNS WITH ANY NA ----
  df_num <- df_num[, colSums(is.na(df_num)) == 0, drop = FALSE]

  # ---- REMOVE EMPTY ROWS AND COLUMNS ----
  df_num <- df_num[rowSums(df_num) > 0, , drop = FALSE]
  df_num <- df_num[, colSums(df_num) > 0, drop = FALSE]

  # ---- CHECK MINIMUM DATA REQUIREMENTS ----
  # Must have at least 2 species and 2 communities
  if (nrow(df_num) < 2) {
    if (verbose) message("Skipping: fewer than 2 species (rows).")
    return(NULL)
  }
  if (ncol(df_num) < 2) {
    if (verbose) message("Skipping: fewer than 2 communities (columns).")
    return(NULL)
  }

  # ---- BUILD ABUNDANCE MATRIX ----
  Abundances <- as.matrix(df_num)

  if (!is.numeric(Abundances)) {
    if (verbose) message("Abundances matrix is not numeric after coercion.")
    return(NULL)
  }

  # ---- ENSURE PROPER ROW AND COLUMN NAMES ----
  if (is.null(rownames(Abundances)) || length(rownames(Abundances)) != nrow(Abundances)) {
    rownames(Abundances) <- paste0("sp_", seq_len(nrow(Abundances)))
  }
  if (is.null(colnames(Abundances)) || length(colnames(Abundances)) != ncol(Abundances)) {
    colnames(Abundances) <- paste0("comm_", seq_len(ncol(Abundances)))
  }

  # ---- COMPUTE COMMUNITY WEIGHTS ----
  # Weights are proportional to total abundance in each community
  Weights <- colSums(Abundances)

  if (any(Weights == 0)) {
    if (verbose) message("Skipping: some communities have zero total abundance.")
    return(NULL)
  }

  Weights <- Weights / sum(Weights)

  # ---- CONSTRUCT METACOMMUNITY OBJECT ----
  mc <- tryCatch({
    entropart::MetaCommunity(Abundances = Abundances, Weights = Weights)
  }, error = function(e) {
    if (verbose) message("MetaCommunity() failed: ", e$message)
    return(NULL)
  })

  # ---- FINAL CHECK ON RESULT ----
  # Ensure MetaCommunity object is valid
  if (!is.null(mc) && !is.numeric(mc$Ns)) {
    if (verbose) message("MetaCommunity object Ns component is not numeric; returning NULL.")
    return(NULL)
  }

  return(mc)
}

#### ---- GENERATE DIVERSITY PLOTS ---- ####
# Period definitions from Pilotto et al. (2022)
rects <- data.frame(
  ystart = c(12000, 9500, 8000, 6500, 4500, 4000, 3500, -500),
  yend   = c(16000,12000, 9500, 8000, 6500, 4500, 4000, 3500),
  col    = factor(
    c("TM 1", "TM 2", "TM 3", "TM 4", "TM 5", "TM 6", "TM 7", "TM 8"),
    levels = c("TM 1", "TM 2", "TM 3", "TM 4", "TM 5", "TM 6", "TM 7", "TM 8")
  )
)

# ---- GENERATE BIAS-CORRECTED DIVERSITY PLOTS ----
# This function takes a named list of sublists, where each sublist contains
# abundance data frames for different time slices. It computes bias-corrected
# Shannon and Simpson diversity using entropart, then generates and saves
# diversity plots with a temporal (y-axis) layout.
#
# Inputs:
#   - list_of_lists: A named list of sublists containing abundance matrices.
#
# Outputs:
#   - A named list of ggplot2 plots, one per top-level list element.
#   - Saves a JPG image for each plot to the analysis/figures directory.

generate_bcDiversity_plots <- function(list_of_lists) {
  plot_list <- list()
  # ---- GET TOP-LEVEL LIST NAMES ----
  top_list_names <- names(list_of_lists)
  if (is.null(top_list_names)) {
    top_list_names <- paste0("List_", seq_along(list_of_lists))
  }

  for (i in seq_along(list_of_lists)) {
    sublist <- list_of_lists[[i]]
    # ---- GET SUBLIST NAMES ----
    sublist_names <- names(sublist)
    if (is.null(sublist_names)) {
      sublist_names <- as.character(seq_along(sublist))
    }
    # ---- INITIALIZE RESULTS DATA FRAME ----
    results <- data.frame(
      SublistIndex = integer(0),
      Metric = character(0),
      Value = numeric(0),
      stringsAsFactors = FALSE
    )
    # ---- CALCULATE DIVERSITY VALUES ----
    for (j in seq_along(sublist)) {
      df <- sublist[[j]]
      mc <- prepare_for_MetaCommunity(df, verbose = FALSE)

      if (!is.null(mc)) {
        shannon <- tryCatch({
          bcShannon(Ns = mc$Ns)
        }, error = function(e) {
          message(sprintf("bcShannon failed for sample %d: %s", j, e$message))
          NA_real_
        })

        simpson <- tryCatch({
          bcSimpson(Ns = mc$Ns)
        }, error = function(e) {
          message(sprintf("bcSimpson failed for sample %d: %s", j, e$message))
          NA_real_
        })

        # Debug output for metric values
        #message(sprintf("Sample %d | Shannon: %s | Simpson: %s", j, shannon, simpson))

        results <- rbind(
          results,
          data.frame(SublistIndex = j, Metric = "Shannon", Value = shannon),
          data.frame(SublistIndex = j, Metric = "Simpson", Value = simpson)
        )
      } else {
        message(sprintf("MetaCommunity object was NULL for sample %d", j))
      }
    }

    if (nrow(results) == 0) {
      message(sprintf("No valid MetaCommunity objects in list %d; no plot created.", i))
      next
    }

    # ---- EXTRACT TIME FROM SUBLIST NAMES ----
    valid_indices <- sort(unique(results$SublistIndex))
    labels <- sublist_names[valid_indices]

    # Try parsing time from names
    numeric_labels <- suppressWarnings(as.numeric(gsub("[^0-9\\.\\-]", "", labels)))

    if (all(!is.na(numeric_labels))) {
      results$Time <- numeric_labels[results$SublistIndex]
    } else {
      results$Time <- valid_indices
    }

    # ---- ORDER RESULTS BY METRIC AND TIME ----
    results <- results[with(results, order(Metric, -Time)), ]

    # Debug output: show top of results
    #cat("\n--- DEBUG: Data for plot", top_list_names[i], "---\n")
    #print(head(results, 10))

    # ---- CALCULATE PLOTTING RANGE (+500 YEAR BUFFER) ----
    min_time <- min(results$Time, na.rm = TRUE)
    max_time <- max(results$Time, na.rm = TRUE)
    plot_min <- min_time - 500
    plot_max <- max_time + 500

    # ---- CROP BACKGROUND RECTANGLES TO TIME RANGE ----
    # Filter rectangles that overlap the plot range (fully or partially)
    rects_filtered <- rects[
      rects$ystart <= plot_max & rects$yend >= plot_min,
    ]

    # Clamp rectangles to visible plot range
    rects_filtered$ystart <- pmin(rects_filtered$ystart, plot_max)
    rects_filtered$yend   <- pmax(rects_filtered$yend, plot_min)

    # ---- DEFINE Y-AXIS TICKS AT 500-YEAR INTERVALS ----
    min_time <- floor(min(results$Time, na.rm = TRUE) / 500) * 500
    max_time <- ceiling(max(results$Time, na.rm = TRUE) / 500) * 500
    time_breaks <- seq(from = max_time, to = min_time, by = -500)  # because axis is reversed

    # ---- GENERATE PLOTS ----
    p <- ggplot() +
      geom_rect(
        data = rects_filtered,
        inherit.aes = FALSE,
        aes(ymin = ystart, ymax = yend, xmin = -Inf, xmax = Inf, fill = col),
        alpha = 0.5
      ) +
      geom_point(
        data = results,
        aes(x = Value, y = Time),
        color = "black", size = 3
      ) +
      geom_path(
        data = results,
        aes(x = Value, y = Time),
        color = "black",
        linetype = "dashed"
      ) +
      facet_wrap(~Metric, scales = "free_x", labeller = labeller(c("Shannon", "Simpson"))) +
      scale_fill_jco(guide = guide_legend(reverse = TRUE)) +
      scale_y_reverse(
        breaks = time_breaks,
        labels = time_breaks,
        expand = expansion(mult = c(0.01, 0.01))
      ) +
      labs(
        title = paste("Bias-corrected diversity indices for", top_list_names[i]),
        x = "Diversity Value",
        y = "Time (Years BP)",
        fill = "Period"
      ) +
      coord_cartesian(ylim = c(plot_max, plot_min)) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 12, colour = "black", angle = 0),
        axis.text.x = element_text(size = 12, colour = "black", angle = 45, vjust = .5),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
        strip.background = element_rect(fill = "#f0f0f0", color = NA),
        strip.text.x = element_text(size = 14, face = "bold", colour = "black"),
        strip.text.y = element_text(size = 12, face = "bold", colour = "black"),
        legend.title = element_text(size = 16, face = "bold", colour = "black"),
        legend.position = "right"
      )

    # ---- SAVE PLOT TO LIST ----
    plot_list[[top_list_names[i]]] <- p

    # ---- SAVE PLOT TO FILE ----
    safe_name <- gsub("[^A-Za-z0-9_]", "_", top_list_names[i])
    file_name <- paste0("003-", safe_name, "-bias-corrected-diversity-regional.jpg")

    ggsave(
      filename = file_name,
      plot = p,
      device = "jpg",
      path = here::here("analysis", "figures"),
      width = 3300,
      height = 4200,
      units = "px",
      dpi = 300
    )

  }

  return(plot_list)
}

# Save images
generate_bcDiversity_plots(Listdf2)

#for testing
#plots <- generate_bcDiversity_plots(Listdf2)
#print(plots[["Iceland"]])  # Replace with actual list name
