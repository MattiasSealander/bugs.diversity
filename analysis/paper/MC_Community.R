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

# List with data split by region
Listdf <- split(raw.mat[,-3], raw.mat$region)

is_all_zero_numeric_df <- function(df) {
  num_data <- df[sapply(df, is.numeric)]

  if (ncol(num_data) == 0) return(TRUE)
  num_vals <- unlist(num_data)
  num_vals <- num_vals[!is.na(num_vals)]
  if (length(num_vals) == 0) return(TRUE)
  all(num_vals == 0)
}

clean_df_remove_empty_communities <- function(df) {
  # Remove columns where the sum is zero (empty communities)
  numeric_cols <- sapply(df, is.numeric)
  numeric_data <- df[, numeric_cols, drop = FALSE]

  # Calculate column sums (communities)
  col_sums <- colSums(numeric_data, na.rm = TRUE)

  # Keep columns where sum > 0
  keep_cols <- col_sums > 0

  # If no columns left, return NULL (discard)
  if (!any(keep_cols)) return(NULL)

  # Keep non-numeric columns and numeric columns with sum > 0
  df_clean <- cbind(
    df[, !numeric_cols, drop = FALSE],    # non-numeric columns (like taxon)
    numeric_data[, keep_cols, drop = FALSE]  # filtered numeric cols
  )

  df_clean
}

Listdf2 <- lapply(Listdf, function(df) {
  splits <- split(df[,-2], df$end.1)

  splits_clean <- Filter(function(subdf) {
    !is_all_zero_numeric_df(subdf)
  }, splits)

  # Clean each subdf by removing empty communities (columns)
  splits_clean <- lapply(splits_clean, function(subdf) {
    cleaned <- clean_df_remove_empty_communities(subdf)
    # If cleaning removed all communities, discard this df
    if (is.null(cleaned)) return(NULL)
    cleaned
  })

  # Remove NULL elements after cleaning
  splits_clean <- Filter(Negate(is.null), splits_clean)

  # Set taxon as rownames if present
  lapply(splits_clean, function(subdf) {
    if ("taxon" %in% names(subdf)) {
      column_to_rownames(subdf, "taxon")
    } else {
      subdf
    }
  })
})

# Your improved prepare_for_MetaCommunity function (from before)
prepare_for_MetaCommunity <- function(df, verbose = TRUE) {
  if ("taxon" %in% names(df)) {
    df <- tibble::column_to_rownames(df, "taxon")
  }

  df_num <- as.data.frame(lapply(df, function(x) as.numeric(as.character(x))))
  df_num <- df_num[, colSums(is.na(df_num)) == 0, drop = FALSE]
  df_num <- df_num[rowSums(df_num) > 0, , drop = FALSE]
  df_num <- df_num[, colSums(df_num) > 0, drop = FALSE]

  if (nrow(df_num) < 2) {
    if (verbose) message("Skipping: fewer than 2 species (rows).")
    return(NULL)
  }
  if (ncol(df_num) < 2) {
    if (verbose) message("Skipping: fewer than 2 communities (columns).")
    return(NULL)
  }

  Abundances <- as.matrix(df_num)

  if (!is.numeric(Abundances)) {
    if (verbose) message("Abundances matrix is not numeric after coercion.")
    return(NULL)
  }

  if (is.null(rownames(Abundances)) || length(rownames(Abundances)) != nrow(Abundances)) {
    rownames(Abundances) <- paste0("sp_", seq_len(nrow(Abundances)))
  }
  if (is.null(colnames(Abundances)) || length(colnames(Abundances)) != ncol(Abundances)) {
    colnames(Abundances) <- paste0("comm_", seq_len(ncol(Abundances)))
  }

  Weights <- colSums(Abundances)
  if (any(Weights == 0)) {
    if (verbose) message("Skipping: some communities have zero total abundance.")
    return(NULL)
  }
  Weights <- Weights / sum(Weights)

  mc <- tryCatch({
    entropart::MetaCommunity(Abundances = Abundances, Weights = Weights)
  }, error = function(e) {
    if (verbose) message("MetaCommunity() failed: ", e$message)
    return(NULL)
  })

  if (!is.null(mc) && !is.numeric(mc$Ns)) {
    if (verbose) message("MetaCommunity object Ns component is not numeric; returning NULL.")
    return(NULL)
  }

  return(mc)
}


generate_bcShannon_plots <- function(list_of_lists) {
  plot_list <- list()

  # Get top-level list names (if available)
  top_list_names <- names(list_of_lists)
  if (is.null(top_list_names)) {
    top_list_names <- paste0("List_", seq_along(list_of_lists))
  }

  for (i in seq_along(list_of_lists)) {
    sublist <- list_of_lists[[i]]

    sublist_names <- names(sublist)
    if (is.null(sublist_names)) {
      sublist_names <- as.character(seq_along(sublist))
    }

    results <- data.frame(
      SublistIndex = integer(0),
      Shannon = numeric(0),
      stringsAsFactors = FALSE
    )

    for (j in seq_along(sublist)) {
      df <- sublist[[j]]
      mc <- prepare_for_MetaCommunity(df, verbose = FALSE)
      if (!is.null(mc)) {
        shannon <- tryCatch({
          bcSimpson(Ns = mc$Ns)
        }, error = function(e) NA_real_)
        results <- rbind(results, data.frame(SublistIndex = j, Shannon = shannon))
      }
    }

    if (nrow(results) == 0) {
      message(sprintf("No valid MetaCommunity objects in list %d; no plot created.", i))
      next
    }

    # Match results indices with names (safe if some sublists were skipped)
    valid_indices <- results$SublistIndex
    labels <- sublist_names[valid_indices]

    p <- ggplot(results, aes(x = SublistIndex, y = Shannon)) +
      geom_rect(
        data = rects,
        inherit.aes = FALSE,
        aes(ymin = xstart, ymax = xend, xmin = -Inf, xmax = Inf, fill = col),
        alpha = 0.2
      ) +
      geom_point(shape = 16, size = 3, show.legend = F) +
      geom_line(group = 1, linetype = "dashed", linewidth = 1) +
      scale_x_reverse(breaks = valid_indices,  # only odd indices (1, 3, 5, ...)
                      labels = labels) +
      coord_flip() +
      labs(title = paste(top_list_names[i]),
           x = "Age (BP)",
           y = "Shannon") +
      theme_minimal() +
      theme(legend.position = "none",
            axis.text.y = element_text(size = 12, colour = "black"),
            axis.text.x = element_text(size = 12, colour = "black", angle = 45, vjust = .5),
            axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
            strip.text.x = element_text(size = 14, face = "bold", colour = "black"),
            strip.text.y = element_text(size = 12, face = "bold", colour = "black"))

    plot_list[[paste0("List_", i)]] <- p
  }

  return(plot_list)
}

plots <-
  generate_bcShannon_plots(Listdf2)

# Show plot for List 1:
print(plots[["List_1"]])


# Period definitions from Pilotto et al. (2022)
rects <- data.frame(
  ystart = c(12000, 9500, 8000, 6500, 4500, 4000, 3500, -500),
  yend   = c(16000,12000, 9500, 8000, 6500, 4500, 4000, 3500),
  col    = factor(
    c("TM 1", "TM 2", "TM 3", "TM 4", "TM 5", "TM 6", "TM 7", "TM 8"),
    levels = c("TM 1", "TM 2", "TM 3", "TM 4", "TM 5", "TM 6", "TM 7", "TM 8")
  )
)


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
      scale_fill_jco() +
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
    file_name <- paste0("003-", safe_name, "-bias-corrected-diversity.jpg")

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

plots <- generate_bcDiversity_plots(Listdf2)

# View a specific plot
print(plots[["British/Irish Isles"]])  # Replace with actual list name

p <- DivPart(q = 1, MC = ListMetaComs[[1]][[1]], Biased = FALSE)
summary(p)
plot(p)

de <- DivEst(q = 1, ListMetaComs[[1]][[1]], Biased = FALSE, Correction = "Best", Simulations = 1000)
plot(de)
summary(de)

bcProfile <- CommunityProfile(bcDiversity, ListMetaComs[[1]][[1]]$Ns)
Profile <- CommunityProfile(Diversity, ListMetaComs[[1]][[1]]$Ps)

plot(bcProfile, type = "l", main = "", xlab = "q", ylab = "Diversity")
lines(y ~ x, data = Profile, lty = 3)
legend("topright", c("Bias Corrected", "Biased"), lty = c(1, 3),
       + inset = 0.02)

ggplot() +
  geom_line(aes(x=bcProfile$x, y=bcProfile$y, color = "Bias corrected"), linetype = "solid", linewidth = 1) +
  geom_line(aes(x=Profile$x, y=Profile$y, color = "Without correction"), linetype = "dashed", linewidth = 1) +
  scale_color_manual(values = c("black","black")) +
  labs(
    title = "γ diversity profile of the the meta-community",
    x = "q",
    y = "Diversity"
  ) +
  theme_bw() +
  theme(legend.position = c(.85,.9),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))
