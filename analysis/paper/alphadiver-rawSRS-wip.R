# ==============================================================
# Alpha diversity – British/Irish Isles (Raw + SRS)
# ==============================================================
pacman::p_load(
  cowplot, data.table, entropart, ggh4x, ggplot2,
  ggsci, tidyverse, IRanges, here, SRS
)

# ---- 1. Import species occurrence data ----
bugs <- fread(
  here::here("analysis/data/raw_data/bugs_europe_extraction_samples_20250612.csv"),
  na.strings = c("", "NA", "NULL"),
  encoding = "UTF-8"
)

# ---- 2. Define temporal periods ----
rects <- data.frame(
  ystart = c(12000, 9500, 8000, 6500, 4500, 4000, 3500, -500),
  yend   = c(16000, 12000, 9500, 8000, 6500, 4500, 4000, 3500),
  col    = factor(
    c("TM 1", "TM 2", "TM 3", "TM 4", "TM 5", "TM 6", "TM 7", "TM 8"),
    levels = c("TM 1", "TM 2", "TM 3", "TM 4", "TM 5", "TM 6", "TM 7", "TM 8")
  )
)

# ---- 3. Prepare temporal data ----
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

# ---- 4. Define 500-year bins ----
range <- data.frame(
  start = c(-500, seq(1, 15501, by = 500)),
  end = seq(0, 16000, by = 500)
)

# ---- 5. Map samples to bins ----
intersection <- findOverlaps(
  query = do.call(IRanges, time.mat),
  subject = do.call(IRanges, range),
  type = "any"
)
hits <- data.frame(time.mat[queryHits(intersection),], range[subjectHits(intersection),]) %>%
  as_tibble(rownames = "sample") %>%
  separate(sample, into = c("sample", "sample_group", "site"), sep = "@") %>%
  inner_join(bugs, by = c("sample", "sample_group")) %>%
  select(-start, -end)

# ---- 6. Raw species abundance matrix (British/Irish Isles only) ----
raw.mat <- hits %>%
  select(site.x, latitude, longitude, sample_group, sample, end.1, taxon, abundance) %>%
  mutate(sample = paste(sample, sample_group, site.x, sep = "@")) %>%
  mutate(
    region = case_when(
      between(latitude, 49.8, 62.6) & between(longitude, -12.6, 1.8) ~ "British/Irish Isles",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(region == "British/Irish Isles") %>%
  select(-latitude, -longitude, -site.x, -sample_group) %>%
  distinct()

# ---- 7. SRS scaling ----
Cmin <- 10
SRS_seed <- 1988
message("Running SRS scaling...")

srs_agg <- raw.mat %>%
  group_by(taxon, sample) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop")

srs_prep <- srs_agg %>%
  pivot_wider(names_from = sample, values_from = abundance, values_fill = 0) %>%
  column_to_rownames("taxon")

srs_prep <- srs_prep[rowSums(srs_prep) > 0, , drop = FALSE]
srs_prep <- srs_prep[, colSums(srs_prep) >= Cmin, drop = FALSE]

taxa_names <- rownames(srs_prep)
set.seed(SRS_seed)
srs_out <- SRS(srs_prep, Cmin = Cmin, set_seed = TRUE, seed = SRS_seed)
rownames(srs_out) <- taxa_names

# Convert SRS to long format
srs_long <- as.data.table(srs_out, keep.rownames = "taxon") %>%
  melt(id.vars = "taxon", variable.name = "sample", value.name = "abundance")

# Restore time info from raw.mat
sample_info <- unique(raw.mat[, c("sample", "end.1")])
srs_long <- merge(srs_long, sample_info, by = "sample", all.x = TRUE)

# ---- 8. Build time-binned community matrices ----
build_Listdf <- function(df) {
  df %>%
    group_by(end.1) %>%
    group_split() %>%
    setNames(map_chr(., ~ as.character(unique(.$end.1)))) %>%
    map(~ {
      # Identify sample columns
      sample_cols <- setdiff(unique(.x$sample), NA_character_)

      # Pivot to wide format with only sample columns
      tab <- pivot_wider(
        .x,
        id_cols = "taxon",
        names_from = sample,
        values_from = abundance,
        values_fill = 0
      )

      # Remove any extra columns that are not samples
      tab <- tab[, c("taxon", sample_cols), drop = FALSE]

      # Convert taxon to rownames
      tab <- column_to_rownames(tab, "taxon")

      # Ensure everything numeric
      tab[] <- lapply(tab, function(x) as.numeric(as.character(x)))
      tab
    })
}

Listdf_raw <- build_Listdf(raw.mat)
Listdf_srs <- build_Listdf(srs_long)

# Check if a matrix/data frame is all zeros
is_all_zero_numeric_df <- function(df) {
  num_data <- df[sapply(df, is.numeric)]
  if (ncol(num_data) == 0) return(TRUE)
  all(unlist(num_data) == 0)
}

# Remove columns (samples) where total abundance is 0
clean_df_remove_empty_communities <- function(df) {
  numeric_cols <- sapply(df, is.numeric)
  numeric_data <- df[, numeric_cols, drop = FALSE]
  keep_cols <- colSums(numeric_data, na.rm = TRUE) > 0
  if (!any(keep_cols)) return(NULL)
  cbind(df[, !numeric_cols, drop = FALSE], numeric_data[, keep_cols, drop = FALSE])
}

Listdf_clean <- lapply(Listdf_raw, function(df) {
  if (is_all_zero_numeric_df(df)) return(NULL)
  clean_df_remove_empty_communities(df)
})

Listdf_clean <- lapply(Listdf_srs, function(df) {
  if (is_all_zero_numeric_df(df)) return(NULL)
  clean_df_remove_empty_communities(df)
})

# Remove NULLs (entirely empty time bins)
Listdf_clean <- Filter(Negate(is.null), Listdf_clean)
Listdf_cleansrs <- Filter(Negate(is.null), Listdf_clean)

# ---- 9. Compute and plot alpha diversity ----
generate_alphaDiversity_BI <- function(Listdf_clean, Listdf_cleansrs, rects,
                                       output_dir = here::here("analysis/figures")) {

  results <- data.frame()

  for (time_name in names(Listdf_raw)) {
    comm_raw <- Listdf_raw[[time_name]]
    comm_srs <- Listdf_srs[[time_name]]

    # Skip empty
    if (is.null(comm_raw) || nrow(comm_raw) == 0) next

    compute_metrics <- function(mat) {
      cov <- tryCatch(entropart::Coverage(mat), error = function(e) NA)
      richness <- tryCatch(sum(rowSums(mat) > 0), error = function(e) NA)
      MC <- tryCatch(entropart::MetaCommunity(Abundances = mat), error = function(e) NULL)
      shannon <- if (!is.null(MC)) tryCatch(as.numeric(entropart::AlphaDiversity(MC, q = 1, Correction = "Best")$Total), error = function(e) NA) else NA
      simpson <- if (!is.null(MC)) tryCatch(as.numeric(entropart::AlphaDiversity(MC, q = 2, Correction = "Best")$Total), error = function(e) NA) else NA
      data.frame(Coverage = cov, `Observed Richness` = richness, Shannon = shannon, `Inverse Simpson` = simpson)
    }

    raw_metrics <- compute_metrics(as.matrix(comm_raw))
    raw_metrics$Time <- time_name
    raw_metrics$Type <- "Raw"

    srs_metrics <- compute_metrics(as.matrix(comm_srs))
    srs_metrics$Time <- time_name
    srs_metrics$Type <- "SRS"

    results <- bind_rows(results, raw_metrics %>% pivot_longer(-c(Time, Type), names_to = "Metric", values_to = "Value"),
                         srs_metrics %>% pivot_longer(-c(Time, Type), names_to = "Metric", values_to = "Value"))
  }

  results$TimeNumeric <- as.numeric(results$Time)
  results$Metric <- factor(results$Metric, levels = c("Coverage", "Observed Richness", "Shannon", "Inverse Simpson"))

  plot_min <- min(results$TimeNumeric, na.rm = TRUE) - 500
  plot_max <- max(results$TimeNumeric, na.rm = TRUE) + 500
  rects_filtered <- rects[rects$ystart <= plot_max & rects$yend >= plot_min, ]
  time_breaks <- pretty(results$TimeNumeric, n = 10)

  # ---- Plot ----
  p <- ggplot() +
    geom_rect(data = rects_filtered, aes(ymin = ystart, ymax = yend, xmin = -Inf, xmax = Inf, fill = col), alpha = 0.5) +
    geom_point(data = subset(results, Type == "SRS"), aes(x = Value, y = TimeNumeric, color = Type, shape = Type), size = 3) +
    geom_path(data = subset(results, Type == "SRS"), aes(x = Value, y = TimeNumeric, group = Metric, color = Type), linetype = "dashed") +
    geom_point(data = subset(results, Type == "Raw"), aes(x = Value, y = TimeNumeric, color = Type, shape = Type), size = 3) +
    geom_path(data = subset(results, Type == "Raw"), aes(x = Value, y = TimeNumeric, group = Metric, color = Type), linetype = "dashed") +
    facet_wrap(~Metric, scales = "free_x") +
    scale_color_manual(values = c(Raw = "black", SRS = "green4")) +
    scale_shape_manual(values = c(Raw = 16, SRS = 16)) +
    scale_fill_jco(guide = guide_legend(reverse = TRUE)) +
    scale_y_reverse(breaks = time_breaks, labels = time_breaks) +
    labs(title = "Alpha diversity – British/Irish Isles", x = "Diversity Value", y = "Time (Years BP)",
         color = "Type", shape = "Type", fill = "Time Periods") +
    coord_cartesian(ylim = c(plot_max, plot_min)) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, vjust = .5),
      strip.text = element_text(size = 14, face = "bold"),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 10)
    )

  ggsave(filename = "003-alpha-diversity-RawSRS-British_Irish_Isles.jpg", plot = p,
         path = output_dir, width = 3300, height = 4200, units = "px", dpi = 300)

  message("✅ Alpha diversity plot generated successfully for British/Irish Isles.")
  list(plot = p, results = results)
}

# ---- 10. Run analysis ----
generate_alphaDiversity_BI(Listdf_clean, Listdf_cleansrs, rects, here::here("analysis", "figures"))
