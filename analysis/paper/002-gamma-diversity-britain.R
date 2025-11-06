# ==============================================================
# Description:
#   Computes gamma diversity metrics (Richness, Shannon, Inverse Simpson,
#   Coverage) for Holocene fossil insect data, comparing raw and
#   SRS-scaled datasets. Both are plotted together:
#     - Raw data = black
#     - SRS scaled = green4
#
#   Fix included: handles single-sample bins gracefully.
# ==============================================================

# ---- 0. Load required packages ----
pacman::p_load(
  cowplot, data.table, entropart, ggh4x, ggplot2,
  ggsci, tidyverse, IRanges, here, SRS
)

# ==============================================================
# 1. Import species occurrence data
# ==============================================================

bugs <- fread(
  here::here("analysis", "data", "raw_data", "bugs_europe_extraction_samples_20250612.csv"),
  na.strings = c("", "NA", "NULL"),
  encoding = "UTF-8"
)

# ==============================================================
# 2. Define temporal periods (TM1–TM8)
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
# ==============================================================

range <- data.frame(
  start = c(-500, seq(1, 15501, by = 500)),
  end = seq(0, 16000, by = 500)
)

# ==============================================================
# 5. Identify overlapping samples with time bins
# ==============================================================

intersection <- findOverlaps(
  query = do.call(IRanges, time.mat),
  subject = do.call(IRanges, range),
  type = "any"
)

# ==============================================================
# 6. Merge sample metadata with bin assignments
# ==============================================================

hits <- data.frame(time.mat[queryHits(intersection),], range[subjectHits(intersection),]) %>%
  as_tibble(rownames = "sample") %>%
  separate(sample, into = c("sample", "sample_group", "site"), sep = "@") %>%
  inner_join(bugs, by = c("sample", "sample_group")) %>%
  select(-start, -end)

# ==============================================================
# 7. Construct raw species abundance matrix
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
# 7b. SRS standardization
# ==============================================================

Cmin <- 10
SRS_seed <- 1988
message("Running SRS standardization...")

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

dt <- as.data.table(srs_out, keep.rownames = "taxon")
srs_long <- melt(dt, id.vars = "taxon",
                 variable.name = "sample", value.name = "abundance")

# ==============================================================
# 7c. Merge region & time info back
# ==============================================================

setDT(srs_long)
srs_long[, c("sample", "sample_group", "site.x", "end.1") :=
           tstrsplit(sample, "@", fixed = TRUE)]

region_info <- as.data.table(raw.mat)[,
                                      c("sample", "sample_group", "site.x", "end.1") :=
                                        tstrsplit(sample, "@", fixed = TRUE)
][, end.1 := as.numeric(end.1)][,
                                .(sample, sample_group, site.x, end.1, region)
]

region_info <- unique(region_info, by = c("sample", "sample_group", "site.x", "end.1"))
srs_long[, end.1 := as.numeric(end.1)]

srs_long <- merge(
  srs_long, region_info,
  by = c("sample", "sample_group", "site.x", "end.1"),
  all.x = TRUE, sort = FALSE
)


# ==============================================================
# 8. Prepare SRS and Raw data for Listdf
# ==============================================================

setDT(srs_long)
srs_long[, sample_id := paste(sample, sample_group, site.x, sep = "@")]
srs_long <- unique(srs_long, by = c("region", "end.1", "taxon", "sample_id"))
srs_long[, end.1 := as.numeric(end.1)]

raw_long <- as.data.table(raw.mat)
raw_long[, c("sample", "sample_group", "site.x", "end.1") :=
           tstrsplit(sample, "@", fixed = TRUE)]
raw_long[, sample_id := paste(sample, sample_group, site.x, sep = "@")]
raw_long <- unique(raw_long, by = c("region", "end.1", "taxon", "sample_id"))
raw_long[, end.1 := as.numeric(end.1)]

# --- Filter to only the British/Irish Isles ---
srs_long <- srs_long[region == "British/Irish Isles"]
raw_long <- raw_long[region == "British/Irish Isles"]


# ==============================================================
# 9. Build Listdf for both datasets
# ==============================================================

create_Listdf <- function(df) {
  df[, {
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
  }, by = region]$region_list
}

Listdf_raw <- create_Listdf(raw_long)
Listdf_SRS <- create_Listdf(srs_long)
names(Listdf_raw) <- names(Listdf_SRS) <- unique(raw_long$region)

# ==============================================================
# 10. Gamma Diversity Calculation and Plotting
# ==============================================================

generate_gammaDiversity_plots <-
  function(Listdf_raw, Listdf_SRS, rects,
           output_dir = here::here("analysis", "figures")) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  plot_list <- list()
  results_list <- list()

  for (region_name in names(Listdf_raw)) {
    region_data_raw <-
      Listdf_raw[[region_name]]
    region_data_srs <-
      Listdf_SRS[[region_name]]
    region_results <-
      data.frame(Time = character(), Metric = character(), Value = numeric(), Type = character())

    for (dataset_type in c("Raw", "SRS")) {
      region_data <-
        if (dataset_type == "Raw") region_data_raw
      else region_data_srs

      for (time_name in names(region_data)) {
        sublist <- region_data[[time_name]]
        if (is.null(sublist) || nrow(sublist) == 0)
          next

        comm_mat <-
          tryCatch(as.matrix(sublist), error = function(e) NULL)
        if (is.null(comm_mat) || ncol(comm_mat) == 0)
          next

        meta <-
          tryCatch(entropart::MetaCommunity(Abundances = comm_mat),
                         error = function(e) NULL)
        if (is.null(meta))
          next

        richness <-
          tryCatch(entropart::GammaDiversity(meta, q = 0), error = function(e) NA_real_)
        shannon  <-
          tryCatch(entropart::GammaDiversity(meta, q = 1), error = function(e) NA_real_)
        simpson  <-
          tryCatch(entropart::GammaDiversity(meta, q = 2), error = function(e) NA_real_)

        region_results <- rbind(
          region_results,
          data.frame(Time = time_name, Metric = "Richness", Value = richness, Type = dataset_type),
          data.frame(Time = time_name, Metric = "Shannon", Value = shannon, Type = dataset_type),
          data.frame(Time = time_name, Metric = "Simpson", Value = simpson, Type = dataset_type)
        )
      }
    }

    region_results$TimeNumeric <- suppressWarnings(as.numeric(region_results$Time))
    region_results <- region_results[!is.na(region_results$TimeNumeric), ]
    region_results$Metric <- factor(region_results$Metric,
                                    levels = c("Richness", "Shannon", "Simpson"))

    # Keep the same aesthetics as before
    plot_min <- min(region_results$TimeNumeric)
    plot_max <- max(region_results$TimeNumeric)
    rects_filtered <- rects[rects$ystart <= plot_max & rects$yend >= plot_min, , drop = FALSE]
    time_breaks <- pretty(region_results$TimeNumeric, n = 10)

    # Plot: identical style, SRS first (behind), Raw on top, both in legend
    p <- ggplot(region_results, aes(x = Value, y = TimeNumeric, color = Type)) +
      geom_rect(data = rects_filtered,
                aes(ymin = ystart, ymax = yend, xmin = -Inf, xmax = Inf, fill = col),
                inherit.aes = FALSE, alpha = 0.5) +

      # SRS points/lines first (behind)
      geom_path(data = subset(region_results, Type == "SRS"),
                aes(group = Metric), linewidth = 1, alpha = 0.8) +
      geom_point(data = subset(region_results, Type == "SRS"),
                 size = 3, alpha = 0.8) +

      # Raw points/lines (on top)
      geom_path(data = subset(region_results, Type == "Raw"),
                aes(group = Metric), linewidth = 1, linetype = "dashed") +
      geom_point(data = subset(region_results, Type == "Raw"),
                 size = 3) +

      facet_wrap(~Metric, scales = "free_x") +
      scale_color_manual(values = c("SRS" = "green4", "Raw" = "black"),
                         name = "Dataset") +
      scale_fill_jco(guide = guide_legend(reverse = TRUE)) +
      scale_y_reverse(breaks = time_breaks, labels = time_breaks) +
      labs(
        title = paste("Gamma Diversity"),
        x = "Hill number",
        y = "Time (Years BP)",
        fill = "Time Periods"
      ) +
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

    safe_name <- gsub("[^A-Za-z0-9_]", "_", region_name)
    ggsave(filename = paste0("004-gamma-diversity-raw-srs-", safe_name, ".jpg"),
           plot = p, path = output_dir, width = 3300, height = 4200, units = "px", dpi = 300)

    plot_list[[region_name]] <- p
    results_list[[region_name]] <- region_results
  }

  message("✅ Gamma diversity plots (Raw vs SRS) generated successfully.")
  return(list(plots = plot_list, results = results_list))
}

# ==============================================================
# 11. Run analysis
# ==============================================================

generate_gammaDiversity_plots(
  Listdf_raw,
  Listdf_SRS,
  rects,
  output_dir = here::here("analysis", "figures")
)
