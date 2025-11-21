# ==============================================================
# Purpose: Summarize fossil insect data for Europe and UK/Ireland during
#          the Late Glacial - Holocene period.
#
# Workflow:
#   1. Load fossil insect occurrence data.
#   2. Filter samples by region, context, and age constraints.
#   3. Assign samples to 500-year temporal bins using IRanges.
#   4. Summarise data at site, sample and abundance level for each bin
#   5. Generate and combine plots for A figure
#   6. Summarise sample age ranges
#   7. Generate age range plot for B figure
#   3. Combine plots using cowplot
#
# Key Notes:
#   - Each unique sample is identified by the concatenation of sample+sample_group+site
#   - Sample age ranges are based on the harmonised dates from SEAD
#
# Output:
#   - Plot: "001-sample-site-abundance-summary.jpg"
# ==============================================================

# ---- 0. Load packages ----
pacman::p_load(
  data.table, tidyverse, cowplot, here, IRanges, viridis
)

# ==============================================================
# 1. Import species occurrence data
# ==============================================================
bugs <- fread(here("analysis", "data", "raw_data", "bugs_europe_extraction_samples_20250612.csv")) %>%
  mutate(sample_id = paste(sample, sample_group, site, sep = "@")) %>%
  as_tibble()

# ==============================================================
# 2. Define temporal periods (TM1–TM8) afte Pilotto et al. (2022)
# ==============================================================
rects <- tibble(
  ystart = c(12000, 9500, 8000, 6500, 4500, 4000, 3500, -500),
  yend   = c(16000, 12000, 9500, 8000, 6500, 4500, 4000, 3500),
  col    = factor(c("TM 1", "TM 2", "TM 3", "TM 4", "TM 5", "TM 6", "TM 7", "TM 8"),
                  levels = c("TM 1", "TM 2", "TM 3", "TM 4", "TM 5", "TM 6", "TM 7", "TM 8"))
)

# ==============================================================
# 3. Filter and prepare temporal range data
# ==============================================================
time_mat <- bugs %>%
  select(country, latitude, longitude, sample_id, sample, sample_group, site, age_older, age_younger, context) %>%
  mutate(age_range = age_older - age_younger,
         mid_age = (age_older + age_younger) / 2,
         region = case_when(
           between(latitude, 49.8, 62.6) & between(longitude, -12.6, 1.8) ~ "British/Irish Isles",
           TRUE ~ "Europe")
  ) %>%
  filter(context == "Stratigraphic sequence",
         sample != "BugsPresence",
         age_range <= 2000,
         between(mid_age, -500, 16000)) %>%
  distinct() %>%
  select(region, site, sample_group, sample_id, age_younger, age_older) %>%
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
  sample_group = time_mat$sample_group[queryHits(intersection)],
  site = time_mat$site[queryHits(intersection)],
  region = time_mat$region[queryHits(intersection)],
  bin_start = range$start[subjectHits(intersection)],
  bin_end   = range$end[subjectHits(intersection)]
)

# ==============================================================
# 7. Summarise hits
# --------------------------------------------------------------
# Summarises the nr. of sites, samples and total abundance for each 500-bin.
# Results are calculated for Europe and the British/Irish Isles.
# ==============================================================

# 7a. Sites per bin
# Europe
sites_all <- hits %>%
  distinct(bin_end, site) %>%
  count(bin_end, name = "sites_all")
# British/Irish Isles
sites_brit <- hits %>%
  filter(region == "British/Irish Isles") %>%
  distinct(bin_end, site) %>%
  count(bin_end, name = "sites_brit")
# Combined
sites_summary <- full_join(
  sites_all, sites_brit, by = "bin_end"
) %>%
  replace_na(list(sites_all = 0L, sites_brit = 0L)) %>%
  arrange(desc(bin_end))

# 7b. Samples per bin (sample_id = sample + sample_group + site)
# Europe
samples_all <- hits %>%
  distinct(bin_end, sample_id) %>%
  count(bin_end, name = "samples_all")
# British/Irish Isles
samples_brit <- hits %>%
  filter(region == "British/Irish Isles") %>%
  distinct(bin_end, sample_id) %>%
  count(bin_end, name = "samples_brit")
# Combined
samples_summary <- full_join(
  samples_all, samples_brit, by = "bin_end"
) %>%
  replace_na(list(samples_all = 0L, samples_brit = 0L)) %>%
  arrange(desc(bin_end))

# 7c. Abundance per bin
# Europe
abundance_all <- hits %>%
  left_join(bugs, by = c("sample_id"), relationship = "many-to-many") %>%
  distinct(bin_end, sample_id, taxon, abundance) %>%
  group_by(bin_end) %>%
  summarise(abundance_all = sum(abundance, na.rm = TRUE), .groups = "drop")
# British/Irish Isles
abundance_brit <- hits %>%
  filter(region == "British/Irish Isles") %>%
  left_join(bugs, by = c("sample_id"), relationship = "many-to-many") %>%
  distinct(bin_end, sample_id, taxon, abundance) %>%
  group_by(bin_end) %>%
  summarise(abundance_brit = sum(abundance, na.rm = TRUE), .groups = "drop")
# Combined
abundance_summary <- full_join(abundance_all, abundance_brit, by = "bin_end") %>%
  replace_na(list(abundance_all = 0L, abundance_brit = 0L)) %>%
  arrange(desc(bin_end))

# 7d. Combine to single tidy table
summary_combined <- sites_summary %>%
  left_join(samples_summary, by = "bin_end") %>%
  left_join(abundance_summary, by = "bin_end") %>%
  arrange(desc(bin_end))

# Pivot longer for plotting convenience
summary_long <- summary_combined %>%
  pivot_longer(
    cols = -bin_end,
    names_to = c("Metric", "Scope"),
    names_sep = "_",
    values_to = "Value",
    values_transform = list(Value = as.integer)
  ) %>%
  mutate(
    # ---- Standardize metric names ----
    Metric = case_when(
      grepl("^sites", Metric) ~ "Sites",
      grepl("^samples", Metric) ~ "Samples",
      grepl("^abundance", Metric) ~ "Abundance",
      TRUE ~ Metric
    ),
    # ---- Standardize region scopes ----
    Scope = case_when(
      Scope %in% c("all", "All") ~ "All_Europe",
      Scope %in% c("brit", "british", "British/Irish Isles") ~ "British/Irish Isles",
      TRUE ~ Scope
    )
  )

# ==============================================================
# 8. Plotting function
# --------------------------------------------------------------
# Create a bar plot for a given metric (Sites, Samples, or Abundance)
# Uses the summary_long table, which includes both 'All_Europe' and 'British/Irish Isles'
# ==============================================================
plot_metric_overlap <- function(data, metric_name, show_x = TRUE, text_x = TRUE) {

  df_metric <- data %>%
    dplyr::filter(Metric == metric_name) %>%
    dplyr::mutate(Scope = factor(Scope, levels = c("All_Europe", "British/Irish Isles")))

  # Base plot
  p <- ggplot(df_metric, aes(x = bin_end, y = Value)) +
    # Background: All_Europe (gray)
    geom_bar(
      data = dplyr::filter(df_metric, Scope == "All_Europe"),
      stat = "identity",
      fill = "gray70",
      color = "black",
      alpha = 0.7,
      position = "identity",
      width = 400,
      just = 0
    ) +
    # Foreground: British/Irish Isles (green)
    geom_bar(
      data = dplyr::filter(df_metric, Scope == "British/Irish Isles"),
      stat = "identity",
      fill = "green4",
      color = "black",
      alpha = 0.9,
      position = "identity",
      width = 400,
      just = 0
    ) +
    scale_x_reverse(
      limits = c(16000, -500),
      breaks = scales::pretty_breaks(n = 10)
    ) +
    labs(
      x = if (show_x) "Time (Years BP)" else NULL,
      y = metric_name
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
      axis.text.y  = element_text(size = 10, colour = "black"),
      axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )

  # Conditionally modify x-axis text
  if (isTRUE(text_x)) {
    p <- p + theme(axis.text.x = element_text(size = 10, colour = "black", angle = 45, vjust = 0.5))
  } else {
    p <- p + theme(axis.text.x = element_blank())
  }

  return(p)
}

# ---- 8b. Generate legend plot ----
legend_plot <- ggplot(summary_long %>% filter(Metric == "Sites"), aes(x = bin_end, y = Value, fill = Scope)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    labels = c("Europe", "UK/Ireland"),
    values = c("All_Europe" = "gray70", "British/Irish Isles" = "green4"),
    name = "Region"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )

# Extract the legend
legend <- get_legend(legend_plot)


# ==============================================================
# 9. Generate individual metric plots
# ==============================================================
fig.sites <- plot_metric_overlap(summary_long, "Sites", show_x = FALSE, text_x = FALSE)
fig.samples <- plot_metric_overlap(summary_long, "Samples", show_x = FALSE, text_x = FALSE)
fig.abundance <- plot_metric_overlap(summary_long, "Abundance", show_x = TRUE, text_x = TRUE)


# ---- 9b. Combine plots with legend ----
fig_combined <- plot_grid(
  fig.sites, fig.samples, fig.abundance,
  ncol = 1, align = "v",
  rel_heights = c(.8, .8, 1),
  label_size = 14
)

# Add legend at the bottom
fig_with_legend <- plot_grid(
  fig_combined,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.05)
)


# ==============================================================
# 10. Filter and prepare sample age ranges
# --------------------------------------------------------------
# Filter samples from UK and Ireland and arrange in order of youngest to oldest
# Add row numbers for figure axis labels.
# ==============================================================
age_ranges <- time_mat %>%
  filter(region == "British/Irish Isles") %>%
  distinct(sample_id, end, start) %>%
  arrange((start)) %>%  # order youngest to oldest
  mutate(sample_row = row_number())  # row number for y-axis


# ==============================================================
# 11. Plot sample age ranges
# ==============================================================
fig_time <-
  ggplot(age_ranges, aes(y = sample_row)) +
  geom_segment(aes(x = end, xend = start, y = sample_row, yend = sample_row, color = end),
               linewidth = 3) +
  scale_x_reverse(limits = c(16000, -500), breaks = scales::pretty_breaks(n = 10)) +
  scale_color_viridis(
    name = "Time (Years BP)",
    direction = -1
  ) +
  labs(
    x = "Time (Years BP)",
    y = "Sample"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none"
  )


# ==============================================================
# 12. Combine final plot
# ==============================================================
fig_final <- plot_grid(
  fig_with_legend,
  fig_time,
  labels = c("A","B"),
  ncol = 1,
  rel_heights = c(1, 0.5)
)


# ==============================================================
# 13. Save figure
# ==============================================================
ggsave(
  filename = "001-sample-site-abundance-summary.jpg",
  plot = fig_final,
  path = here("analysis", "figures"),
  units = "cm",
  width = 30,
  height = 30,
  dpi = 300
)

message("✅ Data has been summarized and figure saved: '001-sample-site-abundance-summary.jpg'")
