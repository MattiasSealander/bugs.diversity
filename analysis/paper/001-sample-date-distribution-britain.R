############################################################
# Purpose: Summarize fossil insect data for Europe during
#          the Late Glacial - Holocene period.
#
# Methods:
#   1. Filter fossil insect data to period and region of study
#   2. Summarize sites, samples and total abundance of insects
#   3. Plot summarized data over time
############################################################

# ---- Load packages ----
pacman::p_load(tidyverse, IRanges, cowplot, here, data.table)

# ==============================================================
# 1. Import species occurrence data
# --------------------------------------------------------------
# The raw dataset includes sample-level fossil insect occurrences across Europe from SEAD.
# ==============================================================
bugs <- fread(
  here::here("analysis/data/raw_data/bugs_europe_extraction_samples_20250612.csv"),
  na.strings = c("", "NA", "NULL"),
  encoding = "UTF-8"
)

# ==============================================================
# 2. Filter and prepare data
# --------------------------------------------------------------
# Filter for samples from natural contexts (Stratigraphic sequences) and valid age ranges.
# Samples are identified by concatenation of site, sample_group, and sample.
# ==============================================================
nat <- bugs %>%
  mutate(
    region = case_when(
      between(latitude, 49.8, 62.6) & between(longitude, -12.6, 1.8) ~ "British/Irish Isles")) %>%
  filter(
    context == "Stratigraphic sequence",
    sample != "BugsPresence",
    between(age_older, -500, 16000),
    between(age_younger, -500, 16000),
    (age_older - age_younger) <= 2000,
    country != "Greenland"
  ) %>%
  distinct(sample, sample_group, site, age_older, age_younger, region) %>%
  mutate(
    start = age_younger,
    end   = age_older
  )


# ==============================================================
# 3. Define 500-year temporal bins
# --------------------------------------------------------------
# Bins are used for assigning samples to time intervals.
# ==============================================================
age_bins <- tibble(
  start = c(-500, seq(1, 15501, by = 500)),
  end   = seq(0, 16000, by = 500)
)


# ==============================================================
# 4. Identify overlapping samples with time bins
# --------------------------------------------------------------
# Uses IRanges to assign each sample to all 500-year bins that overlap its age range.
# ==============================================================
nat_ranges <- IRanges(start = nat$start, end = nat$end)
bin_ranges <- IRanges(start = age_bins$start, end = age_bins$end)

# Find overlaps between sample age ranges and bins
overlaps <- findOverlaps(nat_ranges, bin_ranges, type = "any")


# ==============================================================
# 5. Build overlap hits tibble
# --------------------------------------------------------------
# Produces a table linking each sample to its temporal bin.
# ==============================================================
overlap_hits <- tibble(
  sample  = nat$sample[queryHits(overlaps)],
  sample_group = nat$sample_group[queryHits(overlaps)],
  site    = nat$site[queryHits(overlaps)],
  bin_start = age_bins$start[subjectHits(overlaps)],
  bin_end   = age_bins$end[subjectHits(overlaps)],
  region = nat$region[queryHits(overlaps)]
)


# ==============================================================
# 6. Summarise hits
# --------------------------------------------------------------
# Summarises the nr. of sites, samples and total abundance for each 500-bin.
# Results are calculated for Europe and the British/Irish Isles.
# ==============================================================

# 6a. Sites per bin
# Europe
sites_all <- overlap_hits %>%
  distinct(bin_end, site) %>%
  count(bin_end, name = "sites_all")
# British/Irish Isles
sites_brit <- overlap_hits %>%
  filter(region == "British/Irish Isles") %>%
  distinct(bin_end, site) %>%
  count(bin_end, name = "sites_brit")
# Combined
sites_summary <- full_join(
  sites_all, sites_brit, by = "bin_end"
) %>%
  replace_na(list(sites_all = 0L, sites_brit = 0L)) %>%
  arrange(desc(bin_end))

# 6b. Samples per bin (unique sample identity = sample + sample_group + site)
# Europe
samples_all <- overlap_hits %>%
  distinct(bin_end, sample, sample_group, site) %>%
  mutate(sample_uid = interaction(sample, sample_group, site)) %>%
  count(bin_end, name = "samples_all")
# British/Irish Isles
samples_brit <- overlap_hits %>%
  filter(region == "British/Irish Isles") %>%
  distinct(bin_end, sample, sample_group, site) %>%
  mutate(sample_uid = interaction(sample, sample_group, site)) %>%
  count(bin_end, name = "samples_brit")
# Combined
samples_summary <- full_join(
  samples_all, samples_brit, by = "bin_end"
) %>%
  replace_na(list(samples_all = 0L, samples_brit = 0L)) %>%
  arrange(desc(bin_end))

# 6c. Abundance per bin
# Europe
abundance_all <- overlap_hits %>%
  left_join(bugs, by = c("sample", "sample_group", "site"), relationship = "many-to-many") %>%
  distinct(bin_end, site, sample_group, sample, taxon, abundance) %>%
  group_by(bin_end) %>%
  summarise(abundance_all = sum(abundance, na.rm = TRUE), .groups = "drop")
# British/Irish Isles
abundance_brit <- overlap_hits %>%
  filter(region == "British/Irish Isles") %>%
  left_join(bugs, by = c("sample", "sample_group", "site"), relationship = "many-to-many") %>%
  distinct(bin_end, site, sample_group, sample, taxon, abundance) %>%
  group_by(bin_end) %>%
  summarise(abundance_brit = sum(abundance, na.rm = TRUE), .groups = "drop")
# Combined
abundance_summary <- full_join(abundance_all, abundance_brit, by = "bin_end") %>%
  replace_na(list(abundance_all = 0L, abundance_brit = 0L)) %>%
  arrange(desc(bin_end))

# 6d. Combine to single tidy table
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
# 7. Plotting function
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

# ---- 7b. Generate legend plot ----
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
# 8. Generate individual metric plots
# ==============================================================
fig.sites <- plot_metric_overlap(summary_long, "Sites", show_x = FALSE, text_x = FALSE)
fig.samples <- plot_metric_overlap(summary_long, "Samples", show_x = FALSE, text_x = FALSE)
fig.abundance <- plot_metric_overlap(summary_long, "Abundance", show_x = TRUE, text_x = TRUE)


# ---- 8b. Combine plots with legend ----
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
# 9. Filter and prepare sample age ranges
# --------------------------------------------------------------
# Filter samples from UK and Ireland and arrange in order of youngest to oldest
# Add row numbers for figure axis labels.
# ==============================================================
age_ranges <- nat %>%
  filter(region == "British/Irish Isles") %>%
  distinct(sample, sample_group, site, age_older, age_younger) %>%
  mutate(sample_uid = interaction(site, sample_group, sample, sep = "|")) %>%
  arrange((age_younger)) %>%  # order youngest to oldest
  mutate(sample_row = row_number())  # row number for y-axis


# ==============================================================
# 10. Plot sample age ranges
# ==============================================================
fig_time <-
  ggplot(age_ranges, aes(y = sample_row)) +
  geom_segment(aes(x = age_older, xend = age_younger, y = sample_row, yend = sample_row, color = age_older),
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
# 11. Combine final plot
# ==============================================================
fig_final <- plot_grid(
  fig_with_legend,
  fig_time,
  labels = c("A","B"),
  ncol = 1,
  rel_heights = c(1, 0.5)
)


# ==============================================================
# 12. Save figure
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
