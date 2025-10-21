############################################################
# Purpose: Calculate temporal changes in rank-abundance curves (RAC)
#          across 500-year bins for the British/Irish Isles insect fossil data.
# Methods: Uses IRanges to assign samples to temporal bins,
#          aggregates taxa abundances per bin, and applies
#          codyn::RAC_change() to derive compositional dynamics.
############################################################

# ---- Load packages ----
pacman::p_load(
  codyn, data.table, tidyverse, IRanges, vegan, ggh4x, ggsci, here
)

# ---- 0. Parameters ----
age_min <- -500
age_max <- 16000
bin_width <- 500

# ---- 1. Read and filter data ----
bugs <-
  fread(
  here("analysis/data/raw_data/bugs_europe_extraction_samples_20250612.csv"),
  na.strings = c("", "NA", "NULL"),
  encoding = "UTF-8"
)

bugs_filtered <-
  bugs %>%
  mutate(age_range = age_older - age_younger) %>%
  mutate(
    region =
      case_when(
        between(latitude, 49.8, 62.6) & between(longitude, -12.6, 1.8) ~ "British/Irish Isles",
        TRUE ~ ""
      )
  ) %>%
  filter(
    region == "British/Irish Isles",
    context == "Stratigraphic sequence",
    sample != "BugsPresence",
    between(age_older, age_min, age_max),
    between(age_younger, age_min, age_max),
    age_range <= 2000,
    country != "Greenland",
  )

# ---- 2. Sample metadata ----
sample_meta <-
  bugs_filtered %>%
  distinct(sample, site, sample_group, age_younger, age_older, country) %>%
  transmute(
    sample_id = paste(sample, sample_group, site, sep = "@"),
    start = age_younger,
    end = age_older
  ) %>%
  column_to_rownames("sample_id")

# ---- 3. Temporal bins ----
bins <- tibble(
  start = c(age_min, seq(1, age_max - bin_width + 1, by = bin_width)),
  end   = seq(0, age_max, by = bin_width)
)

# ---- 4. Overlap between samples and bins ----
hits <- findOverlaps(
  query = do.call(IRanges, as.list(sample_meta)),
  subject = do.call(IRanges, bins),
  type = "any"
)

# ---- 5. Combine hits with metadata ----
hits_df <- tibble(
  sample_id = rownames(sample_meta)[queryHits(hits)],
  bin_start = bins$start[subjectHits(hits)],
  bin_end   = bins$end[subjectHits(hits)]
) %>%
  separate(sample_id, into = c("sample", "sample_group", "site"), sep = "@", remove = FALSE) %>%
  inner_join(bugs_filtered, by = c("sample", "sample_group", "site"), relationship = "many-to-many") %>%
  select(sample, sample_group, site, bin_start, bin_end, taxon, abundance)

# ---- 6. Aggregate abundance by bin ----
raw_abund <-
  hits_df %>%
  mutate(end = if_else(bin_end != 0, bin_end * -1, bin_end)) %>%
  group_by(end, taxon) %>%
  summarise(tot_abund = sum(abundance, na.rm = TRUE), .groups = "drop")

# ---- 8. Compute Rank-Abundance Change (RAC) metrics ----
rac <- RAC_change(
  df = raw_abund,
  time.var = "end",
  species.var = "taxon",
  abundance.var = "tot_abund"
) %>%
  # Ensure time variable is numeric and ordered (keep negatives)
  mutate(
    end2 = as.numeric(str_remove(end2, "^-"))
  ) %>%  pivot_longer(
    cols = c(richness_change, evenness_change, rank_change, gains, losses),
    names_to = "metric",
    values_to = "value"
  )

# ---- 9. Factor levels ----
rac$metric <- factor(
  rac$metric,
  levels = c("richness_change", "rank_change", "evenness_change", "gains", "losses")
)

# ---- 10. Background rectangles (Pilotto et al. 2022) ----
rects <- tibble(
  xstart = c(12000, 9500, 8000, 6500, 4500, 4000, 3500, -500),
  xend   = c(16000, 12000, 9500, 8000, 6500, 4500, 4000, 3500),
  col    = factor(paste("TM", 1:8))
)

# ---- 11. Facet labels ----
metric_labels <- c(
  "richness_change" = "Richness change",
  "rank_change"     = "Rank change",
  "evenness_change" = "Evenness change",
  "gains"           = "Species gains",
  "losses"          = "Species losses"
)
metric_labels <- metric_labels[names(metric_labels) %in% unique(rac$metric)]

# ---- 12. Plot ----
fig <- ggplot() +
  geom_rect(
    data = rects,
    aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = col),
    alpha = 0.5
  ) +
  geom_line(
    data = rac,
    aes(x = end2, y = value),
    linetype = "dashed",
    linewidth = 1
  ) +
  geom_point(
    data = rac,
    aes(x = end2, y = value),
    shape = 16,
    size = 3
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "solid",
    color = "black",
    linewidth = 0.5,
    alpha = 0.3
  ) +
  scale_x_reverse(limits = c(16000, -500), breaks = scales::pretty_breaks(n = 10)) +
  scale_fill_jco(name = "Time Periods") +
  coord_flip() +
  facet_grid2(
    ~ metric,
    scales = "free",
    axes = "margins",
    independent = "x",
    labeller = labeller(metric = metric_labels)
  ) +
  labs(x = "Time (Years BP)", y = NULL) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    axis.title.y = element_text(size = 12, face = "bold"),
    strip.text.x = element_text(size = 14, face = "bold", colour = "black")
  )

# ---- 10. Save figure ----
ggsave("004-rank-abundance-curves-britain.jpg",
       fig,
       device = "jpg",
       here::here("analysis", "figures"),
       width=3300,
       height=4200,
       units = "px",
       dpi = 300)

message("✅ RAC analysis complete and figure saved: '0004-rank-abundance-curves-britain.jpg'")
