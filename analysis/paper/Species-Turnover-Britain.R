############################################################
# Purpose: Calculate temporal turnover (total, appearance, disappearance)
#          in fossil insect assemblages across 500-year bins.
# Methods:
#   1. Bin samples into temporal intervals using IRanges
#   2. Aggregate taxon abundances within each time bin
#   3. Applie codyn::turnover() to compute turnover metrics
#   4. Visualizes turnover trajectories with background ecological phases from Pilotto et al. 2022
############################################################

# ---- Load packages ----
pacman::p_load(
  codyn, data.table, ggsci, ggplot2, IRanges, tidyverse, vegan, here
)

# ---- 0. Parameters ----
age_min   <- -500
age_max   <- 16000
bin_width <- 500

# ---- 1. Read fossil insect data ----
bugs <-
  fread(
  here("analysis/data/raw_data/bugs_europe_extraction_samples_20250612.csv"),
  na.strings = c("", "NA", "NULL"),
  encoding = "UTF-8"
)

# ---- 2. Prepare sample-level time ranges ----
time_mat <-
  bugs %>%
  select(country, sample, site, sample_group, age_older, age_younger, context, latitude, longitude) %>%
  mutate(age_range = age_older - age_younger,
         region =
           case_when(
             between(latitude, 49.8, 62.6) & between(longitude, -12.6, 1.8) ~ "British/Irish Isles",
             TRUE ~ ""
           )) %>%
  filter(
    context == "Stratigraphic sequence",
    sample != "BugsPresence",
    between(age_older, age_min, age_max),
    between(age_younger, age_min, age_max),
    age_range <= 2000,
    country != "Greenland",
    region == "British/Irish Isles"
  ) %>%
  distinct() %>%
  select(-age_range, -context, -latitude, -longitude) %>%
  mutate(sample = paste(sample, sample_group, site, sep = "@")) %>%
  column_to_rownames("sample") %>%
  select(-country, -site, -sample_group) %>%
  dplyr::rename(start = age_younger, end = age_older)

# ---- 3. Define temporal bins ----
bins <- tibble(
  start = c(age_min, seq(1, age_max - bin_width + 1, by = bin_width)),
  end   = seq(0, age_max, by = bin_width)
)

# ---- 4. Find overlapping samples per temporal bin ----
hits <- findOverlaps(
  query   = do.call(IRanges, as.list(time_mat)),
  subject = do.call(IRanges, bins),
  type    = "any"
)

# ---- 5. Combine hits with metadata ----
hits_df <- data.frame(
  time_mat[queryHits(hits),],
  bins[subjectHits(hits),]
) %>%
  as_tibble(rownames = "sample") %>%
  separate(sample, into = c("sample", "sample_group", "site"), sep = "@") %>%
  inner_join(bugs, by = c("sample", "sample_group"), relationship = "many-to-many") %>%
  select(-start, -end)

# ---- 6. Aggregate taxon abundance per time bin ----
raw_abund <- hits_df %>%
  select(site.x, sample_group, sample, taxon, end.1, abundance) %>%
  distinct() %>%
  group_by(taxon, end.1) %>%
  summarise(tot_abund = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  # Ensure chronological consistency (negative years before present)
  mutate(end.1 = if_else(end.1 != 0, end.1 * -1, end.1))

# ---- 7. Compute turnover metrics ----
metrics <- c("total", "appearance", "disappearance")

spec_turnover <- map_dfr(metrics, function(met) {
  turnover(df = raw_abund,
           time.var = "end.1",
           species.var = "taxon",
           abundance.var = "tot_abund",
           metric = met) %>%
    mutate(
      end.1 = if_else(end.1 != 0, end.1 * -1, end.1),
      value = .data[[met]],
      type = str_to_title(met)
    ) %>%
    select(end.1, value, type)
})

# ---- 8. Background rectangles (after Pilotto et al. 2022) ----
rects <- tibble(
  xstart = c(12000, 9500, 8000, 6500, 4500, 4000, 3500, -500),
  xend   = c(16000, 12000, 9500, 8000, 6500, 4500, 4000, 3500),
  col    = factor(paste("TM", 1:8))
)

# ---- 9. Plot species turnover ----
fig <- ggplot() +
  # Background ecological periods
  geom_rect(
    data = rects,
    aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = col),
    alpha = 0.5
  ) +
  geom_line(
    data = spec_turnover,
    aes(x = end.1, y = value, linetype = type, color = type),
    linewidth = 1
  ) +
  geom_point(
    data = spec_turnover,
    aes(x = end.1, y = value, color = type),
    size = 2
  ) +

  scale_x_reverse(limits = c(16000, -500), breaks = scales::pretty_breaks(n = 10)) +
  scale_color_manual(
    name = "Metric",
    values = c("green4", "red4", "black")
  ) +
  scale_linetype_manual(
    name = "Metric",
    values = c("longdash", "dotted", "solid")
  ) +
  scale_fill_jco(name = "Time Periods") +
  guides(fill = guide_legend(reverse = TRUE)) +
  coord_flip() +
  labs(
    x = "Age (BP)",
    y = "Relative species turnover"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 12, face = "bold", colour = "black"),
    legend.title = element_text(size = 14, face = "bold", colour = "black"),
    legend.text = element_text(size = 12),
    strip.text.x = element_text(size = 14, face = "bold", colour = "black")
  )

# ---- 10. Save figure ----
ggsave(
  "005-species-turnover-britain.jpg",
  fig,
  device = "jpg",
  path = here("analysis", "figures"),
  width = 2400,
  height = 3600,
  units = "px",
  dpi = 300
)

message("✅ Species turnover analysis complete and figure saved: '005-species-turnover-britain.jpg'")
