############################################################
# Purpose: Perform regional ecological composition analysis
#          of fossil insect data using Bugs ecocodes
# Methods:
#   1. Bin samples into 500-year time intervals
#   2. Calculate proportional ecocode composition per bin
#   3. Normalize using abundance-weighted and non-abundance
#      (presence-based) methods (Buckland 2007)
#   4. Visualizes regionally using stacked bar plots
############################################################

# ---- Load packages ----
pacman::p_load(data.table, ggh4x, ggsci, IRanges, procs, tidyverse, vegan)

# ---- 1. Import site and species data ----
bugs <-
  fread(here::here("analysis/data/raw_data/bugs_europe_extraction_samples_20250612.csv"), na.strings = c("", "NA", "NULL"), encoding = "UTF-8")

# ---- 2. Import ecocodes ----
eco <-
  fread(here::here("analysis/data/raw_data/sead_ecocodes_20250114.csv"), na.strings = c("", "NA", "NULL"), encoding = "UTF-8")

# ---- 3. Filter and structure sample-level metadata ----
# Keep only stratigraphic sequence samples within the target time range,
# create combined sample identifiers, and compute age intervals.
time.mat <- bugs %>%
  select(country, sample, site, sample_group, age_older, age_younger, context) %>%
  mutate(age_range = age_older - age_younger) %>%
  filter(context == "Stratigraphic sequence",
         sample != "BugsPresence",
         between(age_older, -500, 16000),
         between(age_younger, -500, 16000),
         age_range <= 2000,
         country != "Greenland") %>%
  distinct(sample, sample_group, site, .keep_all = TRUE) %>%
  transmute(sample = paste(sample, sample_group, site, sep = "@"),
            start = age_younger,
            end   = age_older)

# ---- 4. Define 500-year temporal bins ----
range <- tibble(start = c(-500, seq(1, 15501, 500)),
                end = seq(0, 16000, 500))

# ---- 5. Compute sample–time bin overlaps ----
# Use IRanges to detect which samples overlap with which time bins.
ir_samples <-
  IRanges(start = time.mat$start, end = time.mat$end)
ir_bins <-
  IRanges(start = range$start, end = range$end)
ov <-
  findOverlaps(ir_samples, ir_bins, type = "any")

# ---- 6. Build table linking samples to time bins ----
# Includes site, sample_group, and bin metadata.
hits <- tibble(
  sample     = time.mat$sample[queryHits(ov)],
  sample_start = time.mat$start[queryHits(ov)],
  sample_end   = time.mat$end[queryHits(ov)],
  bin_start    = range$start[subjectHits(ov)],
  bin_end      = range$end[subjectHits(ov)]
) %>%
  # split the sample id back into components
  separate(col = "sample", into = c("sample", "sample_group", "site"), sep = "@", remove = FALSE) %>%
  # now join original bugs csv (adjust the join keys if needed)
  inner_join(bugs, by = c("sample", "sample_group"), relationship = "many-to-many") %>%
  select(-sample_start, -sample_end)

# ---- 7. Merge ecocode classifications and assign geographic regions ----
eco.species <- eco %>%
  filter(ecocode_system == "Bugs") %>%
  inner_join(hits, by = "taxon") %>%
  transmute(
    sample = paste(sample, sample_group, site.x, sep = "@"),
    ecocode, abundance, taxon, end = bin_end,
    region = case_when(
      between(latitude, 53.5, 74.8) & between(longitude, 2.8, 41) ~ "Scandinavia/Baltic",
      between(latitude, 49.8, 62.6) & between(longitude, -12.6, 1.8) ~ "British/Irish Isles",
      between(latitude, 63, 66.3) & between(longitude, -25, -12) ~ "Iceland",
      latitude < 45 ~ "Meridional",
      TRUE ~ "Continental"
    )
  ) %>%
  distinct()

# ---- 8. Helper function: Prepare ecocode data per region ----
# Creates proportional ecocode composition for each time bin.
# type_label = "Abundance Weighted" (sum) or "No Abundance" (count only).
prep_ecocode_data <- function(data, end_col, value_col, group_cols, type_label) {
  data %>%
    select(end = {{ end_col }}, all_of(group_cols), value = {{ value_col }}) %>%
    group_by(end) %>%
    pivot_wider(
      id_cols = end,
      names_from = all_of(group_cols),
      values_from = value,
      values_fn = function(x) if (type_label == "No Abundance") length(x) else sum(x),
      values_fill = 0
    ) %>%
    ungroup() %>%
    # Normalize each time bin to sum to 100%
    rowwise() %>%
    mutate(across(-end, ~ .x / sum(c_across(-end), na.rm = TRUE) * 100)) %>%
    ungroup() %>%
    pivot_longer(-end, names_to = "name", values_to = "value") %>%
    mutate(type = type_label)
}

# ---- 9. Compute abundance-weighted and non-abundance summaries ----
abund.list <- map(
  split(eco.species, eco.species$region),
  ~ prep_ecocode_data(
    data = .x,
    end_col = end,
    value_col = abundance,
    group_cols = "ecocode",
    type_label = "Abundance Weighted"
  )
)

noabund.list <- map(
  split(eco.species, eco.species$region),
  ~ prep_ecocode_data(
    data = .x %>% distinct(end, taxon, ecocode),
    end_col = end,
    value_col = taxon,
    group_cols = "ecocode",
    type_label = "No Abundance"
  )
)

# ---- 10. Combine and filter by region ----
eco.abund <-
  bind_rows(abund.list, .id = "id") %>% distinct()

eco.noabund <-
  bind_rows(noabund.list, .id = "id") %>% distinct()

# ---- 11. Filter out bins with fewer than 5 ecocodes ----
filter_bins <- function(df) {
  valid_bins <- df %>%
    filter(value > 0) %>%
    group_by(id, end) %>%
    summarise(count = n_distinct(name), .groups = "drop") %>%
    filter(count > 5) %>%
    select(id, end)

  df %>%
    inner_join(valid_bins, by = c("id", "end"))
}

eco.abund  <-
  filter_bins(eco.abund)
eco.noabund <-
  filter_bins(eco.noabund)

# ---- 12. Define ecocode plotting order ----
ecocode_levels <- c(
  "Aquatics", "Indicators: Running water", "Indicators: Standing water",
  "Open wet habitats", "Wetlands/marshes", "Mould beetles", "Halotolerant",
  "Carrion", "Ectoparasite", "General synanthropic", "Stored grain pest",
  "Dung/foul habitats", "Pasture/Dung", "Indicators: Dung", "Disturbed/arable",
  "Sandy/dry disturbed/arable", "Heathland & moorland", "Wood and trees",
  "Indicators: Deciduous", "Indicators: Coniferous", "Dry dead wood",
  "Meadowland"
)

# Apply factor levels to ensure consistent stacking order
eco.abund$name   <-
  factor(eco.abund$name, levels = ecocode_levels)
eco.noabund$name <-
  factor(eco.noabund$name, levels = ecocode_levels)

# ---- 13. Define color palette ----
colors <-
  c("#197EC0FF", "lightskyblue1", "#71D0F5FF", "#709AE1FF", "dodgerblue4", "#370335FF", "#1A9993FF",
    "#8A9197FF", "darkorchid", "#FD8CC1FF", "#D2AF81FF", "#91331FFF", "#FED439FF", "darkred", "salmon1",
    "tomato3", "#D5E4A2FF", "darkslategrey", "#46732EFF", "green4", "brown4", "yellowgreen")

# ---- 14. Shared ggplot settings ----
plot_base <- list(
  geom_bar(position = "stack", stat = "identity", alpha = 0.9),
  guides(fill = guide_legend(nrow = 5, reverse = TRUE)),
  scale_x_reverse(limits = c(16000, -500), breaks = scales::pretty_breaks(10)),
  coord_flip(),
  facet_grid2(~id, scales = "free", axes = "margins", labeller = labeller(id = label_value)),
  theme_bw(),
  theme(
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 14, face = "bold"),
    strip.text.x = element_text(size = 16, face = "bold", colour = "black"),
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12)
  )
)

# ---- 15. Plot abundance-weighted ecocodes ----
fig.abund <- ggplot(eco.abund, aes(fill = name, y = value, x = end)) +
  geom_bar(position = "stack", stat = "identity", alpha = 0.9) +
  scale_fill_manual(values = colors, name = "Habitat") +
  labs(y = "Abundance weighted; %SumRep", x = "Age (BP)") +
  plot_base[-1]   # exclude geom_bar duplication

# ---- 16. Plot presence-based ecocodes ----
fig.noabund <- ggplot(eco.noabund, aes(fill = name, y = value, x = end)) +
  geom_bar(position = "stack", stat = "identity", alpha = 0.9) +
  scale_fill_manual(values = colors, name = "Habitat") +
  labs(y = "No Abundance; %SumRep", x = "Age (BP)") +
  plot_base[-1] # exclude geom_bar duplication

# ---- 17. Save plots ----
ggsave(here::here("analysis", "figures", "004-ecocodes-natural-divided-abund.jpg"),
       fig.abund,
       width=50,
       height=30,
       units="cm",
       dpi=300
)

ggsave(here::here("analysis", "figures", "004-ecocodes-natural-divided-noabund.jpg"),
       fig.noabund,
       width=50,
       height=30,
       units="cm",
       dpi=300
)

message("✅ Ecological composition analysis complete and figures saved: '004-ecocodes-natural-divided-abund.jpg'' and '004-ecocodes-natural-divided-noabund.jpg'")
