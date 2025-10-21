############################################################
# Purpose: Summarize fossil insect data for Europe during
#          the Late Glacial - Holocene period.
#
# Methods:
#   1. Filter fossil insect data to period and sampling context
#   2. Bin samples into 500-year time slices
#   3. Normalize using abundance-weighted and non-abundance
#      (presence-based) methods (Buckland 2007)
#   4. Plot results and save figure
############################################################

# ---- Load packages ----
pacman::p_load(
  data.table, ggh4x, IRanges, tidyverse, here
)

# ---- 1. Import site and species data ----
bugs <- fread(
  here::here("analysis/data/raw_data/bugs_europe_extraction_samples_20250612.csv"),
  na.strings = c("", "NA", "NULL"),
  encoding = "UTF-8"
)

# ---- 2. Import ecocodes ----
eco <- fread(
  here::here("analysis/data/raw_data/sead_ecocodes_20250114.csv"),
  na.strings = c("", "NA", "NULL"),
  encoding = "UTF-8"
)

# ---- 3. Filter and prepare sample metadata ----
time.mat <- bugs %>%
  mutate(age_range = age_older - age_younger,
         region =
           case_when(
             between(latitude, 49.8, 62.6) & between(longitude, -12.6, 1.8) ~ "British/Irish Isles",
             TRUE ~ ""
           )) %>%
  select(sample, sample_group, site, age_range, age_older, age_younger, context, region) %>%
  filter(
    region == "British/Irish Isles",
    context == "Stratigraphic sequence",
    sample != "BugsPresence",
    between(age_older, -500, 16000),
    between(age_younger, -500, 16000),
    age_range <= 2000
  ) %>%
  distinct(sample, sample_group, site, .keep_all = TRUE) %>%
  transmute(
    sample_id = paste(sample, sample_group, site, sep = "@"),
    start = age_younger,
    end   = age_older
  )

# ---- 4. Define 500-year temporal bins ----
age_bins <- tibble(
  start = c(-500, seq(1, 15501, 500)),
  end   = seq(0, 16000, 500)
)

# ---- 5. Compute sample–bin overlaps ----
ir_samples <- IRanges(start = time.mat$start, end = time.mat$end)
ir_bins    <- IRanges(start = age_bins$start, end = age_bins$end)
ov <- findOverlaps(ir_samples, ir_bins, type = "any")

# ---- 6. Build table linking samples to bins ----
hits <- tibble(
  sample_id   = time.mat$sample_id[queryHits(ov)],
  bin_start   = age_bins$start[subjectHits(ov)],
  bin_end     = age_bins$end[subjectHits(ov)]
) %>%
  separate(col = sample_id, into = c("sample", "sample_group", "site"), sep = "@", remove = FALSE) %>%
  inner_join(bugs, by = c("sample", "sample_group"), relationship = "many-to-many")

# ---- 7. Merge ecocode classifications ----
eco.species <- eco %>%
  filter(ecocode_system == "Bugs") %>%
  inner_join(hits, by = "taxon") %>%
  transmute(
    sample_id = paste(sample, sample_group, site.x, sep = "@"),
    ecocode, abundance, taxon, end = bin_end
  ) %>%
  distinct()

# ---- 8. Helper function: Prepare ecocode data ----
prep_ecocode_data <- function(data, end_col, value_col, type_label) {
  data %>%
    select(end = {{ end_col }}, ecocode, value = {{ value_col }}) %>%
    group_by(end) %>%
    pivot_wider(
      names_from = ecocode,
      values_from = value,
      values_fn = function(x) if(type_label == "No Abundance") length(x) else sum(x),
      values_fill = 0
    ) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(across(-end, ~ .x / sum(c_across(-end), na.rm = TRUE) * 100)) %>%
    ungroup() %>%
    pivot_longer(-end, names_to = "ecocode", values_to = "value") %>%
    mutate(type = type_label)
}

# ---- 9. Compute abundance-weighted and presence-based summaries ----
eco.abund   <- prep_ecocode_data(eco.species, end_col = end, value_col = abundance, type_label = "Abundance Weighted")
eco.noabund <- prep_ecocode_data(
  eco.species %>% distinct(end, taxon, ecocode),
  end_col = end, value_col = taxon, type_label = "No Abundance"
)

# ---- 10. Filter out bins with fewer than 5 ecocodes ----
filter_bins <- function(df) {
  valid_bins <- df %>%
    filter(value > 0) %>%
    group_by(end) %>%
    summarise(count = n_distinct(ecocode), .groups = "drop") %>%
    filter(count > 5) %>%
    select(end)

  df %>% inner_join(valid_bins, by = "end")
}

eco.abund   <- filter_bins(eco.abund)
eco.noabund <- filter_bins(eco.noabund)

# ---- 11. Combine abundance and presence-based data ----
eco_combined <- bind_rows(eco.abund, eco.noabund)
eco_combined$type <- factor(eco_combined$type, levels = c("Abundance Weighted", "No Abundance"))

# ---- 12. Define ecocode order and colors ----
ecocode_levels <- c(
  "Aquatics", "Indicators: Running water", "Indicators: Standing water",
  "Open wet habitats", "Wetlands/marshes", "Mould beetles", "Halotolerant",
  "Carrion", "General synanthropic", "Stored grain pest",
  "Dung/foul habitats", "Pasture/Dung", "Indicators: Dung", "Disturbed/arable",
  "Sandy/dry disturbed/arable", "Heathland & moorland", "Wood and trees",
  "Indicators: Deciduous", "Indicators: Coniferous", "Dry dead wood",
  "Meadowland"
)

colors <- c(
  "#197EC0FF", "lightskyblue1", "#71D0F5FF", "#709AE1FF", "dodgerblue4", "#370335FF", "#1A9993FF",
  "#8A9197FF", "#FD8CC1FF", "#D2AF81FF", "#91331FFF", "#FED439FF", "darkred", "salmon1",
  "tomato3", "#D5E4A2FF", "darkslategrey", "#46732EFF", "green4", "brown4", "yellowgreen"
)

eco_combined$ecocode <- factor(eco_combined$ecocode, levels = ecocode_levels)

# ---- 13. Generate combined horizontal figure ----
fig <- ggplot(eco_combined, aes(x = end, y = value, fill = ecocode)) +
  geom_bar(
    stat = "identity",
    position = "stack",
    alpha = 0.9,
    just = 0
    ) +
  scale_fill_manual(values = colors, name = "Habitat") +
  guides(fill = guide_legend(nrow = 5, reverse = TRUE)) +
  scale_x_reverse(limits = c(16000, -500), breaks = scales::pretty_breaks(n = 10)) +
  coord_flip() +
  ggh4x::facet_grid2(~type, scales = "free", axes = "x") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12)
  ) +
  labs(y = "%SumRep", x = "Time (Years BP)")

# ---- 14. Save figure ----
ggsave(
  "005-ecocodes-britain.jpg",
  fig,
  device = "jpg",
  here::here("analysis", "figures"),
  width=40,
  height=50,
  units = "cm",
  dpi = 300)

message("✅ Ecological composition analysis completed and figure saved: 0045ecocodes-britain.jpg")
