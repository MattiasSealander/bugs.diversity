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
  data.table, ggh4x, IRanges, purrr, tidyverse, here
)

# ---- 1. Import site and taxon data ----
bugs <- fread(
  here::here("analysis", "data", "raw_data", "bugs_europe_extraction_samples_20250612.csv"),
  na.strings = c("", "NA", "NULL"),
  encoding = "UTF-8"
)

# ---- 2. Import taxa habitat reference data ----
eco <- fread(
  here::here("analysis", "data", "raw_data", "sead_ecocodes_20250114.csv"),
  na.strings = c("", "NA", "NULL"),
  encoding = "UTF-8"
)

# ---- 3. Filter and prepare sample metadata ----
time_mat <- bugs %>%
  select(country, latitude, longitude, sample, site, sample_group, age_older, age_younger, context) %>%
  mutate(age_range = age_older - age_younger,
         region = case_when(
           between(latitude, 49.8, 62.6) & between(longitude, -12.6, 1.8) ~ "British/Irish Isles",
           TRUE ~ "")
         ) %>%
  filter(
    region == "British/Irish Isles",
    context == "Stratigraphic sequence",
    sample != "BugsPresence",
    age_range <= 2000,
    country != "Greenland"
  ) %>%
  mutate(mid_age = (age_older + age_younger) / 2) %>%
  filter(between(mid_age, -500, 16000)) %>%
  distinct() %>%
  select(-age_range, -context, -mid_age) %>%
  mutate(sample = paste(sample, sample_group, site, sep = "@")) %>%
  column_to_rownames("sample") %>%
  select(age_older, age_younger) %>%
  dplyr::rename(start = age_younger, end = age_older)

# ---- 4. Define 500-year temporal bins ----
range <- data.frame(
  start = c(-500, seq(1, 15501, by = 500)),
  end   = seq(0, 16000, by = 500)
)


# ---- 5. Identify sample–bin overlaps ----
# Samples overlapping each time bin are assigned accordingly.
intersection <- findOverlaps(
  query = do.call(IRanges, time_mat),
  subject = do.call(IRanges, range),
  type = "any"
)


# ---- 6. Merge sample and bin info ----
# Combine sample information with corresponding time bins.
hits <- data.frame(
  time_mat[queryHits(intersection),],
  range[subjectHits(intersection),]
) %>%
  as_tibble(rownames = "sample") %>%
  separate(sample, into = c("sample", "sample_group", "site"), sep = "@") %>%
  inner_join(bugs, by = c("sample", "sample_group"), relationship = "many-to-many") %>%
  select(-start, -end)

# ---- 7. Merge ecocode classifications ----
eco_species <- eco %>%
  filter(ecocode_system == "Bugs") %>%
  inner_join(hits, by = "taxon") %>%
  transmute(
    sample_id = paste(sample, sample_group, site.x, sep = "@"),
    ecocode, abundance, taxon, bin = end.1
  ) %>%
  distinct()

# ---- 8. Helper function: Calculate ecocode reconstruction ----

# Calculate Buckland (2007) diversity measures for ecocode-based habitat reconstructions
#
# This function computes the four ecological diversity measures described by
# Buckland (2007, pp. 109–110) for palaeoecological habitat analysis.
# It accepts taxon habitat (ecocode) data and returns raw and standardized
# (within-sample) habitat representation values using both presence/absence and
# abundance data.
#
# The four methods:
#
# Presence-based
# 1A Raw: Count of taxa representing each habitat class per sample.
# 1B %SumRep: Percent of taxa representing each habitat class per sample.
# %SumRep = number of taxa in habitat class (C) in sample (S) / total taxa across all C in that S × 100
#
# Abundance-weighted:
# 2A Raw: Total individuals representing each habitat class per sample.
# 2B %SumRep: Percent of individuals representing each habitat class per sample.
# %SumRep = total sum of individual of all taxa in habitat class (C) in sample (S) / total individuals across all C in that S x 100
#
# Where:
# S = sample (in this case 500-year bin)
# C = habitat class (ecocode)
calc_ecocodes <- function(data, bin_col, ecocode_col, taxon_col, abundance_col,
                          method = c("1A","1B","2A","2B")) {

  method <- match.arg(method)

  df <- data %>%
    dplyr::select(
      bin       = {{ bin_col }},
      ecocode   = {{ ecocode_col }},
      taxon     = {{ taxon_col }},
      abundance = {{ abundance_col }}
    )

  # Treat each bin as a sample
  if (method == "1A") {
    # Raw presence: count taxon-ecocode pairs
    df_out <- df %>%
      dplyr::mutate(rep = 1) %>%
      dplyr::group_by(bin, ecocode) %>%
      dplyr::summarise(value = sum(rep), .groups = "drop") %>%
      dplyr::mutate(method = "1A Raw")

  } else if (method == "1B") {
    # %SumRep presence: count taxon-ecocode pairs, normalize per bin
    df_out <- df %>%
      distinct(bin, taxon, ecocode) %>%  # unique taxon per bin per ecocode
      group_by(bin, ecocode) %>%
      summarise(rep = n(), .groups = "drop") %>%  # count taxa in each ecocode
      group_by(bin) %>%
      mutate(value = rep / sum(rep) * 100) %>%  # normalize within bin
      ungroup() %>%
      select(bin, ecocode, value) %>%
      mutate(method = "1B %SumRep")

  } else if (method == "2A") {
    # Raw abundance: sum abundances per bin / ecocode
    df_out <- df %>%
      dplyr::group_by(bin, ecocode) %>%
      dplyr::summarise(value = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(method = "2A Raw")

  } else if (method == "2B") {
    # %SumRep abundance: sum abundances per bin / ecocode, normalize per bin
    df_out <- df %>%
      dplyr::group_by(bin, ecocode) %>%
      dplyr::summarise(rep = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
      dplyr::group_by(bin) %>%
      dplyr::mutate(value = rep / sum(rep) * 100) %>%  # normalize per bin
      dplyr::ungroup() %>%
      dplyr::select(bin, ecocode, value) %>%
      dplyr::mutate(method = "2B %SumRep")
  }

  return(df_out)
}

# ---- 9. Compute abundance-weighted and presence-based summaries ----
methods <- c("1B", "2B")

results <- map_df(
  methods,
  ~ calc_ecocodes(
    eco_species,
    bin_col = bin,
    ecocode_col = ecocode,
    taxon_col = taxon,
    abundance_col = abundance,
    method = .x
  )
)

# ---- 10. Filter out bins with fewer than 5 ecocodes ----
results_filtered <- results %>%
  dplyr::filter(value > 0) %>%
  dplyr::group_by(bin, method) %>%
  dplyr::filter(n_distinct(ecocode) >= 5) %>%
  dplyr::ungroup()

# ---- 11 Define ecocode order and colors ----
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

results_filtered$ecocode <- factor(results_filtered$ecocode, levels = ecocode_levels)

# Rename facets
results_filtered$method <- factor(
  results_filtered$method,
  levels = c("2B %SumRep","1B %SumRep"),
  labels = c("Abundance-weighted", "No Abundance")
)

# ---- 12. Generate combined horizontal figure ----
fig <- ggplot(results_filtered, aes(x = bin, y = value, fill = ecocode)) +
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
  ggh4x::facet_grid2(~method, scales = "free", axes = "x") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12)
  ) +
  labs(y = "%SumRep", x = "Time (Years BP)"
       )

# ---- 14. Save figure ----
ggsave(
  "005-ecocodes-britain-rerun.jpg",
  fig,
  device = "jpg",
  here::here("analysis", "figures"),
  width=40,
  height=50,
  units = "cm",
  dpi = 300)

message("✅ Ecological composition analysis completed and figure saved: 005-ecocodes-britain.jpg")
