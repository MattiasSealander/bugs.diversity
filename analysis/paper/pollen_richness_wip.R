# ==============================================================
# Script: Pollen Diversity British Isles
#
# Purpose:
#   Summarise subfossil insect data from the British/Irish Isles
#   to reconstruct past habitat composition using ecocode classifications.
#
# Workflow:
#   1. Load fossil insect occurrence data.
#   2. Filter samples by region, context, and age constraints.
#   3. Assign samples to 500-year temporal bins using IRanges.
#   3. Normalize using abundance-weighted and non-abundance
#      (presence-based) methods (Buckland 2007)
#   4. Plot results and save figure
#
# Output:
#   - Plot: "005-rank-abundance-curves-britain.jpg"
# ==============================================================

# ---- 0. Load packages ----
pacman::p_load(
  ggplot2, janitor, neotoma2, here, geojsonsf, scales, tidyverse, vegan
)

# ==============================================================
# 1. Selection polygon for UK and Ireland
# ==============================================================
sel_polygon <-
  geojsonsf::geojson_sf(
    '{
  "type": "FeatureCollection",
  "features": [
    {
      "type": "Feature",
      "properties": {},
      "geometry": {
        "coordinates": [
          [
            [
              -10.79452894745529,
              50.825889572121724
            ],
            [
              -5.557428185881591,
              49.50597275652875
            ],
            [
              1.2410005177407015,
              50.51965606512752
            ],
            [
              3.2170561943935922,
              53.783830409442956
            ],
            [
              2.6071624670322535,
              57.97385256699957
            ],
            [
              -1.6705187861174693,
              64.36047521584072
            ],
            [
              -12.453439885880101,
              61.09713585644758
            ],
            [
              -10.79452894745529,
              50.825889572121724
            ]
          ]
        ],
        "type": "Polygon"
      }
    }
  ]
}'
  )

# ==============================================================
# 2. Download data using polygon
# ==============================================================
data_selected_downloads <-
  neotoma2::get_datasets(
    loc = sel_polygon
  ) %>%
  neotoma2::filter(
    datasettype == "pollen",
    age_range_old <= 20000
  ) %>%
  # get downloads
  neotoma2::get_downloads()

# ==============================================================
# 3. Extract Sample information
# ==============================================================
data_selected_samples <-
  neotoma2::samples(data_selected_downloads) %>%
  as.data.frame() %>%
  tibble::as_tibble() %>%
  dplyr::mutate(
    datasetid_sampleid = paste0(datasetid, "_", sampleid)
  )

### Get vector of all "pollen" taxa
vec_taxa_pollen <-
  neotoma2::taxa(data_selected_downloads) %>%
  dplyr::filter(element == "pollen") %>%
  purrr::pluck("variablename") %>%
  sort()


### Get pollen counts
data_sample_pollen_counts <-
  data_selected_samples %>%
  dplyr::select("datasetid_sampleid", "value", "variablename") %>%
  # only include selected taxons
  dplyr::filter(
    variablename %in% vec_taxa_pollen
  ) %>%
  # remove duplicates
  dplyr::group_by(
    datasetid_sampleid, variablename
  ) %>%
  dplyr::summarise(
    value = sum(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(variablename) %>%
  # wide format
  tidyr::pivot_wider(
    names_from = "variablename",
    values_from = "value",
    values_fill = 0
  ) %>%
  # clean names
  janitor::clean_names()

# ==============================================================
# 4. Extract Age Information
# ==============================================================
data_sample_age <-
  data_selected_samples %>%
  dplyr::distinct(datasetid_sampleid, depth, age) %>%
  dplyr::arrange(datasetid_sampleid, age)

### Join pollen counts with age information
joined_data <-
  data_sample_pollen_counts %>%
  dplyr::left_join(
    data_sample_age,
    by = "datasetid_sampleid"
  ) %>%
  select(depth, age, everything()) %>%
  dplyr::filter(age <= 16000)

# ==============================================================
# 6. Summarise data by 200-year bins
# ==============================================================
### Filter ages
joined_data_clean <- joined_data %>%
  mutate(age = as.numeric(age)) %>%
  dplyr::filter(!is.na(age), age >= 0, age <= 16000)

### Define 200-year bin edges
edges     <- seq(0, 16000, by = 200)         # 0, 200, ..., 16000
bin_start <- edges[-length(edges)]           # 0, 200, ..., 15800
bin_end   <- edges[-1]                       # 200, 400, ..., 16000

### Assign ages to bins
joined_binned <- joined_data_clean %>%
  mutate(
    bin_idx  = findInterval(age, vec = edges, rightmost.closed = TRUE),  # 1..length(bin_end)
    bin_end  = bin_end[bin_idx],
    bin_start= bin_start[bin_idx]
  )

### Pivot taxa columns to long
id_cols <- c("datasetid_sampleid", "depth", "age", "bin_idx", "bin_start", "bin_end")

taxa_long <- joined_binned %>%
  pivot_longer(
    cols = -all_of(id_cols),
    names_to   = "taxon",
    values_to  = "count"
  ) %>%
  dplyr::filter(!is.na(count), count > 0)

### Sum counts per 200-year bin × taxon
binned_taxa_counts <- taxa_long %>%
  group_by(bin_start, bin_end, taxon) %>%
  summarise(
    total_count = sum(count, na.rm = TRUE),
    n_rows      = n(),             # number of contributing rows (optional)
    .groups     = "drop"
  )

# ==============================================================
# 7. Compute metrics
# ==============================================================
### Pivot to wide matrix
bin_taxa_matrix <- binned_taxa_counts %>%
  select(bin_end, taxon, total_count) %>%
  pivot_wider(names_from = taxon, values_from = total_count, values_fill = 0) %>%
  arrange(bin_end)

### Extract the taxa-only numeric matrix for vegan
comm_mat <- bin_taxa_matrix %>%
  select(-bin_end) %>%
  as.data.frame()

### Compute Shannon H′ and Evenness J per bin
H <- diversity(comm_mat, index = "shannon")            # Shannon H′ (ln base)
S <- specnumber(comm_mat)                              # richness (taxa > 0)
J <- ifelse(S > 0, H / log(S), NA_real_)               # Pielou's evenness J = H′ / ln(S)


### Combine into a tidy frame
bin_metrics <- bin_taxa_matrix %>%
  transmute(
    YearsBP = bin_end,   # using bin end as the time coordinate
    Shannon = H,
    Evenness = J
  ) %>%
  dplyr::filter(!is.na(YearsBP)) %>%
  arrange(YearsBP)

# ==============================================================
# 8. Plot results
# ==============================================================
# Ranges for Shannon (bottom axis) and Evenness (top axis)
Smin <- min(bin_metrics$Shannon, na.rm = TRUE)
Smax <- max(bin_metrics$Shannon, na.rm = TRUE)
Jmin <- min(bin_metrics$Evenness, na.rm = TRUE)
Jmax <- max(bin_metrics$Evenness, na.rm = TRUE)

# Guard against degenerate ranges
if (!is.finite(Smin) || !is.finite(Smax) || Smax <= Smin) {
  stop("Invalid Shannon range: check your data.")
}
if (!is.finite(Jmin) || !is.finite(Jmax) || Jmax <= Jmin) {
  stop("Invalid Evenness range: check your data.")
}

# Linear mapping: x_plot = a + b * J
b <- (Smax - Smin) / (Jmax - Jmin)
a <- Smin - b * Jmin

# 2) Create transformed x for Evenness (so both draw on the same x-axis)
bin_metrics <- bin_metrics %>%
  mutate(
    x_shannon  = Shannon,          # bottom axis
    x_evenness = a + b * Evenness  # transformed for plotting; top axis restored via sec_axis
  )

# 3) Plot both on the same panel with dual x-axes
time_breaks <- seq(0, 16000, by = 1000)

p_pollen <- ggplot(bin_metrics, aes(y = YearsBP)) +
  geom_path(aes(x = x_shannon, linetype = "Shannon H′"),
            colour = "black", linewidth = 1, na.rm = FALSE, alpha = .5) +
  geom_path(aes(x = x_evenness, linetype = "Evenness J"),
            colour = "black", linewidth = 1, na.rm = FALSE) +
  # Keep guide lines within the (padded) range to avoid unintended expansion
  #geom_hline(yintercept = seq(500, 16000, by = 500),  colour = "grey", linetype = "dashed") +
  #geom_hline(yintercept = seq(1000, 16000, by = 1000), colour = "grey") +
  scale_y_reverse(
    limits = c(16000, 0),
    breaks  = time_breaks,                 # should be within 16000–500
    labels  = time_breaks,
    expand  = expansion(add = c(750, 750)),# exact 750-year padding on both ends
    oob     = scales::oob_censor           # drop 200–500 BP data silently
  ) +
scale_x_continuous(
  name    = "Shannon diversity (H′)",
  sec.axis = sec_axis(~ (.-a)/b, name = "Shannon evenness (J)",
                      labels = label_number(accuracy = 0.01))
) +
  scale_linetype_manual(values = c("Shannon H′" = "solid", "Evenness J" = "dashed")) +
  scale_shape_manual(values = c("Shannon H′" = 16, "Evenness J" = 17)) +
  labs(
    title = "Pollen diversity",
    linetype = NULL, shape = NULL,
    y = "Time (Years BP)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x.top      = element_text(size = 11, margin = margin(b = 8)),
    axis.title.x.top     = element_text(size = 11, face = "bold", margin = margin(b = 8)),
    axis.title.x.bottom  = element_text(size = 11, face = "bold", margin = margin(b = 8)),
    axis.text.x          = element_text(size = 11, angle = 45, vjust = 0.5),
    axis.text.y          = element_blank(),
    axis.title.y         = element_blank(),
    plot.title           = element_text(size = 16, face = "bold"),
    legend.position      = "top",
    legend.key.width     = unit(1.6, "lines"),
    legend.key.height    = unit(1.2, "lines")
  )

# ==============================================================
# 9. Plot results
# ==============================================================
save(p_pollen, file = here("analysis", "data", "derived_data", "neotoma_pollen.rds"))
