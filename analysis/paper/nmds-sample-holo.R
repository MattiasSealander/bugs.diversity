# ==============================================================
# NMDS of prehistoric insect communities by region
# Sample-level, 10,000 – -500 BP
# Regions: British Isles, Western Europe, Southern Europe
# ==============================================================

pacman::p_load(
  cowplot, data.table, vegan, ggh4x, ggplot2,
  ggsci, ggnewscale, tidyverse, here
)

# ---- 1. Load fossil insect data ----
bugs <- fread(
  here("analysis", "data", "raw_data", "bugs_europe_extraction_samples_20250612.csv"),
  na.strings = c("", "NA", "NULL"),
  encoding = "UTF-8"
)

# ---- 2. Filter valid samples & compute mean age ----
bugs <- bugs %>%
  filter(
    context == "Stratigraphic sequence",
    sample != "BugsPresence",
    between(age_older, -500, 11700),
    between(age_younger, -500, 11700),
    (age_older - age_younger) <= 2000,
    country != "Greenland"
  ) %>%
  mutate(
    mean_age = (age_older + age_younger) / 2,
    sample_id = paste(sample, sample_group, site, sep = "|")
  ) %>%
  distinct()

# ---- 3. Assign mutually exclusive regions ----
bugs <- bugs %>%
  mutate(region = case_when(
    between(latitude, 49.8, 62.6) & between(longitude, -12.6, 1.8) ~ "British Isles",
    between(latitude, 46.2, 54.2) & between(longitude, -4.8, 8.8) &
      !(between(latitude, 49.8, 62.6) & between(longitude, -12.6, 1.8)) ~ "Western Europe",
    country == "France" & latitude < 46.2 ~ "Southern France",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(region))

# ---- 4. Aggregate abundances per sample × taxon ----
sample_comm <- bugs %>%
  group_by(site, sample_id, taxon, region) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop")

# ---- 5. Pivot to sample × taxa matrix ----
comm <- sample_comm %>%
  pivot_wider(names_from = taxon, values_from = abundance, values_fill = 0)

# ---- 6. Metadata ----
meta <- comm %>% select(site, sample_id, region)
Y <- comm %>% select(-region, -site) %>% column_to_rownames("sample_id")

# ---- 7. Clean & relative abundance ----
Y[is.na(Y)] <- 0
Y <- Y[rowSums(Y) > 0, colSums(Y) > 0]
meta <- meta %>% filter(sample_id %in% rownames(Y))
stopifnot(all(rownames(Y) == meta$sample_id))

Y_rel <- decostand(Y, method = "hellinger")

# ---- 8. Identify and remove dominant outliers (>50% for one taxon) ----
dom <- apply(Y_rel, 1, max)
outliers <- names(dom[dom > 0.5])
if(length(outliers) > 0){
  Y_rel <- Y_rel[!rownames(Y_rel) %in% outliers, ]
  meta  <- meta[!meta$sample_id %in% outliers, ]
}

# ---- 9. Bray–Curtis distance & NMDS ----
dist_mat <- vegdist(Y_rel, method = "bray")
set.seed(123)
mds <- metaMDS(dist_mat, k = 3, autotransform = FALSE, trymax = 50)

# ---- 10. NMDS points with metadata ----
mds_points <- as.data.frame(scores(mds, display = "sites")) %>%
  rownames_to_column("sample_id") %>%
  left_join(meta, by = "sample_id")

# ---- 11. Create site_label for plotting ----
# British Isles grouped, others remain separate
mds_points <- mds_points %>%
  mutate(site_label = if_else(region == "British Isles", "British Isles", site))

# Define factor levels: British Isles first, then Western Europe, then Southern France
site_levels <- c(
  "British Isles",
  mds_points %>% filter(region == "Western Europe") %>% distinct(site_label) %>% pull(site_label),
  mds_points %>% filter(region == "Southern France") %>% distinct(site_label) %>% pull(site_label)
)
mds_points$site_label <- factor(mds_points$site_label, levels = site_levels)

# ---- 12. PERMANOVA ----
permanova_res <- adonis2(Y_rel ~ region, data = meta, method = "bray", permutations = 999)
print(permanova_res)

# ---- 13. Check homogeneity of dispersion ----
bd <- betadisper(dist_mat, meta$region)
anova_res <- anova(bd)

# ---- 14. Convex hulls ----
get_hull <- function(df, group_col, x_col = "NMDS1", y_col = "NMDS2") {
  group_sym <- rlang::sym(group_col)
  df %>%
    dplyr::group_split(!!group_sym) %>%
    purrr::map_dfr(function(g) {
      if (nrow(g) > 2) {
        hull_idx <- chull(g[[x_col]], g[[y_col]])
        g[hull_idx, , drop = FALSE]
      } else {
        g
      }
    })
}

hull_data <- get_hull(mds_points, "region", "NMDS1", "NMDS2")

# ---- 15. Subtitle ----
subtitle_text <- paste0(
  "Stress = ", round(mds$stress, 3),
  "; PERMANOVA F = ", round(permanova_res$F[1], 2),
  ", p = ", permanova_res$`Pr(>F)`[1],
  "; Anova F = ", round(anova_res$F[1], 2),
  ", p = ", round(anova_res$`Pr(>F)`[1])
)
# ---- 16. Define fills/colors/shapes ----
region_hull_fills <- c(
  "British Isles"   = "#2166AC44",
  "Western Europe"  = "#B2182B44",
  "Southern France" = "#4DAF4A44"
)

region_point_colors <- c(
  "British Isles"   = "#08306B",
  "Western Europe"  = "#B2182B",
  "Southern France" = "#00441B"
)

# ---- 17. Safe shape palette (fillable 21–25) ----
safe_shapes <- c(21, 22, 23, 24, 25)

# ---- 18. Assign shapes by region ----
# Southern France unique filled shapes
southern_sites <- mds_points %>%
  filter(region == "Southern France") %>%
  distinct(site) %>%
  arrange(site) %>%
  pull(site)

southern_shapes <- setNames(rep(safe_shapes, length.out = length(southern_sites)), southern_sites)

# Western Europe unique outline shapes
western_sites <- mds_points %>%
  filter(region == "Western Europe") %>%
  distinct(site) %>%
  arrange(site) %>%
  pull(site)

western_shapes <- setNames(rep(safe_shapes, length.out = length(western_sites)), western_sites)

# ---- 19. Assign shape_final, fill_final, color_final ----
mds_points <- mds_points %>%
  mutate(
    shape_final = case_when(
      region == "British Isles"   ~ 21L,  # filled circle
      region == "Southern France" ~ as.integer(southern_shapes[site]),
      region == "Western Europe"  ~ as.integer(western_shapes[site]),
      TRUE                        ~ 16L
    ),
    fill_final = case_when(
      region == "British Isles"   ~ "#08306B",
      region == "Southern France" ~ "#A6D96A",
      TRUE ~ NA_character_
    ),
    color_final = case_when(
      region == "Western Europe"  ~ "#B2182B",
      region == "British Isles"   ~ "black",
      region == "Southern France" ~ "#00441B",
      TRUE ~ "black"
    )
  ) %>%
  mutate(region = factor(region, levels = c("British Isles", "Southern France", "Western Europe"))) %>%
  arrange(region)

# ---- 20. Convex hulls ----
hull_data <- get_hull(mds_points, "region", "NMDS1", "NMDS2")

hull_british <- hull_data %>% filter(region == "British Isles")
hull_western <- hull_data %>% filter(region == "Western Europe")
hull_southern <- hull_data %>% filter(region == "Southern France")

# ---- 21. Legend shape mapping ----
shape_values <- mds_points %>%
  distinct(site_label, shape_final) %>%
  arrange(factor(site_label, levels = unique(mds_points$site_label))) %>%
  deframe()

# ---- 22. Plot NMDS ----
# 2️⃣ Build plot
ggplot() +
  # Hulls (same order as above)
  geom_polygon(
    data = hull_british,
    aes(x = NMDS1, y = NMDS2, group = region),
    fill = region_hull_fills["British Isles"], alpha = 0.4
  ) +
  geom_polygon(
    data = hull_southern,
    aes(x = NMDS1, y = NMDS2, group = region),
    fill = region_hull_fills["Southern France"], alpha = 0.4
  ) +
  geom_polygon(
    data = hull_western,
    aes(x = NMDS1, y = NMDS2, group = region),
    fill = region_hull_fills["Western Europe"], alpha = 0.4
  ) +
  ggnewscale::new_scale_fill() +

  # 3️⃣ Points (same shapes/fills/colors as before, just drawn in new order)
  geom_point(
    data = mds_points,
    aes(
      x = NMDS1, y = NMDS2,
      shape = site_label,
      fill = fill_final,
      color = color_final
    ),
    size = 3.8, stroke = 1.2
  ) +

  # 4️⃣ Keep the shape mapping but hide it from the legend
  scale_shape_manual(
    values = shape_values,
    guide = "none"  # removes site legend
  ) +
  scale_fill_identity(guide = "none") +
  scale_color_identity(guide = "none") +

  # 5️⃣ Add a clean region-color legend
  ggnewscale::new_scale_fill() +
  geom_point(
    data = distinct(mds_points, region),
    aes(x = Inf, y = Inf, fill = region),
    shape = 21, size = 5
  ) +
  scale_fill_manual(
    name = "Region",
    values = c(
      "British Isles"   = alpha("#08306B", .5),
      "Western Europe"  = alpha("#B2182B", .5),
      "Southern France" = alpha("#00441B", .5)
    )
  ) +

  theme_minimal(base_size = 14) +
  labs(
    title = "NMDS of fossil insect communities by region (sample-level)",
    subtitle = subtitle_text,
    x = "NMDS1", y = "NMDS2"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 16)
  )
