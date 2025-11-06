# ==============================================================
# NMDS of prehistoric insect communities by region
# Sample-level, 10,000 – -500 BP
# Regions: British Isles, Western Europe, Southern Europe
# ==============================================================

pacman::p_load(
  cowplot, data.table, vegan, ggh4x, ggplot2,
  ggsci, ggrepel, tidyverse, here
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
    between(latitude, 46.2, 54.2) & between(longitude, -4.8, 8.8) & !(between(latitude, 49.8, 62.6) & between(longitude, -12.6, 1.8)) ~ "Western Europe",
    country == "France" & latitude < 46.2 ~ "Southern France",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(region))

# ---- 4. Aggregate abundances per sample × taxon ----
sample_comm <- bugs %>%
  group_by(sample_id, taxon, region) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop")

# ---- 5. Pivot to sample × taxa matrix ----
comm <- sample_comm %>%
  pivot_wider(names_from = taxon, values_from = abundance, values_fill = 0)

# ---- 6. Metadata ----
meta <- comm %>% select(sample_id, region)
Y <- comm %>% select(-region) %>% column_to_rownames("sample_id")

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

# ---- 11. PERMANOVA ----
permanova_res <- adonis2(Y_rel ~ region, data = meta, method = "bray", permutations = 999)
print(permanova_res)

# ---- 12. Check homogeneity of dispersion ----
bd <- betadisper(dist_mat, meta$region)
anova_res <- anova(bd)

# ---- 13. Convex hulls for NMDS (2D: NMDS1 vs NMDS2) ----
get_hull <- function(df, group_col, x_col = "NMDS1", y_col = "NMDS2") {
  group_sym <- rlang::sym(group_col)
  df %>%
    dplyr::group_split(!!group_sym) %>%
    purrr::map_dfr(function(g) {
      if(nrow(g) > 2) {
        hull_idx <- chull(g[[x_col]], g[[y_col]])
        g[hull_idx, , drop = FALSE]
      } else {
        g
      }
    })
}

hull_data <- get_hull(mds_points, "region")

# ---- 14. subtitle ----
subtitle_text <- paste0(
  "Stress = ", round(mds$stress, 3),
  "; PERMANOVA F = ", round(permanova_res$F[1], 2),
  ", p = ", permanova_res$`Pr(>F)`[1],
  "; Anova F = ", round(anova_res$F[1], 2),
  ", p = ", round(anova_res$`Pr(>F)`[1])
)

# ---- 14. Plot NMDS ----
ggplot(mds_points, aes(x = NMDS1, y = NMDS2, color = region)) +
  geom_polygon(data = hull_data, aes(fill = region, group = region), alpha = 0.3, color = NA) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = sample_id), size = 3, max.overlaps = 10) +
  scale_color_manual(values = c("British Isles" = "#2166AC",
                                "Western Europe" = "#B2182B",
                                "Southern France" = "#4DAF4A")) +
  scale_fill_manual(values = c("British Isles" = "#2166AC",
                               "Western Europe" = "#B2182B",
                               "Southern France" = "#4DAF4A")) +
  theme_minimal(base_size = 14) +
  labs(title = "NMDS of fossil insect communities by region (sample-level)",
       subtitle = subtitle_text,
       x = "NMDS1", y = "NMDS2", color = "Region", fill = "Region") +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 16))

# ---- 15. Save figure ----
ggsave(
  filename = here("analysis", "figures", "NMDS_samples_by_region.jpg"),
  width = 50, height = 50, units = "cm", dpi = 300
)

message("✅ NMDS completed and sample-level regional comparison plot saved.")
