# ==============================================================
# Prehistoric insect community analysis with time-bin assignment
# and PCoA ordination
# ==============================================================

pacman::p_load(
  cowplot, data.table, vegan, ggh4x, ggplot2,
  ggsci, ggrepel, tidyverse, IRanges, here
)

# ---- 1. Load fossil insect data ----
bugs <- fread(
  here("analysis", "data", "raw_data", "bugs_europe_extraction_samples_20250612.csv"),
  na.strings = c("", "NA", "NULL"),
  encoding = "UTF-8"
)

# ---- 2. Filter to region, context & valid age range ----
samples <- bugs %>%
  filter(
    context == "Stratigraphic sequence",
    sample != "BugsPresence",
    between(age_older, -500, 16000),
    between(age_younger, -500, 16000),
    (age_older - age_younger) <= 2000,
    country != "Greenland",
    between(latitude, 49.8, 62.6),
    between(longitude, -12.6, 1.8)
  ) %>%
  mutate(sample_id = paste(sample, sample_group, site, sep = "|")) %>%
  distinct()

# ---- 3. Define 500-yr temporal bins ----
range <- data.frame(
  start = c(-500, seq(1, 15501, by = 500)),
  end   = seq(0, 16000, by = 500)
)

# ---- 4. Extract sample age ranges ----
time.mat <- samples %>%
  select(sample_id, age_older, age_younger) %>%
  mutate(
    start = age_younger,
    end   = age_older
  ) %>%
  select(-age_older, -age_younger)

# ---- 5. Assign samples to time bins ----
intersection <- findOverlaps(
  query   = do.call(IRanges, time.mat[, c("start", "end")]),
  subject = do.call(IRanges, range),
  type = "any"
)

# ---- 6. Build mapping (sample × bin) ----
sample_bin_map <- data.frame(
  time.mat[queryHits(intersection), ],
  range[subjectHits(intersection), ]
) %>%
  as_tibble() %>%
  mutate(timeBin = end.1)

# ---- 7. Join with species abundances ----
sample_bin_map <- sample_bin_map %>%
  inner_join(
    samples %>% select(sample_id, taxon, abundance),
    by = c("sample_id")
  ) %>%
  mutate(sample_bin_id = paste(sample_id, timeBin, sep = "|")) %>%
  select(sample_bin_id, timeBin, taxon, abundance) %>%
  distinct()

# ---- 8. Pivot to community matrix ----
comm <- sample_bin_map %>%
  pivot_wider(names_from = taxon, values_from = abundance, values_fill = 0)

meta <- comm %>% select(sample_bin_id, timeBin) %>% distinct()
meta$timeBin <- factor(meta$timeBin, levels = sort(unique(meta$timeBin)), ordered = TRUE)

Y <- comm %>%
  select(-timeBin) %>%
  column_to_rownames("sample_bin_id")

# ---- 9. Clean & align ----
Y[is.na(Y)] <- 0
Y <- Y[rowSums(Y) > 0, , drop = FALSE]
Y <- Y[, colSums(Y) > 0, drop = FALSE]

meta <- meta %>% filter(sample_bin_id %in% rownames(Y))
stopifnot(all(rownames(Y) == meta$sample_bin_id))

meta <- meta %>%
  mutate(
    timeBin = as.numeric(as.character(timeBin)),
    period = case_when(
      timeBin >= 10500 ~ "Late Glacial",
      timeBin >= 6000 & timeBin < 10500 ~ "Early Holocene",
      timeBin < 6000 ~ "Late Holocene"
    ),
    period = factor(period, levels = c("Late Glacial", "Early Holocene", "Late Holocene"))
  )

# ==============================================================
# ✅ PCoA WORKFLOW (replacing NMDS)
# ==============================================================

# --- Relative abundance standardization ---
Y_rel <- decostand(Y, method = "total")

# Identify & remove dominant-sample outliers (>50% one species)
dom <- apply(Y_rel, 1, max)
outliers <- names(dom[dom > 0.5])
print(outliers)

#Y_rel <- Y_rel[!rownames(Y_rel) %in% outliers, ]
#meta  <- meta[!meta$sample_bin_id %in% outliers, ]

# --- Bray-Curtis distance ---
dist_mat <- vegdist(Y_rel, method = "jaccard")

# ---- Run PCoA (cmdscale) ----
pcoa <- cmdscale(dist_mat, eig = TRUE, k = 3)

# ---- Build score table ----
pcoa_points <- as.data.frame(pcoa$points) %>%
  rownames_to_column("sample_bin_id") %>%
  dplyr::rename(
    PCoA1 = !!colnames(.)[2],
    PCoA2 = !!colnames(.)[3],
    PCoA3 = !!colnames(.)[4]
  ) %>%
  left_join(meta, by = "sample_bin_id")


# ---- Variance explained ----
var_exp <- round(100 * pcoa$eig / sum(pcoa$eig), 2)


# ---- Plot ----
ggplot(pcoa_points, aes(PCoA2, PCoA3, color = period)) +
  geom_point(size = 3, alpha = 0.9) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCoA of fossil insect samples (Bray–Curtis)",
    x = "PCoA 1",
    y = "PCoA 2",
    color = "Period"
  ) +
  scale_color_manual(values = c(
    "Late Glacial"   = "#1f78b4",
    "Early Holocene" = "#33a02c",
    "Late Holocene"  = "#e31a1c"
  )) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# ==============================================================
# PERMANOVA — Do periods differ in community composition?
# ==============================================================

# Test homogeneity of multivariate dispersion (PERMANOVA assumption)
beta_test <- betadisper(dist_mat, meta$period)
anova(beta_test)         # Check if dispersion differs (ns = good)
permutest(beta_test, permutations = 999)

# Run PERMANOVA (Adonis2 recommended)
set.seed(123)
permanova <- adonis2(
  dist_mat ~ period,
  data = meta,
  permutations = 999,
  method = "bray"
)

print(permanova)
