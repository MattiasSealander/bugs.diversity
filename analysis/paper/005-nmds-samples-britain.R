# ==============================================================
# Prehistoric insect community analysis with time-bin assignment
# ==============================================================
# Workflow:
# 1) Load & filter fossil bug samples
# 2) Construct 500-yr temporal bins using IRanges
# 3) Map samples to bins
# 4) Build community matrix (sample_bin × taxa)
# 5) Standardize to relative abundance
# 6) Bray–Curtis distance → NMDS ordination
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
    between(age_older, -500, 10000),
    between(age_younger, -500, 10000),
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

# ---- 5. Assign samples to time bins using IRanges ----
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
  mutate(timeBin = end.1)  # label bin by end-year

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

# Metadata
meta <- comm %>% select(sample_bin_id, timeBin) %>% distinct()
meta$timeBin <- factor(meta$timeBin, levels = sort(unique(meta$timeBin)), ordered = TRUE)

# Community matrix proper
Y <- comm %>%
  select(-timeBin) %>%
  column_to_rownames("sample_bin_id")

# ---- 9. Clean & align ----
Y[is.na(Y)] <- 0
Y <- Y[rowSums(Y) > 0, , drop = FALSE]
Y <- Y[, colSums(Y) > 0, drop = FALSE]

meta <- meta %>% filter(sample_bin_id %in% rownames(Y))

stopifnot(all(rownames(Y) == meta$sample_bin_id))


# ==============================================================
# NMDS WORKFLOW
# ==============================================================

# --- Relative abundance standardization ---
Y_rel <- decostand(Y, method = "total")

# Identify samples where one species >50% abundance
dom <- apply(Y_rel, 1, max)
outliers <- names(dom[dom > 0.5])
print(outliers)

# Optionally remove them:
Y_rel <- Y_rel[!rownames(Y_rel) %in% outliers, ]
meta  <- meta[!meta$sample_bin_id %in% outliers, ]

# --- Bray–Curtis distance ---
dist_mat <- vegdist(Y_rel, method = "bray")

# --- NMDS ---
set.seed(123)
mds <- metaMDS(dist_mat, k = 3, trymax = 100, autotransform = FALSE)

# ---- Scores + metadata ----
mds_points <- as.data.frame(scores(mds, display = "sites")) %>%
  rownames_to_column("sample_bin_id") %>%
  left_join(meta, by = "sample_bin_id")

# ---- Plot ----
ggplot(mds_points, aes(NMDS1, NMDS2, color = timeBin)) +
  geom_point(size = 3) +
  theme_minimal(base_size = 14) +
  labs(
    title = "NMDS of insect community composition (relative abundance)",
    x = "NMDS1", y = "NMDS2", color = "Time Bin"
  )
