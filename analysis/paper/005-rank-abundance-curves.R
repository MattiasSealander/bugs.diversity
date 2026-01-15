# ==============================================================
# Script: Rank-Abundance Curve Change Analysis (RAC)
#
# Purpose:
#   Compute temporal changes in rank-abundance curves (richness,
#   evenness, rank shifts, species gains/losses) across 500-year bins
#   for Lateglacial-Holocene subfossil insect data from UK/Ireland.
#
# Workflow:
#   1. Load fossil insect occurrence data.
#   2. Filter samples by region, context, and age constraints.
#   3. Assign samples to 500-year temporal bins using IRanges.
#   4. Aggregate abundances per bin.
#   5. Compute RAC metrics using codyn::RAC_change().
#   6. Plot results.
#
# Key Notes:
#   - RAC is calculated by making pairwise comparisons between consecutive bins,
#     moving from the oldest to the most recent 500-year bin.
#
# Output:
#   - Plot: "004-rank-abundance-curves.jpg"
# ==============================================================

# ---- 0. Load packages ----
pacman::p_load(
  codyn, data.table, tidyverse, IRanges, ggsci, here
)

# ==============================================================
# 2. Import and filter species occurrence data
# ==============================================================
bugs <- fread(here("analysis", "data", "raw_data", "bugs_europe_extraction_samples_20250612.csv"),
              na.strings = c("", "NA", "NULL"), encoding = "UTF-8") %>%
  mutate(sample_id = paste(sample, sample_group, site, sep = "@"))  # ✅ Construct sample_id early

bugs_filtered <- bugs %>%
  mutate(age_range = age_older - age_younger,
         region = case_when(
           between(latitude, 49.8, 62.6) & between(longitude, -12.6, 1.8) ~ "British/Irish Isles",
           TRUE ~ "")
  ) %>%
  filter(!family_name %in% c("ACANTHOSOMATIDAE","ANTHOCORIDAE", "APHALARIDAE", "APHIDOIDEA", "AUCHENORRHYNCHA", "BIBIONIDAE",
                             "BRACHYCENTRIDAE", "CALOPTERYGIDAE", "CERCOPIDAE", "CHIRONOMIDAE", "CICADELLIDAE", "CICADOMORPHA",
                             "CIXIIDAE", "CORIXIDAE", "CYCLORRHAPHA", "CYDNIDAE", "DELPHACIDAE", "DERMAPTERA", "DIPTERA", "FORFICULIDAE",
                             "FORMICIDAE", "FULGOROMORPHA", "GERRIDAE", "GLOSSOSOMATIDAE", "GOERIDAE", "HEBRIDAE", "HEMIPTERA",
                             "HETEROPTERA", "HOMOPTERA", "HYDROPSYCHIDAE", "HYDROPTILIDAE", "HYMENOPTERA", "LEPIDOPTERA", "LEPIDOSTOMATIDAE",
                             "LEPTOCERIDAE", "LIMNEPHILIDAE", "LYGAEIDAE", "MEMBRACIDAE", "MICROPHYSIDAE", "MICROSPORIDAE", "MIRIDAE",
                             "MOLANNIDAE", "NABIDAE", "NEMATOCERA", "NEMOURIDAE", "ODONATA", "PARASITICA", "PENTATOMIDAE", "PHRYGANEIDAE",
                             "POLYCENTROPIDAE", "PSYCHOMYIIDAE", "PSYLLIDAE", "RAPHIDIIDAE", "SALDIDAE", "SCUTELLERIDAE", "SERICOSTOMATIDAE",
                             "SIALIDAE", "TRICHOPTERA", "THYREOCORIDAE", "TINGIDAE", "TIPULIDAE", "TRIOPSIDAE", "TRIOZIDAE"),
         region == "British/Irish Isles",
         context == "Stratigraphic sequence",
         sample != "BugsPresence",
         age_range <= 2000,
         country != "Greenland") %>%
  mutate(mid_age = (age_older + age_younger) / 2) %>%
  filter(between(mid_age, -500, 16000))

message("✅ Filtered data: ", nrow(bugs_filtered), " records retained.")

# ==============================================================
# 3. Prepare sample metadata (keep sample_id as column)
# ==============================================================
sample_meta <- bugs_filtered %>%
  distinct(sample_id, age_younger, age_older) %>%
  dplyr::rename(start = age_younger, end = age_older)

# ==============================================================
# 4. Define 500-year temporal bins
# ==============================================================
range <- tibble(
  start = c(-500, seq(1, 15501, by = 500)),
  end   = seq(0, 16000, by = 500)
)

# ==============================================================
# 5. Assign samples to bins using IRanges
# ==============================================================
query <- IRanges(start = sample_meta$start, end = sample_meta$end)
subject <- IRanges(start = range$start, end = range$end)

hits <- findOverlaps(query, subject, type = "any")
if (length(hits) == 0) stop("No overlaps found. Check age ranges and bin definitions.")

hits_df <- tibble(
  sample_id = sample_meta$sample_id[queryHits(hits)],
  bin_start = range$start[subjectHits(hits)],
  bin_end   = range$end[subjectHits(hits)]
) %>%
  inner_join(bugs_filtered, by = "sample_id", relationship = "many-to-many") %>%
  select(sample_id, bin_end, taxon, abundance)

message("✅ Assigned samples to bins: ", nrow(hits_df), " records.")

# ==============================================================
# 6. Aggregate abundances per bin
# ==============================================================
raw_abund <- hits_df %>%
  mutate(bin_end = if_else(bin_end != 0, bin_end * -1, bin_end)) %>%
  group_by(bin_end, taxon) %>%
  summarise(tot_abund = sum(abundance, na.rm = TRUE), .groups = "drop")

message("✅ Aggregated abundances across bins.")

# ==============================================================
# 7. Compute RAC metrics using codyn
# ==============================================================
rac <- RAC_change(df = raw_abund,
                  time.var = "bin_end",
                  species.var = "taxon",
                  abundance.var = "tot_abund") %>%
  mutate(bin_numeric = as.numeric(str_remove(bin_end2, "^-"))) %>%
  pivot_longer(cols = c(richness_change, evenness_change, rank_change, gains, losses),
               names_to = "Metric", values_to = "Value")

rac$Metric <- factor(rac$Metric,
                     levels = c("richness_change", "rank_change", "evenness_change", "gains", "losses"))

message("✅ RAC metrics computed for ", length(unique(rac$bin_numeric)), " bins.")

# ==============================================================
# 8. Define periods
# ==============================================================
rects <- tibble(
  ystart = c(11700, 7000, 5000, 0),
  yend   = c(16000, 11700, 7000, 5000),
  col    = factor(c("Lateglacial", "Early Holocene", "Mid-Holocene", "Late Holocene"),
                  levels = c("Lateglacial", "Early Holocene", "Mid-Holocene", "Late Holocene"))
)


metric_labels <- c(
  "richness_change" = "Richness change",
  "rank_change"     = "Rank change",
  "evenness_change" = "Evenness change",
  "gains"           = "Species gains",
  "losses"          = "Species losses"
)

# ==============================================================
# 9. Plot RAC metrics
# ==============================================================
plot_min <- min(rac$bin_numeric)
plot_max <- max(rac$bin_numeric)
rects_filtered <- rects %>% filter(ystart <= plot_max, yend >= plot_min)
time_breaks <- pretty(rac$bin_numeric, n = 10)

p <- ggplot(rac, aes(x = Value, y = bin_numeric)) +
  geom_rect(data = rects_filtered,
            aes(ymin = ystart, ymax = yend, xmin = -Inf, xmax = Inf, fill = col),
            inherit.aes = FALSE, alpha = 0.5) +
  geom_vline(xintercept = 0, color = "black", linetype = "solid", linewidth = 0.5, alpha = 0.3) +
  geom_path(aes(group = Metric), color = "black", linetype = "dashed") +
  geom_point(size = 3, color = "black") +
  facet_wrap(~Metric, nrow = 1, scales = "free_x", labeller = labeller(Metric = metric_labels)) +
  scale_fill_jco(guide = guide_legend(reverse = TRUE)) +
  scale_y_reverse(breaks = seq(16000, -500, by = -2000), labels = seq(16000, -500, by = -2000)) +
  labs(title = "Rank Abundance Curves",
       x = "Value", y = "Time (Years BP)", fill = "Time Periods") +
  coord_cartesian(ylim = c(plot_max, plot_min)) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(size = 12, colour = "black", angle = 45, vjust = .5),
    axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
    axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
    strip.background = element_rect(fill = "#f0f0f0", color = NA),
    strip.text.x = element_text(size = 14, face = "bold", colour = "black"),
    plot.title = element_text(size = 16, face = "bold", colour = "black"),
    legend.title = element_text(size = 16, face = "bold", colour = "black"),
    legend.text = element_text(size = 12, colour = "black"),
    legend.position = "right"
  )

# ==============================================================
# 10. Save figure
# ==============================================================
ggsave(filename = "005-rank-abundance-curves.jpg",
       plot = p, path = here("analysis", "figures"),
       width = 3300, height = 4200, units = "px", dpi = 300)
