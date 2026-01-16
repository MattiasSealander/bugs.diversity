
# ==============================================================
# Script: Gamma Diversity (SRS only) + δ18O (GISP2) — Side-by-side
#
# Purpose:
#   Compute gamma diversity (Shannon, Simpson) using SRS-standardized
#   abundances for Holocene fossil insect data (British/Irish Isles),
#   and compare visually with GISP2 δ18O on aligned time axes.
#
# Output:
#   - Combined plot: "005-gamma-srs-vs-d18o-coleoptera.jpg"
# ==============================================================

# ---- 0. Load packages ----
pacman::p_load(
  data.table, ggsci, tidyverse, here
)

# ==============================================================
# 1. Read data
# ==============================================================

# Read δ18O (GISP2) and harmonize time span
d18o <- fread(here("analysis", "data", "raw_data", "d18o.csv"), na.strings = "NaN") %>%
  as_tibble() %>%
  filter(age_BP_SEAD <= 16000)

# ==============================================================
# 2. Plot δ18O temperature
# ==============================================================

# 2.1 Setup time axis for plotting
time_min <- 0
time_max <- 16000
time_breaks <- seq(time_min, time_max, by = 1000)

# 2.2 Remove NA records and sort by age
d18o_line <- d18o %>%
  dplyr::filter(!is.na(d18O_GISP2_per_mille), !is.na(age_BP_SEAD)) %>%
  dplyr::arrange(age_BP_SEAD)

# 2.3 Generate δ18O plot
p_d18 <- ggplot(d18o_line, aes(x = age_BP_SEAD, y = d18O_GISP2_per_mille)) +
  geom_line(linewidth = .5, color = "black", alpha = 1, aes(group = 1)) +
  geom_vline(xintercept = seq(time_min, time_max, by = 500),  colour = "grey", linetype = "dashed") +
  geom_vline(xintercept = seq(time_min, time_max, by = 1000), colour = "grey") +
  scale_x_reverse(breaks = time_breaks, labels = time_breaks) +
  coord_cartesian(xlim = c(time_max, time_min)) +
  coord_flip() +
  labs(title = "GISP2 δ\u00b9\u2078O (‰)", y = "δ\u00b9\u2078O (‰)", x = "") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y   = element_blank(),
    axis.text.x   = element_text(size = 12, colour = "black", angle = 45, vjust = 0.5),
    axis.title.y  = element_text(size = 12, face = "bold", colour = "black"),
    axis.title.x  = element_text(size = 12, face = "bold", colour = "black"),
    plot.title    = element_text(size = 16, face = "bold", colour = "black"),
    legend.position = "none"
  )

# ==============================================================
# 3. Save plot
# ==============================================================

#Save as RDS
saveRDS(p_d18, file = here::here("analysis", "data", "derived_data", "d18O_plot.rds"))
