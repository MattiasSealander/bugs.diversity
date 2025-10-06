pacman::p_load(codyn, data.table, ggsci, ggplot2, IRanges, tidyverse, vegan)

# Import site and species data
bugs.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "bugs_europe_extraction_samples_20250612.csv"), sep = ";", dec = ".", header = TRUE, na = c("", "NA", "NULL"), encoding = "UTF-8")

# Filter data and pivot insect species with abundance values, grouped by (site + date + sample group + sample)
raw.mat <-
  bugs.csv %>%
  select (country, sample, site, sample_group, age_older, age_younger, context) %>%
  mutate(age_range=age_older - age_younger) %>%
  filter(context == "Stratigraphic sequence", sample != "BugsPresence", between(age_older, -500, 16000) & between(age_younger, -500, 16000), age_range <= 2000,
         country != "Greenland") %>%
  distinct() %>%
  select(-age_range, -context) %>%
  mutate(sample = paste(sample, sample_group, site, sep = "@")) %>%
  column_to_rownames("sample") %>%
  select(-country, -site, -sample_group) %>%
  dplyr::rename(start = age_younger, end = age_older)

# Select age bin(s)
range <- data.frame(
  start = c(-500, seq(1, 15501, by = 500)),
  end = seq(0, 16000, by = 500))

# Find what samples overlap in time with the specified date bin(s)
intersection <-
  findOverlaps(query = do.call(IRanges, raw.mat), subject = do.call(IRanges, range), type = "any")

# Retrieve the hits from the intersection and name of samples
hits <-
  data.frame(raw.mat[queryHits(intersection),], range[subjectHits(intersection),]) %>%
  # retrieve concatenated row name
  as_tibble(rownames = "sample") %>%
  # separate concatenation
  separate(col = sample,
           sep = "@",
           into = c("sample", "sample_group", "site")) %>%
  inner_join(bugs.csv, by = c("sample", "sample_group")) %>%
  select(-start, -end)

raw.abund <-
  hits %>%
  mutate(end.1 = if_else(end.1 != 0, end.1 * -1, end.1)) %>%
  select(end.1, taxon, abundance) %>%
  group_by(end.1, taxon) %>%
  summarise(tot_abund = sum(abundance))

# Perform RAC
rac <-
  RAC_change(df = raw.abund,
             time.var = "end.1",
             species.var = "taxon",
             abundance.var = "tot_abund") %>%
  # Remove negative symbol before year and turn list id into a factor
  mutate(end.1 = as.numeric(sub('.*-', '', end.1))) %>%
  pivot_longer(cols = c(richness_change, evenness_change, rank_change, gains, losses),
               names_to = "metric",
               values_to = "value")

# Breaks for background rectangles, from Pilotto et al. (2022) ecological trait model
rects <-
  data.frame(xstart = c(12000,9500,8000,6500,4500,4000,3500,-500),
             xend = c(16000,12000,9500,8000,6500,4500,4000,3500),
             col = c("TM 1","TM 2","TM 3","TM 4","TM 5","TM 6","TM 7","TM 8"))

# Set background period factor levels for visualisation
rects$col <-
  factor(rects$col, levels = c("TM 1","TM 2","TM 3","TM 4","TM 5","TM 6","TM 7","TM 8"))
# Set metrics factor levels for visualisation
rac$metric <-
  factor(rac$metric, levels = c("richness_change", "rank_change", "evenness_change", "gains", "losses"))

# Labeller for diversity metrics facets
metric.labs <- c(
  'evenness_change' = "Evenness change",
  'rank_change' = "Rank change",
  'richness_change' = "Richness change",
  'gains' = "Species gains",
  'losses' = "Species losses"
)

# Plot rank shifts
fig <-
  ggplot() +
  geom_rect(data = rects,
            aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = col),
            alpha = 0.5) +
  geom_line(data = rac, aes(x=end.1, y=value), linetype = "dashed", linewidth = 1) +
  geom_point(data = rac, aes(y = value, x = end.1), shape = 16, size = 3, show.legend = F) +
  geom_segment(data = rac, aes(y = 0, yend = 0, xend = -500, x = 16000)) +
  scale_x_reverse(limits = c(16000, -500), breaks = scales::pretty_breaks(n = 10)) +
  scale_fill_jco(name = "Time Periods") +
  coord_flip() +
  xlab("Age BP") +
  ylab("") +
  ggh4x::facet_grid2(~ metric, scales = "free", axes = "margins", independent = "x", labeller = labeller(metric = metric.labs)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black", angle = 45, vjust = .5),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
        strip.text.x = element_text(size = 14, face = "bold", colour = "black"),
        strip.text.y = element_text(size = 12, face = "bold", colour = "black"))

#Save figure
ggsave("004-rank-abundance-curves.jpg",
       fig,
       device = "jpg",
       here::here("analysis", "figures"),
       width=3300,
       height=4200,
       units = "px",
       dpi = 300)
