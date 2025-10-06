pacman::p_load(codyn, data.table, ggsci, ggplot2, IRanges, SRS, tidyverse, vegan)

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

# Prepare data for RAC analysis
raw.abund <-
  hits %>%
  # Chronological comparison to previous age bin (i.e., 0-500 vs 500-1000, etc.) requires dates to be negative
  mutate(end.1 = if_else(end.1 != 0, end.1 * -1, end.1)) %>%
  select(end = end.1, taxon, abundance) %>%
  group_by(end, taxon) %>%
  summarise(tot_abund = sum(abundance))

#Interactive SRS shiny app for identifying suitable Cmin in SRS analysis
#if (interactive()) {SRS.shiny.app(data)}

# Prepare data for SRS scaling
srs.prep <-
  hits %>%
  select (site.x, sample_group, sample, end.1, taxon, abundance) %>%
  distinct() %>%
  mutate(sample = paste(sample, sample_group, site.x, end.1, sep = "@")) %>%
  pivot_wider(id_cols = taxon,
              names_from = sample,
              values_from = abundance,
              values_fill = 0) %>%
  column_to_rownames("taxon")

# Bind species row names (dropped during SRS), and the resulting SRS matrix
srs.data <-
  cbind(species = rownames(srs.prep), SRS(srs.prep, 10, set_seed = TRUE, seed = 1988))

# Transpose SRS scaled data and prepare species matrix for diversity analysis
srs.abund <-
  srs.data %>%
  pivot_longer(-species) %>%
  filter(value > 0) %>%
  distinct() %>%
  separate(col = name,
           sep = "@",
           into = c("sample", "sample_group", "site", "end")) %>%
  select(end, taxon = species, abundance = value) %>%
  group_by(end, taxon) %>%
  summarise(tot_abund = sum(abundance)) %>%
  mutate(end = as.numeric(end)) %>%
  mutate(end = if_else(end != 0, end * -1, end)) %>%
  as.data.frame()

# List with data split by region
Listdf <- list(raw.abund, srs.abund)

# Prepare empty lists
rac.list = list()

# Loop through data and generate diversity metrics
for (i in 1:length(Listdf)) {
  rac.list[[i]] = data.frame(
    RAC_change(df = Listdf[[i]],
               time.var = "end",
               species.var = "taxon",
               abundance.var = "tot_abund") %>%
      # Remove negative symbol before year and turn list id into a factor
      mutate(end = as.numeric(sub('.*-', '', end))) %>%
      pivot_longer(cols = c(richness_change, evenness_change, rank_change, gains, losses),
                   names_to = "metric",
                   values_to = "value")
  )
}

rac.raw <- rac.list[[1]]
rac.srs <- rac.list[[2]]

# Breaks for background rectangles, from Pilotto et al. (2022) ecological trait model
rects <-
  data.frame(xstart = c(12000,9500,8000,6500,4500,4000,3500,-500),
             xend = c(16000,12000,9500,8000,6500,4500,4000,3500),
             col = c("TM 1","TM 2","TM 3","TM 4","TM 5","TM 6","TM 7","TM 8"))

# Set background period factor levels for visualisation
rects$col <-
  factor(rects$col, levels = c("TM 1","TM 2","TM 3","TM 4","TM 5","TM 6","TM 7","TM 8"))
# Set metrics factor levels for visualisation
rac.raw$metric <-
  factor(rac.raw$metric, levels = c("richness_change", "rank_change", "evenness_change", "gains", "losses"))
# Set metrics factor levels for visualisation
rac.srs$metric <-
  factor(rac.srs$metric, levels = c("richness_change", "rank_change", "evenness_change", "gains", "losses"))

# Labeller for diversity metrics facets
metric.labs <- c(
  'evenness_change' = "Evenness change",
  'rank_change' = "Rank change",
  'richness_change' = "Richness change",
  'gains' = "Species gains",
  'losses' = "Species losses"
)

# Plot rank shifts
fig.raw <-
  ggplot() +
  geom_rect(data = rects,
            aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = col),
            alpha = 0.5) +
  geom_line(data = rac.raw, aes(x=end, y=value), linetype = "dashed", linewidth = 1) +
  geom_point(data = rac.raw, aes(y = value, x = end), shape = 16, size = 3, show.legend = F) +
  annotate(geom = "segment",x = -500, xend = 16000, y = 0, yend = 0) +
  scale_x_reverse(limits = c(16000, -500), breaks = scales::pretty_breaks(n = 10)) +
  scale_fill_jco(name = "Time Periods") +
  #coord_flip() +
  xlab("Age BP") +
  ylab("") +
  ggh4x::facet_grid2(metric ~ ., scales = "free", axes = "margins", independent = "x", labeller = labeller(metric = metric.labs)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black", angle = 45, vjust = .5),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
        strip.text.x = element_text(size = 14, face = "bold", colour = "black"),
        strip.text.y = element_text(size = 12, face = "bold", colour = "black"))

# Plot rank shifts
fig.srs <-
  ggplot() +
  geom_rect(data = rects,
            aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = col),
            alpha = 0.5) +
  geom_line(data = rac.srs, aes(x=end, y=value), linetype = "dashed", linewidth = 1) +
  geom_point(data = rac.srs, aes(y = value, x = end), shape = 16, size = 3, show.legend = F) +
  annotate(geom = "segment",x = -500, xend = 16000, y = 0, yend = 0) +
  scale_x_reverse(limits = c(16000, -500), breaks = scales::pretty_breaks(n = 10)) +
  scale_fill_jco(name = "Time Periods") +
  #coord_flip() +
  xlab("Age BP") +
  ylab("") +
  ggh4x::facet_grid2(metric ~ ., scales = "free", axes = "margins", independent = "x", labeller = labeller(metric = metric.labs)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black", angle = 45, vjust = .5),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
        strip.text.x = element_text(size = 14, face = "bold", colour = "black"),
        strip.text.y = element_text(size = 12, face = "bold", colour = "black"))

fig <-
  cowplot::plot_grid(fig.raw, fig.srs,
                   align = "h",
                   nrow = 1,
                   labels = c("A", "B"))

#Save figure
ggsave("004-rank-abundance-curves.jpg",
       fig,
       device = "jpg",
       here::here("analysis", "figures"),
       width=3300,
       height=4200,
       units = "px",
       dpi = 300)
