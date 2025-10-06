pacman::p_load(ggh4x, ggsci, ggplot2, tidyverse, vegan, procs, data.table, IRanges, SRS)

#Import site and species data
bugs.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "bugs_europe_extraction_samples_20250612.csv"), sep = ";", dec = ".", header = TRUE, na = c("", "NA", "NULL"), encoding = "UTF-8")

# Filter data to period of interest from natural deposits (context == "Stratigraphic sequence")
# Exclude samples with age range > 2000 years and samples that are temporary presence only records ("BugsPresence")
bugs <-
  bugs.csv %>%
  mutate(age_range=age_older - age_younger) %>%
  filter(context == "Stratigraphic sequence", sample != "BugsPresence", between(age_older, -500, 16000) & between(age_younger, -500, 16000), age_range <= 2000,
         country != "Greenland") %>%
  distinct()

# Filter data and pivot insect species with abundance values, grouped by (site + date + sample group + sample)
time.mat <-
  bugs %>%
  select (sample, site, sample_group, age_older, age_younger) %>%
  mutate(sample = paste(sample, sample_group, site, sep = "@")) %>%
  distinct() %>%
  column_to_rownames("sample") %>%
  select(-site, -sample_group) %>%
  dplyr::rename(start = age_younger, end = age_older)

# Select date bin(s)
range <- data.frame(
  start = c(-500, seq(1, 15501, by = 500)),
  end = seq(0, 16000, by = 500))

#find what samples overlap in time with the specified date bin(s)
intersection <-
  findOverlaps(query = do.call(IRanges, time.mat), subject = do.call(IRanges, range), type = "any")

#Retrieve the hits from the intersection and name of samples
hits <-
  data.frame(time.mat[queryHits(intersection),], range[subjectHits(intersection),]) %>%
  #retrieve concatenated row name
  as_tibble(rownames = "sample") %>%
  #separate concatenation
  separate(col = sample,
           sep = "@",
           into = c("sample", "sample_group", "site")) %>%
  inner_join(bugs.csv, by = c("sample", "sample_group"), relationship = "many-to-many") %>%
  select(-start, -end)

#filter data and pivot insect species with abundance values, grouped by (site + date + sample group + sample)
raw.mat <-
  hits %>%
  select (site = site.x, sample_group, sample, end = end.1, taxon, abundance) %>%
  distinct() %>%
  group_by(site, sample_group, sample, end) %>%
  pivot_wider(id_cols = c(site, sample_group, sample, end),
              names_from = taxon,
              values_from = abundance,
              values_fill = 0) %>%
  mutate(sample = paste(site, sample_group, sample, sep = "@")) %>%
  column_to_rownames("sample") %>%
  select(-site, -sample_group)

# Calculate species richness per sample and then summarise mean and standard error per age bin
raw.rich <-
  cbind(end = raw.mat[,1], richness = specnumber(raw.mat[,2:2326])) %>%
  as.data.frame() %>%
  group_by(end) %>%
  summarise(err = sd(richness, na.rm = TRUE)/sqrt(length(richness)),
            mean_rich = round(mean(richness, na.rm = TRUE), 2)) %>%
  mutate(end = as.numeric(end))

# Filter lowest sample abundance to rarefaction threshold (10)
rare.mat <-
  raw.mat %>%
  filter(rowSums(.[,2:2326]) >= 10)

# Calculate rarefied species richness per sample and then summarise mean and standard error per age bin
rare.rich <-
  cbind(end = rare.mat[,1], richness = rarefy(rare.mat[,2:2326], sample = 10)) %>%
  as.data.frame() %>%
  group_by(end) %>%
  summarise(mean_rich = round(mean(richness), 2),
            err = sd(richness)/sqrt(length(richness))) %>%
  mutate(end = as.numeric(end))

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
srs.mat <-
  srs.data %>%
  pivot_longer(-species) %>%
  pivot_wider(names_from = species,
              values_from = value) %>%
  separate(col = name,
           sep = "@",
           into = c("sample", "sample_group", "site", "end")) %>%
  as.data.frame() %>%
  distinct() %>%
  mutate(sample = paste(site, sample_group, sample, sep = "@")) %>%
  column_to_rownames("sample") %>%
  select(-site, -sample_group)

# Calculate richness on SRS scaled values
srs.rich <-
  cbind(end = srs.mat[,1],
        srs.rich = specnumber(srs.mat[,2:2326])) %>%
  as.data.frame() %>%
  mutate(srs.rich = as.numeric(srs.rich),
         end = as.numeric(end)) %>%
  group_by(end) %>%
  # Calculate mean and standard error of richness by age bin
  summarize(mean_rich = round(mean(srs.rich, na.rm = T), 2),
            err = sd(srs.rich)/sqrt(length(srs.rich)))

# Bind to one data frame
rich.mat <-
  rbind((raw.rich %>%
           mutate(type = "Raw")),
        (rare.rich %>%
           mutate(type = "Rarefaction")),
        (srs.rich %>%
           mutate(type = "SRS")))

#Breaks for background rectangles, from Pilotto et al, for natural context data
rects.n <-
  data.frame(xstart = c(12000,9500,8000,6500,4500,4000,3500,-500),
             xend = c(16000,12000,9500,8000,6500,4500,4000,3500),
             col = c("TM 1","TM 2","TM 3","TM 4","TM 5","TM 6","TM 7","TM 8"))
# Set factor levels
rects.n$col <-
  factor(rects.n$col, levels = c("TM 1","TM 2","TM 3","TM 4","TM 5","TM 6","TM 7","TM 8"))
rich.mat$type <-
  factor(rich.mat$type, levels = c("Raw","Rarefaction","SRS"))

# Plot mean richness for each type by 500 year bins with error as ribbon
fig <-
  ggplot() +
  geom_rect(data = rects.n,
            aes(ymin = xstart, ymax = xend, xmin = -Inf, xmax = Inf, fill = col),
            alpha = 0.5) +
  geom_ribbon(data = rich.mat, aes(y = end, xmin = mean_rich - err, xmax = mean_rich + err), alpha = 0.5) +
  #geom_errorbar(data = rare.n, aes(y = end, xmin = mean_rich - err, xmax = mean_rich + err), lwd = 0.5) +
  geom_path(data = rich.mat, aes(y = end, x = mean_rich), linewidth = 1, linetype = "dashed") +
  geom_point(data = rich.mat, aes(y = end, x = mean_rich), shape = 21, size = 3, fill = "green4", col = "black") +
  scale_y_reverse(limits = c(16000, -500), breaks = scales::pretty_breaks(n = 10)) +
  scale_fill_jco(name = "Time Periods") +
  facet_grid2(~ type, scales = "free", axes = "margins", independent = "x", labeller = labeller(type = rich.mat$type)) +
  labs(y="Age (BP)",
       x="Mean richness",
       color=NULL) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 14, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black", hjust = 0.5),
        strip.text.x = element_text(size = 14, face = "bold", colour = "black"),
        strip.text.y = element_text(size = 12, face = "bold", colour = "black"))

# Save figure
ggsave("002-richness-comparison.jpg",
       fig,
       device = "jpg",
       here::here("analysis", "figures"),
       width=1800,
       height=2600,
       units = "px",
       dpi = 300)
