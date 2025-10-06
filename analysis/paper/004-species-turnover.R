pacman::p_load(codyn, ggsci, ggplot2, IRanges, tidyverse, vegan)

# Import site and species data
bugs.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "bugs_europe_extraction_samples_20250612.csv"), sep = ";", dec = ".", header = TRUE, na = c("", "NA", "NULL"), encoding = "UTF-8")

# Filter data and pivot insect species with abundance values, grouped by (site + date + sample group + sample)
time.mat <-
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
  findOverlaps(query = do.call(IRanges, time.mat), subject = do.call(IRanges, range), type = "any")

# Retrieve the hits from the intersection and name of samples
hits <-
  data.frame(time.mat[queryHits(intersection),], range[subjectHits(intersection),]) %>%
  # retrieve concatenated row name
  as_tibble(rownames = "sample") %>%
  # separate concatenation
  separate(col = sample,
           sep = "@",
           into = c("sample", "sample_group", "site")) %>%
  inner_join(bugs.csv, by = c("sample", "sample_group"), relationship = "many-to-many") %>%
  select(-start, -end)

# summarise taxon per age bin
raw.mat <-
  hits %>%
  select(site.x, sample_group, sample, taxon, end.1, abundance) %>%
  distinct() %>%
  group_by(taxon, end.1) %>%
  summarise(tot_abund = sum(abundance)) %>%
  # To ensure the turnover function makes chronological comparisons correctly values before 0 BP need to be changed to negative values
  mutate(end.1 = if_else(end.1 != 0, end.1 * -1, end.1))

# Calculate relative total turnover within replicates
total.res <- turnover(df = raw.mat,
                      time.var = "end.1",
                      species.var = "taxon",
                      abundance.var = "tot_abund",
                      metric = "total") %>%
  # Change age bin values back to positive for visualisaton purposes
  mutate(end.1 = if_else(end.1 != 0, end.1 * -1, end.1)) %>%
  dplyr::rename(value = total) %>%
  mutate(type = "Total turnover")

# Calculate relative total turnover within replicates
appear.res <- turnover(df = raw.mat,
                      time.var = "end.1",
                      species.var = "taxon",
                      abundance.var = "tot_abund",
                      metric = "appearance") %>%
  # Change age bin values back to positive for visualisaton purposes
  mutate(end.1 = if_else(end.1 != 0, end.1 * -1, end.1)) %>%
  dplyr::rename(value = appearance) %>%
  mutate(type = "Appearance")

# Calculate relative total turnover within replicates
disappear.res <- turnover(df = raw.mat,
                      time.var = "end.1",
                      species.var = "taxon",
                      abundance.var = "tot_abund",
                      metric = "disappearance") %>%
  # Change age bin values back to positive for visualisaton purposes
  mutate(end.1 = if_else(end.1 != 0, end.1 * -1, end.1)) %>%
  dplyr::rename(value = disappearance) %>%
  mutate(type = "Disappearance")

# Merged dataframes
spec.turnover <-
  bind_rows(total.res, appear.res, disappear.res)

# Breaks for background rectangles, from Pilotto et al, for natural context data
rects <-
  data.frame(xstart = c(12000,9500,8000,6500,4500,4000,3500,-500),
             xend = c(16000,12000,9500,8000,6500,4500,4000,3500),
             col = c("TM 1","TM 2","TM 3","TM 4","TM 5","TM 6","TM 7","TM 8"))

# Set factor levels for visualisation
rects$col <-
  factor(rects$col, levels = c("TM 1","TM 2","TM 3","TM 4","TM 5","TM 6","TM 7","TM 8"))

fig <-
  ggplot() +
  geom_rect(data = rects,
            aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = col),
            alpha = 0.5) +
  geom_line(data = spec.turnover, aes(x=end.1, y=value, linetype = type, color = type), linewidth = 1) +
  geom_point(data = spec.turnover, aes(x=end.1, y=value, color = type), size = 2) +
  scale_x_reverse(limits = c(16000, -500), breaks = scales::pretty_breaks(n = 10)) +
  scale_color_manual(name = "Metric",,
                     values = c("green4", "red4", "black")) +
  scale_linetype_manual(name = "Metric",
                        values=c("longdash", "dotted", "solid"))+
  scale_fill_jco(name = "Time Periods") +
  coord_flip() +
  xlab("Age (BP)") +
  ylab("Relative species turnover") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
        axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
        legend.title = element_text(size = 14, face = "bold", colour = "black"),
        legend.text = element_text(size = 12))

#Save figure
ggsave("004-species-turnover.jpg",
       fig,
       device = "jpg",
       here::here("analysis", "figures"),
       width=2400,
       height=3600,
       units = "px",
       dpi = 300)
