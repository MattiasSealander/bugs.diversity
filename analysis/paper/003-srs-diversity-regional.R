pacman::p_load(ggplot2, ggsci, tidyverse, cowplot, IRanges, SRS, vegan)

# Import site and species data
bugs.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "bugs_europe_extraction_samples_20250612.csv"), sep = ";", dec = ".", header = TRUE, na = c("", "NA", "NULL"), encoding = "UTF-8")

# Prepare data for age bin intersection analysis
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

# Find what samples overlap in time with the specified age bin(s)
intersection <-
  findOverlaps(query = do.call(IRanges, time.mat), subject = do.call(IRanges, range), type = "any")

# Retrieve the hits from the intersection and name of samples
hits <-
  data.frame(time.mat[queryHits(intersection),], range[subjectHits(intersection),]) %>%
  # Retrieve concatenated row name
  as_tibble(rownames = "sample") %>%
  # Separate concatenated sample name
  separate(col = sample,
           sep = "@",
           into = c("sample", "sample_group", "site")) %>%
  # Join age bins to fossil insect data
  inner_join(bugs.csv, by = c("sample", "sample_group")) %>%
  select(-start, -end)

# Filter data to period of interest from natural deposits (context == "Stratigraphic sequence")
# Exclude samples with age range > 2000 years and samples that are temporary presence only records ("BugsPresence")
raw.mat <-
  hits %>%
  select (site.x, latitude, longitude, sample_group, sample, end.1, taxon, abundance) %>%
  distinct() %>%
  group_by(site.x, end.1, sample_group, sample, latitude, longitude) %>%
  mutate(place = case_when(
    between(latitude, 53.5, 74.8) & between(longitude, 2.8, 41) ~ "Scandinavia/Baltic",
    between(latitude, 49.8, 62.6) & between(longitude, -12.6, 1.8) ~ "British/Irish Isles",
    between(latitude, 63, 66.3) & between(longitude, -25, -12) ~ "Iceland",
    latitude < 45 ~ "Meridional")) %>%
  mutate(region = ifelse(is.na(place),
                         "Continental",
                         place)) %>%
  select(-latitude, -longitude, -place) %>%
  distinct() %>%
  mutate(sample = paste(sample, sample_group, site.x, sep = "@"))

# Prepare empty list for diversity metrics
srs.list = list()
div.list = list()

# Add species matrix to list to use in for loop
Listdf <- list(
  (raw.mat %>%
     filter(region == "Scandinavia/Baltic") %>%
     select(-region) %>%
     pivot_wider(id_cols = taxon,
                 names_from = sample,
                 values_from = abundance,
                 values_fill = 0) %>%
     column_to_rownames("taxon")),
  (raw.mat %>%
     filter(region == "British/Irish Isles") %>%
     select(-region) %>%
     pivot_wider(id_cols = taxon,
                 names_from = sample,
                 values_from = abundance,
                 values_fill = 0) %>%
     column_to_rownames("taxon")),
  (raw.mat %>%
     filter(region == "Iceland") %>%
     select(-region) %>%
     pivot_wider(id_cols = taxon,
                 names_from = sample,
                 values_from = abundance,
                 values_fill = 0) %>%
     column_to_rownames("taxon")),
  (raw.mat %>%
     filter(region == "Continental") %>%
     select(-region) %>%
     pivot_wider(id_cols = taxon,
                 names_from = sample,
                 values_from = abundance,
                 values_fill = 0) %>%
     column_to_rownames("taxon")),
  (raw.mat %>%
      filter(region == "Meridional") %>%
      select(-region) %>%
      pivot_wider(id_cols = taxon,
                  names_from = sample,
                  values_from = abundance,
                  values_fill = 0) %>%
      column_to_rownames("taxon"))
)

# Loop through data and generate diversity metrics
for (i in 1:length(Listdf)) {
    srs.list[[i]] = data.frame(cbind(
      species = rownames(Listdf[[i]]), SRS(Listdf[[i]], 10, set_seed = TRUE, seed = 1988)),
      check.names = FALSE
      ) %>%
      pivot_longer(-species) %>%
      pivot_wider(names_from = species,
                  values_from = value) %>%
      separate(col = name,
               sep = "@",
               into = c("sample", "sample_group", "site")) %>%
      as.data.frame() %>%
      # Join with date intersection data to get age bins per sample
      inner_join(hits %>%
                   select(site.x, sample_group, sample, end.1), by = join_by("site" == "site.x", "sample_group", "sample")) %>%
      distinct() %>%
      mutate(sample = paste(site, sample_group, sample, sep = "@")) %>%
      column_to_rownames("sample") %>%
      select(-site, -sample_group) %>%
      select(end.1, everything())

    div.list[[i]] = data.frame(cbind(
        end = srs.list[[i]][,1],
        "raw.shan" = diversity(srs.list[[i]][,2:ncol(srs.list[[i]])], index = "shannon"), # Shannon
        "raw.simps" = diversity(srs.list[[i]][,2:ncol(srs.list[[i]])], index = "invsimpson"), # Simpson
        "raw.rich" = specnumber(srs.list[[i]][,2:ncol(srs.list[[i]])])),
        check.names = FALSE) %>%
  rownames_to_column("sample") %>%
  mutate(raw.even = raw.shan/log(raw.rich)) %>%
  group_by(end) %>%
  summarize(shan_mean = round(mean(raw.shan, na.rm = T), 2),
            #mean_even = round(mean(raw.even, na.rm = T), 2),
            rich_mean = round(mean(raw.rich, na.rm = T), 2),
            simps_mean = round(mean(raw.simps, na.rm = T), 2),
            shan_err = sd(raw.shan)/sqrt(length(raw.shan)),
            rich_err = sd(raw.rich)/sqrt(length(raw.rich)),
            simps_err = sd(raw.simps)/sqrt(length(raw.simps)))
}

# Breaks for background rectangles, from Pilotto et al. (2022) paper, for natural context data
# Pilotto et al. link: https://pmc.ncbi.nlm.nih.gov/articles/PMC9233931/#s2
rects <-
  data.frame(xstart = c(12000,9500,8000,6500,4500,4000,3500,-500),
             xend = c(16000,12000,9500,8000,6500,4500,4000,3500),
             col = c("TM 1","TM 2","TM 3","TM 4","TM 5","TM 6","TM 7","TM 8"))
# Set factor levels for visualisation
rects$TM <-
  factor(rects$col, levels = c("TM 1","TM 2","TM 3","TM 4","TM 5","TM 6","TM 7","TM 8"))

# Prepare table for facetted plot, "id" indicates which list object results come from, also rename age bin column to "end"
dresults <- lapply(div.list, as.data.frame) %>%
  bind_rows(.id = "id") %>%
  group_by(id) %>%
  select(-rich_err, -shan_err, -simps_err) %>%
  pivot_longer(cols = c(shan_mean, rich_mean, simps_mean),
               names_to = c("type", ".value"),
               names_pattern = "(.*)_(.*)",
               values_to = c("value")) %>%
  ungroup() %>%
  inner_join(
    (lapply(div.list, as.data.frame) %>%
       bind_rows(.id = "id") %>%
       group_by(id) %>%
       pivot_longer(cols = c(shan_err, rich_err, simps_err),
                    names_to = c("type", ".value"),
                    names_pattern = "(.*)_(.*)",
                    values_to = c("value")) %>%
       select(id, end, type, err))
    , by = c("id", "end", "type"))

# Labeller for region facets
region.labs <- c(
  '1' = "Scandinavia/Baltic",
  '2' = "British/Irish Isles",
  '3' = "Iceland",
  '4' = "Continental",
  '5' = "Meridional")

# Labeller for diversity metrics facets
metric.labs <- c(
  'rich' = "Richness",
  'shan' = "Shannon",
  'simps' = "Simpson")

# Plot the figure
fig <-
  ggplot() +
  geom_rect(data = rects,
            aes(ymin = xstart, ymax = xend, xmin = -Inf, xmax = Inf, fill = col),
            alpha = 0.5) +
  geom_ribbon(data = dresults, aes(y = end, xmin = mean - err, xmax = mean + err), alpha = 0.5, show.legend = F) +
  #geom_errorbar(data = dresults, aes(y = end, xmin = mean - err, xmax = mean + err), lwd = 0.5, show.legend = F) +
  geom_point(data = dresults, aes(y = end, x = mean), shape = 21, size = 3, fill = "green4", col = "black", show.legend = F) +
  scale_y_reverse(limits = c(16500, -500), breaks = scales::pretty_breaks(n = 10)) +
  ggsci::scale_fill_jco(name = "Time Periods") +
  labs(y="Age (BP)",
       x="",
       color=NULL) +
  #facet_grid(id ~ type, scales = "free_x", labeller = labeller(id = region.labs, type = metric.labs)) +
  ggh4x::facet_grid2(id ~ type, scales = "free", axes = "margins", independent = "x", labeller = labeller(id = region.labs, type = metric.labs)) +
  #facet_wrap(, ncol = 3, labeller = labeller(id = region.labs, type = metric.labs)) +
  theme_bw() +
  theme(legend.position = "right",
        legend.title = element_text(size = 16, face = "bold", colour = "black"),
        legend.text = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
        #axis.title.x = element_text(size = 12, face = "bold", colour = "black")
        plot.title = element_text(size = 12, face = "bold", colour = "black", hjust = -0.1),
        strip.text.x = element_text(size = 14, face = "bold", colour = "black"),
        strip.text.y = element_text(size = 12, face = "bold", colour = "black"))

#Save figure
ggsave("003-srs-diversity-natural-divided-inv-fixed.jpeg",
       fig,
       device = "jpeg",
       here::here("analysis", "figures"),
       scale = 1,
       width=30,
       height=50,
       units = "cm",
       dpi = 300)
