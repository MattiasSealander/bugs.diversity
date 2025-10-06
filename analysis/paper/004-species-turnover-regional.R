pacman::p_load(codyn, ggh4x, ggsci, ggplot2, IRanges, tidyverse, vegan)
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
  select(-start, -end) %>%
  mutate(place = case_when(
    between(latitude, 53.5, 74.8) & between(longitude, 2.8, 41) ~ "Scandinavia/Baltic",
    between(latitude, 49.8, 62.6) & between(longitude, -12.6, 1.8) ~ "British/Irish Isles",
    between(latitude, 63, 66.3) & between(longitude, -25, -12) ~ "Iceland",
    latitude < 45 ~ "Meridional")) %>%
  mutate(region = ifelse(is.na(place),
                         "Continental",
                         place))

# summarise taxon per region and age bin
raw.mat <-
  hits %>%
  select(region, site.x, sample_group, sample, taxon, end = end.1, abundance) %>%
  distinct() %>%
  group_by(region, taxon, end) %>%
  summarise(tot_abund = sum(abundance)) %>%
  # To ensure the turnover function makes chronological comparisons correctly values before 0 BP need to be changed to negative values
  mutate(end = if_else(end != 0, end * -1, end))

# List with data split by region
Listdf <- split(raw.mat[,-1], raw.mat$region)

# Prepare empty lists
app.list = list()
disapp.list = list()
total.list = list()

# Loop through data and generate diversity metrics
for (i in 1:length(Listdf)) {
  app.list[[i]] = data.frame(
    turnover(df = Listdf[[i]],
             time.var = "end",
             species.var = "taxon",
             abundance.var = "tot_abund",
             metric = "appearance") %>%
      # Change age bin values back to positive for visualisaton purposes
      mutate(end = if_else(end != 0, end * -1, end)) %>%
      dplyr::rename(value = appearance) %>%
      mutate(type = "Appearance")
  )

  disapp.list[[i]] = data.frame(
    turnover(df = Listdf[[i]],
             time.var = "end",
             species.var = "taxon",
             abundance.var = "tot_abund",
             metric = "disappearance") %>%
      # Change age bin values back to positive for visualisaton purposes
      mutate(end = if_else(end != 0, end * -1, end)) %>%
      dplyr::rename(value = disappearance) %>%
      mutate(type = "Disappearance")
  )
  total.list[[i]] = data.frame(
    turnover(df = Listdf[[i]],
             time.var = "end",
             species.var = "taxon",
             abundance.var = "tot_abund",
             metric = "total") %>%
      # Change age bin values back to positive for visualisaton purposes
      mutate(end = if_else(end != 0, end * -1, end)) %>%
      dplyr::rename(value = total) %>%
      mutate(type = "Total turnover")
  )
}

# Prepare table for facetted plot, "id" indicates which list object results come from, also rename age bin column to "end"
app.df <-
  lapply(app.list, as.data.frame) %>%
  bind_rows(.id = "id")
disapp.df <-
  lapply(disapp.list, as.data.frame) %>%
  bind_rows(.id = "id")
total.df <-
  lapply(total.list, as.data.frame) %>%
  bind_rows(.id = "id")

# Merged dataframes
spec.turnover <-
  bind_rows(app.df, disapp.df, total.df) %>%
  mutate(id = factor(id, levels = c("5", "3", "1", "2", "4"),
                     labels = c("Scandinavia/Baltic", "Iceland", "British/Irish Isles", "Continental", "Meridional")))

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
  geom_line(data = spec.turnover, aes(x=end, y=value, linetype = type, color = type), linewidth = 1) +
  geom_point(data = spec.turnover, aes(x=end, y=value, color = type), size = 2) +
  scale_x_reverse(limits = c(16000, -500), breaks = scales::pretty_breaks(n = 10)) +
  scale_color_manual(name = "Metric",,
                     values = c("green4", "red4", "black")) +
  scale_linetype_manual(name = "Metric",
                        values=c("longdash", "dotted", "solid"))+
  scale_fill_jco(name = "Time Periods") +
  facet_grid2(id ~ ., scales = "free", axes = "margins", labeller = labeller(type = spec.turnover$type)) +
  #coord_flip() +
  xlab("Age (BP)") +
  ylab("Relative species turnover") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
        axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
        strip.text.x = element_text(size = 16, face = "bold", colour = "black"),
        strip.text.y = element_text(size = 16, face = "bold", colour = "black"),
        legend.title = element_text(size = 14, face = "bold", colour = "black"),
        legend.text = element_text(size = 12))

#Save figure
ggsave("004-relative-total-turnover-divided-fixed.jpg",
       fig,
       device = "jpg",
       here::here("analysis", "figures"),
       width=40,
       height=30,
       units = "cm",
       dpi = 300)
