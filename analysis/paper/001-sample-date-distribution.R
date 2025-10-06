pacman::p_load(ggplot2, tidyverse, IRanges)

#Import site and species data
bugs.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "bugs_europe_extraction_samples_20250612.csv"), sep = ";", dec = ".", header = TRUE, na = c("", "NA", "NULL"), encoding = "UTF-8")

nat <-
  bugs.csv %>%
  mutate(age_range=age_older - age_younger) %>%
  filter(context == "Stratigraphic sequence", sample != "BugsPresence", between(age_older, -500, 16000) & between(age_younger, -500, 16000),
         age_range <= 2000, !country == "Greenland") %>%
  select (sample, site, sample_group, age_older, age_younger, context) %>%
  distinct() %>%
  mutate(sample = paste(sample, sample_group, site, sep = "@")) %>%
  column_to_rownames("sample") %>%
  select(-site, -sample_group) %>%
  dplyr::rename(start = age_younger, end = age_older)

#prepare age bin(s)
range <- data.frame(
  start = c(-500, seq(1, 15501, by = 500)),
  end = seq(0, 16000, by = 500))

#find what samples overlap in time with the specified age bin(s)
intersect <-
  findOverlaps(query = do.call(IRanges, nat), subject = do.call(IRanges, range), type = "any")

#Retrieve the hits from the intersection and summarise the nr of samples per age bin
hits.samples <-
  data.frame(nat[queryHits(intersect),], range[subjectHits(intersect),]) %>%
  #retrieve concatenated row name
  as_tibble(rownames = "sample") %>%
  group_by(end.1) %>%
  summarise(samples = n_distinct(sample))

#Retrieve the hits from the intersection and name of samples
hits.sites <-
  data.frame(nat[queryHits(intersect),], range[subjectHits(intersect),]) %>%
  #retrieve concatenated row name
  as_tibble(rownames = "sample") %>%
  #separate concatenation
  separate(col = sample,
           sep = "@",
           into = c("sample", "sample_group", "site")) %>%
  group_by(end.1) %>%
  summarise(sites = n_distinct(site))

#Retrieve the hits from the intersection and summarise the nr of samples per age bin
hits.abund <-
  data.frame(nat[queryHits(intersect),], range[subjectHits(intersect),]) %>%
  #retrieve concatenated row name
  as_tibble(rownames = "sample") %>%
  # separate concatenation
  separate(col = sample,
           sep = "@",
           into = c("sample", "sample_group", "site")) %>%
  inner_join(bugs.csv, by = c("sample", "sample_group")) %>%
  select(-start, -end) %>%
  select(end.1, sample, sample_group, site.x, abundance) %>%
  group_by(end.1) %>%
  summarise(abundance = sum(abundance))

#plot abundance
fig.sites <-
  ggplot() +
  geom_bar(data = hits.sites, aes(y=sites, x=end.1),
           position="stack", stat="identity", just = 0,
           fill = "green4",
           color = "black", alpha = 0.8) +
  labs(y="Sites",
       x="",
       color=NULL) +
  #coord_flip() +
  scale_x_reverse(limits = c(16000, -500), breaks = scales::pretty_breaks(n = 10)) +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 10, face = "bold", colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.x = element_text(size = 12, face = "bold", colour = "black"))

#plot abundance
fig.samples <-
  ggplot() +
  geom_bar(data = hits.samples, aes(y=samples, x=end.1),
           position="stack", stat="identity", just = 0,
           fill = "green4",
           color = "black", alpha = 0.8) +
  labs(y="Samples",
       x="Age bin (BP)",
       color=NULL) +
  #coord_flip() +
  scale_x_reverse(limits = c(16000, -500), breaks = scales::pretty_breaks(n = 10)) +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 10, face = "bold", colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.x = element_text(size = 12, face = "bold", colour = "black"))

#plot abundance
fig.abundance <-
  ggplot() +
  geom_bar(data = hits.abund, aes(y=abundance, x=end.1),
           position="stack", stat="identity", just = 0,
           fill = "green4",
           color = "black", alpha = 0.8) +
  labs(y="Abundance",
       x="Age bin (BP)",
       color=NULL) +
  #coord_flip() +
  #scale_y_continuous(trans='log10') +
  scale_x_reverse(limits = c(16000, -500), breaks = scales::pretty_breaks(n = 10)) +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 10, face = "bold", colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.x = element_text(size = 12, face = "bold", colour = "black"))

fig <-
  cowplot::plot_grid(fig.sites, fig.samples, fig.abundance,
                     align = "v",
                     labels = c('A', 'B'),
                     ncol = 1,
                     label_size = 12)

#Save figure
ggsave("001-sample-date-distribution2.jpg",
       fig,
       device = "jpg",
       here::here("analysis", "figures"),
       width=2400,
       height=1800,
       units = "px",
       dpi = 300)
