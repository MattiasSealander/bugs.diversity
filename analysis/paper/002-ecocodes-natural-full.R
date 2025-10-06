pacman::p_load(ggh4x, ggsci, ggplot2, tidyverse, procs, data.table, IRanges)

#Import site and species data
bugs.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "bugs_europe_extraction_samples_20250612.csv"), sep = ";", dec = ".", header = TRUE, na = c("", "NA", "NULL"), encoding = "UTF-8")

#Import ecocodes
eco.csv <-
  read.csv2(here::here("analysis", "data", "raw_data", "sead_ecocodes_20250114.csv"), sep = ",", dec = ".", header = TRUE, na = c("", "NA", "NULL"), encoding = "UTF-8")

#filter data and pivot insect species with abundance values, grouped by (site + date + sample group + sample)
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

#select date bin(s)
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

#join taxa with ecocodes
eco.species <-
  eco.csv %>%
  select(-author_name, -family_name, -genus_name, -species) %>%
  inner_join(hits, by = "taxon", relationship = "many-to-many")

#filter data to sites from "natural" contexts and pivot to wide format with counts per species
abund.mat <-
  hits %>%
  select(site.x, sample_group, sample, end.1) %>%
  distinct() %>%
  inner_join(eco.species, by = c("site.x", "sample_group", "sample")) %>%
  filter(ecocode_system == "Bugs") %>%
  #left_join(sample.filt, by = c("sample_group", "sample")) %>%
  select (end = end.1.x, abundance, ecocode, sample_group, sample) %>%
  distinct() %>%
  group_by(end) %>%
  pivot_wider(id_cols = end,
              names_from = ecocode,
              values_from = c(abundance),
              values_fn = function(x) sum(x),
              values_fill = 0) %>%
  ungroup() %>%
  #standardise the counts according to Buckland (2007) https://urn.kb.se/resolve?urn=urn:nbn:se:umu:diva-1105
  #this is "Abundance weighted; %SumRep"
  #code based on https://stackoverflow.com/questions/66710455/how-to-divide-a-columns-values-by-the-sum-of-multiple-column-values-by-row
  mutate(test = select(., 2:23)/select(., 2:23) %>%
           rowSums(na.rm = TRUE) * 100) %>%
  select(-2:-23) %>%
  unnest(cols = c(test)) %>%
  pivot_longer(2:23) %>%
  mutate(type = "Abundance Weighted")


#join filtered data with ecocodes, using the Bugs system, and leave out narrower indicator codes
#summarise and standardise environmental reconstruction
noabu.mat <-
  hits %>%
  select(site.x, sample_group, sample, end.1) %>%
  distinct() %>%
  inner_join(eco.species, by = c("site.x", "sample_group", "sample")) %>%
  filter(ecocode_system == "Bugs") %>%
  select (end = end.1.x, abundance, taxon, ecocode) %>%
  distinct() %>%
  group_by(end) %>%
  pivot_wider(id_cols = end,
              names_from = ecocode,
              values_from = taxon,
              values_fn = list(taxon = length)) %>%
  ungroup() %>%
  #standardise the counts according to Buckland (2007) https://urn.kb.se/resolve?urn=urn:nbn:se:umu:diva-1105
  #this is "No abundance; %SumRep"
  #below code based on https://stackoverflow.com/questions/66710455/how-to-divide-a-columns-values-by-the-sum-of-multiple-column-values-by-row
  mutate(test = select(., 2:23)/select(., 2:23) %>%
           rowSums(na.rm = TRUE) * 100) %>%
  select(-2:-23) %>%
  unnest(cols = c(test)) %>%
  pivot_longer(2:23) %>%
  mutate(type = "No Abundance")

# Bind to one data frame
ecocodes <-
  rbind((abund.mat %>%
           mutate(type = "Abundance Weighted")),
        (noabu.mat %>%
           mutate(type = "No Abundance")))

#set factor levels for ecocodes
ecocodes$name <-
  factor(ecocodes$name, levels = c("Aquatics", "Indicators: Running water", "Indicators: Standing water", "Open wet habitats", "Wetlands/marshes", "Mould beetles", "Halotolerant",
                                   "Carrion", "Ectoparasite", "General synanthropic", "Stored grain pest", "Dung/foul habitats","Pasture/Dung", "Indicators: Dung", "Disturbed/arable",
                                   "Sandy/dry disturbed/arable", "Heathland & moorland", "Wood and trees", "Indicators: Deciduous", "Indicators: Coniferous", "Dry dead wood",
                                   "Meadowland"))

#set factor levels for ecocodes
noabu.mat$name <-
  factor(noabu.mat$name, levels = c("Aquatics", "Indicators: Running water", "Indicators: Standing water", "Open wet habitats", "Wetlands/marshes", "Mould beetles", "Halotolerant",
                                    "Carrion", "Ectoparasite", "General synanthropic", "Stored grain pest", "Dung/foul habitats","Pasture/Dung", "Indicators: Dung", "Disturbed/arable",
                                    "Sandy/dry disturbed/arable", "Heathland & moorland", "Wood and trees", "Indicators: Deciduous", "Indicators: Coniferous", "Dry dead wood",
                                    "Meadowland"))

habitat <-
  c("Aquatics", "Indicators: Running water", "Indicators: Standing water", "Open wet habitats", "Wetlands/marshes", "Mould beetles", "Halotolerant",
    "Carrion", "Ectoparasite", "General synanthropic", "Stored grain pest", "Dung/foul habitats","Pasture/Dung", "Indicators: Dung", "Disturbed/arable",
    "Sandy/dry disturbed/arable", "Heathland & moorland", "Wood and trees", "Indicators: Deciduous", "Indicators: Coniferous", "Dry dead wood",
    "Meadowland")


colors <-
  c("#197EC0FF", "lightskyblue1", "#71D0F5FF", "#709AE1FF", "dodgerblue4", "#370335FF", "#1A9993FF",
    "#8A9197FF", "darkorchid", "#FD8CC1FF", "#D2AF81FF", "#91331FFF", "#FED439FF", "darkred", "salmon1",
    "tomato3", "#D5E4A2FF", "darkslategrey", "#46732EFF", "green4", "brown4", "yellowgreen")

# Stacked No abundance, sumrep
fig <-
  ecocodes %>%
  ggplot(aes(fill=name, y=value, x=end)) +
  geom_bar(position="stack", stat="identity", just = 0, alpha = 0.9) +
  xlab("Age (BP)") +
  ylab("%SumRep") +
  ggtitle("") +
  scale_fill_manual(values = colors,
                    name = "Habitat") +
  guides(fill = guide_legend(nrow = 5,
                             reverse = T)) +
  scale_x_reverse(limits = c(16000, -500), breaks = scales::pretty_breaks(n = 10)) +
  coord_flip() +
  facet_grid2(~ type, scales = "free", axes = "margins", labeller = labeller(type = ecocodes$type)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 16, face = "bold", colour = "black"),
        #axis.text.x = element_text(angle = 90),
        axis.title.y = element_text(size = 16, face = "bold", colour = "black"),
        legend.title = element_text(size = 16, face = "bold", colour = "black"),
        strip.text.x = element_text(size = 16, face = "bold", colour = "black"),
        strip.text.y = element_text(size = 16, face = "bold", colour = "black"),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.text = element_text(size = 12),
        legend.position = "bottom")

#Save figure
ggsave("002-ecocodes-natural-fixed.jpg",
       fig,
       device = "jpg",
       here::here("analysis", "figures"),
       width=40,
       height=50,
       units = "cm",
       dpi = 300)
