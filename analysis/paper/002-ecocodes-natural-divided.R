pacman::p_load(ggh4x, ggsci, ggplot2, tidyverse, vegan, procs, data.table, IRanges)

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
  inner_join(hits, by = "taxon", relationship = "many-to-many") %>%
  filter(ecocode_system == "Bugs") %>%
  select(site.x, latitude, longitude, sample_group, sample, end.1, ecocode, abundance, taxon) %>%
  mutate(sample = paste(sample, sample_group, site.x, sep = "@")) %>%
  distinct() %>%
  mutate(place = case_when(
    between(latitude, 53.5, 74.8) & between(longitude, 2.8, 41) ~ "Scandinavia/Baltic",
    between(latitude, 49.8, 62.6) & between(longitude, -12.6, 1.8) ~ "British/Irish Isles",
    between(latitude, 63, 66.3) & between(longitude, -25, -12) ~ "Iceland",
    latitude < 45 ~ "Meridional")) %>%
  mutate(region = ifelse(is.na(place),
                         "Continental",
                         place)) %>%
  select(-site.x, -sample_group, -latitude, -longitude)


# Species matrix of raw values
# Filter data and pivot insect species with abundance values, grouped by (site + date + sample group + sample)
abund.mat <-
  eco.species %>%
  select (end = end.1, abundance, ecocode) %>%
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

# List with data split by region
Listdf <- split(eco.species[,-7], eco.species$region)

# Prepare empty list
abund.list = list()
noabund.list = list()

# Loop through data and generate diversity metrics
for (i in 1:length(Listdf)) {
  abund.list[[i]] = data.frame(
    Listdf[[i]] %>%
      select (end = end.1, abundance, ecocode) %>%
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
      mutate(test = select(., 2:ncol(.))/select(., 2:ncol(.)) %>%
               rowSums(na.rm = TRUE) * 100) %>%
      select(-2:-(ncol(.)-1)) %>%
      unnest(cols = c(test)) %>%
      pivot_longer(2:ncol(.)) %>%
      mutate(type = "Abundance Weighted")
  )

  noabund.list[[i]] = data.frame(
    Listdf[[i]] %>%
      select (end = end.1, abundance, taxon, ecocode) %>%
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
      mutate(test = select(., 2:ncol(.))/select(., 2:ncol(.)) %>%
               rowSums(na.rm = TRUE) * 100) %>%
      select(-2:-(ncol(.)-1)) %>%
      unnest(cols = c(test)) %>%
      pivot_longer(2:ncol(.)) %>%
      mutate(type = "No Abundance")
  )
}

# Prepare table for facetted plot, "id" indicates which list object results come from, also rename age bin column to "end"
eco.abund <-
  lapply(abund.list, as.data.frame) %>%
  bind_rows(.id = "id")

eco.noabund <-
  lapply(noabund.list, as.data.frame) %>%
  bind_rows(.id = "id")

# Count of ecocodes per period and region
count.abu <-
  eco.abund %>%
  filter(value > 0) %>%
  select(id, end, name) %>%
  group_by(id, end) %>%
  summarise(count = n_distinct(name)) %>%
  ungroup() %>%
  # Remove age bins with less than 5 ecocodes present
  filter(count > 5) %>%
  select(id, end)

# Join to main dataframes to remove age bins with less than 5 ecocodes present
eco.abund <-
  eco.abund %>%
  inner_join(count.abu, by = c("id", "end"))

# Count of ecocodes per period and region
count.noabund <-
  eco.noabund %>%
  filter(value > 0) %>%
  select(id, end, name) %>%
  group_by(id, end) %>%
  summarise(count = n_distinct(name)) %>%
  ungroup() %>%
  # Remove age bins with less than 5 ecocodes present
  filter(count > 5) %>%
  select(id, end)

# Join to main dataframes to remove age bins with less than 5 ecocodes present
eco.noabund <-
  eco.noabund %>%
  inner_join(count.noabund, by = c("id", "end"))

#set factor levels for ecocodes
eco.abund$name <-
  factor(eco.abund$name, levels = c("Aquatics", "Indicators: Running water", "Indicators: Standing water", "Open wet habitats", "Wetlands/marshes", "Mould beetles", "Halotolerant",
                                    "Carrion", "Ectoparasite", "General synanthropic", "Stored grain pest", "Dung/foul habitats","Pasture/Dung", "Indicators: Dung", "Disturbed/arable",
                                    "Sandy/dry disturbed/arable", "Heathland & moorland", "Wood and trees", "Indicators: Deciduous", "Indicators: Coniferous", "Dry dead wood",
                                    "Meadowland"))

#set factor levels for ecocodes
eco.noabund$name <-
  factor(eco.noabund$name, levels = c("Aquatics", "Indicators: Running water", "Indicators: Standing water", "Open wet habitats", "Wetlands/marshes", "Mould beetles", "Halotolerant",
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

# Set region factor levels for visualisation
eco.abund$id <-
  factor(eco.abund$id, levels = c("5", "3", "1", "2", "4"))
eco.noabund$id <-
  factor(eco.noabund$id, levels = c("5", "3", "1", "2", "4"))

# Labeller for region facets
region.labs <- c(
  '1' = "British/Irish Isles",
  '2' = "Continental",
  '3' = "Iceland",
  '4' = "Meridional",
  '5' = "Scandinavia/Baltic")

# Stacked No abundance, sumrep
fig.abund <-
  eco.abund %>%
  ggplot(aes(fill=name, y=value, x=end)) +
  geom_bar(position="stack", stat="identity", just = 0, alpha = 0.9) +
  labs(y="Abundance weighted; %SumRep",
       x="Age (BP)",
       color=NULL) +
  scale_fill_manual(values = colors,
                    name = "Habitat") +
  guides(fill = guide_legend(nrow = 5,
                             reverse = T)) +
  scale_x_reverse(limits = c(16000, -500), breaks = scales::pretty_breaks(n = 10)) +
  coord_flip() +
  facet_grid2(~ id, scales = "free", axes = "margins", labeller = labeller(id = region.labs)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 14, face = "bold"),
        #axis.text.x = element_text(angle = 90),
        #axis.title.y = element_text(size = 14, face = "bold", colour = "black"),
        strip.text.x = element_text(size = 16, face = "bold", colour = "black"),
        strip.text.y = element_text(size = 16, face = "bold", colour = "black"),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.title = element_text(size = 12, face = "bold", colour = "black"),
        legend.text = element_text(size = 12),
        legend.position = "bottom")

# Stacked No abundance, sumrep
fig.noabund <-
  eco.noabund %>%
  ggplot(aes(fill=name, y=value, x=end)) +
  geom_bar(position="stack", stat="identity", just = 0, alpha = 0.9) +
  labs(y="No Abundance; %SumRep",
       x="Age (BP)",
       color=NULL) +
  scale_fill_manual(values = colors,
                    name = "Habitat") +
  guides(fill = guide_legend(nrow = 5,
                             reverse = T)) +
  scale_x_reverse(limits = c(16000, -500), breaks = scales::pretty_breaks(n = 10)) +
  coord_flip() +
  facet_grid2(~ id, scales = "free", axes = "margins", labeller = labeller(id = region.labs)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 14, face = "bold"),
        #axis.text.x = element_text(angle = 90),
        #axis.title.y = element_text(size = 14, face = "bold", colour = "black"),
        strip.text.x = element_text(size = 16, face = "bold", colour = "black"),
        strip.text.y = element_text(size = 16, face = "bold", colour = "black"),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.title = element_text(size = 12, face = "bold", colour = "black"),
        legend.text = element_text(size = 12),
        legend.position = "bottom")

#Save figure
ggsave("004-ecocodes-natural-divided-abund-fixed2.jpg",
       fig.abund,
       device = "jpg",
       here::here("analysis", "figures"),
       width=50,
       height=30,
       units = "cm",
       dpi = 300)

#Save figure
ggsave("004-ecocodes-natural-divided-noabund-fixed2.jpg",
       fig.noabund,
       device = "jpg",
       here::here("analysis", "figures"),
       width=50,
       height=30,
       units = "cm",
       dpi = 300)
