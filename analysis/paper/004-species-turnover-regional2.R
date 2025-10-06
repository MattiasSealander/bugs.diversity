pacman::p_load(codyn, ggh4x, ggsci, ggplot2, IRanges, SRS, tidyverse, vegan)
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

#Interactive SRS shiny app for identifying suitable Cmin in SRS analysis
#if (interactive()) {SRS.shiny.app(data)}

# Prepare data for SRS scaling
srs.prep <-
  hits %>%
  select (site.x, sample_group, sample, taxon, end = end.1, abundance) %>%
  distinct() %>%
  mutate(sample = paste(sample, sample_group, site.x, end, sep = "@")) %>%
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
  inner_join((hits %>% select(region, site.x, sample_group, sample)),
             by = c("sample", "sample_group", "site" = "site.x"), relationship = "many-to-many") %>%
  distinct() %>%
  select(region, end, taxon = species, abundance = value) %>%
  group_by(region, end, taxon) %>%
  summarise(tot_abund = sum(abundance)) %>%
  mutate(end = as.numeric(end)) %>%
  mutate(end = if_else(end != 0, end * -1, end)) %>%
  as.data.frame()

# List with data split by region
Listdf <- list(split(raw.mat[,-1], raw.mat$region), (split(srs.abund[,-1], srs.abund$region)))

# Initialize the output lists with the same structure as Listdf
app.list <- vector("list", length(Listdf))
disapp.list <- vector("list", length(Listdf))
total.list <- vector("list", length(Listdf))

# Loop through data and generate diversity metrics
for (i in seq_along(Listdf)) {
  app.list[[i]] <- vector("list", length(Listdf[[i]]))

  for (j in seq_along(Listdf[[i]])) {
    # Only process if the current element is a data frame
    if (is.data.frame(Listdf[[i]][[j]])) {
      result <- turnover(
        df = Listdf[[i]][[j]],
        time.var = "end",
        species.var = "taxon",
        abundance.var = "tot_abund",
        metric = "appearance"
      ) %>%
        mutate(
          end = if_else(.data$end != 0, .data$end * -1, .data$end),
          value = .data$appearance,
          type = "Appearance"
        ) %>%
        select(end, value, type)  # Rearrange columns if needed

      # Save result
      app.list[[i]][[j]] <- result
    }
  }
  disapp.list[[i]] <- vector("list", length(Listdf[[i]]))

  for (j in seq_along(Listdf[[i]])) {
    # Only process if the current element is a data frame
    if (is.data.frame(Listdf[[i]][[j]])) {
      result <- turnover(
        df = Listdf[[i]][[j]],
        time.var = "end",
        species.var = "taxon",
        abundance.var = "tot_abund",
        metric = "disappearance"
      ) %>%
        mutate(
          end = if_else(.data$end != 0, .data$end * -1, .data$end),
          value = .data$disappearance,
          type = "Disappearance"
        ) %>%
        select(end, value, type)  # Rearrange columns if needed

      # Save result
      disapp.list[[i]][[j]] <- result
    }
  }
  total.list[[i]] <- vector("list", length(Listdf[[i]]))

  for (j in seq_along(Listdf[[i]])) {
    # Only process if the current element is a data frame
    if (is.data.frame(Listdf[[i]][[j]])) {
      result <- turnover(
        df = Listdf[[i]][[j]],
        time.var = "end",
        species.var = "taxon",
        abundance.var = "tot_abund",
        metric = "total"
      ) %>%
        mutate(
          end = if_else(.data$end != 0, .data$end * -1, .data$end),
          value = .data$total,
          type = "Total"
        ) %>%
        select(end, value, type)  # Rearrange columns if needed

      # Save result
      total.list[[i]][[j]] <- result
    }
  }
}

# Prepare table for facetted plot, "id" indicates which list object results come from
app.raw.df <-
  lapply(app.list[[1]], as.data.frame) %>%
  bind_rows(.id = "id")
app.srs.df <-
  lapply(app.list[[2]], as.data.frame) %>%
  bind_rows(.id = "id")
disapp.raw.df <-
  lapply(disapp.list[[1]], as.data.frame) %>%
  bind_rows(.id = "id")
disapp.srs.df <-
  lapply(disapp.list[[2]], as.data.frame) %>%
  bind_rows(.id = "id")
total.raw.df <-
  lapply(total.list[[1]], as.data.frame) %>%
  bind_rows(.id = "id")
total.srs.df <-
  lapply(total.list[[2]], as.data.frame) %>%
  bind_rows(.id = "id")

# Merged dataframes
spec.raw.turnover <-
  bind_rows(app.raw.df, disapp.raw.df, total.raw.df) %>%
  mutate(id = factor(id, levels = c("5", "3", "1", "2", "4"),
                     labels = c("Scandinavia/Baltic", "Iceland", "British/Irish Isles", "Continental", "Meridional")))
# Merged dataframes
spec.srs.turnover <-
  bind_rows(app.srs.df, disapp.srs.df, total.srs.df) %>%
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

fig.raw <-
  ggplot() +
  geom_rect(data = rects,
            aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = col),
            alpha = 0.5) +
  geom_line(data = spec.raw.turnover, aes(x=end, y=value, linetype = type, color = type), linewidth = 1) +
  geom_point(data = spec.raw.turnover, aes(x=end, y=value, color = type), size = 2) +
  scale_x_reverse(limits = c(16000, -500), breaks = scales::pretty_breaks(n = 10)) +
  scale_color_manual(name = "Metric",,
                     values = c("green4", "red4", "black")) +
  scale_linetype_manual(name = "Metric",
                        values=c("longdash", "dotted", "solid"))+
  scale_fill_jco(name = "Time Periods") +
  facet_grid2(id ~ ., scales = "free", axes = "margins", labeller = labeller(type = spec.raw.turnover$type)) +
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


fig.srs <-
  ggplot() +
  geom_rect(data = rects,
            aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = col),
            alpha = 0.5) +
  geom_line(data = spec.srs.turnover, aes(x=end, y=value, linetype = type, color = type), linewidth = 1) +
  geom_point(data = spec.srs.turnover, aes(x=end, y=value, color = type), size = 2) +
  scale_x_reverse(limits = c(16000, -500), breaks = scales::pretty_breaks(n = 10)) +
  scale_color_manual(name = "Metric",,
                     values = c("green4", "red4", "black")) +
  scale_linetype_manual(name = "Metric",
                        values=c("longdash", "dotted", "solid"))+
  scale_fill_jco(name = "Time Periods") +
  facet_grid2(id ~ ., scales = "free", axes = "margins", labeller = labeller(type = spec.srs.turnover$type)) +
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
ggsave("004-relative-total-turnover-divided-raw.jpg",
       fig.raw,
       device = "jpg",
       here::here("analysis", "figures"),
       width=40,
       height=30,
       units = "cm",
       dpi = 300)

#Save figure
ggsave("004-relative-total-turnover-divided-srs.jpg",
       fig.srs,
       device = "jpg",
       here::here("analysis", "figures"),
       width=40,
       height=30,
       units = "cm",
       dpi = 300)
