##### Analysis and Figure production for Chapter 3

### Load necessary packages
packages <- c("tidyverse", "dplyr", "tidyr", "readxl", "phyloseq", 
              "ggtext", "ggh4x", "cowplot", "patchwork", "ggpmisc",
              "scales", "gridExtra", "broom", "Hmisc",  "readxl", 
              "geosphere", "vegan")


funlist <-  lapply(packages, function(x) {
  if (x %in% rownames(installed.packages())) {
    require(x, character.only = TRUE)
  } else {
    install.packages(x, character.only = TRUE)
    require(x, character.only = TRUE)
  }
})

options(max.print = 100) #Personal preference

### Setup date for saving files
todays_date <- format(Sys.Date(), "%m%d%Y")

### Source the metadata and ASV QC scripts
source("analysis/asv_data_process-2024.r")
source("analysis/metadata_process-2024.r") 
### OR
### Load in data
## Find most recent file
pattern = "_processed-"
relevant_files <- list.files("data", pattern, recursive = T, full.names = T)
recent_files <- do.call("rbind", lapply(relevant_files, file.info)) %>%
  rownames_to_column(var = "file") %>%
  mutate(group = lapply(regmatches(file,
                            regexec("data.\\s*(.*?)\\s*.fjord",
                                    file)),
                            function(x)x[2]),
         t_diff = difftime(Sys.Date(), ctime, units = "days")) %>%
  group_by(group) %>%
  reframe(file = file[which.min(t_diff)]) %>%
  pull(file)
lapply(recent_files, load)

### Output = diversity_df, asv_table_rare, metadata_filter

## Directories
figure_files <- c("~/Documents/Alaina/FjordPhyto/Chapter_3/Figures/")

## Data processing -----
### Filter out data to include only same samples across all dfs
metadata_filter <- metadata_w_samples %>%
  #drop_na(secchi_depth) %>%
  rename(sample = Genetics_18sv9_Sample_ID) %>%
  filter(!sample %in% samples_removed)

### Find the ASVs that are present in more than 5 samples and have more than 5 
### reads across all samples
asv_filter <- asv_table_rare %>%
  drop_na(phytogroups) %>%
  filter(sample %in% unique(metadata_filter$sample)) %>%
  group_by(Feature.ID) %>%
  reframe(total_reads = sum(rare_reads, na.rm = T),
          total_samples = length(unique(sample[rare_reads > 0]))) %>%
  filter(total_reads > 5 & total_samples > 5) %>%
  pull(Feature.ID) %>%
  unique()

### Only include the above ASVs
asv_table_filter <- asv_table_rare %>%
  filter(rare_reads > 0) %>%
  drop_na(phytogroups) %>%
  filter(sample %in% unique(metadata_filter$sample)) %>%
  filter(Feature.ID %in% asv_filter)

diversity_filter <- diversity_df %>%
  filter(sample %in% unique(metadata_filter$sample))

#### Plotting matter ----

## Set up grouped color scale
group <- "phytogroups"
subgroup <- "Species"

taxatable <- asv_table_filter %>% 
  select(Feature.ID, Kingdom, Supergroup, Division, Class, Order, Family, Genus, Species, phytogroups) %>%
  distinct() %>%
  as.data.frame()
in_biom_tax <- tax_table(taxatable %>% select(-Feature.ID) %>% as.matrix()) 
taxa_names(in_biom_tax) <- taxatable$Feature.ID

categories <- aggregate(as.formula(paste(subgroup, group, sep="~" )),
                        in_biom_tax, function(x) length(unique(x)))
categories_pal <- c("Red-Blue", "YlGn", "Teal", "Purp","OrYel", "Burg", "Light Grays")
names(categories_pal) <- names(group_colors)

sp_colors <- unlist(lapply(1:nrow(categories),
                           function(i){
                             colorspace::sequential_hcl(
                               n = categories$Species[i],
                               pal = c(
                                 categories_pal[categories$phytogroups[i]]))}))

sp_colors_df <- asv_table_filter %>%
  distinct(phytogroups, Species) %>%
  arrange(phytogroups, Species) %>%
  mutate(colors = sp_colors) %>%
  group_by(phytogroups) %>%
  group_modify(~add_row(.x, .before = 0)) %>% 
  ungroup() %>%
  mutate(
    phytogroups = ifelse(phytogroups == "Greenalgae", "Green algae",
                         phytogroups),
    colors = ifelse(is.na(Species), "white", colors),
    sublabel = ifelse(is.na(Species),
                      sprintf("**%s**", 
                              as.character(phytogroups)), 
                      as.character(Species)),
    shapes = case_when(phytogroups == "Diatoms" ~ 22,
                       phytogroups == "Dinoflagellates" ~ 23,
                       phytogroups == "MAST" ~ 24,
                       TRUE ~ 21)) %>%
  as.data.frame()


group_colors <- c("#F06400", "#00F064","#008CF0","grey10",
                  "gold2","#6400F0","#F0008C")
group_colors_df <- data.frame(phytogroups = categories$phytogroups,
                              colors = group_colors,
                              labels = sp_colors_df$sublabel[is.na(sp_colors_df$Species)])

phyto_labels <- c(
  richness_Cryptophytes = "Cryptophytes",
  richness_Diatoms = "Diatoms", 
  richness_Dinoflagellates = "Dinoflagellates",
  richness_Haptophytes = "Haptophytes",
  richness_MAST = "MAST", 
  richness_Rhodophytes = "Rhodophytes", 
  richness_Greenalgae = "Greenalgae")



my_theme = theme_linedraw() + theme(text = element_text(size = 14), strip.background = element_blank(), strip.text = element_text(face = "bold", color = "black"))

### Set up common variables for quick actions
env_vars <- c("Surface_Temperature_", "Surface_Salinity_", "days_since")
p_groups <- names(group_colors)

full_df <- left_join(metadata_filter, asv_table_filter) %>%
  mutate(prop_reads = rare_reads/7000) %>%
  left_join(diversity_filter)

#### Set up voxel labels and ids for plotting
voxel_ids <- expand_grid(
  days_cuts = factor(c("Early", "Mid", "Late")),
  sal_cuts = factor(c("High", "Mid", "Low")),
  temp_cuts = factor(c("Low", "Mid", "High"))) %>%
  mutate(
    group_id = 1:nrow(.),
    group_box = factor(1:nrow(.),
                       levels = 1:nrow(.),
                       labels = paste0(
                         days_cuts, ": ",
                         sal_cuts, " Sal, ",
                         temp_cuts, " Temp")),
    x_val = rep(c(-1.2, 0.2, 1.2), 9),
    y_val = rep(rep(c(34.4, 33.4, 32.4), each = 3), 3))

#### Two options: Maximum for each species over all 
#### (which temp, sal, day of season, any season)
bins_df <- full_df %>%
  drop_na(phytogroups) %>%
  filter(prop_reads == max(prop_reads, na.rm = T),
         .by = c("phytogroups", "Genus", "Species")) %>%
  select(phytogroups, Genus, Species, prop_reads, seasonyear, 
         Surface_Temperature_, Surface_Salinity_, days_since) %>%
  mutate(
    days_cuts = factor(case_when(days_since <= 50 ~ "Early",
                                 days_since <= 100 &
                                   days_since > 50 ~ "Mid",
                                 TRUE ~ "Late"),
                       levels = c("Early", "Mid", "Late")),
    sal_cuts = factor(case_when(Surface_Salinity_ <= 32.5 ~ "Low",
                                Surface_Salinity_ <= 33.5 &
                                  Surface_Salinity_ > 32.5 ~ "Mid",
                                TRUE ~ "High")),
    temp_cuts = factor(case_when(Surface_Temperature_ <= 0 ~ "Low",
                                 Surface_Temperature_ <= 1 &
                                   Surface_Temperature_ > 0 ~ "Mid",
                                 TRUE ~ "High"))) %>%
  left_join(., voxel_ids)

bins_df_season <- full_df %>%
  drop_na(phytogroups) %>%
  filter(prop_reads == max(prop_reads, na.rm = T),
         .by = c("phytogroups", "Genus", "Species", "seasonyear")) %>%
  select(phytogroups, Genus, Species, prop_reads, seasonyear, 
         Surface_Temperature_, Surface_Salinity_, days_since) %>%
  mutate(
    days_cuts = factor(case_when(days_since <= 50 ~ "Early",
                                 days_since <= 100 &
                                   days_since > 50 ~ "Mid",
                                 TRUE ~ "Late"),
                       levels = c("Early", "Mid", "Late")),
    sal_cuts = factor(case_when(Surface_Salinity_ <= 32.5 ~ "Low",
                                Surface_Salinity_ <= 33.5 &
                                  Surface_Salinity_ > 32.5 ~ "Mid",
                                TRUE ~ "High")),
    temp_cuts = factor(case_when(Surface_Temperature_ <= 0 ~ "Low",
                                 Surface_Temperature_ <= 1 &
                                   Surface_Temperature_ > 0 ~ "Mid",
                                 TRUE ~ "High"))) %>%
  left_join(., voxel_ids)

write.csv(bins_df_season, paste0(figure_files, "bins_df_season.csv"))
write.csv(bins_df, paste0(figure_files, "bins_df.csv"))
          
## Figure: Environmental variables over time -----

plot_a <- metadata_filter %>%
  filter(seasonyear!="2022-2023") %>%
  ggplot() +
  geom_point(aes(x = days_since, y = Surface_Temperature_)) +
  geom_smooth(aes(x = days_since, y = Surface_Temperature_),
              method = "gam", color = "black") +
  geom_point(aes(x = days_since, y = Surface_Salinity_-32), color = "red") +
  geom_smooth(aes(x = days_since, y = Surface_Salinity_-32), color = "red",
              method = "gam") +
  theme_bw() +
  my_theme +
  scale_x_continuous(name = "Days since Nov. 1", expand = c(1E-2,1E-2)) +
  scale_y_continuous(name = "Temperature",
                     sec.axis = sec_axis(name = "Salinity", trans = ~.+32)) +
  theme(axis.text.y.right = element_text(color = "red"),
        axis.title.y.right = element_text(color = "red"),
        axis.title.x = element_blank())

plot_b <- metadata_filter %>%
  filter(seasonyear!="2022-2023") %>%
  ggplot() +
  geom_point(aes(x = days_since, y = Surface_Temperature_)) +
  geom_smooth(aes(x = days_since, y = Surface_Temperature_),
              method = "gam", color = "black") +
  geom_point(aes(x = days_since, y = Surface_Salinity_-32), color = "red") +
  geom_smooth(aes(x = days_since, y = Surface_Salinity_-32), color = "red",
              method = "gam") +
  theme_bw() +
  my_theme +
  scale_x_continuous(name = "Days since Nov. 1", expand = c(1E-2,1E-2)) +
  scale_y_continuous(name = "Temperature",
                     sec.axis = sec_axis(name = "Salinity", trans = ~.+32)) +
  theme(axis.text.y.right = element_text(color = "red"),
        axis.title.y.right = element_text(color = "red")) +
  facet_wrap(~seasonyear, nrow = 1)

plot_c <- plot_a / plot_b

ggsave(paste0(figure_files, "ts-time-all-season.png"),
       plot = plot_c,
       width = 10, height = 5, dpi = 300)


##### Figure: Diversity over time ----
#### Linear regression of the relationship between days since, salinity, and temperature versuse the richness of each phytoplankton group

 lm_vals <- do.call("rbind", full_df %>%
  select(-grep("_all", colnames(.))) %>% 
  select(-contains("_phytogroups")) %>%
  pivot_longer(cols = contains("richness_"), names_to = "group",
               values_to = "richness") %>%
  drop_na(richness) %>%
  pivot_longer(cols = c(
    "days_since", "Surface_Temperature_", "Surface_Salinity_"), 
    names_to = "var",
    values_to = "values") %>%
  filter(richness > 0) %>%
  drop_na(richness) %>%
  select(group, richness, values, var) %>%
  distinct() %>%
  group_by(group) %>%
  group_map(~{do.call("rbind",
    .x %>%
      group_by(var) %>%
      group_map(~{
        .x %>%
      lm(richness~values, data = .) %>%
      summary %>%
          glance %>%
          as.data.frame() %>%
          mutate(
            group = unique(.x$group),
            var = unique(.x$var)
          )
      }, .keep = T))}, .keep = T))

write.csv(lm_vals, paste0(figure_files, "lm_vals.csv"))


plot_a <- full_df %>%
  select(-grep("_all", colnames(.))) %>% 
  select(-contains("_phytogroups")) %>%
  pivot_longer(cols = contains("richness_"), names_to = "group",
               values_to = "richness") %>%
  drop_na(richness) %>%
  pivot_longer(cols = c(
    "days_since", "Surface_Temperature_", "Surface_Salinity_"), 
    names_to = "var",
    values_to = "values") %>%
  filter(richness > 0) %>%
  mutate(group = factor(group, levels = names(phyto_labels),
                        labels = gsub("richness_", "", names(phyto_labels)))) %>%
  ggplot() +
  geom_point(aes(x = values, y = richness, color = group)) +
  stat_poly_line(aes(x = values, y = richness, color = group),
                 method = "lm", show.legend = F) +
  theme_bw() +
  my_theme +
  facet_grid(~var, scales = "free",
             labeller = labeller(var = c(
               days_since = "Days since Nov. 1",
               Surface_Salinity_ = "Surface Salinity (PSU)",
               Surface_Temperature_ = "Surface Temperature (*C)"
             )), switch = "x") +
  scale_color_manual(name = "", values = group_colors,
                     labels = c("Greenalgae" = "Green algae")) +
  theme(strip.placement = "outside",
        axis.title.x = element_blank()) +
  labs(y = "# of species") +
  coord_cartesian(expand = 0)

ggsave(
  filename = paste0(figure_files, "all_richness_vars.png"),
  plot_a,
  width = 10,
  height = 5,
  units = "in",
  dpi = 300
)

#### Figure: Correlation between species and variables ----
bio_species <- full_df %>%
  filter(!is.na(phytogroups)) %>%
  select(Feature.ID, Species, phytogroups, rare_reads, sample, days_since,
         Surface_Salinity_, Surface_Temperature_) %>%
  pivot_wider(id_cols = c("sample", "days_since", "Surface_Temperature_", "Surface_Salinity_"), names_from = "Species", values_from = "rare_reads", values_fn = function(r)sum(r,na.rm = T), values_fill = 0) %>% 
  column_to_rownames(var = "sample") 

corr_plot_r <- rcorr(as.matrix(bio_species), type = "pearson")$r
corr_plot_p <- rcorr(as.matrix(bio_species), type = "pearson")$P

plot_a <- as.data.frame(corr_plot_r) %>%
  rownames_to_column("var2") %>%
  pivot_longer(cols = -var2, names_to = "var1", values_to = "r") %>%
  mutate(P = as.vector(corr_plot_p)) %>%
  filter(var1 != var2) %>%
  filter(var2 %in% env_vars & !var1 %in% env_vars ) %>%
  mutate(spec = factor(var1, levels = rev(sp_colors_df$sublabel))) %>%
  filter(P < 0.05) %>%
  ggplot() +
  geom_point(aes(x = factor(var2), y = spec, color = r, size = P)) +
  scale_color_gradient2(low = "blue", high = "red") +
  theme_bw() +
  labs(x = "", y = "", fill = "Correlation") +
  my_theme +
  theme(axis.text.y = element_markdown(size = 8)) +
  scale_y_discrete(drop = F) + 
  scale_x_discrete(labels = c("Days since Nov. 1", 
                              "Salinity",
                              "Temperature")) +
  scale_size_continuous(range = c(2, 5), transform = "reverse")

ggsave(
  filename = paste0(figure_files, "species_corr.png"),
  plot_a,
  width = 6,
  height = 9,
  units = "in",
  dpi = 300
)

##### Figure: Maximum for each Species -----
max_sub <- full_df %>%
  group_by(Species) %>%
  reframe(maxr = max(prop_reads, na.rm = T),
          temp = Surface_Temperature_[which.max(prop_reads)],
          salinity = Surface_Salinity_[which.max(prop_reads)])

plot_a <- full_df %>%
  filter(prop_reads > 0) %>%
  ggplot() +
  geom_point(aes(y = Surface_Temperature_, 
                 x = Surface_Salinity_,
                 color = prop_reads)) +
  geom_point(data = max_sub, aes(y = temp, x = salinity),
             color = "red", shape = 21, size = 3, stroke = 1) +
  facet_wrap(~Species) +
  theme_bw() +
  my_theme +
  scale_y_continuous(name = "Surface Temperature") +
  scale_x_continuous(name = "Surface Salinity") +
  scale_color_viridis_c(name = "Prop. reads") +
  theme(strip.text = element_text(size = 9))

ggsave(paste0(figure_files, "species_TS-max_red.png"),
       plot = plot_a,
       width = 18, height = 12, dpi = 300)

#### Figure: Voxels ----
plot_a <- bins_df %>%
  mutate(Species = factor(Species, levels = sp_colors_df$sublabel)) %>%
  ggplot() +
  geom_vline(xintercept = c(32.5, 33.5),
             linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = c(0, 1) , linetype = "dashed", color = "grey40") +
  geom_jitter(aes(y = Surface_Temperature_,
                  x = Surface_Salinity_,
                  fill = Species), shape = 21,
              color = "black", size = 5,
              width = 0.05, height = 0.25) +
  scale_fill_manual(name = "", values = sp_colors_df$colors,
                    labels = sp_colors_df$sublabel, drop = F) +
  theme_bw() +
  my_theme +
  theme(legend.text = element_markdown(size = 9)) +
  labs(y = bquote(Surface~Temperature~(degree*C)),
       x = "Surface Salinity (PSU)") +
  scale_shape_discrete(name = "")

ggsave(paste0(figure_files, "dotplot_species-voxel-all-jitter2.png"),
       plot = plot_a,
       width = 14, height = 6, dpi = 300)

#### Figure: Voxels by season ----

plot_a <- bins_df_season %>%
  mutate(Species = factor(Species, levels = sp_colors_df$sublabel)) %>%
  ggplot() +
  geom_vline(xintercept = c(32.5, 33.5),
             linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = c(0, 1) , linetype = "dashed", color = "grey40") +
  geom_jitter(aes(y = Surface_Temperature_,
                  x = Surface_Salinity_,
                  fill = Species), shape = 21,
              color = "black", size = 3,
              width = 0.05, height = 0.25) +
  scale_fill_manual(name = "", values = sp_colors_df$colors,
                    labels = sp_colors_df$sublabel, drop = F) +
  theme_bw() +
  my_theme +
  theme(legend.text = element_markdown(size = 9)) +
  labs(y = bquote(Surface~Temperature~(degree*C)),
       x = "Surface Salinity (PSU)") +
  scale_shape_discrete(name = "") +
  facet_wrap(~seasonyear)

ggsave(paste0(figure_files, "dotplot_species-voxel-all-jitter2.png"),
       plot = plot_a,
       width = 14, height = 6, dpi = 300)

#### Figure: Change over a season ----
gg_test_season1 <- merge(
  bins_df_season,
  bins_df_season %>% 
    filter(seasonyear == "2017-2018") %>% 
    group_by(Species) %>% 
    reframe(
      x_1 = mean(Surface_Salinity_),
      y_1 = mean(Surface_Temperature_),
      d_1 = mean(days_since)), 
  by = "Species")

gg_test_groups <- merge(
  bins_df_season,
  bins_df_season %>% 
    filter(seasonyear == "2017-2018") %>% 
    group_by(phytogroups) %>% 
    reframe(
      x_1 = mean(Surface_Salinity_),
      y_1 = mean(Surface_Temperature_),
      d_1 = mean(days_since)), 
  by = "phytogroups")

change_species <- ggplot(gg_test_season1, aes(
  Surface_Salinity_,
  Surface_Temperature_,
  color = seasonyear)) +
  geom_point(size = 2) +
  geom_segment(aes(
    x = x_1,
    y = y_1,
    xend = Surface_Salinity_,
    yend=Surface_Temperature_)) +
  facet_wrap(~Species) +
  theme_bw() +
  my_theme +
  labs(x = "Surface Salinity (PSU)", 
       y = bquote(Surface~Temperature~(degree*C)),
       color = "") +
  theme(strip.text = element_text(size = 7))

ggsave(paste0(figure_files, "change_max_over_season-species.png"),
       plot = change_species,
       width = 15, height = 10, dpi = 300)

test_text <- data.frame(
  Species = "",
  x = c(-37.5, 50, -1.25, 1, -1.75, 1),
  label = c("Earlier", "Later",
            "Saltier", "Fresher",
            "Warmer", "Colder"),
  variable = c("diff_days", "diff_days",
               "diff_sal", "diff_sal",
               "diff_temp", "diff_temp"))

gg_shift_groups <- gg_test_groups %>%
  group_by(phytogroups) %>%
  filter(seasonyear != "2017-2018") %>%
  reframe(
    diff_sal = mean(x_1 - Surface_Salinity_),
    sd_sal = sd(x_1 - Surface_Salinity_),
    diff_temp = mean(Surface_Temperature_ - y_1),
    sd_temp = sd(Surface_Temperature_ - y_1),
    diff_days = mean(days_since - d_1),
    sd_days = sd(days_since - d_1)) %>%
  pivot_longer(cols = c(diff_sal, diff_temp, diff_days),
               names_to = "variable",
               values_to = "value") %>%
  pivot_longer(cols = c(sd_sal, sd_temp, sd_days),
               names_to = "variable2",
               values_to = "sd") %>%
  mutate(Species = paste0("**", phytogroups, "**"),
         Species = case_when(Species == "**Greenalgae**" ~ "**Green algae**",
                             TRUE ~ Species)) %>%
  filter(gsub("diff_","", variable) == gsub("sd_","", variable2))


avg_change <- gg_test_season1 %>%
  group_by(Species) %>%
  filter(seasonyear != "2017-2018") %>%
  reframe(
    diff_sal = mean(x_1 - Surface_Salinity_),
    diff_temp = mean(Surface_Temperature_ - y_1),
    diff_days = mean(days_since - d_1)) %>%
  pivot_longer(cols = c(diff_sal, diff_temp, diff_days),
               names_to = "variable",
               values_to = "value") %>%
  mutate(Species = factor(Species, 
                          levels = c(rev(sp_colors_df$sublabel), ""))) %>%
  ggplot() +
  geom_point(aes(y = Species, x = value)) +
  geom_point(data = gg_shift_groups, aes(y = Species, x = value,
                                         color = Species),
             show.legend = F, size = 3) +
  geom_segment(data = gg_shift_groups, aes(
    y = Species, x = value - sd, xend = value + sd,
    color = Species), show.legend = F) +
  facet_wrap(~variable, scales = "free_x", switch = "x",
             labeller = labeller(variable = c(
               "diff_days" = "Avg. shift in peak day\nfrom 2017-2018",
               "diff_sal" = "Avg. shift in peak salinity\nfrom 2017-2018",
               "diff_temp" = "Avg. shift in peak temperature\nfrom 2017-2018"))) +
  scale_y_discrete(name = "", drop = F) +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  geom_label(data = test_text, aes(
    y = Species, x = x, label = label, hjust = 0.5, vjust = 0.5)) +
  scale_color_manual(name = "", values = group_color_df$colors,
                     labels = group_colors_df$labels) +
  theme_bw() +
  my_theme +
  theme(axis.text.y = element_markdown(size = 9),
        axis.title.x = element_blank(),
        strip.placement = "outside")
  
ggsave(paste0(figure_files, "avg_change_since2017.png"),
       plot = avg_change,
       width = 9, height = 12, dpi = 300)



