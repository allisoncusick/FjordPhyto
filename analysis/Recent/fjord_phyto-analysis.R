#### Load packages ----
packages <- c("tidyverse", "vegan", "lubridate", "dplyr", "tidyr", "magrittr",
              "readr", "phyloseq")


funlist <-  lapply(packages, function(x) {
  if (x %in% rownames(installed.packages())) {
    require(x, character.only = T)
  }else{
    install.packages(x, character.only = T); require(x, character.only = T)
  }
})

options(max.print = 100)

### Set directories ----
data_home <- "~/Documents/GitHub/FjordPhyto/data/"
figures_local <- "~/Documents/Fjord_Phyto/Figures/"
### Source the clean_metadata file OR load the QCd files ----
source("~/Documents/GitHub/FjordPhyto/analysis/Recent/clean_metadata-taxonomy.R")
## -> This provides: metadata and asv_table_rare


#### Analysis ----
##### Set up rules and filters for this specific analysis ----

## 1) Remove ASVs that has a total number of reads less than 5
asv_table_filter <- asv_table_rare %>%
  group_by(Feature.ID) %>%
  mutate(n_reads = sum(reads,na.rm = T)) %>%
  ungroup() %>%
  filter(n_reads > 5)
## We are left with this list of unique ASVs
unique_asvs <- unique(asv_table_filter$Feature.ID)

## 2) Keep only the phytoplankton ASVs
phyto_asvs <- unique(taxa_table_split$Feature.ID[!is.na(taxa_table_split$phytogroups)])

## 3) Remove the April sample from the 2022 samples
april_sample <- metadata$samples[metadata$month == 4]
metadata_filter <- metadata %>% filter(samples %in% samples_keep &
                                  samples != april_sample)


##### Set up phyloseq ----
## We need the OTU table (i.e. asv_table), the Sample table (i.e. metadata),
## and the Taxonomy table (i.e. taxa_table_split)
###### OTU Table ----
## Filter ASV table based off of the rules above and put into right format
subset_bio <- subset(asv_table_filter,
                     samples %in% unique(metadata_filter$samples) &
                       Feature.ID %in% phyto_asvs)

asvs_wide <- subset_bio %>%
  dplyr::select(Feature.ID, samples, reads) %>%
  pivot_wider(., values_from = reads, id_cols = Feature.ID, names_from = samples)
asvs_wide_df <- as.data.frame(asvs_wide[,-1]) #gets rid of Feature.ID column

rownames(asvs_wide_df) <- asvs_wide$Feature.ID #make FeatureID the rownames

in_biom <- otu_table(asvs_wide_df, taxa_are_rows = T)
taxa_names(in_biom) <- asvs_wide$Feature.ID

## Now were are left with "in_biom" where the rows are the Feature IDs and the columns are the samples

###### Sample Table ----
in_biom_metadata <- sample_data(metadata_filter)
sample_names(in_biom_metadata) <- in_biom_metadata$samples

###### Taxonomy Table -----
## Keep only the phytoplankton
taxatable <- taxa_table_split %>%
  filter(Feature.ID %in% phyto_asvs)
## Reorder the columns so that "phytogroups" is higher than "Genus" and "Species"
in_biom_tax <- tax_table(as.matrix(taxatable[,c(2:7, 11, 8:9)]))
taxa_names(in_biom_tax) <- taxatable$Feature.ID
colnames(in_biom_tax) <- colnames(taxatable[,c(2:7, 11, 8:9)])

## Quick check that the dimensions of the tables max sense:
## They should all have the same number of rows 

phylo_fjordphyto <- merge_phyloseq(in_biom, in_biom_tax, in_biom_metadata)


##### Define diversity metrics from rarefied reads -----

## Define list of all the groups/categories to caluclate diversity across
diversity_group <- c("phytogroups", unique(taxatable$phytogroups))

for (i in 1:length(diversity_group)){
  
  if (diversity_group[i] == "phytogroups"){
    taxa_pull <- taxatable
  }else{
    taxa_pull <- taxatable %>% filter(phytogroups == diversity_group[i])
  }
  
  
  piv_all <- subset_bio %>%
    filter(Feature.ID %in% taxa_pull$Feature.ID) %>%
    dplyr::select(c(Feature.ID, samples, reads))
  

  piv_all <- piv_all %>% group_by(samples) %>%
    mutate(prop_reads = reads/sum(reads)) %>% ungroup() #all reads in a sample, proportion of each taxa in each sample
  piv_all$reads <- NULL
  piv_all$raw_reads <- NULL
  
  piv_all <- piv_all %>%
    pivot_wider(names_from = "Feature.ID", values_from = "prop_reads",
                values_fill = 0)
  
  
  richness <- specnumber(piv_all[,-1])
  shannon <- diversity(piv_all[,-1], MARGIN = 1, index = "shannon")
  evenness <- shannon/log(richness)
  simpson <- diversity(piv_all[,-1], MARGIN = 1, index = "simpson")
  
  shannon_name <- paste0("shannon_", diversity_group[i])
  even_name <- paste0("evenness_", diversity_group[i])
  simpson_name <- paste0("simpson_", diversity_group[i])
  richness_name <- paste0("richness_", diversity_group[i])
  
  df_output <- data.frame(samples = piv_all[,1],
                          shannon = shannon,
                          evenness = evenness,
                          simpson = simpson,
                          richness = richness)
  colnames(df_output)[2:ncol(df_output)] <- str_replace(c(shannon_name, even_name, simpson_name, richness_name), " ", "_")
  
  if (i == 1){
    df_return <- df_output 
  }else{
    df_return <- left_join(df_return, df_output, by = "samples")
  }
}

##### Compile the full data set ----
full_df <- left_join(subset_bio, metadata_filter) %>%
  left_join(., taxatable) %>%
  left_join(., df_return)



## Calculate the days since variable 
days_since <- yday(full_df$Date) - yday("2000-11-01")
days_since[days_since < 0] <- days_since[days_since < 0] + 365
full_df$days_since <- days_since

species.filter <- full_df %>%
  group_by(Species, season) %>%
  reframe(sum = length(unique(samples[reads > 0]))) %>%
  filter(sum > 5) 

full_df_filter <- full_df %>%
  filter(reads > 0 & region == "middle")

### Frequencys --
full_df %>%
  group_by(Site_Name_2, season) %>%
  reframe(t_s = length(unique(samples))) %>%
  ggplot() +
  geom_density(aes(x = t_s, y = ..count.., color = season, group = season)) 



#### Run the assemblages ----
bio <- full_df %>%
  select(Feature.ID, Species, reads, samples) %>%
  pivot_wider(.,
  id_cols = "samples", names_from = "Species", values_from = "reads", values_fn = function(r)sum(r,na.rm = T), values_fill = 0) %>%
  column_to_rownames(var = "samples") 

correlation_matrix <- cor(bio, method = "spearman")
hierarchical_result <- hclust(dist(1 - correlation_matrix), method = "ward.D2")

plot(hierarchical_result, cex = 0.6, hang = -1)
rect.hclust(hierarchical_result, k = 3, border = 2:10)
groups <- cutree(hierarchical_result, k = 3)

myannotation <- as.data.frame(groups)
names(myannotation)[1] = "Cluster" 
myannotation$Cluster <- factor(myannotation$Cluster, levels= 1:3, 
                               labels=1:3)

pheatmap::pheatmap(correlation_matrix,
                   main = "Species Co-occurrence Matrix (Spearman's Correlation)",
                   fontsize_row = 6,
                   fontsize_col = 6, clustering_method = "ward.D2",
                   cutree_rows = 3, cutree_cols = 3,
                   annotation_col = myannotation)


bio_NMS1 <-  metaMDS(dist(1 - correlation_matrix),
                     k = 2,
                     maxit = 200, 
                     trymax = 100,
                     try = 50,
                     wascores= TRUE,
                     weakties = T,
                     expand = TRUE,
                     previous.best = T,
                     autotransform = FALSE,
                     na.rm = T)


data_scores_1 <- as.data.frame(scores(bio_NMS1$points)) %>%
  rownames_to_column(var = "Species") 

nmds_clust <- ggplot(data_scores_1) + 
  geom_point(aes(x = MDS1, y = MDS2,
                 color = factor(
                   Species, levels = names(groups), labels = groups))) +
  stat_ellipse(aes(x = MDS1, y = MDS2,
                   color = factor(
                     Species, levels = names(groups), labels = groups))) +
  theme_bw() + 
  theme(legend.position = c(0.23,0.94),
        legend.direction = "horizontal",
        legend.background = element_rect(color = "black", fill = "white"),
        axis.text = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 12),
        panel.background = element_rect(color = "black"),
        panel.grid = element_blank())
nmds_clust



#### Idk yet ----

model.list.names <- full_df %>%
  filter(reads > 0 & region == "middle") %>%
  group_by(phytogroups, Genus, Species, season) %>% group_keys() %>% mutate(list.id = paste(phytogroups, Genus, Species,season,sep = ":"))

species.model.list <- full_df %>%
  filter(reads > 0 & region == "middle") %>%
  group_by(phytogroups, Genus, Species, season) %>%
  group_map({~
      if(length(unique(.x$samples[.x$reads >0])) < 5 | nrow(.x) == 0){
        NULL
      }else{
      .x %>%
      loess(log10(reads) ~ days_since, ., span = 0.75, degree = 2)
      }
  })
names(species.model.list) <- model.list.names$list.id


loess.df <- do.call("cbind", lapply(species.model.list, function(dr){
        if(is.null(dr)){
          vector(length = length(lo.x))
        }else{
        predict(dr, lo.x)
        }
      })) %>%
  as_tibble() %>%
  select_if(~ !all(. == 0)) %>%
  mutate(days_since = 1:140) %>%
  pivot_longer(cols = !days_since, names_to = c("phytogroups", "Genus", "Species", "season"), names_sep = ":", values_to = "loess.fit")

loess.df$loess.fit[loess.df$loess.fit < 0] <- 0

day_since_df <- loess.df %>%
  group_by(phytogroups, Genus, Species, season) %>%
  reframe(sd = sqrt(
    sum(loess.fit*(days_since - days_since[which.max(loess.fit)])^2, na.rm = T)/
      (((length(loess.fit>0)-1)/length(loess.fit))*sum(
        loess.fit, na.rm = T))),
    day_max = days_since[which.max(loess.fit)])

days_order <- day_since_df %>%
  filter(season == "2017-2018")  %>%
  arrange(desc(day_max)) %>% pull(Species)

day_since_df %>%
  mutate(Species = factor(Species, levels = c(unique(day_since_df$Species)[!unique(day_since_df$Species) %in% days_order], unique(days_order)))) %>%
  ggplot() +
  geom_point(aes(y = Species, x = day_max, color = season)) +
  geom_errorbarh(aes(y = Species, xmin = day_max-sd, xmax = day_max+sd, color = season)) +
  labs(x = "Days since start of season") +
  scale_fill_manual(name = "", values = c("black", "blue", "green4", "red")) +
  my_theme +
  theme(axis.text.y = element_text(size = 6))


combined_df <- left_join(day_since_df, loess.df) %>%
  mutate(period = factor(case_when(day_max < 50 ~ "Early",
                            day_max >= 51 & day_max < 100 ~ "Mid",
                            day_max >= 100 ~ "Late"), levels = c("Early", "Mid", "Late")))

plot.periods <- data.frame(xmin = c(0, 51, 100), xmax = c(50, 100, 140), ymax = Inf, ymin = -Inf, period = c("Early","Mid","Late"))

day_change_plot <- combined_df %>%
  mutate(Species = factor(Species,
                          levels = c
                          (unique(day_since_df$Species)[
                            !unique(day_since_df$Species) %in% days_order],
                            unique(days_order))),
         phytogroups = factor(phytogroups,
                              levels = unique(combined_df$phytogroups),
                              labels = c("Cryptophytes", "Diatoms",
                                         "Dinoflagellates",
                                         "Green algae", "Haptophytes",
                                         "MAST", "Rhodophytes"))) %>%
  ggplot() +
  geom_point(aes(y = Species, x = day_max)) +
  geom_errorbarh(aes(y = Species, xmin = day_max-sd, xmax = day_max+sd)) +
  geom_rect(data = plot.periods, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = period), alpha = 1/5, show.legend = F) +
  labs(x = "Days since start of season") +
  scale_fill_manual(name = "", values = c("blue", "green4", "red")) +
  my_theme +
  theme(axis.text.y = element_text(size = 6),
        strip.text.y = element_text(
          size = 10, angle = 0, hjust = 0),
        strip.clip = "off") +
  facet_grid(phytogroups~season, scales = "free_y", space = "free_y")


ggsave(
  paste0(figures_local, "day_change_plot.png"),
  dpi = 300,
  width = 7,
  height = 7,
  unit = "in")

all_sp <- unique(day_since_df$Species)
### Day change organize by numbers 
days_order <- day_since_df %>%
  pivot_wider(id_cols = c(phytogroups, Genus, Species), names_from = season,
              values_from = day_max) %>%
  group_by(phytogroups) %>%
  arrange(desc(`2017-2018`), .by_group = T) %>% 
  mutate(sp_num = 1:n()) %>%
  pull(sp_num, Species)
  

plot.periods <- data.frame(xmin = c(0, 51, 100), xmax = c(50, 100, 140), ymax = Inf, ymin = -Inf, period = c("Early","Mid","Late"))
  
phytogroup_labels <- c("Cryptophytes", "Diatoms",
                       "Dinoflagellates",
                       "Green algae", "Haptophytes",
                       "MAST", "Rhodophytes")
names(phytogroup_labels) <- unique(combined_df$phytogroups)


mini_combined <- combined_df %>%
     group_by(phytogroups, Genus, Species, season) %>%
    distinct(day_max, sd)
  
sp_order_plot <- day_since_df %>%
  pivot_wider(id_cols = c(phytogroups, Genus, Species), names_from = season,
              values_from = day_max) %>%
  pivot_longer(cols = !c(phytogroups, Genus, Species),
               names_to = "season", values_to = "day_max") %>%
  mutate(day_max = replace_na(day_max, 141)) %>%
  pivot_wider(id_cols = c(phytogroups, Genus, Species), names_from = season,
              values_from = day_max) %>%
  mutate(phytogroups = factor(phytogroups,
                       levels = rev(unique(combined_df$phytogroups)))) %>%
  group_by(phytogroups) %>%
  arrange(desc(`2017-2018`), .by_group = T) %>% 
  mutate(Species = factor(Species, levels = Species, labels = Species)) %>%
  ungroup() %>%
  mutate(sp_num = rev(1:n())) %>%
  pivot_longer(cols = !c(phytogroups, Genus, Species, sp_num),
               names_to = "season", values_to = "day_max") %>%
  left_join(., mini_combined) %>%
  mutate(sd = replace_na(sd, 0)) %>%
  group_by(season) %>%
  group_map(~{
   .x %>%
      group_by(season, phytogroups) %>%
      arrange(desc(day_max),.by_group = T) %>%
      mutate(Species = factor(Species, levels = Species, labels = sp_num)) %>%
      ggplot() +
      geom_point(aes(y = Species, x = day_max,
                     color = ifelse(day_max == 141, NA, "black"))) +
      geom_errorbarh(aes(y = Species, xmin = day_max-sd, xmax = day_max+sd,
                         color = ifelse(day_max == 141, NA, "black"))) +
      geom_rect(data = plot.periods, aes(
        xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = period),
        alpha = 1/5, show.legend = F) +
      scale_fill_manual(name = "", values = c("blue", "green4", "red")) +
      scale_color_identity() +
      scale_x_continuous(name = "Days since start of season",
                         limits = c(-73,214),
                        breaks = c(0, 100, 200)) +
      {if(unique(.x$season)=="2017-2018"){
        scale_y_discrete(name = "",
                         labels = function(y){
                           paste0("***",
                                  rev(unique(.x$Species))[as.numeric(y)],
                                  "***",
                                  " - ", y)})
      }}+
      my_theme +
      theme(axis.text.y = ggtext::element_markdown(size = 6),
            strip.text.y = element_text(
              size = 10, angle = 0, hjust = 0),
            strip.clip = "off") +
      facet_grid(phytogroups~season, scales = "free_y", space = "free_y",
                 labeller = labeller(phytogroups = phytogroup_labels)) +
      {if(unique(.x$season)%in%unique(day_since_df$season)[1:3]){
        theme(strip.text.y = element_blank())}}+
      {if(unique(.x$season)%in%unique(day_since_df$season)[2:4]){
        theme(axis.title.y = element_blank())}}
  }, .keep = T) %>%
  patchwork::wrap_plots(., ncol = 4, axes = "collect_x", axis_titles = "collect_x")


ggsave(
  paste0(figures_local, "day_change_plot-ordered.png"),
  sp_order_plot,
  dpi = 300,
  width = 11,
  height = 7,
  unit = "in")

library(ggtext)
# arrange(desc(!! rlang::sym(unique(combined_df$season)[i])),.by_group = T) %>% 


full_df_clust <- left_join(left_join(full_df_filter, myannotation %>% rownames_to_column(var = "Species")), day_since_df)

mini_combined <- left_join(day_since_df,  myannotation %>%
                             rownames_to_column(var = "Species")) %>%
  group_by(Cluster, phytogroups, Genus, Species, season) %>%
  distinct(day_max, sd)

sp_order_plot <- left_join(day_since_df,  myannotation %>%
                             rownames_to_column(var = "Species")) %>%
  pivot_wider(id_cols = c(Cluster, phytogroups, Genus, Species),
              names_from = season,
              values_from = day_max) %>%
  pivot_longer(cols = !c(Cluster, phytogroups, Genus, Species),
               names_to = "season", values_to = "day_max") %>%
  mutate(day_max = replace_na(day_max, 141)) %>%
  pivot_wider(id_cols = c(Cluster, phytogroups, Genus, Species),
              names_from = season,
              values_from = day_max) %>%
  group_by(Cluster) %>%
  arrange(desc(`2017-2018`), .by_group = T) %>% 
  mutate(Species = factor(Species, levels = Species, labels = Species)) %>%
  ungroup() %>%
  mutate(sp_num = rev(1:n())) %>%
  pivot_longer(cols = !c(Cluster, phytogroups, Genus, Species, sp_num),
               names_to = "season", values_to = "day_max") %>%
  left_join(., mini_combined) %>%
  mutate(sd = replace_na(sd, 0)) %>%
  group_by(season) %>%
  group_map(~{
    .x %>%
      group_by(season, Cluster) %>%
      arrange(desc(day_max),.by_group = T) %>%
      mutate(Species = factor(Species, levels = Species, labels = sp_num)) %>%
      ggplot() +
      geom_point(aes(y = Species, x = day_max,
                     color = ifelse(day_max == 141, NA, "black"))) +
      geom_errorbarh(aes(y = Species, xmin = day_max-sd, xmax = day_max+sd,
                         color = ifelse(day_max == 141, NA, "black"))) +
      geom_rect(data = plot.periods, aes(
        xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = period),
        alpha = 1/5, show.legend = F) +
      scale_fill_manual(name = "", values = c("blue", "green4", "red")) +
      scale_color_identity() +
      scale_x_continuous(name = "Days since start of season",
                         limits = c(-73,214),
                         breaks = c(0, 100, 200)) +
      {if(unique(.x$season)=="2017-2018"){
        scale_y_discrete(name = "",
                         labels = function(y){
                           paste0("***",
                                  rev(unique(.x$Species))[as.numeric(y)],
                                  "***",
                                  " - ", y)})
      }}+
      my_theme +
      theme(axis.text.y = ggtext::element_markdown(size = 6),
            strip.text.y = element_text(
              size = 10, angle = 0, hjust = 0),
            strip.clip = "off") +
      facet_grid(Cluster~season, scales = "free_y", space = "free_y") +
      {if(unique(.x$season)%in%unique(day_since_df$season)[1:3]){
        theme(strip.text.y = element_blank())}}+
      {if(unique(.x$season)%in%unique(day_since_df$season)[2:4]){
        theme(axis.title.y = element_blank())}}
  }, .keep = T) %>%
  patchwork::wrap_plots(., ncol = 4, axes = "collect_x", axis_titles = "collect_x")




loess_fits_genus_plots <- combined_df %>%
  group_by(Genus) %>%
  group_map(~{
    .x %>%
  ggplot() +
  geom_line(aes(x = days_since, y = loess.fit, color = Species), linewidth = 1) +
  geom_segment(data = .x %>%
              group_by(Species, season) %>%
              reframe(
                day.max = unique(day_max),
                fit.max = loess.fit[day.max]),
              aes(x = day.max, xend = day.max,
                  y = -Inf, yend = fit.max,
                  color = Species),
              group = 1,
             linetype = "dashed", show.legend = F) +
  scale_color_brewer(palette = "Set2") +
  ggtitle(unique(.x$Genus)) +
  labs(x = "Days since start of season", y = bquote(log[10]~reads~(LOESS~fit))) +
  facet_wrap(~season) +
  my_theme
  }, .keep = T)

ggsave(
  filename = paste0(figures_local, "loess_genus_plots.pdf"), 
  plot = gridExtra::marrangeGrob(loess_fits_genus_plots, nrow=1, ncol=1), 
  width = 7, height = 5, units = "in")


periods_by_season <- combined_df %>%
  ggplot() +
  geom_line(aes(x = days_since, y = loess.fit, color = period, group = Species)) +
  facet_wrap(~season) +
  my_theme +
  labs(x = "Days since start of season",
       y = bquote(log[10]~reads~(LOESS~fit)), color = "")

ggsave(
  paste0(figures_local, "periods_by_season.png"),
  dpi = 300,
  width = 5,
  height = 5,
  unit = "in")
  
day_since_df <- full_df %>%
  filter(reads > 0 & region == "middle") %>%
  group_by(samples) %>%
  mutate(total_reads = sum(reads)) %>%
  group_by(Species, samples) %>%
  mutate(prop_reads = sum(reads)/total_reads) %>%
  group_by(Species, season) %>%
  mutate(sd = sqrt(
    sum(prop_reads*(days_since - days_since[which.max(prop_reads)])^2)/
                      (((length(prop_reads>0)-1)/length(prop_reads))*sum(
                        prop_reads))),
          day_max = days_since[which.max(prop_reads)]) 

#### Figures ----
#### 
#### The communities are similar but they do not succeed one to another.
#### Behrenfeld: Not resource partitioning but the food web in the ocean is so efficient that grazing shapes the community 
#### B

my_theme = theme_linedraw() + theme(text = element_text(size = 14), strip.background = element_blank(), strip.text = element_text(face = "bold", color = "black"))


days_order <- day_since_df %>%
  filter(season == "2017-2018")  %>%
  arrange(desc(day_max)) %>% pull(Species)

day_since_df %>%
  mutate(Species = factor(Species, levels = c(unique(day_since_df$Species)[!unique(day_since_df$Species) %in% days_order], unique(days_order)))) %>%
  ggplot() +
  geom_point(aes(y = Species, x = day_max, color = season)) +
  geom_errorbarh(aes(y = Species, xmin = day_max-sd, xmax = day_max+sd, color = season)) +
  labs(x = "Days since start of season") +
  scale_color_manual(name = "", values = c("black", "blue", "green4", "red")) +
  theme(axis.text.y = element_text(size = 6))


group_reads_loess <- full_df %>%
  filter(reads > 0 & phytogroups == "Diatoms") %>%
  group_by(phytogroups) %>%
  group_map(~{
    .x %>%
      group_by(Genus) %>%
      group_map(~{
        if (length(unique(.x$Species)) > 1){
          .x %>%
            ggplot() +
            geom_point(aes(x = days_since, y = log10(reads), color = Species),
                       alpha = 1/5, show.legend = T) +
            geom_smooth(aes(x = days_since, y = log10(reads),
                            group = Species, color = Species),
                        method = "loess",
                        show.legend = F, se = F) +
            facet_wrap(~Genus) +
            scale_x_continuous(name = "Days since start of season",
                               limits = c(0,140)) +
            scale_y_continuous(name = bquote(log[10]~reads)) +
            scale_color_discrete(name = unique(.x$Genus)) +
            theme_bw() +
            theme(strip.background = element_blank(),
                  strip.text = element_text(face = "bold"),
                  legend.position = c(0.2,0.8),
                  legend.text = element_text(size = 8),
                  legend.margin = margin(0,0,0,0, unit = "pt")) +
            guides(color = guide_legend(override.aes = list(alpha = 1)))
        }else{
          .x %>%
            ggplot() +
            geom_point(aes(x = days_since, y = log10(reads)), color = "black",
                       show.legend = T, alpha = 1/5) +
            geom_smooth(aes(x = days_since, y = log10(reads),
                            group = Species, color = Species), color = "black",
                        method = "loess",
                        show.legend = F, se = F) +
            facet_wrap(~Genus) +
            scale_x_continuous(name = "Days since start of season",
                               limits = c(0,140)) +
            scale_y_continuous(name = bquote(log[10]~reads)) +
            scale_color_discrete(name = "") +
            theme_bw() +
            theme(strip.background = element_blank(),
                  strip.text = element_text(face = "bold"),
                  legend.position = c(0.2,0.8),
                  legend.text = element_text(size = 8),
                  legend.margin = margin(0,0,0,0, unit = "pt"))
        }
      }, .keep = T) %>%
      patchwork::wrap_plots(., guides = "collect")}, .keep = T) 

ggsave(
  filename = paste0(figures_local, "group_reads_loess.pdf"), 
  plot = gridExtra::marrangeGrob(group_reads_loess, nrow=7, ncol=1), 
  width = 10, height = 10, units = "in")





### Who is in "Always early", "Mostly early", etc....




for (sp in 1:length(unique(full_df$Species))){
  
}
species_list <- unique(full_df$Species)

gam.1 <-full_df %>%
  filter(reads > 0 & Species == species_list[[1]] & season == "2017-2018") %>%
  mgcv::gam(reads ~ poly(days_since, 1), data = .) 

gam.2 <- full_df %>%
  filter(reads > 0 & Species == species_list[[60]] & season == "2017-2018") %>%
  stats::loess(reads ~ days_since, data = .) 
summary(gam.2)
chisq.test(predict(gam.2))
anova(gam.1, gam.2)



full_df %>%
  filter(reads > 0 & region == "middle" & Species == "Dino-Group-I-Clade-1_X_sp." & season == "2021-2022") %>%
  ggplot() +
  geom_point(aes(x = days_since, y = reads)) +
  stat_smooth(aes(x = days_since, y = reads), method = "loess", span = 0.75, method.args = list(degree = 2)) +
  my_theme

day_since_df %>%
  filter(reads > 0 & Genus == "Chaetoceros") %>%
  ggplot() +
  geom_point(aes(x = days_since, y = log10(reads), color = Species),
             show.legend = T, alpha = 1/5) +
  geom_smooth(aes(x = days_since, y = log10(reads),
                  group = Species, color = Species),
              method = "loess", span = 0.75, 
              method.args = list(degree = 2),
              show.legend = F, se = F) +
  facet_wrap(~season)
  
### Gam for reagion, gam for seaosn??

day_since_df %>%
  ggplot() +
  geom_point(aes(y = Site_Name_2, x = day_max, color = season)) +
  # geom_errorbarh(aes(y = Site_Name_2, xmin = day_max-sd, xmax = day_max+sd, color = season)) +
  labs(x = "Days since start of season") +
  scale_color_manual(name = "", values = c("black", "blue", "green4", "red")) +
  my_theme



 