packages <- c(
  "tidyverse", "ragg", "patchwork", "ggmap", "vegan",
  "RColorBrewer", "geosphere", "lubridate", "dplyr",
  "ggplot2", "tidyr", "dendextend", "factoextra",
  "cluster", "NbClust", "pheatmap", "fantaxtic", "phyloseq",
  "tidyverse", "magrittr", "ggnested", "knitr", "kableExtra",
  "gridExtra", "ggrepel", "sf", "ggh4x", "ggmagnify",
  "cowplot", "rstatix", "ggpubr", "mapproj")


funlist <-  lapply(packages, function(x) {
  if (x %in% rownames(installed.packages())) {
    require(x, character.only = T)
  }else{
    install.packages(x, character.only = T); require(x, character.only = T)
  }
})

options(max.print = 100)

#where to deposit figures made
figures_home <- "~/Documents/GitHub/FjordPhyto/figures/"
data_home <- "~/Documents/GitHub/FjordPhyto/data/"

## Colors and labels ---
region_colors <- c("red4", "green4", "blue", "purple")
names(region_colors) <- c("southern", "middle", "northern", "shetlands")
region_labels <- c("Southern", "Middle", "Northern", "Shetlands")
names(region_labels) <-c("southern", "middle", "northern", "shetlands")

season_colors <- c("#c64575", "#5fb651", "#9e5cc5", "#c58942")
names(season_colors) <- c("2017-2018", "2018-2019", "2019-2020", "2021-2022")

month_colors <- c("#4bb092", "#ca5740", "#6587cd", "#b4b542", "#c979ad", "#6b7f38")
names(month_colors) <- c("11", "12", "1", "2", "3", "4'")

partition_colors <- c("#E5B2F9", "#D587FA", "#47195C")
names(partition_colors) <- c("A", "B", "C")

cluster_colors <- c("#30C5FF", "#555B6E", "#F46036")
names(cluster_colors) <- c(1:3)

group_colors <- c("#F06400", "#00F064","#008CF0","gold2","#6400F0","#F0008C","grey10")
names(group_colors) <- c("Cryptophytes", "Diatoms", "Dinoflagellates", "Haptophytes",
                         "MAST", "Rhodophytes", "Greenalgae")


# load in data from data_processing_w_cleanmeta_rerunqiime.R
load("~/Documents/GitHub/FjordPhyto/data/WAP_set_meta_rerun.Rdata") #see the README file for how this was re-created
#or 
#source("Data/data_processing_w_cleanmeta_rerunqiime.R")

#### Metadata ----
metadata$year[metadata$year == "2020" & metadata$month == "12" & metadata$day == "19"] <- "2019"
metadata$Date <-  as.Date(paste(metadata$day,metadata$month,metadata$year,sep = '/'),format = '%d/%m/%Y')

#make SEASONS brackets 
metadata$month <- factor(as.character(metadata$month), levels = c("11","12","1","2","3","4"))
#metadata$Date <- mdy(metadata$Date)
metadata$season <- "2017-2018"
metadata$season[which(mdy("06-01-2019") > metadata$Date & metadata$Date > mdy("09-01-2018"))] <- "2018-2019"
metadata$season[which(mdy("06-01-2020") > metadata$Date & metadata$Date > mdy("09-01-2019"))] <- "2019-2020"
metadata$season[which(mdy("06-01-2021") > metadata$Date & metadata$Date > mdy("09-01-2020"))] <- "2020-2021" #word of note on 2020-2021 code, there should not be any seasons since it was PANDEMIC year no travel
metadata$season[which(mdy("06-01-2022") > metadata$Date & metadata$Date > mdy("09-01-2021"))] <- "2021-2022"

#make REGIONS brackets "northern" (anything north of -63) "middle" (anything between -63 and -65) "southern" (anything south of -65)
metadata$region <-"shetlands" #assign everything something, start with northern then parse/ID other labels
#metadata$region[which( latitude less than -63 == northern, latitude 63 - 65 = middle, latitude >65 southern)]
metadata$region[which(metadata$Lat < -63 & metadata$Lat > -66)] <-"middle"
metadata$region[which(metadata$Lat < -66 & metadata$Lat > -73)] <-"southern"
metadata$region[which(metadata$Long < -50 & metadata$Long > -58)] <-"northern"

metadata <-  metadata %>%
  mutate(Site_Name_2 = ifelse(
    grepl("Transect", Site_Name) & grepl("Neko", Site_Name),
    "Neko Harbour",
    Site_Name),
    Site_Name_2 = ifelse(Site_Name_2 == "Mikkelson Harbour", "Mikkelsen Harbour", Site_Name_2),
    Site_Name_2 = ifelse(
      grepl("Ple", Site_Name_2),
      "Pleneau Bay",
      Site_Name_2)) %>%
  drop_na(Long, Lat) %>%
  filter(Site_Name_2 != "")

metadata_raw <- metadata

### Fixing order of site names
fix_names <- read_csv(paste0(data_home, "site_names_edit_manual.csv"), skip = 0)

### Final metadata config
metadata <- metadata %>%
  filter(Site_Name_2 != "Whalers Bay") %>%
  mutate(Site_Name_2 = ifelse(Site_Name_2 == "Andvord Bay Transect Station 1 Useful",
                                            "Useful Island", Site_Name_2),
         Site_Name_2 = ifelse(Site_Name_2 == "Andvord Bay Transect Station 2 Errera",
                              "Errera Channel", Site_Name_2),
         Site_Name_2 = ifelse(Site_Name_2 == "Andvord Bay Transect Station 3 Middle",
                              "Mid Andvord Bay", Site_Name_2),
         Site_Name_2 = ifelse(Site_Name_2 == "Andvord Bay Transect Station 5 Interior",
                              "Bagshawe Glacier", Site_Name_2),
         Site_Name_2 = ifelse(Site_Name_2 == "Horshoe Island",
                              "Horseshoe Island", Site_Name_2),
    Site_Name_2 = factor(Site_Name_2, levels = fix_names$Rename),
    site_ids = factor(Site_Name_2, levels = fix_names$Rename, labels = fix_names$site_num)) %>%
  group_by(Site_Name_2) %>%
  mutate(lon = mean(Long),
         lat = mean(Lat),
         region = factor(region, levels = c("southern", "middle", "northern", "shetlands")))


### Figure 1 - Sampling effort --- 
map <- map_data("world")

sta_map <- ggplot(metadata %>%
                    mutate(region = factor(region, 
                                           levels = c("shetlands","northern", "middle", "southern")))) + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black", linewidth = 0.25) +
  coord_map(projection = "ortho",xlim = c(-70,-55), ylim = c(-68,-61.8),orientation = c(-100,-80,-12.5)) +
  geom_point(aes(x = lon, y = lat), size = 1.5) +
  geom_label_repel(data = metadata %>%
                     mutate(region = factor(region, 
                                            levels = c("shetlands","northern", "middle", "southern"))) %>%
                     group_by(Site_Name_2, site_ids, lon, lat) %>%
                     reframe(region = unique(region)),
                   aes(x = lon, y = lat, label = site_ids, fill = region),
                   color = "white",
                   segment.color="black",
                   min.segment.length = 0, max.overlaps = 500,
                   label.padding = unit(0.1, "lines"),
                   size = 6
  ) +
  labs(x = "Longitude", y = "Latitude", fill = "Region") +
  scale_fill_manual(values = region_colors, labels = region_labels) +
  theme(text = element_text(size = 18),
        axis.text = element_text(size = 15, color = "black"),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.key = element_blank(),
        strip.background = element_blank(),
        legend.position = c(0.8,0.25)) + 
  guides(
    fill = guide_legend(
      override.aes = aes(label = "")
    )
  )

ant_inset <- ggplot(data = map %>% filter(region == "Antarctica")) +
  geom_polygon(aes(x = long,y = lat, group=group),
               color=NA, fill="grey60") +
  geom_rect(aes(xmin = -70, xmax = -55, ymin = -68, ymax = -61.8), color = "red", fill = NA) +
  coord_polar(start = 4*pi/3, clip = "on") +
  theme_map() +
  theme(
    panel.background = element_rect(color = "black", linewidth = 1)
  )

sta_map_inset <- ggdraw(sta_map, clip = "on") +
  draw_plot(ant_inset, x = 0.1, y = 0.7, width = 0.3, height = 0.3)


axis_color <- c(rep("red4",3), rep("green4",21), rep("blue",5), rep("purple", 3))
sta_samples <- metadata %>%
  filter(month != 4) %>%
  group_by(region, site_ids, season, month) %>%
  summarize(
    total_samples = length(unique(samples))) %>%
  ggplot() + 
  geom_tile(aes(x = month, y = site_ids, fill = total_samples)) +
  coord_cartesian(expand = 0) +
  theme_bw() +
  labs(x = "", y = "") +
  scale_x_discrete(labels = c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr")) +
  scale_fill_gradient("Samples", low = "lightblue", high = "blue4") +
  theme(
    axis.text.x = element_text(size = 15, color = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 13, color = axis_color),
    strip.text = element_text(face = "bold", size = 16),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 15)
  ) +
  facet_grid(~season)

station_map_sampling <- plot_grid(sta_map_inset, sta_samples, align = "h", axis = "t",
                                  labels = c("(A)", "(B)"), nrow = 1, rel_widths = c(1,1.5))


ggsave(
  filename = paste0(figures_home, "map_sampling-effort.jpg"),
  station_map_sampling,
  width = 16,
  height = 5.5,
  units = "in",
  dpi = 300
)

### Fix taxa names ----
# split taxa and assign phytogroups - can always go back and reassign phytogroups (e.g., pennates vs centrics, mixotrophs vs autotrophs)
taxa_table_split <- separate(taxa_table, pr2_Taxon, sep = ";",
                             into = c("Kingdom",
                                      "Supergroup",
                                      "Division",
                                      "Class",
                                      "Order",
                                      "Family",
                                      "Genus",
                                      "Species"))

#to remove NAs from the taxa table, try to attach _unknown if its unkown, but i think this messed something up downstream
for (i in 1:length(unique(taxa_table_split$Feature.ID))){
  taxa_vec <- subset(taxa_table_split, Feature.ID == unique(taxa_table_split$Feature.ID)[i])
  taxa_na <- which(is.na(taxa_vec))
  
  if (length(taxa_na) > 0){
    number_of_X <- 1:length(taxa_na)
    first_taxa <- colnames(taxa_table_split)[taxa_na[1]-1]
    if (grepl("_X", taxa_vec[first_taxa])){
      number_of_X <- number_of_X + 1
    }
    last_taxa <- colnames(taxa_table_split)[taxa_na[length(taxa_na)]]
    X_vector <- unlist(lapply(number_of_X, function(r){paste0(rep("X",r), collapse = "")}))
    new_list <- paste(gsub("_X", "", taxa_vec[first_taxa]), X_vector, sep = "_")
    names(new_list) <- colnames(taxa_table_split)[taxa_na]
    
    taxa_vec[names(new_list)] <- new_list
    taxa_table_split[
      taxa_table_split$Feature.ID == unique(taxa_table_split$Feature.ID)[i],
      names(new_list)] <- new_list
    
    if (last_taxa == "Species"){
      taxa_table_split[
        taxa_table_split$Feature.ID == unique(taxa_table_split$Feature.ID)[i],
        "Species"]  <- paste0(taxa_table_split[
          taxa_table_split$Feature.ID == unique(taxa_table_split$Feature.ID)[i],
          "Genus"], "_sp.")
    } 
    
  }
  
}

#### Define phytoplankton groups ----
#create my phytoplankton groups - justify how i did this - is this comprehensive of "phytoplankton" if it also includes rhodophytes? 
taxa_table_split <- taxa_table_split %>%
  mutate(
    phytogroups = case_when(
      Class == "Bacillariophyta" ~ "Diatoms",
      Division == "Dinoflagellata" ~ "Dinoflagellates",
      Division %in% c("Chlorophyta", "Prasinodermophyta") ~ "Greenalgae", #note: i just deleted the space in between, see later code
      Division == "Haptophyta" ~ "Haptophytes",
      Division == "Rhodophyta" ~ "Rhodophytes",
      Division == "Cryptophyta"~ "Cryptophytes",
      grepl("MAST-", Class) ~ "MAST"))




asv_table_wide <- pivot_wider(asv_table,
                              id_cols = "Feature.ID",
                              names_from = "samples",
                              values_from = "reads",
                              values_fill = 0)


asv_table_fix <- pivot_longer(
  asv_table_wide,
  cols = colnames(asv_table_wide)[-1],
  names_to = "samples",
  values_to = "reads") %>%
  group_by(Feature.ID) %>%
  mutate(total_samples = length(unique(samples[reads > 0])),
         n_reads = sum(reads)) 


## ASVs that were present in other samples and run with those samples but not present in WAP
not_present <- asv_table_fix %>% filter(total_samples == 0) %>% pull(Feature.ID) %>% unique()

asv_table_2 <- asv_table_fix %>%
  filter(total_samples > 0) # Can change minimum number of samples present 
asv_table_2 <- asv_table_2[,-ncol(asv_table_2)]
asv_table_raw <- asv_table_2
 ## Config final raw ASV file


##### Rarefy data set ------
asv_table_wide_2 <- pivot_wider(
  asv_table_raw, id_cols = "samples", names_from = "Feature.ID", values_from = "reads")

raremax <- 7000 #can change this cut off to something else if justified 
samples_keep <- asv_table_wide_2[rowSums(asv_table_wide_2[,-1]) >= raremax,1]$samples
samples_removed <- asv_table_wide_2[rowSums(asv_table_wide_2[,-1]) < raremax,1]$samples
metadata <- metadata %>%
  filter(samples %in% samples_keep)

Srare <- rrarefy(asv_table_wide_2[rowSums(asv_table_wide_2[,-1]) >= raremax,-1], raremax)
df_rare <- Srare %>%
  as.data.frame() %>%
  mutate(samples = samples_keep) %>%
  pivot_longer(cols = !samples, names_to = "Feature.ID", values_to = "reads")
colnames(asv_table_2)[3] <- "raw_reads"
asv_table_rare <- left_join(asv_table_2, df_rare, by = c("Feature.ID", "samples"))
## Final rarefied asv table config
asv_rarefy <- asv_table_raw %>%
  left_join(., metadata, by = "samples") %>%
  group_by(region, Feature.ID) %>%
  drop_na(region) %>%
  reframe(sum_reads = sum(reads, na.rm = T)) %>%
  pivot_wider(id_cols = "region", names_from = "Feature.ID", values_from = "sum_reads") %>%
  select(-region) %>%
  rarecurve(., step = 100, sample = raremax)
asv_rarefy <- as.vector(asv_table_raw %>%
                          group_by(Feature.ID) %>%
                          reframe(sum_reads = sum(reads, na.rm = T)) %>%
                          pull(sum_reads))
asv_rarefy[is.na(asv_rarefy)] <- 0

#### Supplementary Figure 1 -----
out <- rarecurve(asv_table_wide_2[rowSums(asv_table_wide_2[,-1]) >= raremax,-1],
          step = 20, sample = raremax)
names(out) <- as.vector(asv_table_wide_2[rowSums(asv_table_wide_2[,-1]) >= raremax,1]$samples)

protox <- mapply(FUN = function(x, y) {
  mydf <- as.data.frame(x)
  colnames(mydf) <- "Species"
  mydf$samples <- y
  mydf$samplesize <- attr(x, "Subsample")
  mydf
}, x = out, y = as.list(names(out)), SIMPLIFY = FALSE)

xy <- do.call(rbind, protox)
rownames(xy) <- NULL 

region_rare <- left_join(xy, metadata, by = "samples") %>%
  drop_na(season, region) %>%
  mutate(region = factor(region, levels = c("shetlands", "northern", "middle", "southern"))) %>%
ggplot() +
  geom_line(aes(x = samplesize, y = Species, color = region, group = samples)) +
  theme_bw() +
  geom_vline(aes(xintercept = raremax)) +
  scale_color_manual(name = "Region",
                     values = region_colors,
                     labels = region_labels) +
  labs(x = "Sample Size", y = "Number of species") +
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = "black"),
        legend.position = c(0.8,0.2),
        legend.box.background = element_rect(color = "black"),
        text = element_text(size = 14, color = "black"))
  
season_rare <- left_join(xy, metadata, by = "samples") %>%
  drop_na(season, region) %>%
  ggplot() +
  geom_line(aes(x = samplesize, y = Species, color = season, group = samples)) +
  theme_bw() +
  geom_vline(aes(xintercept = raremax)) +
  scale_color_manual(name = "Season",
                     values = season_colors) +
  labs(x = "Sample Size", y = "Number of species") +
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = "black"),
        legend.position = c(0.8,0.2),
        legend.box.background = element_rect(color = "black"),
        text = element_text(size = 14, color = "black"))

month_rare <- left_join(xy, metadata, by = "samples") %>%
  drop_na(season, region) %>%
  ggplot() +
  geom_line(aes(x = samplesize, y = Species, color = month, group = samples)) +
  theme_bw() +
  geom_vline(aes(xintercept = raremax)) +
  scale_color_manual(name = "Month",
                     values = month_colors) +
  labs(x = "Sample Size", y = "Number of species") +
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = "black"),
        legend.position = c(0.8,0.2),
        legend.box.background = element_rect(color = "black"),
        text = element_text(size = 14, color = "black"))

rare_curves_plot <- region_rare + season_rare + month_rare &
  theme(legend.position = c(0.8, 0.25), panel.grid = element_blank()) &
  plot_annotation(tag_levels = 'A', tag_suffix = ")", tag_prefix = "(")

ggsave(
  filename = paste0(figures_home, "Fig_S1-rare_curves.jpg"),
  rare_curves_plot,
  width = 13,
  height = 5.5,
  units = "in",
  dpi = 300
)

### Filter out data based on total number of reads ----
asv_table_filter <- asv_table_rare %>%
  group_by(Feature.ID) %>%
  mutate(n_reads = sum(reads,na.rm = T)) %>%
  filter(n_reads > 5)


## Within the Western Antarctica Peninsula,
## across four seasons (2017-2018, 2018-2019, 2019-2020, 2021-2022) 
## 271 number of samples were collected using XX methods.
## samples missing metadata (e.g., GPS etc) were removed
metadata_raw %>% pull(samples) %>% unique() %>% length()
## The full 18S dataset includes XX ASVs and Xx total reads.
asv_table_raw %>% pull(Feature.ID) %>% unique() %>% length() 
asv_table_raw %>% pull(reads) %>% sum()

## To normalize library sizes, we rarefied each sample down to 7000 reads.
## XX samples had less than 7000 total reads and were excluded from further analysis
length(unique(samples_removed))
## leaving XX number of samples.
length(unique(samples_keep))
## For statistical significance, we removed ASVs that had less than 5
## reads over the full rarefied dataset which accounted for about 58% of ASVs 
ordered_reads <- asv_table_rare %>%
  group_by(Feature.ID) %>%
  mutate(n_reads = sum(reads,na.rm = T)) %>%
  arrange(n_reads) %>%
  pull(n_reads)

#So from the start, how many reads have we cut? 
## We categorized the remaining XX ASVs and xxx reads into 7 distinct phytoplankton groups.
#CHECK THAT NUMBER 521, cuz after phytogroup categorization 192 and NAs 332 that adds up to 529
asv_table_filter %>% filter(reads > 0) %>% pull(Feature.ID) %>% unique() %>% length()
## Diatoms include XX asvs totaling XX reads within the Class Baciliariphyta which includes XX Genus and XX species,
## Dinos include XX asvs totaling XX reads within the Class Baciliariphyta which includes XX Genus and XX species,
## Green algae include XX asvs totaling XX reads within the Class Baciliariphyta which includes XX Genus and XX species,
## Haptos include XX asvs totaling XX reads within the Class Baciliariphyta which includes XX Genus and XX species,
## Cryptos include XX asvs totaling XX reads within the Class Baciliariphyta which includes XX Genus and XX species,
## MAST include XX asvs totaling XX reads within the Class Baciliariphyta which includes XX Genus and XX species,
## Rhodos include XX asvs totaling XX reads within the Class Baciliariphyta which includes XX Genus and XX species,
## We removed XX ASVs that are not classified within one of these phytoplankton groups.
left_join(asv_table_filter, taxa_table_split, by = "Feature.ID") %>%
  group_by(phytogroups) %>%
  reframe("# ASVs" = length(unique(Feature.ID[reads > 0 & !is.na(reads)])),
          "# reads" = sum(reads, na.rm = T),
          "# Genera" = length(unique(Genus[reads > 0 & !is.na(reads)])),
          "# Species" = length(unique(Species[reads > 0 & !is.na(reads)]))) %>%
  mutate(phytogroups = if_else(is.na(phytogroups), 'Other', phytogroups),
         phytogroups = if_else(phytogroups == "Greenalgae", "Green algae", phytogroups),
          phytogroups = factor(phytogroups,
                              levels = c("Diatoms", "Dinoflagellates",
                                         "Cryptophytes", "MAST","Haptophytes",
                                         "Rhodophytes", "Green algae", "Other"))) %>%
  arrange(phytogroups) %>%
  mutate(" " = phytogroups) %>%
  select(-phytogroups) %>%
  relocate(" ") %>%
  kbl(caption = "Table 2. Total number of ASVs, reads, Genera, and Species within each phytoplankton group") %>%
  kable_classic(full_width = F, html_font = "Cambria")



#### Define diversity metrics from rarefied reads -----

taxa_table_split_2 <- subset(taxa_table_split,
                             Feature.ID %in% unique(asv_table_filter$Feature.ID))
diversity_group <- c("all",
                     "phytogroups",
                     unique(taxa_table_split_2$phytogroups)[
                       !is.na(unique(taxa_table_split_2$phytogroups))])

for (i in 1:length(diversity_group)){
  
  if (diversity_group[i] == "all"){
    taxa_pull <- taxa_table_split_2
  } else if (diversity_group[i] == "phytogroups"){
    taxa_pull <- taxa_table_split_2 %>% filter(!is.na(phytogroups))
  }else{
    taxa_pull <- taxa_table_split_2 %>% filter(phytogroups == diversity_group[i])
  }
  
  
  piv_all <- asv_table_filter %>%
    filter(Feature.ID %in% taxa_pull$Feature.ID) %>%
    dplyr::select(-c(raw_reads, total_samples, n_reads))
  
  
  #piv_all$Feature.ID <- paste0("X",piv_all$Feature.ID)
  
  piv_all <- piv_all %>% group_by(samples) %>%
    mutate(prop_reads = reads/sum(reads)) %>% ungroup() #all reads in a sample, proportion of each taxa in each sample
  piv_all$reads <- NULL
  piv_all$raw_reads <- NULL
  
  piv_all <- piv_all %>%
    pivot_wider(names_from = "Feature.ID", values_from = "prop_reads", values_fill = 0)
  
  
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
    df_return <- df_output #cbind(df_output_2, piv_all[,-1])
  }else{
    df_return <- left_join(df_return, df_output, by = "samples")
  }
}

### metadata = metadata; asv_table_filter = filtered, rarefied asv table,
### taxa_table_split_2 = taxa table; df_retun = diversity metrics
full_df <- left_join(asv_table_filter, taxa_table_split_2, by ="Feature.ID") %>%
  filter(!is.na(phytogroups)) %>%
  left_join(., metadata, by = c("samples")) %>%
  drop_na(Site_Name_2)

## Fig. 2
full_df %>%
  filter(!is.na(phytogroups)) %>%
  mutate(phyto2 = factor(phytogroups, levels = c(
    "Diatoms", "Dinoflagellates", "Cryptophytes","MAST",
    "Haptophytes", "Rhodophytes", "Greenalgae"))) %>%
  group_by(phyto2, Species) %>%
  reframe(total = sum(reads, na.rm = T)) %>%
  group_by(phyto2) %>%
  arrange(desc(total), .by_group = T) %>%
  pull(Species)

species_order <- full_df %>%
  filter(!is.na(phytogroups)) %>%
  mutate(phyto2 = factor(phytogroups, levels = c(
    "Diatoms", "Dinoflagellates", "Cryptophytes","MAST",
    "Haptophytes", "Rhodophytes", "Greenalgae"))) %>%
  group_by(phyto2, Species) %>%
  reframe(total = sum(reads, na.rm = T)) %>%
  group_by(phyto2) %>%
  arrange(phyto2, desc(total)) %>%
  pull(Species)

species_color <- full_df %>%
  filter(!is.na(phytogroups)) %>%
  mutate(phyto2 = factor(phytogroups, levels = c(
    "Greenalgae", "Rhodophytes", "Haptophytes", "MAST",
    "Cryptophytes", "Dinoflagellates", "Diatoms"))) %>%
  group_by(phyto2, Species) %>%
  reframe(total = sum(reads, na.rm = T)) %>%
  group_by(phyto2) %>%
  arrange(phyto2, total) %>%
  mutate(col = case_when(
    phyto2 == "Diatoms" ~ "#00F064",
    phyto2 == "Dinoflagellates" ~ "#008CF0",
    phyto2 == "Greenalgae" ~ "grey10",
    phyto2 == "MAST" ~ "#6400F0",
    phyto2 == "Cryptophytes" ~ "#F06400",
    phyto2 ==  "Haptophytes" ~ "gold2", 
    phyto2 == "Rhodophytes" ~ "#F0008C"))  %>%
  mutate(Species = gsub("__", "_", gsub("X","", Species)))

### Species order
full_df %>%
  group_by(site_ids) %>%
  mutate(max = sum(reads, na.rm = T)) %>%
  group_by(phytogroups, Species, region, site_ids) %>%
  reframe(rel = log10(sum(reads, na.rm = T)/unique(max)[1] * 100)) %>%
  mutate(Species =  gsub("__", "_", gsub("X","", factor(Species, levels = species_color$Species))))
fig_2 <- full_df %>%
  group_by(site_ids) %>%
  mutate(max = sum(reads, na.rm = T)) %>%
  group_by(phytogroups, Species, region, site_ids) %>%
  reframe(rel = log10(sum(reads, na.rm = T)/unique(max)[1]*100)) %>%
  mutate(Species =  factor(gsub("__", "_", gsub("X","", Species)), levels = species_color$Species)) %>%
  ggplot() +
  geom_tile(aes(x = site_ids, y = Species, fill = rel))+
  scale_x_discrete(name = "", position = "top") +
  scale_fill_distiller("Relative\nabundance",
                       palette = "RdYlBu",
                      na.value = "white",
                      limits = c(-1, 2),
                      breaks = c(-1,0,1,2),
                      labels = c("<0.1%", "1%", "10%", "100%"),
                      oob=scales::squish) +
  facet_wrap(~region, scales = "free_x", ncol = 4,
             labeller = labeller(region = c(
               "southern" = "Southern",
               "middle" = "Middle",
               "northern" = "Northern",
               "shetlands" = "Shetlands"
             ))) +
  force_panelsizes(cols = c(1,5,1.5,1)) +
  coord_cartesian(expand = 0) +
  theme(strip.placement = "outside",
        axis.text.y = element_text(color = species_color$col, size = 8),
        axis.text.x = element_text(color = "black", size = 10),
        strip.background = element_blank(),
        strip.text = element_text(size = 11, face = "bold"))


ggsave(
  filename = paste0(figures_home, "heatmap_species-region_relabun.jpg"),
  fig_2,
  width = 9.5,
  height = 10,
  units = "in",
  dpi = 300
)

fig_3 <- full_df %>%
  filter(month != 4) %>%
  group_by(season, month) %>%
  mutate(max = sum(reads, na.rm = T)) %>%
  group_by(phytogroups, Species, season, month) %>%
  reframe(rel = log10(sum(reads, na.rm = T)/unique(max)[1] * 100)) %>%
  mutate(Species =  factor(gsub("__", "_", gsub("X","", Species)), levels = species_color$Species)) %>%
  ggplot() +
  geom_tile(aes(x = month, y = Species, fill = rel))+
  scale_x_discrete(name = "", position = "top") +
  scale_fill_distiller("Relative\nabundance",
                       palette = "RdYlBu",
                       na.value = "white",
                       limits = c(-1, 2),
                       breaks = c(-1,0,1,2),
                       labels = c("<0.1%", "1%", "10%", "100%"),
                       oob=scales::squish) +
  facet_wrap(~season, scales = "free_x", ncol = 4) +
  force_panelsizes(cols = c(1,1,1,1)) +
  coord_cartesian(expand = 0) +
  theme(strip.placement = "outside",
        axis.text.y = element_text(color = species_color$col, size = 8),
        axis.text.x = element_text(color = "black", size = 10),
        strip.background = element_blank(),
        strip.text = element_text(size = 11, face = "bold"))

ggsave(
  filename = paste0(figures_home, "heatmap_species-season-month_relabun.jpg"),
  fig_3,
  width = 9.5,
  height = 9,
  units = "in",
  dpi = 300
)

full_df %>%
  group_by(month) %>%
  mutate(max = sum(reads, na.rm = T)) %>%
  group_by(phytogroups, Species, month) %>%
  reframe(rel = log10(sum(reads, na.rm = T)/unique(max)[1] * 100)) %>%
  mutate(Species = factor(Species, levels = species_color$Species)) %>%
  ggplot() +
  geom_tile(aes(x = month, y = Species, fill = rel))+
  scale_x_discrete(name = "", position = "top") +
  scale_fill_gradient("Log10\n% Rel.\nabundance",
                      low = "lightblue1", high = "purple4",
                      na.value = "white",
                      limits=c(-3, 2),
                      labels = ) +
  coord_cartesian(expand = 0) +
  theme(strip.placement = "outside",
        axis.text.y = element_text(color = species_color$col, size = 7),
        axis.text.x = element_text(color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size = 11, face = "bold"))

# Fig. 3
richness <- specnumber(piv_all[,-1])
shannon <- diversity(piv_all[,-1], MARGIN = 1, index = "shannon")
evenness <- shannon/log(richness)
simpson <- diversity(piv_all[,-1], MARGIN = 1, index = "simpson")
plot_list <- left_join(full_df, df_return, by = "samples") %>%
  filter(month != 4) %>%
  group_by(Site_Name_2, lat, lon) %>%
  reframe(gamma = max(richness_phytogroups),
          alpha = mean(richness_phytogroups),
          evenness = mean(evenness_phytogroups),
          shannon = mean(shannon_phytogroups)) %>%
  pivot_longer(c(gamma, alpha, evenness, shannon)) %>%
  mutate(name = factor(name,
                       levels = c("gamma", "alpha", "evenness", "shannon"),
                       labels = c("Total species richness (# ASVs)",
                                  "Mean species richness (# ASVs)",
                                  "Mean species evenness",
                                  "Mean Shannon diversity index"))) %>%
  group_by(name, .add = T) %>%
  group_split() %>% 
  map(
    ~{
      gg <- ggplot(., aes(x = lon, y = lat, color = value)) + 
        geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black", linewidth = 0.25) +
        #coord_map(projection = "ortho",xlim = c(-70,-55), ylim = c(-68,-61.8),orientation = c(-100,-80,-12.5)) +
        geom_point(alpha = 1, size = 3) +
        scale_color_viridis_c() +
        labs(x = "Longitude", y = "Latitude", color = .x$name) +
        theme(text = element_text(size = 15),
              panel.grid = element_line(color = "grey"),
              panel.background = element_blank(),
              panel.border = element_rect(fill = NA, color = "black"),
              legend.key = element_blank(),
              strip.background = element_blank(),
              legend.position = "bottom") +
        guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = 20))
      
      from <- c(xmin = -64.2, xmax = -62, ymin = -66, ymax = -64)
      to <- c(-60, -55, -68.5, -65)
      
      gg_out <- gg +
        geom_magnify(from = from, to = to, shape = "ellipse",
                     inset.linetype = 0, proj.combine = F, proj = "corresponding") +
        coord_map(projection = "ortho",
                  xlim = c(-70,-55),
                  ylim = c(-68.5,-61.8), orientation = c(-100,-80,-12.5))
      return(gg_out)
    }) %>%
  plot_grid(plotlist = .,
            align = 'hv',
            ncol = 2,
            labels = c("(A)", "(B)", "(C)", "(D)"))
plot_list

ggsave(
  filename = paste0(figures_home, "diversity-map.jpg"),
  plot_list,
  width = 9.8,
  height = 7,
  units = "in",
  bg = "white",
  dpi = 300
)

left_join(full_df, df_return, by = "samples") %>%
  filter(month != 4) %>%
  group_by(Site_Name_2, lat, lon) %>%
  reframe(gamma = max(richness_phytogroups),
          alpha = mean(richness_phytogroups),
          evenness = mean(evenness_phytogroups),
          shannon = mean(shannon_phytogroups)) %>%
  pivot_longer(c(gamma, alpha, evenness, shannon)) %>%
  mutate(name = factor(name,
                       levels = c("gamma", "alpha", "evenness", "shannon"),
                       labels = c("Total species richness (# ASVs)",
                                  "Mean species richness (# ASVs)",
                                  "Mean species evenness",
                                  "Mean Shannon diversity index"))) %>%
  group_by(name, .add = T) %>%
  group_split() %>%
  map(~{
    .x %>%
      arrange(value)
  })


# Fig. 4 
region_div_avg <- left_join(full_df, df_return, by = "samples") %>%
  filter(month != 4) %>%
  pivot_longer(c(richness_phytogroups, evenness_phytogroups, shannon_phytogroups)) %>%
  mutate(name = factor(name,
                       levels = c("richness_phytogroups", "evenness_phytogroups",
                                  "shannon_phytogroups")),
         region = factor(region,
                         levels = c("shetlands", "northern", "middle", "southern"))) %>%
group_by(name, region, month) %>%
    reframe(mean = mean(value,na.rm = T))


stat_test <-  left_join(full_df, df_return, by = "samples") %>%
  filter(month != 4) %>%
  pivot_longer(c(richness_phytogroups, evenness_phytogroups, shannon_phytogroups)) %>%
  mutate(name = factor(name,
                       levels = c("richness_phytogroups", "evenness_phytogroups",
                                  "shannon_phytogroups")),
         region = factor(region,
                         levels = c("shetlands", "northern", "middle", "southern"))) %>%
  group_by(name, month) %>%
  tukey_hsd(value ~ region) %>%
  add_xy_position(x = "month", group = "region",scales = "free_y")

stat_test_season <-  left_join(full_df, df_return, by = "samples") %>%
  filter(month != 4) %>%
  pivot_longer(c(richness_phytogroups, evenness_phytogroups, shannon_phytogroups)) %>%
  mutate(name = factor(name,
                       levels = c("richness_phytogroups", "evenness_phytogroups",
                                  "shannon_phytogroups")),
         region = factor(region,
                         levels = c("shetlands", "northern", "middle", "southern"))) %>%
  group_by(name, region) %>%
  tukey_hsd(value ~ month)



fig_4 <- left_join(full_df, df_return, by = "samples") %>%
  filter(month != 4) %>%
  pivot_longer(c(richness_phytogroups, evenness_phytogroups, shannon_phytogroups)) %>%
  mutate(name = factor(name,
                       levels = c("richness_phytogroups", "evenness_phytogroups",
                                  "shannon_phytogroups")),
         region = factor(region,
                         levels = c("shetlands", "northern", "middle", "southern"))) %>%
  ggplot() +
  geom_boxplot(aes(x = month, y = value, fill = region), position = position_dodge2(preserve = "single")) +
  stat_pvalue_manual(stat_test, tip.length = 0, step.increase = 0) +
  facet_wrap(name~., scales = "free_y", ncol = 1, switch = "y",
             labeller = labeller(name = c(
               "richness_phytogroups" = "Species richness",
               "evenness_phytogroups" = "Species evenness",
               "shannon_phytogroups" = "Shannon diversity index"))) +
  scale_fill_manual(values = region_colors, labels = region_labels) +
  theme(
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    panel.grid = element_blank(),
    axis.text = element_text(color = "black", size = 10),
    strip.placement = "outside",
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_blank()) +
  labs(x = "Month", y = "", fill = "Region")

ggsave(
  paste0(figures_home, "fig_4-signifiance.jpg"),
  fig_4,
  width = 7,
  height = 6,
  units = "in",
  dpi = 300
)


## Fig. 5
phyto_asvs <- unique(taxa_table_split_2$Feature.ID[!is.na(taxa_table_split_2$phytogroups)])
april_sample <- metadata$samples[metadata$month == 4]

subset_bio <- subset(asv_table_filter,
                     samples %in% samples_keep &
                       samples != april_sample &
                       Feature.ID %in% phyto_asvs)

taxa_ids <- taxa_table_split_2 %>%
  select(Species,Feature.ID) %>% distinct()

taxa_names_cluster <- setNames(as.character(taxa_ids$Feature.ID), taxa_ids$Species)

subset_bio$Species <- factor(subset_bio$Feature.ID, levels = taxa_names_cluster, labels = names(taxa_names_cluster))

bio <- pivot_wider(
  subset_bio[colnames(subset_bio)[colnames(subset_bio) != 'raw_reads']],
  id_cols = "samples", names_from = "Species", values_from = "reads", values_fn = function(r)sum(r,na.rm = T)) %>%
  column_to_rownames(var = "samples")

bio_remove <-  bio[,apply(bio, 2, function(r){!all(r == 0)})]
#colnames(bio_remove) <- gsub("__", "_", gsub("X","", colnames(bio_remove)))

#bio_distmat1 <- vegdist(t(bio_remove), binary=FALSE, method = "bray")

spearman.cor <- bio_remove  %>% 
  cor(use="pairwise.complete.obs", method = "spearman")

spearman.cor <- as.dist(1 - spearman.cor)

bio_NMS1 <-  metaMDS(spearman.cor,
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

hc1 <- hclust(spearman.cor, method = "ward.D2")

plot(hc1, cex = 0.6, hang = -1)
rect.hclust(hc1, k = 3, border = 2:10)
groups <- cutree(hc1, k = 3)

nmds_clust <- ggplot(data_scores_1) + 
  geom_point(aes(x = MDS1, y = MDS2,
                 color = factor(
                   Species, levels = names(groups), labels = groups))) +
  stat_ellipse(aes(x = MDS1, y = MDS2,
                   color = factor(
                     Species, levels = names(groups), labels = groups))) +
  theme_bw() + 
  scale_color_manual(name = "Cluster",
                     values = cluster_colors) +
  theme(legend.position = c(0.23,0.94),
        legend.direction = "horizontal",
        legend.background = element_rect(color = "black", fill = "white"),
        axis.text = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 12),
        panel.background = element_rect(color = "black"),
        panel.grid = element_blank())
nmds_clust


myannotation <- as.data.frame(cutree(hc1, k = 3))
names(myannotation)[1] = "Cluster" 
myannotation$Cluster <- factor(myannotation$Cluster, levels= 1:3, 
                               labels=1:3)

cluster_df <- as.data.frame(myannotation)
cluster_df$Species <- rownames(cluster_df)

sample_cluster_df <- subset_bio %>%
  left_join(., metadata, by = "samples", keep = F) %>%
  left_join(.,taxa_table_split_2, by = c("Feature.ID", "Species")) %>%
  left_join(., df_return , by = "samples") %>%
  drop_na(Site_Name_2, phytogroups) %>%
  left_join(., cluster_df, by = "Species") %>%
  filter(reads > 0)


ann_colors <- list(Cluster = cluster_colors)

fig_5 <- pheatmap(bio_remove  %>% 
                    cor(use="pairwise.complete.obs", method = "spearman"),
                  main = "Species Co-occurrence Matrix (Spearman's Correlation)",
                  fontsize_row = 6, fontsize_col = 6,clustering_method = "ward.D2",
                  cutree_rows = 3, cutree_cols = 3,
                  legend_breaks =seq(-1, 1, length.out=5),
                  annotation_col=myannotation, annotation_colors=ann_colors)

ggsave(
  filename = paste0(figures_home, "spearman-rank.jpg"),
  fig_5,
  width = 10.1,
  height = 9.1,
  units = "in",
  dpi = 300)





#### Rethinking partitions -- related them back to samples -----

### Days since "start of the season" -> Nov. 1
days_since <- yday(sample_cluster_df$Date) - yday("2000-11-01")
days_since[days_since < 0] <- days_since[days_since < 0] + 365

sample_cluster_df$days_since <- days_since
#"Porosira_sp.","Geminigera_cryophila","Dino-Group-I-Clade-1_X_sp."
day_since_df <- sample_cluster_df %>%
  group_by(samples) %>%
  mutate(total_reads = sum(reads)) %>%
  group_by(samples, Cluster) %>%
  mutate(prop_reads = sum(reads)/total_reads) %>%
  #filter(Species %in% c("Phaeocystis_sp.", "Porosira_sp.","Geminigera_cryophila","Dino-Group-I-Clade-1_X_sp.")) %>% 
  group_by(Cluster, Species, season, region) %>%
  reframe(sd = sqrt(sum(prop_reads*(days_since - days_since[which.max(prop_reads)])^2)/
                      (((length(prop_reads>0)-1)/length(prop_reads))*sum(prop_reads))),
          day_max = days_since[which.max(prop_reads)]) %>%
  filter(region == "middle")

days_order <- day_since_df %>%
  filter(season == "2017-2018")  %>%
  arrange(desc(day_max)) %>% pull(Species)
day_since_df %>%
  mutate(Species = factor(Species, levels = c(unique(day_since_df$Species)[!unique(day_since_df$Species) %in% days_order], unique(days_order)))) %>%
ggplot() +
  geom_point(aes(y = Species, x = day_max, color = season)) +
  geom_errorbarh(aes(y = Species, xmin = day_max-sd, xmax = day_max+sd, color = season)) +
  facet_wrap(Cluster~region) +
  scale_color_manual(name = "", values = c("black", "blue", "green4", "red")) +
  theme(axis.text.y = element_text(size = 6))
  

sample_cluster_df %>%
  group_by(phytogroups, season, days_since) %>%
  reframe(mean_reads = mean(reads,na.rm = T),
          sd_reads = sd(reads,na.rm = T)) %>%
  ggplot() +
  geom_point(aes(x = days_since, y = mean_reads, color = phytogroups)) +
  geom_smooth(aes(x = days_since, y = mean_reads, color = phytogroups, group = phytogroups), method = "gam", se = F) +
  facet_wrap(~season) +
  scale_y_log10() +
  scale_color_manual(values = group_colors)

sample_cluster_df %>%
  group_by(phytogroups, season, days_since) %>%
  reframe(mean_days = mean(days_since,na.rm = T)) %>%
  ggplot() +
  geom_point(aes(x = days_since, y = mean_days, color = phytogroups)) +
  geom_smooth(aes(x = days_since, y = mean_days, color = phytogroups, group = phytogroups), method = "gam", se = F) +
  facet_wrap(~season) +
  scale_y_log10() +
  scale_color_manual(values = group_colors)


sample_cluster_df %>%
  filter(phytogroups == "Diatoms") %>%
  mutate(succesion_period = case_when())
  ggplot() +
  geom_point(aes(x = days_since, y = reads)) +
  geom_smooth(aes(x = days_since, y = reads, color = Genus, group = Genus), method = "loess", show.legend = F, se = F) 


### Do the cluster/ groups be the same? idk why not?
left_join(partitions_df, cluster_groups, by = "Species")
cluster_groups <- as.data.frame(groups)
cluster_groups$Species <- rownames(cluster_groups)


# prop_species <- sample_cluster_df %>%
#   filter(month != 4) %>% 
#   group_by(clusters) %>%
#   mutate(total_reads = sum(reads)) %>%
#   group_by(clusters, Species) %>%
#   reframe(species_reads = sum(reads),
#           prop_reads = sum(reads)/total_reads) %>%
#   distinct() %>%
#   filter(species_reads > 0) %>%
#   group_by(clusters) %>%
#   slice_min(prop_reads, n = 10) %>%
#   ggplot() + 
#   geom_bar(aes(x = clusters, y = prop_reads, fill = Species), stat = "identity") +
#   theme_bw() +
#   labs(x = "Cluster", y = "Proportion of reads") +
#   coord_cartesian(expand = 0)



### Phyloseq
library(phyloseq)
phyto_asvs <- unique(taxa_table_split_2$Feature.ID[!is.na(taxa_table_split_2$phytogroups)])
april_sample <- metadata$samples[metadata$month == 4]

subset_bio <- subset(asv_table_filter,
                     samples %in% samples_keep &
                       samples != april_sample &
                       Feature.ID %in% phyto_asvs)

asvs_wide <- subset_bio %>%
  dplyr::select(Feature.ID, samples, reads) %>%
  pivot_wider(., values_from = reads, id_cols = Feature.ID, names_from = samples)
asvs_wide_df <- as.data.frame(asvs_wide[,-1]) #gets rid of Feature.ID column
asvs_wide_df[is.na(asvs_wide_df)] <- 0 #if there are NA reads, make 0 (this may have been in my processing stage when i was merging the sequencing runs - MUST CHECK WHY)

rownames(asvs_wide_df) <- asvs_wide$Feature.ID #make FeatureID the rownames
asvs_wide_df[is.na(asvs_wide_df)]

in_biom <- otu_table(asvs_wide_df, taxa_are_rows = T)
taxa_names(in_biom) <- asvs_wide$Feature.ID
#Make sample_data
# in_biom_metadata <- sample_data(left_join(data_scores_1, metadata, by = "samples") %>%
#                                   drop_na(month) %>%
#                                   mutate(clusters = factor(
#                                     samples, levels = names(groups), labels = groups)))
# in_biom_metadata <- sample_data(metadata)
sample_names(in_biom_metadata) <- in_biom_metadata$samples


#Make tax_table
taxatable <- left_join(taxa_table_split_2, cluster_df, by = "Species") %>% drop_na(Cluster)
in_biom_tax <- tax_table(taxatable[,c(2:7, 11:12, 8:9)])
taxa_names(in_biom_tax) <- taxatable$Feature.ID
colnames(in_biom_tax) <- colnames(taxatable[,c(2:7, 11:12, 8:9)])

phylo_fjordphyto <- merge_phyloseq(in_biom, in_biom_tax, in_biom_metadata)
top_nested <- nested_top_taxa(phylo_fjordphyto,
                              top_tax_level = "Clusters",
                              nested_tax_level = "Species",
                              n_top_taxa = 7, 
                              n_nested_taxa = 5)


top_asv <- top_taxa(phylo_fjordphyto, n_taxa = 10, grouping = "clusters",
                    by_proportion = T, tax_level = "Species", FUN = sum)
ps_obj <- top_asv$ps_obj
pal <- taxon_colours(ps_obj,
                     tax_level = "phytogroups",
                     merged_label = "Other", 
                     merged_clr = "grey90",
                     palette = NULL,
                     base_clr = "#008CF0")
ps_tmp <- ps_obj %>% name_na_taxa(include_rank = T, 
                                  na_label = "<tax> (<rank>)")
#View(tax_table(psdf))

# ps_tmp <- ps_tmp %>% label_duplicate_taxa(tax_level = "Species",
#                                           asv_as_id = F,
#                                           duplicate_label = "<tax> <id>")
psdf <- psmelt(ps_tmp)
psdf <- left_join(psdf, top_asv$top_taxa) %>%
  mutate(
    phytogroups = case_when(is.na(tax_rank) ~ "Other", TRUE ~ phytogroups),
    Genus = case_when(is.na(tax_rank) ~ "Other", TRUE ~ Genus),
    Species = case_when(is.na(tax_rank) ~ "Other", TRUE ~ Species))
psdf <- move_label(psdf = psdf,
                   col_name = "phytogroups",
                   label = "Other", 
                   pos = 0)
psdf <- move_nested_labels(psdf,
                           top_level = "phytogroups", 
                           nested_level = "Species",
                           top_merged_label = "Other", 
                           nested_label = gsub("<tax>", "", "Other <tax>"), 
                           pos = Inf)
psdf_rel <- top_asv$top_taxa %>%
  group_by(clusters) %>%
  mutate(Abundance = abundance / sum(abundance)) %>%
  ungroup() %>%
  mutate(Species = gsub("__", "_", gsub("X","", Species)))


top_10_taxa <- ggnested(psdf_rel, aes_string(main_group = "phytogroups",
                                            sub_group = "Species", 
                               x = "clusters", y = "Abundance"),
         main_palette = group_colors) + 
  scale_y_continuous(name = "Normalized relative abundance",
                     expand = c(0,0)) +
  scale_x_discrete(name = "Clusters", expand = c(1E-3, 1E-3)) +
  theme_nested(theme_light) + 
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(size = 12)) +
  geom_col(position = position_fill())

fig_6 <- nmds_clust + top_10_taxa +
  plot_annotation(tag_levels = 'A', tag_suffix = ")", tag_prefix = "(") &
  theme(axis.text = element_text(size = 10, color = "black")) 
fig_6  


ggsave(filename = paste0(figures_home, "fig_6-2.jpg"),
       fig_6,
       width = 15,
       height = 5,
       units = "in",
       dpi = 300)


### Figure 5

fig_5 <- ggplot(sample_cluster_df)  +
  geom_point(aes(x = Date, y = Site_Name_2, color = clusters)) +
  facet_grid(~season, scales = "free_x") +
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 12, color = axis_color),
    strip.text = element_text(face = "bold", size = 14),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white", color = "black"),
    panel.grid = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  ) +
  scale_x_date(name = "") +
  scale_y_discrete(name = "") +
  scale_color_manual(name = "Cluster", values = cluster_colors)

ggsave(
  filename = paste0(figures_home, "cluster-space-time.png"),
  fig_5,
  width = 8.5,
  height = 6,
  units = "in")




## Fig. 7
region_div_avg <- sample_cluster_df %>%
  pivot_longer(c(richness_phytogroups, evenness_phytogroups, shannon_phytogroups)) %>%
  mutate(name = factor(name,
                       levels = c("richness_phytogroups", "evenness_phytogroups",
                                  "shannon_phytogroups"),
                       labels = c("Species richness",
                                  "Species evenness",
                                  "Shannon diversity index"))) %>%
  pivot_longer(c(month, region, season), names_to = "facet", values_to = "facet_values") %>%
  mutate(facet_values = factor(facet_values,
                               levels = c(
                                 "southern", "middle", "northern", "shetlands",
                                 "2017-2018",
                                 "2018-2019",
                                 "2019-2020",
                                 "2021-2022",
                                 "11",
                                 "12",
                                 "1",
                                 "2",
                                 "3"),
                               labels = c("Southern",
                                          "Middle",
                                          "Northern",
                                          "Shetlands",
                                          "2017-2018",
                                          "2018-2019",
                                          "2019-2020",
                                          "2021-2022",
                                          "11",
                                          "12",
                                          "1",
                                          "2",
                                          "3")),
         facet = factor(facet, levels = c("region", "season", "month"),
                        labels = c("Region", "Season", "Month"))) %>%
  group_by(name, facet, facet_values, clusters) %>%
  reframe(mean = mean(value,na.rm = T)) %>% View()


stat_test <-  sample_cluster_df %>%
  pivot_longer(c(richness_phytogroups, evenness_phytogroups, shannon_phytogroups)) %>%
  mutate(name = factor(name,
                       levels = c("richness_phytogroups", "evenness_phytogroups",
                                  "shannon_phytogroups"),
                       labels = c("Species richness",
                                  "Species evenness",
                                  "Shannon diversity index"))) %>%
  pivot_longer(c(month, region, season), names_to = "facet", values_to = "facet_values") %>%
  mutate(facet_values = factor(facet_values,
                               levels = c(
                                 "southern", "middle", "northern", "shetlands",
                                 "2017-2018",
                                 "2018-2019",
                                 "2019-2020",
                                 "2021-2022",
                                 "11",
                                 "12",
                                 "1",
                                 "2",
                                 "3"),
                               labels = c("Southern",
                                          "Middle",
                                          "Northern",
                                          "Shetlands",
                                          "2017-2018",
                                          "2018-2019",
                                          "2019-2020",
                                          "2021-2022",
                                          "11",
                                          "12",
                                          "1",
                                          "2",
                                          "3")),
         facet = factor(facet, levels = c("region", "season", "month"),
                        labels = c("Region", "Season", "Month"))) %>%
  filter(facet_values != "Shetlands") %>%
  group_by(name, facet_values) %>%
  tukey_hsd(value ~ clusters) %>%
  add_xy_position(x = "facet_values", group = "name", scales = "free_y")
new_row <- stat_test[7,]
new_row$facet_values <- "Shetlands"

stat_test <- rbind(stat_test, new_row)
stat_test_season <-  left_join(full_df, df_return, by = "samples") %>%
  filter(month != 4) %>%
  pivot_longer(c(richness_phytogroups, evenness_phytogroups, shannon_phytogroups)) %>%
  mutate(name = factor(name,
                       levels = c("richness_phytogroups", "evenness_phytogroups",
                                  "shannon_phytogroups")),
         region = factor(region,
                         levels = c("shetlands", "northern", "middle", "southern"))) %>%
  group_by(name, region) %>%
  tukey_hsd(value ~ month)
  

fig_7 <- sample_cluster_df %>%
  pivot_longer(c(richness_phytogroups, evenness_phytogroups, shannon_phytogroups)) %>%
  mutate(name = factor(name,
                       levels = c("richness_phytogroups", "evenness_phytogroups",
                                  "shannon_phytogroups"),
                       labels = c("Species richness",
                                  "Species evenness",
                                  "Shannon diversity index"))) %>%
  pivot_longer(c(month, region, season), names_to = "facet", values_to = "facet_values") %>%
  mutate(facet_values = factor(facet_values,
                         levels = c(
                           "southern", "middle", "northern", "shetlands",
                           "2017-2018",
                           "2018-2019",
                           "2019-2020",
                           "2021-2022",
                           "11",
                           "12",
                           "1",
                           "2",
                           "3"),
                         labels = c("Southern",
                                    "Middle",
                                    "Northern",
                                    "Shetlands",
                                    "2017-2018",
                                    "2018-2019",
                                    "2019-2020",
                                    "2021-2022",
                                    "11",
                                    "12",
                                    "1",
                                    "2",
                                    "3")),
  facet = factor(facet, levels = c("region", "season", "month"),
                        labels = c("Region", "Season", "Month"))) %>%
  ggplot() +
  stat_boxplot(aes(x = facet_values, y = value, fill = clusters),
               position = position_dodge2(preserve = "single")) +
  facet_grid(name~facet, scales = "free", switch = "both") +
  scale_fill_manual(name = "Cluster", values = cluster_colors) +
  theme(
    legend.text = element_text(color = "black", size = 12),
    legend.title = element_text(color = "black", size = 13),
    axis.text = element_text(color = "black", size = 12),
    strip.placement = "outside",
    strip.text = element_text(face = "bold", size = 13),
    strip.background = element_blank(),
    panel.background = element_rect(color = "black", fill = "white"),
    panel.grid = element_blank()) +
  labs(x = "", y = "")

ggsave(
  filename = paste0(figures_home, "cluster-diversity-all.svg"),
  fig_7,
  width = 13,
  height = 7.5,
  units = "in",
  dpi = 300)




## Supplemental Figures
library(ggpubr)

sample_cluster_df %>%
  filter(phytogroups == "Cryptophytes") %>%
  group_by(samples) %>%
  reframe(reads = sum(reads,na.rm = T),
          shan = mean(shannon_Cryptophytes))
  distinct(samples, shannon_Cryptophytes)

fig_SY <- sample_cluster_df %>%
  pivot_longer(c(shannon_Cryptophytes,
                 shannon_Diatoms,
                 shannon_MAST,
                 shannon_Dinoflagellates,
                 shannon_Greenalgae,
                 shannon_Haptophytes,
                 shannon_Rhodophytes)) %>%
  mutate(name = factor(name, levels = c("shannon_Diatoms",
                                        "shannon_Dinoflagellates",
                                        "shannon_Cryptophytes",
                                        "shannon_MAST",
                                        "shannon_Haptophytes",
                                        "shannon_Rhodophytes",
                                        "shannon_Greenalgae"),
                       labels = c("Diatoms", "Dinoflagellates", "Cryptophytes",
                                  "MAST", "Haptophytes", "Rhodophytes", "Green algae"))) %>%
  ggplot(data = ., aes(x = value, y = shannon_phytogroups)) +
  geom_point(aes(x = value, y = shannon_phytogroups)) +
  scale_y_continuous(name = "Shannon diversity of all phytogplankton", limits = c(0, 4.5)) +
  scale_x_continuous(name = "Shannon diversity of specific phytoplankton groups") +
  stat_smooth(aes(x = value, y = shannon_phytogroups), method = "glm") +
  stat_cor(label.x.npc = "left", label.y.npc = "top", method = "pearson", label.sep = ",",
           aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~','~"))) +
  facet_wrap(~name) +
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 11),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 11))

ggsave(
  filename = paste0(figures_home, "fig_SY.png"),
  fig_SY,
  width = 8,
  height = 8,
  units = "in")



## comparison to other datasets
load("~/Desktop/18sv9_tara_polar.Rdata")
# polar_dat; polar_tax
microscope <- readxl::read_xlsx("~/Desktop/TaxaList_FjordPhyto_Microscopy_Mascioni.xlsx")
microscope <- microscope %>%
  mutate(species = if_else(is.na(specificEpithet), NA,
                           paste0(genus,"_",specificEpithet))) %>%
  relocate(species, .before = "specificEpithet") %>%
  as.data.frame()

for (i in 1:nrow(microscope)){
  taxa_vec <- microscope[i,3:9]
  taxa_na <- which(is.na(taxa_vec))
    
  if (length(taxa_na) > 0){
    number_of_X <- 1:length(taxa_na)
    first_taxa <- colnames(taxa_vec)[taxa_na[1]-1]
    if (grepl("_X", taxa_vec[first_taxa])){
      number_of_X <- number_of_X + 1
    }
    last_taxa <- colnames(taxa_vec)[taxa_na[length(taxa_na)]]
    X_vector <- unlist(lapply(number_of_X, function(r){paste0(rep("X",r), collapse = "")}))
    new_list <- paste(gsub("_X", "", taxa_vec[first_taxa]), X_vector, sep = "_")
    names(new_list) <- colnames(taxa_vec)[taxa_na]
    
    taxa_vec[names(new_list)] <- new_list
    microscope[i, names(new_list)] <- as.vector(new_list)
    
    if (last_taxa == "species"){
      microscope[i,"species"]  <- paste0(microscope[i,"genus"], "_sp.")
    } 
  }
}
polar_tax <- separate(polar_tax, Taxon, sep = ";",
                             into = c("Kingdom",
                                      "Supergroup",
                                      "Division",
                                      "Class",
                                      "Order",
                                      "Family",
                                      "Genus",
                                      "Species")) %>%
  filter(Confidence >= 0.97)

for (i in 1:length(unique(polar_tax$Feature.ID))){
  taxa_vec <- subset(polar_tax, Feature.ID == unique(polar_tax$Feature.ID)[i])
  taxa_na <- which(is.na(taxa_vec))
  
  if (length(taxa_na) > 0){
    number_of_X <- 1:length(taxa_na)
    first_taxa <- colnames(polar_tax)[taxa_na[1]-1]
    if (grepl("_X", taxa_vec[first_taxa])){
      number_of_X <- number_of_X + 1
    }
    last_taxa <- colnames(polar_tax)[taxa_na[length(taxa_na)]]
    X_vector <- unlist(lapply(number_of_X, function(r){paste0(rep("X",r), collapse = "")}))
    new_list <- paste(gsub("_X", "", taxa_vec[first_taxa]), X_vector, sep = "_")
    names(new_list) <- colnames(polar_tax)[taxa_na]
    
    taxa_vec[names(new_list)] <- new_list
    polar_tax[
      polar_tax$Feature.ID == unique(polar_tax$Feature.ID)[i],
      names(new_list)] <- new_list
    
    if (last_taxa == "Species"){
      polar_tax[
        polar_tax$Feature.ID == unique(polar_tax$Feature.ID)[i],
        "Species"]  <- paste0(polar_tax[
          polar_tax$Feature.ID == unique(polar_tax$Feature.ID)[i],
          "Genus"], "_sp.")
    } 
  }
}
polar_tax_fix <- polar_tax


polar_tax_fix <- polar_tax_fix %>% mutate(
  phytogroups = case_when(
    Class == "Bacillariophyta" ~ "Diatoms",
    Division == "Dinoflagellata" ~ "Dinoflagellates",
    Division %in% c("Chlorophyta", "Prasinodermophyta") ~ "Greenalgae", #note: i just deleted the space in between, see later code
    Division == "Haptophyta" ~ "Haptophytes",
    Division == "Rhodophyta" ~ "Rhodophytes",
    Division == "Cryptophyta"~ "Cryptophytes",
    grepl("MAST-", Class) ~ "MAST")) %>%
  filter(!is.na(phytogroups))

tara_meta_1 <- read_delim("~/Downloads/PRJEB6610_protist_amplicon.txt", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)
tara_meta_1$accession <- sub('.*_', '', str_sub(tara_meta_1$sample_title, start= -12))
accession_want <- sub(".*_", "", c("TARA_N000001042",
                    "TARA_N000001034",
                    "TARA_N000001026",
                    "TARA_N000001030",
                    "TARA_N000001040",
                    "TARA_N000001032",
                    "TARA_N000001024",
                    "TARA_N000001028",
                    "TARA_N000001036",
                    "TARA_N000001362",
                    "TARA_N000001440",
                    "TARA_N000001006",
                    "TARA_N000001438",
                    "TARA_N000001442"))

run_accessions <- tara_meta_1 %>%
  filter(accession %in% accession_want) %>%
  pull(run_accession)

polar_dat %>%
  rownames_to_column(var = "run_accession") %>%
  pull(run_accession)

x <- list("Tara Polar" = unique(polar_tax_fix$Species),
          "FjordPhyto" = unique(taxa_table_split_2$Species),
          "Microscopy" = unique(microscope$species))

plyr::ldply(x, rbind)
all_dfs <- data.frame(tara = sort(x[[1]]))
all_dfs <- left_join(all_dfs, sort(c(x[[2]], rep(NA,length(x[[1]]) - length(x[[2]])))))
           fjordphyto = sort(c(x[[2]], rep(NA,length(x[[1]]) - length(x[[2]])))),
           microscopy = arrange(c(x[[3]], rep(NA,length(x[[1]]) - length(x[[3]])))))
all_dfs %>%
  arrange(tara)

ggVennDiagram(x)

#compare species taxa list between Tara Polar, FjordPhyto and Microscopy Mascioni
all_dfs <- data.frame(names = sort(unique(unlist(x))))
write.csv(all_dfs %>%
  mutate(tara = all_dfs$name %in% x[[1]],
         fjordphyto = all_dfs$name %in% x[[2]],
         microscopy = all_dfs$name %in% x[[3]]) %>%
    arrange(desc(microscopy)),
  file = "~/Desktop/tara_fjord_microscopy-compare.csv")



## asv read 'common' rare dominant etc. 

