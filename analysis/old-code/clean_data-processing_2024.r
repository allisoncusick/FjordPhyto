
require(vegan)
require(pheatmap)
require(R.utils)
require(data.table)
##### Data Processing -----


# Load the data
#### Rarefy data set ------
asv_table_wide_2 <- full_asv %>%
  select(Feature.ID, sample, reads) %>%
pivot_wider(., id_cols = "sample", names_from = "Feature.ID", values_from = "reads") %>%
  replace(is.na(.), 0)

raremax <- 7000 #can change this cut off to something else if justified 
## Keep samples (rows) where the sum of all the reads is greater than raremax
samples_keep <- asv_table_wide_2[rowSums(asv_table_wide_2[,-1]) >= raremax,1]$sample
## Remove samples (rows) where the sum of all the reads is less than raremax
samples_removed <- asv_table_wide_2[rowSums(asv_table_wide_2[,-1],) < raremax,1]$sample

## Rarefy using vegan package
Srare <- vegan::rrarefy(asv_table_wide_2[rowSums(asv_table_wide_2[,-1]) >= raremax,-1], raremax)
df_rare <- Srare %>%
  as.data.frame() %>%
  mutate(sample = samples_keep) %>%
  pivot_longer(cols = !sample,
               names_to = "Feature.ID", values_to = "reads")%>%
  replace(is.na(.), 0)

colnames(full_asv)[colnames(full_asv) == "reads"] <- "raw_reads"

not_present <- df_rare %>%
  group_by(Feature.ID) %>%
  reframe(total_reads = sum(reads)) %>%
  filter(total_reads < 5) %>%
  pull(Feature.ID) %>%
  unique()

asv_table_rare <- left_join(full_asv, df_rare,
                            by = c("Feature.ID", "sample")) %>%
  filter(!Feature.ID %in% not_present & !sample %in% samples_removed)
  

diversity_group <- c("phytogroups", unique(taxa_table_split$phytogroups))
diversity_group <- diversity_group[!is.na(diversity_group)]

for (i in 1:length(diversity_group)){
  
  if (diversity_group[i] == "phytogroups"){
    taxa_pull <- taxa_table_split
  }else{
    taxa_pull <- taxa_table_split %>% filter(phytogroups == diversity_group[i])
  }
  
  piv_all <- asv_table_rare %>%
    filter(Feature.ID %in% taxa_pull$Feature.ID) %>%
    dplyr::select(c(Feature.ID, sample, reads))
  
  
  # piv_all <- piv_all %>% group_by(sample) %>%
  #   mutate(prop_reads = reads/sum(reads)) %>% ungroup() #all reads in a sample, proportion of each taxa in each sample
  #piv_all$reads <- NULL
  #piv_all$raw_reads <- NULL
  
  piv_all <- piv_all %>%
    pivot_wider(names_from = "Feature.ID", values_from = "reads", values_fn = mean, values_fill = 0)
  
  
  richness <- specnumber(piv_all[,-1])
  shannon <- diversity(piv_all[,-1], MARGIN = 1, index = "shannon")
  evenness <- shannon/log(richness)
  simpson <- diversity(piv_all[,-1], MARGIN = 1, index = "simpson")
  rel_abun <- rowSums(piv_all[,-1])/raremax

  shannon_name <- paste0("shannon_", diversity_group[i])
  even_name <- paste0("evenness_", diversity_group[i])
  simpson_name <- paste0("simpson_", diversity_group[i])
  richness_name <- paste0("richness_", diversity_group[i])
  rel_abun_name <- paste0("rel_abun_", diversity_group[i])
  
  df_output <- data.frame(sample = piv_all[,1],
                          shannon = shannon,
                          evenness = evenness,
                          simpson = simpson,
                          richness = richness,
                          rel_abun = rel_abun)
  colnames(df_output)[2:ncol(df_output)] <- str_replace(c(shannon_name, even_name, simpson_name, richness_name, rel_abun_name), " ", "_")
  
  if (i == 1){
    df_return <- df_output 
  }else{
    df_return <- left_join(df_return, df_output, by = "sample")
  }
}






require("ggnewscale")
###Evidence for the relationship between warmer+fresherwater = crypots
metadata_div %>%
  filter(Project_Name == "WAP" &
           Genetics_18sv9_Sample_ID %in% samples_keep &
           Surface_Salinity_ > 30) %>%
  ggplot(aes(
    x = as.numeric(Surface_Temperature_),
    y = as.numeric(Surface_Salinity_),
    color = as.numeric(chlor_a))) +
  geom_point() +
  labs(
    x = "Surface Temperature (C)",
    y = "Surface Salinity (PSU)",
    color = "Chlorophyll-a (mg/m^3)"
  )

metadata_div %>%
  filter(Project_Name == "WAP" &
           Genetics_18sv9_Sample_ID %in% samples_keep &
           Surface_Salinity_ > 30) %>%
  ggplot(aes(
    x = as.numeric(Surface_Temperature_),
    y = as.numeric(Surface_Salinity_),
    color = as.numeric(richness_Diatoms))) +
  geom_point() +
  labs(
    x = "Surface Temperature (C)",
    y = "Surface Salinity (PSU)",
    color = "Richness (Diatoms)"
  )
  
metadata_div %>%
  filter(Project_Name == "WAP" &
           Genetics_18sv9_Sample_ID %in% samples_keep &
           Surface_Salinity_ > 30) %>%
  ggplot(aes(
    x = as.numeric(Surface_Temperature_),
    y = as.numeric(Surface_Salinity_),
    color = as.numeric(MLD_1))) +
  geom_point() +
  labs(
    x = "Surface Temperature (C)",
    y = "Surface Salinity (PSU)",
    color = "Richness (Diatoms)"
  )

metadata_div %>%
  filter(Project_Name == "WAP" &
           Genetics_18sv9_Sample_ID %in% samples_keep &
           Surface_Salinity_ > 30) %>%
  ggplot(aes(
    x = as.numeric(Surface_Temperature_),
    y = as.numeric(Surface_Salinity_),
    color = as.numeric(shannon_Diatoms))) +
  geom_point() +
  labs(
    x = "Surface Temperature (C)",
    y = "Surface Salinity (PSU)",
    color = "Shannon (Diatoms)"
  )



#bar or lineplot

metadata_div %>%
  filter(Project_Name == "WAP" &
           Genetics_18sv9_Sample_ID %in% samples_keep &
           Surface_Salinity_ > 30) %>%
  pivot_longer(cols = c("rel_abun_Cryptophytes", "rel_abun_Dinoflagellates", "rel_abun_Diatoms", "rel_abun_Haptophytes", "rel_abun_Greenalgae", "rel_abun_MAST"),
               names_to = "group", values_to = "rel_abun") %>%
  ggplot(aes(x = Surface_Temperature_,
             #y = rel_abun,
             fill = group)) +
  geom_bar()

metadata_div %>%
  filter(Project_Name == "WAP" &
           Genetics_18sv9_Sample_ID %in% samples_keep &
           Surface_Salinity_ > 30) %>%
  ggplot() +
  geom_point(aes(x = longitude, y = latitude, color = month))

metadata_div %>%
  filter(month == "01" & year == "2022")


#### ASV Clustering thing ----
bio <- asv_table_rare %>%
  filter(!is.na(phytogroups)) %>%
  select(Species, rare_reads, sample) %>%
  filter(sample %in% unique(metadata_div$Genetics_18sv9_Sample_ID)) %>%
  replace_na() %>%
  pivot_wider(.,
              id_cols = "sample", names_from = "Species", values_from = "rare_reads", values_fn = function(r)sum(r,na.rm = T), values_fill = 0) %>%
  column_to_rownames(var = "sample")
bio <- bio[,-(which(colSums(bio)==0))]

bio_distmat1 <- vegdist(bio, binary=FALSE, method = "bray")

correlation_matrix <- cor(t(bio), method = "spearman")
hierarchical_result <- hclust(bio_distmat1, method = "ward.D2")
plot(hierarchical_result, cex = 0.6, hang = -1)
rect.hclust(hierarchical_result, k = 3, border = 2:10)
groups <- cutree(hierarchical_result, k = 3)

bio_NMS1 <- metaMDS(bio_distmat1,
                    distance = "bray",
                    k = 2,
                    maxit = 999, 
                    trymax = 500,
                    wascores= TRUE, expand = TRUE, autotransform = FALSE)

bio.spp_scrs1 <-   sppscores(bio_NMS1) <- bio
#Para agregar las especies significativas al dataframe
bio.spp.fit <- envfit(bio_NMS1, bio, permutations = 999) # this fits species vectors
spp.scrs <- as.data.frame(scores(bio.spp.fit, display = "vectors", arrow.mul=2.5)) #save species intrinsic values into dataframe #arrow.mul es un vector que multiplica el largo de las flechas
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs)) #add species names to dataframe
spp.scrs <- cbind(spp.scrs, pval = bio.spp.fit$vectors$pvals) #add pvalues to dataframe so you can select species which are significant
spp.scrs<- cbind(spp.scrs, abrev = abbreviate(spp.scrs$Species, minlength = 6)) #abbreviate species names
sig.spp.scrs <- subset(spp.scrs, pval<=0.05) #subset data to show species significant at 0.05

subset_meta <- metadata_div %>%
  select(-grep("_all", colnames(.))) %>%
  select(Genetics_18sv9_Sample_ID, days_since, month, seasonyear, region, Surface_Salinity_, Surface_Temperature_, chlor_a, MLD_5, Meltwater_Fraction_, secchi_depth) %>%
  rename(sample = "Genetics_18sv9_Sample_ID") %>%
  filter(sample %in% rownames(bio))


env_fit <- envfit(bio_NMS1, subset_meta, permutations = 999, na.rm = T)
env.scrs <- as.data.frame(scores(env_fit, display = "vectors", arrow.mul=2.5))
env.scrs <- cbind(env.scrs, variable = rownames(env.scrs)) 
env.scrs <- cbind(env.scrs, pval = env_fit$vectors$pvals)

data.scores1 <- as.data.frame(scores(bio_NMS1$points)) 
data.scores1 <- data.frame(data.scores1, type = groups[rownames(data.scores1)])
ggplot() +
  geom_point(data = data.scores1, aes(x = MDS1, y = MDS2, color = factor(type))) +
geom_segment(data = sig.spp.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = sig.spp.scrs, aes(x=NMDS1, y=NMDS2, label = Species), cex = 3, direction = "both")+ #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
  labs(title = "Ordination with phytogroup vectors")

ggplot() +
  geom_point(data = data.scores1, aes(x = MDS1, y = MDS2, color = factor(type))) +
  geom_segment(data = env.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = env.scrs, aes(x=NMDS1, y=NMDS2, label = variable), cex = 3, direction = "both")+ #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
  labs(title = "Ordination with environmental vectors")


data_scores_1 <- as.data.frame(scores(bio_NMS1$points)) %>%
  rownames_to_column(var = "sample") %>%
  mutate(clust = factor(sample, levels = names(groups), labels = groups)) %>%
  left_join(subset_meta, by = "sample") %>%
  filter(!is.na(month)) %>%
  filter(Surface_Salinity_ > 30)

ggplot() +
  geom_point(data = data_scores_1,
             aes(x = days_since,
                 y = Meltwater_Fraction_,
                 color = clust))


ano <- anosim(bio_distmat1, data.scores_2$days_since, distance = "bray", permutations = 9999)
ano <- anosim(bio_distmat1, data.scores_2$type, distance = "bray", permutations = 9999)

plot(ano)

data_scores_1 <- as.data.frame(scores(bio_NMS1$points)) %>%
  rownames_to_column(var = "sample") %>%
  mutate(clust = factor(sample, levels = names(groups), labels = groups)) %>%
  left_join(subset_meta, by = "sample") %>%
  filter(!is.na(month)) %>%
  filter(Surface_Salinity_ > 30)



data.scores_2 <- left_join(data.scores1 %>% rownames_to_column("sample"), subset_meta)

correlation_matrix <- cor(bio, method = "spearman")
hierarchical_result <- hclust(dist(1 - correlation_matrix), method = "ward.D2")

plot(hierarchical_result, cex = 0.6, hang = -1)
rect.hclust(hierarchical_result, k = 3, border = 2:10)
groups <- cutree(hierarchical_result, k = 3)

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

subset_meta <- metadata_div %>%
  select(Genetics_18sv9_Sample_ID, month, seasonyear, region, latitude, Surface_Salinity_, Surface_Temperature_, chlor_a, MLD_5, Meltwater_Fraction_, secchi_depth, days_since) %>%
  rename(sample = "Genetics_18sv9_Sample_ID") %>%
  filter(sample %in% rownames(bio))

#Get the vectors for bioenv.fit
bio_fit <- envfit(bio_NMS1, bio, permutations = 999)

env_fit <- envfit(bio_NMS1, subset_meta, permutations = 999, na.rm = T)

df_biofit<-scores(bio_fit,display=c("vectors"))
df_biofit<-df_biofit*vegan:::ordiArrowMul(df_biofit)
df_biofit<-as.data.frame(df_biofit)

#Get the vectors for env.fit
df_envfit<-scores(env_fit,display=c("vectors"))
df_envfit<-df_envfit*vegan:::ordiArrowMul(df_envfit)
df_envfit<-as.data.frame(df_envfit)

df<-scores(bio_NMS1,display=c("sites")) 
df <- data.frame(df, type = groups[rownames(df)])

ggplot() +
  geom_point(data = df, aes(x = NMDS1, y = NMDS2, color = factor(type))) +
  geom_segment(data=df_envfit, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")),color="#4C005C",alpha=0.5) +
  geom_text(data=as.data.frame(df_envfit*1.1),aes(NMDS1, NMDS2, label = rownames(df_envfit)),color="#808080",alpha=0.5)


data_scores_1 <- as.data.frame(scores(bio_NMS1$points)) %>%
  rownames_to_column(var = "sample") %>%
  mutate(clust = factor(sample, levels = names(groups), labels = groups)) %>%
  left_join(subset_meta, by = "sample") %>%
  filter(!is.na(month)) %>%
  filter(Surface_Salinity_ > 30)




# Pairwise correlation between samples (columns)
cols.cor <- cor(bio, use = "pairwise.complete.obs", method = "pearson")
# Pairwise correlation between rows (genes)
rows.cor <- cor(t(bio), use = "pairwise.complete.obs", method = "pearson")

pheatmap::pheatmap(
  bio, scale = "row", 
  clustering_distance_cols = as.dist(1 - cols.cor),
  clustering_distance_rows = as.dist(1 - rows.cor),
  clustering_method = "ward.D2",
  fontsize = 5
)


correlation_matrix <- cor(t(bio), method = "spearman")
hierarchical_result <- hclust(dist(1 - correlation_matrix), method = "ward.D2")

plot(hierarchical_result, cex = 0.6, hang = -1)
rect.hclust(hierarchical_result, k = 3, border = 2:10)
groups <- cutree(hierarchical_result, k = 3)

myannotation <- as.data.frame(groups)
names(myannotation)[1] = "Cluster" 
myannotation$Cluster <- factor(myannotation$Cluster, levels= 1:3, 
                               labels=1:3)

pheatmap(correlation_matrix,
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


ordi<-ordisurf(bio_NMS1, subset_meta$chlor_a,
               plot = FALSE)
ordi.grid <- ordi$grid #extracts the ordisurf object
str(ordi.grid) #it's a list though - cannot be plotted as is
ordi.mite <- expand.grid(x = ordi.grid$x, y = ordi.grid$y) #get x and ys
ordi.mite$z <- as.vector(ordi.grid$z) #unravel the matrix for the z scores
ordi.mite.na <- data.frame(na.omit(ordi.mite)) #gets rid of the nas
NMDS=data.frame(x=bio_NMS1$point[,1],y=bio_NMS1$point[,2],Type=groups)

ggplot()+
  stat_contour(data = ordi.mite.na, aes(x = x, y = y, z = z, color = ..level..), stat = "identity") +
  geom_point(data=NMDS,aes(x,y,fill=factor(Type)),pch=21,size=3) +
  scale_colour_continuous(high = "darkgreen", low = "darkolivegreen1") +
  labs(color = "Temperature") +
  theme_bw()


ggplot(data_scores_1) +
  geom_point(aes(x = MDS1, y = MDS2, fill = shannon_Diatoms),
             pch = 21, color = "transparent", size = 3) +
  stat_ellipse(aes(x = MDS1, y = MDS2,
                   color = clust)) +
  theme_bw() + 
  theme(
        legend.background = element_blank(),
        axis.text = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 12),
        panel.background = element_rect(color = "black"),
        panel.grid = element_blank()) +
  labs(color = "Sample\nCluster", fill = "Surface\nTemp.")



subset_meta <- metadata_div %>%
  select(Genetics_18sv9_Sample_ID, month, seasonyear, region, latitude, Surface_Salinity_, Surface_Temperature_, chlor_a, MLD_5, Meltwater_Fraction_, secchi_depth, days_since, grep("shannon_", colnames(.)), grep("evenness_", colnames(.)),grep("richness_", colnames(.))) %>%
  rename(sample = "Genetics_18sv9_Sample_ID") %>%
  filter(sample %in% rownames(bio))


data_scores_1 <- as.data.frame(scores(bio_NMS1$points)) %>%
  rownames_to_column(var = "sample") %>%
  mutate(clust = factor(sample, levels = names(groups), labels = groups)) %>%
  left_join(subset_meta, by = "sample") %>%
  filter(!is.na(month)) %>%
filter(Surface_Salinity_ > 30)

nmds_clust <- ggplot(data_scores_1) + 
  geom_point(aes(x = MDS1, y = MDS2,
                 color = clust)) +
  stat_ellipse(aes(x = MDS1, y = MDS2,
                   color = clust)) +
  theme_bw() + 
  theme(legend.position = c(0.23,0.94),
        legend.direction = "horizontal",
        legend.background = element_rect(color = "black", fill = "white"),
        axis.text = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 12),
        panel.background = element_rect(color = "black"),
        panel.grid = element_blank()) +
  scale_color_discrete(name = "Sample\nCluster")
nmds_clust

bio_spp.fit <- envfit(bio_NMS1, bio, permutations = 999) # this fits species vectors
subset_fit <- data_scores_1[complete.cases(data_scores_1),]

bio.env.fit <- envfit(bio_NMS1, subset_fit[,4:ncol(subset_fit)], permutations = 999, na.rm = T) 
spp.scrs <- as.data.frame(scores(bio.spp.fit, display = "vectors", arrow.mul=2.5)) 
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs)) #add species names to dataframe
spp.scrs <- cbind(spp.scrs, pval = bio.spp.fit$vectors$pvals) #add pvalues to dataframe so you can select species which are significant
#spp.scrs<- cbind(spp.scrs, abrev = abbreviate(spp.scrs$Species, minlength = 6)) #abbreviate species names
sig.spp.scrs <- subset(spp.scrs, pval<=0.05) #subset data to show species 

bio_distmat1 <- vegdist(bio, binary=FALSE, method = "bray")


#Get the vectors for bioenv.fit
df_biofit<-scores(,display=c("vectors"))
df_biofit<-df_biofit*vegan:::ordiArrowMul(df_biofit)
df_biofit<-as.data.frame(df_biofit)

#Get the vectors for env.fit
df_envfit<-scores(env.fit,display=c("vectors"))
df_envfit<-df_envfit*vegan:::ordiArrowMul(df_envfit)
df_envfit<-as.data.frame(df_envfit)



testing_meta <- subset_meta %>%
  filter(sample%in%rownames(bio)) %>%
  pull(days_since)

ano <- anosim(bio_distmat1, testing_meta, distance = "bray", permutations = 9999)
plot(ano)


### Generic exploratory figs

ggplot(data_scores_1) + 
  geom_point(aes(x = days_since, y = richness_phytogroups, color = clust))

ggplot(data_scores_1) + 
  geom_point(aes(x = days_since, y = richness_Diatoms, color = clust))

ggplot(data_scores_1) + 
  geom_point(
    aes(x = Surface_Salinity_, y = Surface_Temperature_,
        color = clust, size = Meltwater_Fraction_))

ggplot(data_scores_1) +
  geom_point(aes(x = days_since, y = MLD_5, color = clust)) +
  facet_wrap(~seasonyear)

ggplot(data_scores_1) +
  geom_point(aes(x = days_since, y = Meltwater_Fraction_, color = clust)) +
  facet_wrap(~seasonyear)

ggplot(data_scores_1) +
  geom_point(aes(x = days_since, y = Meltwater_Fraction_)) +
  geom_smooth(aes(x = days_since, y = Meltwater_Fraction_), method = "glm") +
  facet_wrap(~seasonyear) +
  theme_bw() + 
  theme(legend.position = c(0.23,0.94),
        legend.direction = "horizontal",
        legend.background = element_rect(color = "black", fill = "white"),
        axis.text = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 12),
        panel.background = element_rect(color = "black"),
        panel.grid = element_blank()) +
  scale_color_discrete(name = "Sample\nCluster") 
  labs(x = "Days since Oct. 1", y = "Meltwater Frac.")

ggplot(data_scores_1) +
  geom_point(aes(x = days_since, y = Surface_Salinity_)) +
  geom_smooth(aes(x = days_since, y = Surface_Salinity_), method = "glm") +
  facet_wrap(~seasonyear) +
  theme_bw() +
  my_theme +
  labs(x = "Days since Oct. 1", y = "Surface Salinity")

generic_facs <- c("Surface_Temperature_", "Surface_Salinity_", "chlor_a", "MLD_5", "Meltwater_Fraction_", "secchi_depth")

lapply(generic_facs, function(x){
  ggplot(data_scores_1) +
    geom_point(aes(x = days_since, y = !!sym(x))) +
    geom_smooth(aes(x = days_since, y = !!sym(x)), method = "glm") +
    facet_wrap(~seasonyear) +
    theme_bw() +
    my_theme +
    labs(x = "Days since start of season", y = gsub("_", " ", x))
  
  saveas <- paste0("figures/exploratory_figs/", x,"days-since_plot.png")
  ggsave(filename = saveas, plot = last_plot(), width = 7, height = 5, units = "in", dpi = 300)
})

lapply(generic_facs, function(x){
  ggplot(data_scores_1) +
    geom_point(aes(x = days_since, y = !!sym(x), color = clust)) +
    geom_smooth(aes(x = days_since, y = !!sym(x)),
                  group = 1, method = "glm") +
    facet_wrap(~seasonyear) +
    theme_bw() + 
    theme(
          legend.background = element_rect(color = "black", fill = "white"),
          axis.text = element_text(size = 11, color = "black"),
          axis.title = element_text(size = 12),
          panel.background = element_rect(color = "black"),
          panel.grid = element_blank()) +
    scale_color_discrete(name = "Sample\nCluster") +
    labs(x = "Days since start of season", y = gsub("_", " ", x))
  
  saveas <- paste0("figures/exploratory_figs/", x,"days-since_plot-clusters.png")
  ggsave(filename = saveas, plot = last_plot(), width = 7, height = 5, units = "in", dpi = 300)
})


div_facs <- c("shannon_phytogroups", "evenness_phytogroups", "richness_phytogroups", "shannon_Diatoms", "evenness_Diatoms", "richness_Diatoms", "shannon_Cryptophytes", "evenness_Cryptophytes", "richness_Cryptophytes", "shannon_Dinoflagellates", "evenness_Dinoflagellates", "richness_Dinoflagellates", "shannon_Haptophytes", "evenness_Haptophytes", "richness_Haptophytes", "shannon_Greenalgae", "evenness_Greenalgae", "richness_Greenalgae", "shannon_MAST", "evenness_MAST", "richness_MAST")

lapply(div_facs, function(x) {
  data_scores_1 %>%
  pivot_longer(cols = c("Surface_Temperature_", "Surface_Salinity_", "chlor_a", "MLD_5", "Meltwater_Fraction_", "secchi_depth"),
               names_to = "variable", values_to = "value") %>%
ggplot() +
  geom_point(aes(x = value, y = !!sym(x),
                 color = clust)) +
  facet_wrap(~variable, scales = "free_x",
             labeller = labeller(variable = c(
               Surface_Temperature_ = "Surface Temperature",
               Surface_Salinity_ = "Surface Salinity",
               chlor_a = "Chlorophyll-a",
               MLD_5 = "MLD (5 m)",
               Meltwater_Fraction_ = "Meltwater Fraction",
               secchi_depth = "Secchi Depth"
             )), switch = "x") +
  scale_color_discrete(name = "Sample\nCluster") +
  theme_bw() +
  my_theme +
  theme(strip.placement = "outside") +
  labs(y = gsub("_", " ", x), x = "")
  
  saveas <- paste0("figures/exploratory_figs/", x,"_envi-vars.png")
  ggsave(filename = saveas, plot = last_plot(), width = 8, height = 5, units = "in", dpi = 300)
})


data_scores_1 %>%
  pivot_longer(cols = c("Surface_Temperature_", "Surface_Salinity_", "chlor_a", "MLD_5", "Meltwater_Fraction_", "secchi_depth"),
               names_to = "variable", values_to = "value") %>%
ggplot() +
  geom_histogram(aes(x = value, fill = clust)) +
  facet_wrap(~variable, scales = "free_x",
             labeller = labeller(variable = c(
               Surface_Temperature_ = "Surface Temperature",
               Surface_Salinity_ = "Surface Salinity",
               chlor_a = "Chlorophyll-a",
               MLD_5 = "MLD (5 m)",
               Meltwater_Fraction_ = "Meltwater Fraction",
               secchi_depth = "Secchi Depth"
             )), switch = "x") +
  scale_fill_discrete(name = "Sample\nCluster") +
  theme_bw() +
  my_theme +
  theme(strip.placement = "outside") +
  labs(y = "# samples", x = "")

  ggsave(plot = last_plot(),
    filename = "figures/exploratory_figs/envi-vars_hist-clusters.png",
    width = 9, height = 6, units = "in", dpi = 300)

  

lapply(c("richness_", "evenness_", "shannon_"), function(x){
  data_scores_1 %>%
  pivot_longer(cols = grep(x, colnames(.)),
               names_to = "variable", values_to = "value") %>%
  ggplot() +
  geom_histogram(aes(x = value, fill = clust)) +
  facet_wrap(~variable, scales = "free", switch = "x",
             labeller = as_labeller(function(v){
               paste0(capitalize(gsub("_", " (", v)), ")")})) +
  scale_fill_discrete(name = "Sample\nCluster") +
  theme_bw() +
  my_theme +
  theme(strip.placement = "outside") +
  labs(y = "# samples", x = "")
  
  saveas <- paste0("figures/exploratory_figs/", x,"hist-clusters.png")
  ggsave(filename = saveas, plot = last_plot(), width = 9, height = 6, units = "in", dpi = 300)
})


data_scores_1 %>%
  pivot_longer(cols = grep("richness_", colnames(.)),
               names_to = "variable", values_to = "value") %>%
  ggplot() +
  geom_histogram(aes(x = value, fill = clust)) +
  facet_wrap(~variable, scales = "free", switch = "x") +
  scale_color_discrete(name = "Sample\nCluster") +
  theme_bw() +
  my_theme +
  theme(strip.placement = "outside") +
  labs(y = "# samples", x = "")


data_scores_1 %>%
  pivot_longer(cols = grep("richness_", colnames(.)),
               names_to = "variable", values_to = "value") %>%
  ggplot() +
  geom_point(aes(x = Surface_Salinity_, y = value, color = variable)) +
  geom_smooth(aes(x = Surface_Salinity_, y = value, color = variable),
              method = "loess")


correlation_matrix <- cor(data_scores_1[,8:ncol(data_scores_1)], method = "spearman")


data_scores_1 %>%
  pivot_longer(cols = grep("shannon_", colnames(.)),
               names_to = "variable", values_to = "value") %>%
  filter(variable %in% c("shannon_Diatoms", "shannon_Dinoflagellates")) %>%
  ggplot() +
  geom_point(aes(x = Surface_Salinity_, y = Surface_Temperature_,
                 color = variable, size = value)) +
  scale_color_manual(values = c("hotpink","black"))

# Heatmap showing the significant Pearson’s correlations (P≤0.05) between environmental factors and relative amplicon frequency of phytoplankton taxa. Pairwise correlations without an assigned colored dot represent correlations that are not significant.

subset_meta <- metadata_div %>%
  select(Genetics_18sv9_Sample_ID, month, latitude, Surface_Salinity_, Surface_Temperature_, chlor_a, MLD_5, Meltwater_Fraction_, secchi_depth, days_since) %>%
  rename(sample = "Genetics_18sv9_Sample_ID")


bio2 <- asv_table_rare %>%
  filter(!is.na(phytogroups)) %>%
  select(phytogroups, rare_reads, sample) %>%
  filter(sample %in% unique(metadata_div$Genetics_18sv9_Sample_ID)) %>%
  pivot_wider(.,
              id_cols = "sample", names_from = "phytogroups", values_from = "rare_reads", values_fn = function(r)sum(r,na.rm = T), values_fill = 0) %>%
  left_join(subset_meta) %>%
  column_to_rownames(var = "sample") %>%
  mutate_all(as.numeric) %>%
  drop_na()


correlation_matrix <- cor(bio2, method = "spearman", use = "pairwise.complete.obs")

pheatmap::pheatmap(
  correlation_matrix)

