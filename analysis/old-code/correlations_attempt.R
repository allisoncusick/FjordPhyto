
bio <- asv_table_rare %>%
  filter(!is.na(phytogroups)) %>%
  select(Species, rare_reads, sample) %>%
  filter(sample %in% unique(metadata_div$Genetics_18sv9_Sample_ID)) %>%
  replace_na() %>%
  pivot_wider(.,
              id_cols = "sample", names_from = "Species", values_from = "rare_reads", values_fn = function(r)sum(r,na.rm = T), values_fill = 0) %>%
  column_to_rownames(var = "sample")
bio <- bio[,-(which(colSums(bio)==0))]

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
data_scores_cluster <- as.data.frame(scores(bio_NMS1)) %>%
  rownames_to_column("Species") %>%
  mutate(groups = factor(Species, levels = names(groups), labels = groups))

ggplot(data_scores_cluster) +
  geom_point(aes(x=NMDS1, y = NMDS2, color = groups))


bio_distmat1 <- vegdist(bio, binary=FALSE, method = "bray")

bio_NMS1 <- metaMDS(bio_distmat1,
                    distance = "bray",
                    k = 2,
                    maxit = 999, 
                    trymax = 500,
                    wascores= TRUE, expand = TRUE, autotransform = FALSE)

subset_meta <- metadata_div %>%
  select(Genetics_18sv9_Sample_ID, month, seasonyear, region, latitude, longitude, Surface_Salinity_, Surface_Temperature_, chlor_a, MLD_5, Meltwater_Fraction_, secchi_depth, days_since, Density__Derived) %>%
  rename(sample = "Genetics_18sv9_Sample_ID") %>%
  filter(sample %in% rownames(bio)) %>%
  group_by(sample) %>%
  mutate(den_surf = unlist(Density__Derived)[1]) 

cluster_sample <- asv_table_rare %>%
  filter(!is.na(phytogroups)) %>%
  select(Species, rare_reads, sample) %>%
  mutate(clust = factor(Species, levels = names(groups), labels = groups)) %>%
  drop_na(clust) %>%
  filter(Species != "Porosira_sp.") %>%
  group_by(sample) %>%
  mutate(total_reads = sum(rare_reads)) %>%
  group_by(sample, clust) %>%
  reframe(prop_reads = sum(rare_reads)/unique(total_reads))


data_scores_1 <- as.data.frame(scores(bio_NMS1$points)) %>%
  rownames_to_column(var = "sample") %>%
  left_join(subset_meta, by = "sample") %>%
  filter(!is.na(month)) %>%
  filter(Surface_Salinity_ > 30) %>%
  left_join(cluster_sample) %>%
  left_join(metadata_div %>%
              select(Genetics_18sv9_Sample_ID, all_of(div_facs)),
            by = c("sample" = "Genetics_18sv9_Sample_ID"))


ggplot(data_scores_1) +
  geom_point(aes(x = Surface_Salinity_, y = Surface_Temperature_, color = prop_reads, size = MLD_5)) +
  facet_wrap(~clust, labeller = labeller(
    clust = c("1" = "Mixed diatoms", "2" = "Porosira dom.", "3" = "Mixed cryptos")
  )) +
  scale_color_gradientn(colours = c("purple", "orange")) +
  labs(x = "Surface Salinity (PSU)", y = "Surface Temperature (C)",
       color = "Prop. reads", size = "MLD (m)",
       title = "Without Porosira") +
  theme_bw() +
  my_theme

ggplot(data_scores_1) +
  geom_point(aes(x = days_since, y = Meltwater_Fraction_, color = prop_reads, size = MLD_5)) +
  facet_wrap(~clust, labeller = labeller(
    clust = c("1" = "Mixed diatoms", "2" = "Porosira dom.", "3" = "Mixed cryptos")
  )) +
  scale_color_gradientn(colours = c("purple", "orange")) +
  labs(x = "Days since start of season", y = "Melwater Fraction",
       color = "Prop. reads", size = "MLD (m)",
       title = "With Porosira") +
  theme_bw() +
  my_theme

ggsave(filename = "figures/exploratory_figs/cluster_time-wPorosira.png", plot = last_plot(), width = 8, height = 5, units = "in", dpi = 300)

install.packages("ggoce")
library(ggoce)

data_scores_1 %>%
  distinct(sample, .keep_all = T) %>%
  select(Surface_Salinity_, Surface_Temperature_, den_surf)

size_var <- c("MLD_5", "Meltwater_Fraction_", "month", "latitude", "region")
lapply(size_var, function(d) {
data_scores_1 %>%
  group_by(sample) %>%
  mutate(dom = case_when(
    prop_reads > 0.65 ~ "dominant", 
    TRUE ~ "minor")) %>%
  group_by(sample) %>% 
  mutate(dom = ifelse(
    n_distinct(dom) == 1, "mixed", as.character(clust[grepl("dominant", dom)]))) %>%
  filter(clust == dom | dom == "mixed") %>%
  distinct(sample, dom, .keep_all = T) %>%
  ggplot() +
  geom_isopycnal(aes(x = Surface_Salinity_, y = Surface_Temperature_), 
                 n_breaks = 10, label_placer = label_placer_isopycnal("t")) +
  geom_point(aes(x = Surface_Salinity_, y = Surface_Temperature_, color = dom, size = !!sym(d))) +
  scale_color_manual(values = c(cluster_colors, "mixed" = "grey")) +
  labs(x = "Surface Salinity (PSU)", y = "Surface Temperature (C)",
       color = "Dominant\nCluster", size = gsub("_", " ", capitalize(d)),
       title = "") +
  theme_bw() +
  my_theme +
    facet_wrap(~seasonyear)
  
  saveas <- paste0("figures/exploratory_figs/", d,"-T-S_season.png")
  ggsave(filename = saveas, plot = last_plot(), width = 7, height = 7, units = "in", dpi = 300)
})


size_var <- c("MLD_5", "region")
lapply(size_var, function(d) {
  data_scores_1 %>%
    group_by(sample) %>%
    mutate(dom = case_when(
      prop_reads > 0.65 ~ "dominant", 
      TRUE ~ "minor")) %>%
    group_by(sample) %>% 
    mutate(dom = ifelse(
      n_distinct(dom) == 1, "mixed", as.character(clust[grepl("dominant", dom)]))) %>%
    filter(clust == dom | dom == "mixed") %>%
    distinct(sample, dom, .keep_all = T) %>%
    ggplot() +
    geom_point(aes(x = days_since, y = Meltwater_Fraction_, color = dom, size = !!sym(d))) +
    scale_color_manual(values = c(cluster_colors, "mixed" = "grey")) +
    labs(x = "Days since start of season", y = "Meltwater Fraction",
         color = "Dominant\nCluster", size = gsub("_", " ", capitalize(d)),
         title = "") +
    theme_bw() +
    my_theme +
    facet_wrap(~seasonyear)
  
  saveas <- paste0("figures/exploratory_figs/", d,"-days_meltwater-season.png")
  ggsave(filename = saveas, plot = last_plot(), width = 7, height = 7, units = "in", dpi = 300)
})

data_scores_1 %>%
  group_by(sample) %>%
  mutate(dom = case_when(
    prop_reads > 0.65 ~ "dominant", 
    TRUE ~ "minor")) %>%
  group_by(sample) %>% 
  mutate(dom = ifelse(
    n_distinct(dom) == 1, "mixed", as.character(clust[grepl("dominant", dom)]))) %>%
  filter(clust == dom | dom == "mixed") %>%
  distinct(sample, dom, .keep_all = T) %>%
  ggplot() +
  geom_isopycnal(aes(x = Surface_Salinity_, y = Surface_Temperature_), 
                 n_breaks = 10, label_placer = label_placer_isopycnal("t")) +
  geom_point(aes(x = Surface_Salinity_, y = Surface_Temperature_, color = dom, size = Meltwater_Fraction_)) +
  scale_color_manual(values = c(cluster_colors, "mixed" = "grey")) +
  labs(x = "Days since start of season", y = "Melwater Fraction",
       color = "Prop. reads", size = "MLD (m)",
       title = "With Porosira") +
  theme_bw() +
  my_theme +
  facet_wrap(~seasonyear)





data_scores_1 %>%
  group_by(sample) %>%
  mutate(dom = case_when(
    prop_reads > 0.65 ~ "dominant", 
    TRUE ~ "minor")) %>%
  group_by(sample) %>% 
  mutate(dom = ifelse(
    n_distinct(dom) == 1, "mixed", as.character(clust[grepl("dominant", dom)]))) %>%
  filter(clust == dom | dom == "mixed") %>%
  distinct(sample, dom, .keep_all = T) %>% pull(den_surf)
ggplot() +
  stat_contour(aes(x = Surface_Salinity_, y = Surface_Temperature_, z = den_surf)) +
  geom_point(aes(x = Surface_Salinity_, y = Surface_Temperature_,
                     color = dom, size = MLD_5)) +
  scale_color_manual(values = c(cluster_colors, "mixed" = "grey")) +
  labs(x = "Surface Salinity (PSU)", y = "Surface Temperature (C)",
       color = "Prop. reads", size = "MLD (m)",
       title = "Without Porosira") +
  theme_bw() +
  my_theme


generic_facs <- c("Surface_Temperature_", "Surface_Salinity_", "MLD_5", "Meltwater_Fraction_", "secchi_depth")

lapply(generic_facs, function(r){
data_scores_1 %>%
    pivot_wider(names_from = clust, values_from = prop_reads) %>%
    pivot_longer(cols = c("1", "2", "3", div_facs)) %>%
    select(name, value, !!sym(r)) %>%
    rename("variable" = !!sym(r)) %>%
  group_by(name) %>%
  group_map({~
    .x %>%
     glm(formula = reformulate(paste0(generic_facs, collapse = "+"), "value"), data = .) %>%
      summary()
  })
  })


test_glm <- lapply(generic_facs, function(r){
  data_scores_1 %>%
    pivot_wider(names_from = clust, values_from = prop_reads) %>%
    pivot_longer(cols = c("1", "2", "3", div_facs)) %>%
    select(name, value, !!sym(r)) %>% 
    rename("variable" = !!sym(r)) %>%
    group_by(name) %>%
    group_map({~
        .x %>%
        glm(formula = value ~ variable, data = .) %>%
        summary()
    })
})


lapply(seq_along(generic_facs), function(r){
  lapply(seq_along(c("1","2","3", div_facs)), function(t){
    data.frame(variable = generic_facs[r],
                group = c("1","2","3", div_facs)[t],
                p_val = test_glm[[r]][[t]]$coefficients[2,4],
               aic = test_glm[[r]][[t]]$aic)
    
  }) %>%
    bind_rows()
}) %>% bind_rows() %>%
  mutate(signif = if_else(p_val < 0.05, "Significant", "Not significant")) %>%
  ggplot() +
  geom_point(aes(x = group, y = variable, color = signif, size = aic)) +
  theme_bw() +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(name = "", values = c("Significant" = "red", "Not significant" = "grey80")) +
  labs(x = "Grouping", y = "Envi. Var", color = "", size = "AIC") +
  scale_y_discrete(labels = function(x) gsub("_", " ", capitalize(x)))

ggsave(filename = "figures/exploratory_figs/GLM-AIC-compare.png", plot = last_plot(), width = 7, height = 5, units = "in", dpi = 300)


test_glm[[1]][[1]]$coefficients[2,4]
chisq.test(test_glm[[1]][[3]]$data$prop_reads,
           test_glm[[1]][[3]]$fitted.values)
summary(test_glm[[1]][[3]])

ggplot(data_scores_1) +
  geom_point(aes(x = MDS1, y = MDS2, size=prop_reads, color = clust)) +
  facet_wrap(~clust)

require("Hmisc")
require("corrplot")

groups


bio <- asv_table_rare %>%
  filter(!is.na(phytogroups)) %>%
  select(Species, phytogroups, rare_reads, sample) %>%
  filter(sample %in% unique(metadata_div$Genetics_18sv9_Sample_ID)) %>%
  replace_na() %>%
  mutate(cluster = factor(Species,
                          levels = names(groups), labels = groups)) %>%
  drop_na(cluster) %>%
  pivot_wider(.,
              id_cols = "sample", names_from = "phytogroups", values_from = "rare_reads", values_fn = function(r)sum(r,na.rm = T), values_fill = 0) %>% 
  left_join(metadata_div %>%
              select(Genetics_18sv9_Sample_ID, all_of(generic_facs)),
            by = c("sample" = "Genetics_18sv9_Sample_ID")) %>%
  column_to_rownames(var = "sample")
  
corr_plot_r <- rcorr(as.matrix(bio), type = "pearson")$r
corr_plot_p <- rcorr(as.matrix(bio), type = "pearson")$P

as.data.frame(corr_plot_r) %>%
  rownames_to_column("var2") %>%
pivot_longer(cols = -var2, names_to = "var1", values_to = "r") %>%
  mutate(P = as.vector(corr_plot_p)) %>%
  filter(var1 != var2) %>%
  filter(var1 %in% c("Surface_Temperature_", "Surface_Salinity_", "MLD_5", "Meltwater_Fraction_", "secchi_depth")) %>%
  filter(var2 %in% c("Diatoms", "Cryptophytes", "Dinoflagellates", "Greenalgae", "Haptophytes", "MAST")) %>%
  ggplot() +
  geom_point(aes(x = var1, y = var2, color = r, size = P)) +
  scale_color_gradient2(low = "blue", high = "red") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "", color = "Correlation") +
  my_theme

corr_plot <- rcorr(as.matrix(bio), type = "pearson")
diag(corr_plot$P) <- 0

file_path= "figures/exploratory_figs/corr_plot-phytogroups.png"
png(height=5, width=7, units = "in",res = 300, file=file_path, type = "cairo")

corrplot(corr_plot$r, type = "upper", tl.col = "black", tl.srt = 45, addCoef.col = "grey30", number.cex = 0.5, tl.cex = 0.5,
                     p.mat = corr_plot$P, sig.level = 0.05, insig = "blank")
dev.off()
ggplot(data_scores_1) +
  geom_point(aes(x = MDS1, y = MDS2, size=prop_reads, color = clust)) +
  facet_wrap(~clust)
