
full_df %>%
  ggplot() +
  geom_point(aes(x = days_since, y = Surface_Salinity_)) +
  geom_smooth(aes(x = days_since, y = Surface_Salinity_),
              method = "gam", formula = y~poly(x,2))

full_df %>%
  mutate(month = factor(month, levels = c("11", "12", "01", "02", "03"))) %>%
  ggplot() +
  geom_boxplot(aes(x = month, y = Surface_Salinity_)) 

full_df %>%
  ggplot() +
  geom_point(aes(x = Surface_Temperature_, y = Surface_Salinity_)) +
  geom_smooth(aes(x = Surface_Temperature_, y = Surface_Salinity_),
              method = "gam", formula = y~poly(x,2)) +
  facet_wrap(~day)


group_color_df <- data.frame(colors = as.character(group_colors),
                             names = names(phyto_labels),
                             labels = phyto_labels)


library(ggpmisc)

plot_a <- left_join(full_df, diversity_filter) %>%
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
  stat_poly_eq(aes(x = values, y = richness, color = group),
               formula = y ~ x, parse = T, size = 3) +
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
  coord_cartesian(ylim = c(0, 35), expand = 0)
  

ggsave(
  filename = "~/Documents/Alaina/richness_all-group-r2.png",
  plot_a,
  width = 9,
  height = 4,
  units = "in",
  dpi = 300
)

phyto_labels <- c(
  richness_Cryptophytes = "Cryptophytes",
  richness_Diatoms = "Diatoms", 
  richness_Dinoflagellates = "Dinoflagellates",
  richness_Haptophytes = "Haptophytes",
  richness_MAST = "MAST", 
  richness_Rhodophytes = "Rhodophytes", 
  richness_Greenalgae = "Green algae")

plot_a <- left_join(full_df, diversity_filter) %>%
  select(-grep("_all", colnames(.))) %>% 
  select(-contains("_phytogroups")) %>%
  pivot_longer(cols = contains("richness_"), names_to = "group",
               values_to = "richness") %>%
  drop_na(richness) %>%
  filter(richness>0) %>%
  ggplot() +
  geom_point(aes(x = Surface_Temperature_, y = richness)) +
  facet_wrap(~group, scales = "free_y",
             labeller = labeller(group = phyto_labels)) +
  theme_bw() +
  my_theme +
  labs(x = "Temperature", y = "# species")

ggsave(
  filename = "~/Documents/Alaina/group_richness_temp.png",
  plot_a,
  width = 8,
  height = 4,
  units = "in",
  dpi = 300
)

plot_a <- left_join(full_df, diversity_filter) %>%
  pivot_longer(cols = c("days_since", "Surface_Temperature_", "Surface_Salinity_"), names_to = "var",
               values_to = "values") %>%
  ggplot() +
  geom_point(aes(x = values, y = richness_phytogroups)) +
  facet_wrap(~var, ncol = 1, scales = "free_x",
             switch = "x", labeller = labeller(var = c(
               days_since = "Days since start of Oct. 1",
               Surface_Salinity_ = "Surface Salinity",
               Surface_Temperature_ = "Surface Temperature"
             ))) +
  theme_bw() +
  my_theme +
  labs(x = "", y = "# Species") +
  theme(
    strip.placement = "outside",
    strip.text = element_text(face = "plain")
  )

ggsave(
  filename = "~/Documents/Alaina/all_richness_vars.png",
  plot_a,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)

left_join(full_df, diversity_filter) %>%
  select(-c("richness_all", "richness_phytogroups")) %>%
  pivot_longer(cols = contains("richness_"), names_to = "var",
               values_to = "values") %>%
  ggplot() +
  geom_bar(aes(x = days_since, y = values, fill = var), stat = "identity") +
  theme_bw() +
  my_theme +
  labs(x = "", y = "# Species")

plot_a <- metadata_filter %>%
  ggplot() +
  geom_point(aes(x = days_since, y = Surface_Temperature_)) +
  geom_smooth(aes(x = days_since, y = Surface_Temperature_),
              method = "gam", color = "black") +
  geom_point(aes(x = days_since, y = Surface_Salinity_-32), color = "red") +
  geom_smooth(aes(x = days_since, y = Surface_Salinity_-32), color = "red",
              method = "gam") +
  theme_bw() +
  my_theme +
  labs(title = "A.") +
  scale_x_continuous(name = "Days since Oct. 1", expand = c(1E-2,1E-2)) +
  scale_y_continuous(name = "Temperature",
                     sec.axis = sec_axis(name = "Salinity", trans = ~.+32)) +
  theme(axis.text.y.right = element_text(color = "red"),
        axis.title.y.right = element_text(color = "red"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

plot_b <- left_join(full_df, diversity_filter) %>%
  select(-c("richness_all", "richness_phytogroups")) %>%
  pivot_longer(cols = contains("richness_"), names_to = "var",
               values_to = "values") %>%
  ggplot() +
  geom_bar(aes(x = days_since, y = values, fill = var), stat = "identity") +
  scale_x_continuous(limits = c(0, 140), expand = c(1E-2,1E-2)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  my_theme +
  labs(x = "Days since Oct. 1", y = "# Species", title = "B.") +
  scale_fill_manual(name = "", values = as.character(group_colors),
                     labels = phyto_labels) +
  theme(legend.position = c(0.1, 0.7),
        legend.background = element_blank())

plot_c <- plot_a / plot_b


ggsave(
  filename = "~/Documents/Alaina/temp_sal_richness_bar.png",
  plot_c,
  width = 9,
  height = 7,
  units = "in",
  dpi = 300
)

metadata_filter %>%
  ggplot() +
  geom_point(aes(x = days_since, y = Surface_Temperature_)) +
  geom_smooth(aes(x = days_since, y = Surface_Temperature_),
              method = "gam", color = "black") +
  geom_point(aes(x = days_since, y = Surface_Salinity_-32), color = "red") +
  geom_smooth(aes(x = days_since, y = Surface_Salinity_-32), color = "red",
              method = "gam") +
  theme_bw() +
  my_theme +
  labs(title = "A.") +
  scale_x_continuous(name = "Days since Oct. 1", expand = c(1E-2,1E-2)) +
  scale_y_continuous(name = "Temperature",
                     sec.axis = sec_axis(name = "Salinity", trans = ~.+32)) +
  theme(axis.text.y.right = element_text(color = "red"),
        axis.title.y.right = element_text(color = "red"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())





left_join(full_df, diversity_filter) %>%
  select(-c("richness_all", "richness_phytogroups")) %>%
  pivot_longer(cols = contains("richness_"), names_to = "var",
               values_to = "values") %>%
  ggplot() +
  geom_point(aes(x = days_since, y = values, color = var)) +
  theme_bw() +
  my_theme +
  labs(x = "", y = "# Species")

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

bins_df <- full_df %>%
  drop_na(phytogroups) %>%
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
  filter(prop_reads == max(prop_reads, na.rm = T),
         .by = c("phytogroups", "Genus", "Species")) %>%
  left_join(., voxel_ids)


bins_df %>%
  ggplot() +
  geom_point(aes(x = days_since, y = Surface_Salinity_)) +
  facet_wrap(~days_cuts, scales = "free_x")



bins_df %>%
  mutate(Species = factor(Species, levels = sp_colors_df$sublabel)) %>%
  ggplot() +
  geom_text(data = voxel_ids, aes(x = x_val, y = y_val, label = group_id),
            fontface = "bold") +
  # geom_text(data = bins_missing, aes(x = x_val+0.5, y = y_val-0.3),
  #           label = "n.d.", fontface = "bold") +
  geom_hline(yintercept = c(32.5, 33.5),
             linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(0, 1) , linetype = "dashed", color = "grey40") +
  geom_jitter(aes(x = Surface_Temperature_,
                  y = Surface_Salinity_,
                  fill = Species), shape = 21,
              color = "black", size = 3,
              width = 0.25, height = 0.1) +
  scale_fill_manual(name = "", values = sp_colors_df$colors,
                    labels = sp_colors_df$sublabel, drop = F) +
  facet_wrap(~factor(days_cuts)) +
  theme_bw() +
  my_theme +
  theme(legend.text = element_markdown(size = 9)) +
  labs(x = bquote(Surface~Temperature~(degree*C)),
       y = "Salinity") +
  scale_shape_discrete(name = "")


env_facs <- c("days_since", "Surface_Temperature_", "Surface_Salinity_")

bio_id <- full_df %>%
  filter(!is.na(phytogroups)) %>%
  select(Feature.ID, Species, phytogroups, rare_reads, sample, days_since,
         Surface_Salinity_, Surface_Temperature_) %>%
  pivot_wider(id_cols = c("sample", "days_since", "Surface_Temperature_", "Surface_Salinity_"), names_from = "Feature.ID", values_from = "rare_reads", values_fn = function(r)sum(r,na.rm = T), values_fill = 0) %>% 
  column_to_rownames(var = "sample") 

bio_species <- full_df %>%
  filter(!is.na(phytogroups)) %>%
  select(Feature.ID, Species, phytogroups, rare_reads, sample, days_since,
         Surface_Salinity_, Surface_Temperature_) %>%
  pivot_wider(id_cols = c("sample", "days_since", "Surface_Temperature_", "Surface_Salinity_"), names_from = "Species", values_from = "rare_reads", values_fn = function(r)sum(r,na.rm = T), values_fill = 0) %>% 
  column_to_rownames(var = "sample") 

bio_group <- full_df %>%
  filter(!is.na(phytogroups)) %>%
  select(Feature.ID, Species, phytogroups, rare_reads, sample, days_since,
         Surface_Salinity_, Surface_Temperature_) %>%
  pivot_wider(id_cols = c("sample", "days_since", "Surface_Temperature_", "Surface_Salinity_"), names_from = "phytogroups", values_from = "rare_reads", values_fn = function(r)sum(r,na.rm = T), values_fill = 0) %>% 
  column_to_rownames(var = "sample") 

bio_all <- left_join(left_join(bio_id, bio_species), bio_group)

corr_plot_r <- rcorr(as.matrix(bio_species), type = "pearson")$r
corr_plot_p <- rcorr(as.matrix(bio_species), type = "pearson")$P

plot_a <- as.data.frame(corr_plot_r) %>%
  rownames_to_column("var2") %>%
  pivot_longer(cols = -var2, names_to = "var1", values_to = "r") %>%
  mutate(P = as.vector(corr_plot_p)) %>%
  filter(var1 != var2) %>%
  filter(var2 %in% env_facs & !var1 %in% env_facs ) %>%
  filter(P < 0.05) %>%
  mutate(spec = factor(var1, levels = rev(sp_colors_df$sublabel))) %>%
  # left_join(full_df %>% select(Feature.ID, Species, phytogroups) %>% 
  #             distinct(), by = c("var1" = "Feature.ID")) %>%
  #filter(v == "Diatoms") %>%
  ggplot() +
  geom_tile(aes(x = var2, y = spec, fill = r), size = 3) +
  scale_fill_gradient2(low = "blue", high = "red") +
  theme_bw() +
  labs(x = "", y = "", fill = "Correlation") +
  my_theme +
  theme(axis.text.y = element_markdown(size = 8)) +
  scale_y_discrete(drop = F) + 
  scale_x_discrete(labels = c("Days since Oct. 1", 
                              "Salinity",
                              "Temperature"),
                   expand = c(0,0))

ggsave(
  filename = "~/Documents/Alaina/species_corr.png",
  plot_a,
  width = 6,
  height = 9,
  units = "in",
  dpi = 300
)






