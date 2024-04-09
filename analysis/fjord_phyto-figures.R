##### Code to run figures ----
##### Load packages -----
library(ggrepel)
library(cowplot)
library("ggmagnify")
library(fantaxtic)
library(ggnested)
library(patchwork)
##### Load files from analysis ----

### Set up colors -----
region_colors <- c("red4", "green4", "blue", "purple")
names(region_colors) <- c("southern", "middle", "northern", "shetlands")
region_labels <- c("Southern", "Middle", "Northern", "Shetlands")
names(region_labels) <-c("southern", "middle", "northern", "shetlands")

season_colors <- c("#c64575", "#5fb651", "#9e5cc5", "#c58942")
names(season_colors) <- c("2017-2018", "2018-2019", "2019-2020", "2021-2022")

month_colors <- c("#4bb092", "#ca5740", "#6587cd", "#b4b542", "#c979ad", "#6b7f38")
names(month_colors) <- c("11", "12", "1", "2", "3", "4'")

cluster_colors <- c("#30C5FF", "#555B6E", "#F46036")
names(cluster_colors) <- c(1:3)

group_colors <- c("#F06400", "#00F064","#008CF0","gold2","#6400F0","#F0008C","grey10")
names(group_colors) <- c("Cryptophytes", "Diatoms", "Dinoflagellates", "Haptophytes",
                         "MAST", "Rhodophytes", "Greenalgae")

my_theme = theme_linedraw() + theme(text = element_text(size = 14), strip.background = element_blank(), strip.text = element_text(face = "bold", color = "black"))

##### Figure 1-----
### Figure 1 - Sampling effort --- 
map <- map_data("world")

sta_map <- ggplot(metadata %>%
                    mutate(region = factor(region, 
                                           levels = c(
                                             "shetlands","northern", "middle", "southern")))) + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group),
               fill = "grey", color = "black", linewidth = 0.25) +
  coord_map(projection = "ortho",
            xlim = c(-70,-55), ylim = c(-68,-61.8),orientation = c(-100,-80,-12.5)) +
  geom_point(aes(x = lon, y = lat), size = 1.5) +
  geom_label_repel(data = metadata %>%
                     mutate(region = factor(region, 
                                            levels = c(
                                              "shetlands","northern", "middle", "southern"))) %>%
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

axis_color <- c(rep(region_colors["southern"],3),
                rep(region_colors["middle"],21),
                rep(region_colors["northern"],5),
                rep(region_colors["shetlands"], 3))
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
  filename = paste0(figures_local, "map_sampling-effort.jpg"),
  station_map_sampling,
  width = 16,
  height = 5.5,
  units = "in",
  dpi = 300
)


#### Heatmap of rel abun by region -----
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

fig_2 <- full_df %>%
  group_by(site_ids) %>%
  mutate(max = sum(reads, na.rm = T)) %>%
  group_by(phytogroups, Species, region, site_ids) %>%
  reframe(rel = log10(sum(reads, na.rm = T)/unique(max)[1]*100)) %>%
  mutate(Species =  factor(gsub("__", "_", gsub("X","", Species)),
                           levels = species_color$Species),
         phytogroups = factor(phytogroups,
                              levels = c("Diatoms", "Dinoflagellates",
                                         "MAST", "Greenalgae",
                                         "Haptophytes", "Rhodophytes",
                                         "Cryptophytes"), 
                              labels = c("Diatoms", "Dinoflagellates",
                                         "MAST", "Green algae",
                                         "Haptophytes", "Rhodophytes",
                                         "Cryptophytes"))) %>%
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
  facet_grid(phytogroups~region, scales = "free", space = "free",
             labeller = labeller(region = c(
               "southern" = "Southern",
               "middle" = "Middle",
               "northern" = "Northern",
               "shetlands" = "Shetlands"
             ))) +
  coord_cartesian(expand = 0) +
  theme(strip.placement = "outside",
        axis.text.y = element_text(color = "black", size = 7, face = "italic"),
        axis.text.x = element_text(color = "black", size = 10),
        strip.background = element_blank(),
        strip.text = element_text(size = 11, face = "bold"),
        strip.text.y = element_text(angle = 0, hjust = 0),
        strip.clip = "off")

ggsave(
  filename = paste0(figures_local, "heatmap_species-region_relabun.jpg"),
  fig_2,
  width = 11,
  height = 11,
  units = "in",
  dpi = 300
)

#####  Heatmap of species over time (just middle regions) ----
heatmap_time <- full_df %>%
  filter(reads > 0) %>%
  group_by(season, month) %>%
  mutate(max = sum(reads, na.rm = T)) %>%
  group_by(phytogroups, Species, season, month) %>%
  reframe(rel = log10(sum(reads, na.rm = T)/unique(max)[1]*100)) %>%
  mutate(Species =  factor(gsub("__", "_", gsub("X","", Species)),
                           levels = species_color$Species),
         phytogroups = factor(phytogroups,
                              levels = c("Diatoms", "Dinoflagellates",
                                         "MAST", "Greenalgae",
                                         "Haptophytes", "Rhodophytes",
                                         "Cryptophytes"), 
                              labels = c("Diatoms", "Dinoflagellates",
                                         "MAST", "Green algae",
                                         "Haptophytes", "Rhodophytes",
                                         "Cryptophytes"))) %>%
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
  facet_grid(phytogroups~season, scales = "free", space = "free") +
  coord_cartesian(expand = 0) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        strip.placement = "outside",
        axis.text.y = element_text(color = "black", size = 7, face = "italic"),
        axis.text.x = element_text(color = "black", size = 10),
        strip.background = element_blank(),
        strip.text = element_text(size = 11, face = "bold"),
        strip.text.y = element_text(angle = 0, hjust = 0),
        strip.clip = "off")

ggsave(
  filename = paste0(figures_local, "heatmap_species-month_relabun-all.jpg"),
  heatmap_time,
  width = 11,
  height = 11,
  units = "in",
  dpi = 300
)
#####  Diversity map by region -----

plot_list <- full_df %>%
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
        geom_polygon(data = map, aes(x=long, y = lat, group = group),
                     fill = "grey", color = "black", linewidth = 0.25) +
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


######  Diversity over time -----

month_day_range <- full_df_filter %>%
  group_by(month) %>%
  reframe(min = range(days_since)[1],
          max =  range(days_since)[2]) %>%
  mutate(min = ifelse(month == 11, 0, min),
         mid = (min + max)/2,
         label = c("November", "December", "January", "February", "March"),
         max = ifelse(month == 3, NA, max))


div_time <- full_df_filter %>%
  pivot_longer(c(richness_phytogroups, evenness_phytogroups, shannon_phytogroups)) %>%
  mutate(name = factor(name,
                       levels = c("richness_phytogroups", "evenness_phytogroups",
                                  "shannon_phytogroups"),
                       labels = c("Species richness (# ASVs)",
                                  "Species evenness",
                                  "Shannon diversity index"))) %>%
  distinct(samples, season, month, name, value) %>%
  ggplot() +
  geom_boxplot(aes(x = month, y = value, fill = season)) +
  # geom_hline(data = month_day_range %>%
  #              mutate(name = factor("Shannon diversity index",
  #                                   levels = c("Species richness (# ASVs)",
  #                                              "Species evenness","Shannon diversity index"))),
  #            aes(yintercept = 0)) +
  # geom_label(data = month_day_range %>%
  #              mutate(name = factor("Shannon diversity index",
  #                                   levels = c("Species richness (# ASVs)",
  #                                              "Species evenness","Shannon diversity index"))),
  #            aes(x = mid, y = 0, label = label),
  #            label.size = NA, color = "grey40") +
  #geom_vline(data = month_day_range, aes(xintercept = max), color = "grey70") +
  # geom_point(aes(x = month, y = value, color = season)) +
  # stat_smooth(aes(x = month, y = value, color = season, group = season),
  #             show.legend = F, alpha = 1/5) +
  facet_grid(name~., scales = "free_y", switch = "y") +
  scale_fill_manual(name = "", 
                    values = season_colors) +
  my_theme +
  theme(strip.placement = "outside",
        legend.spacing.x = unit(0, "lines")) +
  labs(x = "Month", y = "")
  
ggsave(
  filename = paste0(figures_local, "div_time-e.jpg"),
  div_time,
  width = 8,
  height = 6.5,
  units = "in",
  dpi = 300
)


#### Fig. 6 - In between plots -----
pull.genus <- full_df_filter %>%
  filter(Species != "Porosira_sp.") %>%
  group_by(phytogroups, Genus, Species) %>%
  reframe(total = sum(reads, na.rm = T)) %>%
  ungroup() %>%
  slice_max(order_by = total, by = phytogroups) %>% pull(Species)

full_df_filter %>%
  #filter(Genus == "MAST-12B_XX") %>%
  group_by(phytogroups) %>%
  group_map(~{
    .x %>%
      group_by(Genus) %>%
      {if(length(unique(days_since[reads > 0])) > 10)
        .x %>%
     ggplot() +
      geom_point(aes(x = days_since, y = log10(reads)),
                 alpha = 1/5, show.legend = T) +
      # geom_smooth(aes(x = days_since, y = log10(reads)),
      #             method = "loess",
      #             show.legend = F, se = F) +
      #facet_wrap(~phytogroups) +
      scale_x_continuous(name = "Days since start of season",
                         limits = c(0,140)) +
      scale_y_continuous(name = bquote(log[10]~reads)) +
      my_theme
    }else{NULL}
  })

        
species.to.fit <- full_df_filter %>%
  group_by(phytogroups, Genus, Species) %>%
  reframe(pts = length(unique(days_since[reads >0]))) %>%
  filter(pts >= 10) %>% pull(Species)

full_df_filter %>%
  filter(phytogroups == "Diatoms" & Species %in% species.to.fit) %>%
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
      wrap_plots(., guides = "collect")}, .keep = T)

genus.to.plot <- c("Polar-centric-Mediophyceae_XX", "Odontella", )       


full_df_filter %>%
  group_by(samples) %>%
  mutate(sample_tot = sum(reads)) %>%
  group_by(phytogroups, Genus, Species, season, samples) %>%
  reframe(days_since = unique(days_since),
          sp_rel = sum(reads)/sample_tot) %>%
  filter(Genus == "Porosira") %>%
  ggplot() +
  geom_point(aes(x = days_since, y = sp_rel, color = season)) +
  stat_smooth(aes(x = days_since, y = sp_rel, color = season),
              method = "loess", alpha = 1/5)
  

model.list.names <- full_df_filter %>%
  filter(Species %in% species.filter) %>%
  group_by(phytogroups, Genus, Species, season) %>% group_keys() %>% mutate(list.id = paste(phytogroups, Genus, Species,season,sep = ":"))

lo.x <- 1:140

species.model.list <- full_df_filter %>%
  filter(Species %in% unique(species.filter$Species)) %>%
  group_by(samples) %>%
  mutate(sample_tot = sum(reads)) %>%
  group_by(phytogroups, Genus, Species, season, samples) %>%
  reframe(days_since = unique(days_since),
          sp_rel = sum(reads)/sample_tot) %>%
  select(-samples) %>%
  nest(data = c(days_since, sp_rel)) %>%
  mutate(
    model = map(data, ~ {
      if(length(unique(.$days_since[.$sp_rel >0])) < 5 | nrow(.) == 0){
      NULL
    }else{
      .x %>%
        loess(sp_rel ~ days_since, ., span = 0.75, degree = 2)
    }}),
    fit = map(model, ~{
      if(is.null(.)){
      vector(length = length(lo.x))
    }else{
      predict(., lo.x)
    }
    }%>% as_tibble() %>% mutate(days_since = lo.x, value = ifelse(value<0,0,value))),
    day_max = map(fit, ~{.$days_since[which.max(.$value)[1]]}),
    day_sd = map(fit, ~{
      sqrt(
        sum(.$value*(.$days_since - .$days_since[which.max(.$value)[1]])^2, na.rm = T)/
          (((length(.$value>0)-1)/length(.$value))*sum(.$value, na.rm = T)))
    })
  )
  

genus.plot <- species.model.list %>%
  mutate(total_samples = map(data, ~{
    length(unique(.$days_since[.$sp_rel>0 & !is.na(.$sp_rel)]))})) %>%
 unnest(total_samples)%>%
  group_by(phytogroups) %>%
  filter(Genus != "Porosira") %>%
  slice_max(total_samples) %>% pull(Genus)


tag_facet2 <-  function(p, open="(", close = ")",
                        tag_pool = letters,
                        x = 0, y = 0.5,
                        hjust = 0, vjust = 0.5, 
                        fontface = 2, ...){
  
  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  nm <- names(gb$layout$facet$params$rows)
  
  tags <- paste0(open,tag_pool[unique(lay$PANEL)],close)
  
  tl <- lapply(tags, grid::textGrob, x=x, y=y,
               hjust=hjust, vjust=vjust, gp=grid::gpar(fontface=fontface))
  
  g <- ggplot_gtable(gb)
  g <- gtable::gtable_add_rows(g, grid::unit(1,"line"), pos = 0)
  lm <- unique(g$layout[grepl("strip",g$layout$name), "l"])
  tm <- unique(g$layout[grepl("strip",g$layout$name), "t"])
  g <- gtable::gtable_add_grob(g, grobs = tl, t=1, l=lm)
  grid::grid.newpage()
  grid::grid.draw(g)
}

tag_facet2(example_max_plot)

example_max_plot <- species.model.list %>%
  select(-c(model)) %>%
  filter(Genus %in% "Porosira") %>%
  mutate(data = map(data, ~{left_join(data.frame(days_since = 1:140),.) %>%
      group_by(days_since) %>%
      reframe(sp_rel = mean(sp_rel, na.rm = T))})) %>%
 unnest(c(data, fit, day_max, day_sd), names_sep = ".") %>%
  ggplot() +
  geom_rect(aes(xmin = day_max - day_sd, xmax = day_max+day_sd, ymin = -Inf, ymax = Inf),
            fill = "red", alpha = 1/5, stat = "unique") +
  geom_point(aes(x = data.days_since, y = data.sp_rel)) +
  geom_line(aes( x = fit.days_since, y = fit.value, group = 1)) +
  geom_vline(aes(xintercept = day_max), color = "red", linetype = "dashed") +
  facet_wrap(~season) +
  my_theme +
  labs(y = "Relative proportion of reads", x = "Days since start of season") 
  

ggsave(
  filename = paste0(figures_local, "example_porosira.png"), 
  plot = example_max_plot, 
  width = 5, height = 5, units = "in")





#### Spearman cor ----
bio <- full_df %>%
  filter(reads > 0) %>%
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
rownames(myannotation) <- gsub("__", "_", gsub("X","", rownames(myannotation)))
names(myannotation)[1] = "Cluster" 
myannotation$Cluster <- factor(myannotation$Cluster, levels= 1:3, 
                               labels=1:3)
ann_colors <- list(Cluster = cluster_colors)


colnames(correlation_matrix) <- gsub("__", "_", gsub("X","", colnames(correlation_matrix)))
rownames(correlation_matrix) <- gsub("__", "_", gsub("X","", rownames(correlation_matrix)))


cluster_cor_plot <- pheatmap::pheatmap(correlation_matrix,
                   main = "Species Co-Occurrence Matrix (Spearman's Correlation)",
                   fontsize_row = 6,
                   fontsize_col = 6, clustering_method = "ward.D2",
                   cutree_rows = 3, cutree_cols = 3,
                   annotation_col = myannotation, annotation_colors = ann_colors)

ggsave(
  filename = paste0(figures_local, "cluster_correlation-2.png"), 
  plot = cluster_cor_plot, 
  width = 9, height = 9, units = "in" )


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
                   Species, levels = rownames(myannotation), labels = myannotation$Cluster))) +
  stat_ellipse(aes(x = MDS1, y = MDS2,
                   color = factor(
                     Species, levels = rownames(myannotation), labels = myannotation$Cluster))) +
  my_theme +
  theme(legend.position = c(0.23,0.94),
        legend.direction = "horizontal",
        legend.background = element_rect(color = "black", fill = "white"),
        axis.text = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 12),
        panel.background = element_rect(color = "black"),
        panel.grid = element_blank()) +
  scale_color_manual(name = "Cluster", values = cluster_colors)
nmds_clust

add_cluster <- left_join(full_df, as.data.frame(groups) %>%
                           rownames_to_column(var = "Species"),
                         by = "Species")


names(group_colors)[7]<- "Green algae"

library("ggpattern")


my_pattern_function <- function(pat_var, group_var, subgroup_var) {
  # Define patterns based on unique values of group and subgroup variables
  patterns <- c("dots", "circles", "squares", "stars")
  
  # Create a data frame with patterns assigned to group and subgroup combinations
  patterns_df <- expand.grid(group = unique(group_var), subgroup = unique(subgroup_var))
  patterns_df$pattern <- rep(c(group_patterns, subgroup_patterns), length.out = nrow(patterns_df))
  
  # Return a vector of patterns corresponding to each combination of group and subgroup
  return(patterns_df$pattern)
}


cluster_bars <- add_cluster %>%
  drop_na(groups) %>%
  group_by(groups) %>%
  mutate(clust_sum = sum(reads, na.rm = T)) %>%
  group_by(Species, groups, phytogroups) %>%
  reframe(sp_rel = sum(reads, na.rm = T)/clust_sum) %>%
  distinct(Species, groups, phytogroups, .keep_all = T) %>%
  group_by(groups) %>%
  slice_max(sp_rel, n = 5) %>%
  ungroup() %>%
  mutate(
        phytogroups = factor(phytogroups,
                              levels = c("Diatoms", "Dinoflagellates",
                                        "MAST", "Greenalgae",
                                        "Haptophytes", "Rhodophytes",
                                        "Cryptophytes"), 
                              labels = c("Diatoms", "Dinoflagellates",
                                         "MAST", "Green algae",
                                         "Haptophytes", "Rhodophytes",
                                         "Cryptophytes")),
         pat = case_when(phytogroups == "Diatoms" & groups == "2"~"wave",
                         phytogroups == "Diatoms" & groups == "3"~"stripe",
                         TRUE~"none"),
         pat = if_else(Species == "Porosira_sp.", "circle", pat)) %>%
  ggnested_pattern(data = ., aes_string(main_group = "phytogroups",
                                sub_group = "Species", 
                                x = "groups", y = "sp_rel", pattern = "pat"),
           main_palette = group_colors) + 
  geom_col_pattern(pattern_fill="white", pattern_color="white", pattern_key_scale_factor=0.25,
                   pattern_size = 0.1) +
  scale_y_continuous(name = "Relative abundance",
                     expand = c(0,0)) +
  scale_x_continuous(name = "Cluster", expand = c(1E-3,1E-3), breaks = c(1,2,3)) +
  theme_nested(theme_linedraw) + 
  theme(
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(size = 12),
    legend.title = element_blank(),
    panel.border = element_rect(color = "black")) +
  guides(fill=guide_legend(ncol = 1))


cluster_id_plot <- nmds_clust + cluster_bars +
  plot_annotation(tag_levels = 'A', tag_suffix = ")", tag_prefix = "(") +
  theme(text = element_text(size = 13))


ggsave(filename = paste0(figures_local, "cluster_ids-pattern.jpg"),
       cluster_id_plot,
       width = 15,
       height = 6,
       units = "in",
       dpi = 300)


dom_cluster_region <- add_cluster %>%
  drop_na(groups) %>%
  group_by(season, days_since, Site_Name_2) %>%
  mutate(sample_sum = sum(reads, na.rm = T)) %>%
  group_by(season, days_since, Site_Name_2, groups) %>%
  reframe(clust_rel = sum(reads, na.rm = T)/sample_sum) %>%
  distinct() %>%
  group_by(season, days_since, Site_Name_2) %>%
  slice_max(clust_rel) %>%
  ggplot() + 
  geom_point(aes(x = days_since, y = Site_Name_2, color = factor(groups))) +
  facet_grid(~season) +
  my_theme +
  labs(x = "Days since start of season", y = "") +
  scale_color_manual(name = "Dominant\ncluster", values = cluster_colors) +
  theme(
    axis.text.y = element_text(color = axis_color)
  )

ggsave(filename = paste0(figures_local, "cluster_regions.jpg"),
       dom_cluster_region,
       width = 10,
       height = 7,
       units = "in",
       dpi = 300)


cluster_time <- add_cluster %>%
  drop_na(groups) %>%
  group_by(season, days_since, Site_Name_2) %>%
  mutate(sample_sum = sum(reads, na.rm = T)) %>%
  group_by(season, days_since, Site_Name_2, groups) %>%
  reframe(clust_rel = sum(reads, na.rm = T)/sample_sum) %>%
  distinct() %>%
  ggplot() + 
  geom_point(aes(x = days_since, y = as.numeric(clust_rel), color = factor(groups)), alpha = 1/2) +
  geom_smooth(aes(x = days_since, y = as.numeric(clust_rel), color = factor(groups)),
             show.legend = F, method = "loess", alpha = 1/5) +
  facet_grid(~season) +
  my_theme +
  theme(legend.position = c(0.12, 0.95),
        legend.direction = "horizontal",
        legend.background = element_rect(color = "black")) +
  labs(x = "Days since start of season") +
  scale_y_continuous(name = "Proportion of reads", limits = c(0,1.15)) +
  scale_color_manual(name = "Cluster", values = cluster_colors,
                     guide = guide_legend(override.aes = list(alpha = 1)))


ggsave(filename = paste0(figures_local, "cluster_time.jpg"),
       cluster_time,
       width = 10,
       height = 5,
       units = "in",
       dpi = 300)




#sp_order_plot
plot.periods <- data.frame(xmin = c(0, 51, 100), xmax = c(50, 100, 140), ymax = Inf, ymin = -Inf, period = c("Early","Mid","Late"))

phytogroup_labels <- c("Cryptophytes", "Diatoms",
                       "Dinoflagellates",
                       "Green algae", "Haptophytes",
                       "MAST", "Rhodophytes")
names(phytogroup_labels) <- unique(species.model.list$phytogroups)


species.model.list %>%
  unnest(c(day_max,day_sd)) %>%
  select(-c(data,model,fit)) %>%
  pivot_wider(id_cols = c(phytogroups, Genus, Species), names_from = season,
              values_from = day_max) %>%
  pivot_longer(cols = !c(phytogroups, Genus, Species),
               names_to = "season", values_to = "day_max") %>%
  mutate(day_max = replace_na(day_max, 141)) %>%
  pivot_wider(id_cols = c(phytogroups, Genus, Species), names_from = season,
              values_from = day_max) %>%
  mutate(phytogroups = factor(phytogroups,
                              levels = rev(unique(species.model.list$phytogroups)))) %>%
  group_by(phytogroups) %>%
  arrange(desc(`2017-2018`), .by_group = T) %>% 
  mutate(Species = factor(Species, levels = Species, labels = Species)) %>%
  ungroup() %>%
  mutate(sp_num = rev(1:n())) %>%
  pivot_longer(cols = !c(phytogroups, Genus, Species, sp_num),
               names_to = "season", values_to = "day_max") %>%
  left_join(., species.model.list %>%
              unnest(c(day_max,day_sd)) %>%
              select(-c(data,model,fit, day_max)))  %>%
  mutate(sd = replace_na(day_sd, 0)) %>%
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
      {if(unique(.x$season)%in%unique(species.model.list$season)[1:3]){
        theme(strip.text.y = element_blank())}}+
      {if(unique(.x$season)%in%unique(species.model.list$season)[2:4]){
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

species.model.list %>%
  unnest(c(day_max,day_sd)) %>%
  select(-c(data,model,fit)) %>%
  left_join(., as.data.frame(groups) %>%
              rownames_to_column(var = "Species")) %>%
  ggplot() +
  geom_boxplot(aes(x = factor(groups, levels = c(3,2,1)), y = day_max, fill = factor(groups))) +
  coord_flip() +
  labs(y = "Days since start of season", x = "Clusters") +
  scale_fill_manual(name = "", values = cluster_colors) +
  my_theme +
  facet_grid(~season)


ggsave(
  paste0(figures_local, "day_change_plot.png"),
  dpi = 300,
  width = 7,
  height = 7,
  unit = "in")



