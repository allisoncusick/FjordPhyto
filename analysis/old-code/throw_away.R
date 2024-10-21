##### Junk code

### Add diversity values to full metadata
metadata_div <- left_join(metadata_filter, df_return, by = c("Genetics_18sv9_Sample_ID" = "sample"))

metadata_div[
  grepl("Andvord", metadata_div$site_name_nickname) |
    grepl("Transect", metadata_div$site_name_nickname),] %>%
  mutate_if(is.list, ~paste(unlist(.), collapse = '|')) %>%
  mutate_if(str_detect(., "\\|"), ~str_split(., "\\|")) %>% pull(Temperature_)
#if the character string has a | in it, then convert to a list and separate values by "|" 
transect_sub <- transect_sub %>%
  mutate_if(str_detect(., "\\|"), ~str_split(., "\\|")) %>%
  unnest(cols = everything())

asv_table_rare %>%
  filter(sample %in% transect_sub$Genetics_18sv9_Sample_ID) %>%
  write.csv('transect_subset-asv.csv', row.names = FALSE)




metadata_ch2 <- metadata_w_samples %>%
  mutate(
    days_since = yday(date_local) - yday("2000-11-01"),
    days_since = ifelse(days_since < 0, days_since + 365, days_since))

ggplot(metadata_ch2) + 
  geom_point(aes(x = date_local, y = as.numeric(MLD_5), color = seasonyear)) +
  geom_smooth(aes(x = date_local, y = as.numeric(MLD_5), color = seasonyear), method = "gam", se = FALSE) +
  facet_wrap(~seasonyear, scales = "free_x")

ggplot(metadata_ch2) + 
  geom_point(aes(x = longitude, y = latitude, color = as.numeric(chlor_a))) +
  facet_wrap(~seasonyear)

ggplot(metadata_ch2, 
       aes(
         x = as.numeric(Surface_Temperature_),
         y = as.numeric(chlor_a),
         color = seasonyear)) + 
  geom_point() +
  geom_smooth(method = "gam", se = FALSE)

ggplot(metadata_ch2, 
       aes(
         x = days_since,
         y = as.numeric(chlor_a),
         color = as.numeric(Surface_Temperature_))) + 
  geom_point() +
  geom_smooth(method = "gam", se = FALSE)

colnames(metadata_ch2)

map <- map_data("world") %>%
  filter(region %in% c("Antarctica", "Falkland Islands", "South Georgia and the South Sandwich Islands"))

fix_names <- read_csv(paste0(data_home, "site_names_edit_manual.csv"), skip = 0)

require("geosphere")
dist_vec <- matrix(ncol = 2, nrow = length(unique(metadata_w_samples$location)))
for (r in 1:length(unique(metadata_w_samples$location))) {
  lon1 <- unique(metadata_w_samples$longitude[metadata_w_samples$location == "Stonington Island"])
  lat1 <- unique(metadata_w_samples$latitude[metadata_w_samples$location == "Stonington Island"])
  lon2 <- mean(unique(metadata_w_samples$longitude[metadata_w_samples$location == unique(metadata_w_samples$location)[r]]))
  lat2 <- mean(unique(metadata_w_samples$latitude[metadata_w_samples$location == unique(metadata_w_samples$location)[r]]))
  dist_vec[r,2] <- distm(c(lon1, lat1), c(lon2, lat2), fun = distHaversine)
  dist_vec[r,1] <- unique(metadata_w_samples$location)[r]
}

dist_df <- arrange(as.data.frame(dist_vec), V2) %>%
  mutate(site_num = 1:nrow(dist_vec))
colnames(dist_df) <- c("location", "distance_to_stonington", "site_num")


metadata_w_samples$site_ids <- factor(metadata_w_samples$location, levels = dist_df$location, labels = dist_df$site_num)

metadata_w_samples$region <-"shetlands" #assign everything something, start with northern then parse/ID other labels
#metadata$region[which( latitude less than -63 == northern, latitude 63 - 65 = middle, latitude >65 southern)]
metadata_w_samples$region[which(metadata_w_samples$latitude < -63 & metadata_w_samples$latitude > -66)] <-"middle"
metadata_w_samples$region[which(metadata_w_samples$latitude < -66 & metadata_w_samples$latitude > -73)] <-"southern"
metadata_w_samples$region[which(metadata_w_samples$longitude < -50 & metadata_w_samples$longitude > -58)] <-"northern"
metadata_w_samples_2 <- metadata_w_samples %>%
  filter(Genetics_18sv9_Sample_ID %in% rownames(bio))

sta_map <- ggplot(metadata_w_samples_2 %>%
                    mutate(region = factor(region, 
                                           levels = c(
                                             "shetlands","northern", "middle", "southern")))) + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group),
               fill = "grey", color = "black", linewidth = 0.25) +
  coord_map(projection = "ortho",
            xlim = c(-70,-55), ylim = c(-68,-61.8),orientation = c(-100,-80,-12.5)) +
  geom_point(aes(x = longitude, y = latitude), size = 1.5) +
  geom_label_repel(data = metadata_w_samples %>%
                     mutate(region = factor(region, 
                                            levels = c(
                                              "shetlands","northern", "middle", "southern"))) %>%
                     group_by(location, site_ids, longitude, latitude) %>%
                     reframe(region = unique(region)),
                   aes(x = longitude, y = latitude, label = site_ids, fill = region),
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
  theme(
    panel.grid = element_line(color = "grey"),
    axis.title = element_blank(),
    panel.background = element_rect(color = "black", fill = "transparent", linewidth = 0.5))
#require("cowplot")
sta_map_inset <- ggdraw(sta_map, clip = "on") +
  draw_plot(ant_inset, x = 0.12, y = 0.65, width = 0.3, height = 0.3)

axis_color <- c(rep(region_colors["southern"],2),
                rep(region_colors["middle"],16),
                rep(region_colors["northern"],5),
                rep(region_colors["shetlands"], 3))
sta_samples <- metadata_w_samples_2 %>%
  filter(!is.na(site_ids)) %>%
  group_by(region, site_ids, seasonyear, month) %>%
  reframe(
    total_samples = length(unique(Genetics_18sv9_Sample_ID))) %>%
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
  facet_grid(~seasonyear)

station_map_sampling <- plot_grid(sta_map, sta_samples, align = "h", axis = "t",
                                  labels = c("(A)", "(B)"), nrow = 1, rel_widths = c(1,1.5))

ggsave(
  filename = "figures/map_sampling-effort-part2.jpg",
  station_map_sampling,
  width = 16,
  height = 5.5,
  units = "in",
  dpi = 300
)



data_scores_1 %>%
  filter(location %in% grepl("Transect", location))

##### Analysis for Ch. 2 (now 3)

#source()
my_theme = theme_linedraw() + theme(text = element_text(size = 14), strip.background = element_blank(), strip.text = element_text(face = "bold", color = "black"))

# diversity_df, asv_table_rare, metadata_filter
metadata_filter <- metadata_w_samples %>%
  drop_na(secchi_depth) %>%
  rename(sample = Genetics_18sv9_Sample_ID) %>%
  filter(!sample %in% samples_removed)

asv_filter <- asv_table_rare %>%
  drop_na(phytogroups) %>%
  filter(sample %in% unique(metadata_filter$sample)) %>%
  group_by(Feature.ID) %>%
  reframe(total_reads = sum(rare_reads, na.rm = T),
          total_samples = length(unique(sample[rare_reads > 0]))) %>%
  filter(total_reads > 5 & total_samples > 5) %>%
  pull(Feature.ID) %>%
  unique()

asv_table_filter <- asv_table_rare %>%
  filter(rare_reads > 0) %>%
  drop_na(phytogroups) %>%
  filter(sample %in% unique(metadata_filter$sample)) %>%
  filter(Feature.ID %in% asv_filter)

diversity_filter <- diversity_df %>%
  filter(sample %in% unique(metadata_filter$sample))

#### Driving questions: How does salinity and temperature influence the microbial community?

## Pull the 9 most sampled locations and filter the data

stations <- c("Cierva Cove","Cuverville Island", "Danco Island",
              "Neko Harbour", "Orne Harbour", "Paradise Harbour", "Petermann Island", "Wilhelmina Bay")
env_vars <- c("Surface_Temperature_", "Surface_Salinity_", "month", "Meltwater_Fraction_", "days_since", "MLD_20", "secchi_depth")
p_groups <- unique(asv_table_filter$phytogroups)

## Raw values 
lapply(stations, function(s){
  lapply(env_vars, function(e){
    #lapply(p_groups, function(r){
    plot_a <- left_join(metadata_filter, asv_table_filter) %>%
      filter(location == s) %>%
      mutate(value = !!sym(e)) %>%
      ggplot() +
      geom_point(aes(x = value, y = log10(rare_reads/7000))) +
      stat_smooth(aes(x = value, y = log10(rare_reads/7000)),
                  method = "loess") +
      theme_bw() +
      my_theme +
      theme(strip.text = element_text(size = 8)) +
      facet_wrap(~phytogroups) +
      labs(y = "log10(Proportion of reads)",
           x = capitalize(gsub("_", " ", e))) +
      ggtitle(paste0("Station = ", s))
    
    plot_name <- paste0("figures/_loess-log_over_", e, "at-station_", s, ".png")
    
    ggsave(plot = plot_a,
           filename = plot_name,
           width = 14, height = 10, dpi = 300)
  })
})
})
## Max bins
lapply(stations, function(s){
  lapply(env_vars, function(e){
    lapply(p_groups, function(r){
      plot_a <- left_join(metadata_filter, asv_table_filter) %>%
        filter(phytogroups == r) %>%
        mutate(bins =  as.numeric(as.character(cut(!!sym(e),
                                                   breaks = seq(min(!!sym(e)),
                                                                max(!!sym(e)),
                                                                length.out = 30),
                                                   labels = seq(min(!!sym(e)),
                                                                max(!!sym(e)),
                                                                length.out = 30)[-30],
                                                   include.lowest = T)))) %>%
        group_by(across(all_of(c("phytogroups", "Species", "bins")))) %>%
        reframe(reads = sum(rare_reads, na.rm = T),
                prop_reads = reads/7000,
                max_prop = max(prop_reads)) %>%
        ggplot() +
        geom_point(aes(x = bins, y = log10(max_prop))) +
        stat_smooth(aes(x = bins, y = log10(max_prop)), method = "loess") +
        theme_bw() +
        my_theme +
        theme(strip.text = element_text(size = 8)) +
        facet_wrap(~Species) +
        labs(y = "log10(Max reads)", x = capitalize(gsub("_", " ", e))) +
        ggtitle(paste0("Phytoplankton group: ", r))
      
      plot_name <- paste0("figures/phyto_group_", r,
                          "_loess-log-max_over_", e, ".png")
      
      ggsave(plot = plot_a,
             filename = plot_name,
             width = 14, height = 10, dpi = 300)
    })
  })
})

group_plot <- "All_groups"
pca_meta <- metadata_filter %>%
  mutate(month = as.numeric(month)) %>%
  select(sample, all_of(env_vars)) %>%
  select(-c(Meltwater_Fraction_, month))

pca_asv <- asv_table_filter %>%
  select(sample, phytogroups, rare_reads) %>%
  pivot_wider(id_cols = "sample", names_from = "phytogroups", values_from = "rare_reads", values_fn = function(r){log10(sum(r, na.rm = T))}) 

pca_df <- left_join(pca_meta, pca_asv) %>%
  mutate_at(vars(-sample), as.numeric) %>%
  column_to_rownames("sample")
pca_df <- scale(pca_df)
# require("missMDA")
# require("factoextra")
nb <- estim_ncpPCA(pca_df, ncp.max = 5)
pca_res <- imputePCA(pca_df, ncp = nb$ncp)
res.pca <- prcomp(pca_res$completeObs)
mvdf_pca <- res.pca
table_name <- paste0("figures/pca_summary-", group_plot, ".csv")
write.csv(summary(mvdf_pca)$importance, file = table_name)
plot_a <- fviz_contrib(mvdf_pca, choice = "var", axes = 1:2) +
  ggtitle(paste0("Phytoplankton group: ", group_plot))
plot_name <- paste0("figures/contrib_plot_", group_plot, ".png")
ggsave(filename = plot_name,
       plot = plot_a,
       dpi = 300, 
       width = 5,
       height = 5,
       units = "in")
plot_b <- fviz_pca_var(mvdf_pca, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = T) +
  ggtitle(paste0("Phytoplankton group: ", group_plot))
plot_name <- paste0("figures/pca_biplot_", group_plot, ".png")
ggsave(filename = plot_name,
       plot = plot_b,
       dpi = 300, 
       width = 7,
       height = 7,
       units = "in")


lapply(stations, function(e){
  pca_meta <- metadata_filter %>%
    filter(location == e) %>% 
    mutate(month = as.numeric(month)) %>%
    select(sample, all_of(env_vars)) %>%
    drop_na() %>%
    column_to_rownames("sample")
  
  pca_out <- prcomp(pca_meta, scale. = T)
  summary(pca_out)
  
  pca_arrows <- as.data.frame(pca_out$rotation) %>%
    rownames_to_column("var")
  pca_arrows$x_start <- 0
  pca_arrows$y_start <- 0
  
  pca_add <- cbind(pca_meta, pca_out$x) %>%
    rownames_to_column("sample")
  
  left_join(pca_add, asv_table_filter) %>%
    ggplot() +
    geom_point(aes(x = PC1, y = PC2))+
    geom_segment(data = pca_arrows,
                 aes(x = x_start, y = y_start, xend = PC1*4, yend = PC2*4,
                     color = var),
                 arrow = arrow()) +
    ggtitle(paste0("Location = ", e))
})

library(vegan)

cca_test <- cca(pca_asv, pca_meta)
summary(cca_test)

veg_1 = as.data.frame(cca_test$CCA$biplot)
veg_1["env"] = row.names(veg_1)

veg_2 = as.data.frame(cca_test$CCA$v)
veg_2["species"] = row.names(veg_2)

ggplot() +
  geom_point(data = veg_1, aes(x = CCA1, y = CCA2), color = "red") +
  geom_point(data =
               veg_2, aes(x = CCA1, y = CCA2), color = "blue") +
  geom_text_repel(data = veg_1,
                  aes(x = CCA1, y = CCA2, label = veg_1$env),
                  nudge_y = -0.05)

pca_asv <- asv_table_filter %>%
  select(sample, Species, rare_reads) %>%
  pivot_wider(id_cols = "sample", names_from = "Species", values_from = "rare_reads", values_fn = sum, values_fill = 0) %>%
  filter(sample %in% rownames(pca_meta)) %>%
  column_to_rownames("sample")

pca_out_asv <- prcomp(pca_asv, scale. = T)
summary(pca_out_asv)

pca_arrows <- as.data.frame(pca_out$rotation) %>%
  rownames_to_column("var")
pca_arrows$x_start <- 0
pca_arrows$y_start <- 0


pca_add <- cbind(pca_meta, pca_out$x) %>%
  rownames_to_column("sample")

left_join(pca_add, asv_table_filter) %>%
  ggplot() +
  geom_point(aes(x = PC1, y = PC2))+
  geom_segment(data = pca_arrows,
               aes(x = x_start, y = y_start, xend = PC1*4, yend = PC2*4,
                   color = var),
               arrow = arrow())


install.packages("ca")
require(ca)
left_join(metadata_filter, asv_table_filter) %>%
  select(sample, Feature.ID, rare_reads, all_of(env_vars)) %>%
  pivot_wider(id_cols = c(Feature.ID, rare_reads))

mytable <- with(mydata, table(metadata_filter,asv_table_filter))


max_e_df_reduce <-Reduce(function(x, y) merge(x, y, by = c("phytogroups","Genus", "Species")), max_e_df)


max_e_df <- lapply(env_vars, function(e){
  do.call("rbind", lapply(p_groups, function(r){
    left_join(metadata_filter, asv_table_filter) %>%
      mutate(prop_reads = rare_reads/7000) %>%
      filter(phytogroups == r) %>%
      mutate(bins =  as.numeric(
        as.character(cut(as.numeric(!!sym(e)),
                         breaks = seq(min(as.numeric(!!sym(e)), na.rm = T),
                                      max(as.numeric(!!sym(e)), na.rm = T),
                                      length.out = 30),
                         labels = seq(min(as.numeric(!!sym(e)), na.rm = T),
                                      max(as.numeric(!!sym(e)), na.rm = T),
                                      length.out = 30)[-30],
                         include.lowest = T)))) %>%
      group_by(across(all_of(c(
        "phytogroups", "Genus", "Species", "bins")))) %>%
      reframe(max_prop = max(prop_reads)) %>%
      group_by(phytogroups, Species) %>%
      filter(max_prop == max(max_prop)) %>%
      rename_with(~paste0("max_reads_", e), "max_prop") %>%
      rename_with(~paste0("bins_", e), "bins") %>% View()
    select(-c(prop_reads))
  }))
})


max_e_df <- left_join(metadata_filter, asv_table_filter) %>%
  filter(location %in% stations) %>%
  mutate(prop_reads = rare_reads/7000) %>%
  mutate_at(env_vars, ~as.numeric(
    as.character(cut(as.numeric(.),
                     breaks = seq(min(as.numeric(.), na.rm = T),
                                  max(as.numeric(.), na.rm = T),
                                  length.out = 50),
                     labels = seq(min(as.numeric(.), na.rm = T),
                                  max(as.numeric(.), na.rm = T),
                                  length.out = 50)[-50],
                     include.lowest = T)))) %>%
  group_by(across(all_of(c(
    "phytogroups", "Genus", "Species", "location", env_vars)))) %>%
  reframe(max_prop = max(prop_reads)) %>%
  group_by(phytogroups, Species, location) %>%
  filter(max_prop == max(max_prop))

lapply(env_vars[-4], function(e){
  lapply(p_groups, function(r){
    plot_a <- max_e_df %>%
      filter(phytogroups == r) %>%
      mutate(value = !!sym(e)) %>%
      ggplot(data = ., aes(x = days_since,
                           y = value,
                           color = Genus)) +
      geom_point() +
      theme_bw() +
      my_theme +
      theme(strip.text = element_text(size = 8)) +
      scale_x_continuous(name = "Days since start of season") +
      scale_y_continuous(name = capitalize(gsub("_", " ", e))) +
      labs(color = "Genus") +
      ggtitle(paste0("Phytoplankton group: ", r))
    
    ggsave(plot = plot_a,
           filename = paste0("figures/max_bin_",e,"over_days-color_",r,".png"),
           width = 10, height = 6, dpi = 300)
  })
})

max_e_df %>%
  ggplot(data = ., aes(x = days_since,
                       y = ,
                       color = log10(max_prop))) +
  geom_point() +
  theme_bw() +
  my_theme +
  theme(strip.text = element_text(size = 8)) +
  scale_x_continuous(name = "Days since start of season")

left_join(metadata_filter, asv_table_filter) %>%
  group_by(across(all_of(c("sample", "Species", env_vars)))) %>%
  reframe(reads = sum(rare_reads, na.rm = T),
          prop_reads = reads/7000) %>%
  filter(Species %in% unique(top_15$Species)) %>%
  mutate(across(c(env_vars),
                ~as.numeric(as.character(
                  cut(.x,
                      breaks = seq(min(.x),
                                   max(.x),
                                   length.out = 15),
                      labels = seq(min(.x),
                                   max(.x),
                                   length.out = 15)[-15],
                      include.lowest = T))),
                .names = "{.col}bin")) %>% 
  group_by(across(all_of(c("Species", paste0(env_vars, "bin"))))) %>%
  reframe(max_reads = max(reads)) %>%
  mutate(days_since_bin = case_when(days_sincebin <= 50 ~ "Early", 
                                    days_sincebin <= 100 &
                                      days_sincebin > 50 ~ "Mid",
                                    TRUE ~ "Late")) %>%
  ggplot(data = ., aes(x = Surface_Temperature_bin,
                       y = Surface_Salinity_bin,
                       color = days_since_bin,
                       size = max_reads)) +
  geom_point() +
  theme_bw() +
  my_theme +
  theme(strip.text = element_text(size = 8)) +
  facet_wrap(~Species) +
  labs(x = "Temperature", y = "Salinity", color = "Season phase",
       size = "Max reads") +
  scale_color_manual(values = c(
    "Early" = "#A0DDFF", "Mid" = "#EF8354", "Late" = "#624CAB"),
    labels = c("Early (< 50 days)", "Mid (51 - 100 days)",
               "Late (> 100 days)"))



### Top 15 species
top_15 <- asv_table_filter %>%
  select(Species, rare_reads) %>%
  mutate(total_reads = sum(rare_reads, na.rm = T)) %>%
  group_by(Species) %>%
  reframe(sp_rel = sum(rare_reads, na.rm = T)/total_reads) %>%
  distinct(Species, .keep_all = T) %>%
  slice_max(sp_rel, n = 15)

#### Fit a density estimate to each species within each location
r = "Surface_Salinity_"

full_bins <- left_join(metadata_filter, asv_table_filter) %>%
  filter(location %in% stations) %>%
  mutate(bins = cut(!!sym(r),
                    breaks = seq(min(!!sym(r)),
                                 max(!!sym(r)),
                                 by = 0.01),
                    labels = seq(min(!!sym(r)),
                                 max(!!sym(r)),
                                 by = 0.01)[-15],
                    include.lowest = T)) %>%
  pull(bins) %>% levels() %>% as.character() %>% as.numeric()

species.loess.df <- left_join(metadata_filter, asv_table_filter) %>%
  filter(location %in% stations) %>%
  mutate(bins = as.numeric(
    as.character(cut(!!sym(r),
                     breaks = seq(min(!!sym(r)),
                                  max(!!sym(r)),
                                  length.out = 15),
                     labels = seq(min(!!sym(r)),
                                  max(!!sym(r)),
                                  length.out = 15)[-15],
                     include.lowest = T)))) %>%
  group_by(location) %>%
  mutate(loc_sum = sum(rare_reads)) %>%
  group_by(phytogroups, Genus, Species, location, bins) %>%
  reframe(sp_rel = sum(rare_reads)/loc_sum) %>%
  nest(data = c(bins, sp_rel)) %>%
  mutate(
    model = map(data, ~ {
      if(length(unique(.$bins[.$sp_rel >0])) < 5 | nrow(.) == 0){
        NULL
      }else{
        .x %>%
          loess(sp_rel ~ bins, ., span = 0.75, degree = 2)
      }}),
    fit = map(model, ~{
      if(is.null(.)){
        vector(length = length(full_bins))
      }else{
        predict(., full_bins)
      }
    }%>% as_tibble() %>% mutate(bins = full_bins, value = ifelse(value<0,0,value))),
    bin_max = map(fit, ~{.$bins[which.max(.$value)[1]]}),
    bin_sd = map(fit, ~{
      sqrt(
        sum(.$value*(.$bins - .$bins[which.max(.$value)[1]])^2, na.rm = T)/
          (((length(.$value>0)-1)/length(.$value))*sum(.$value, na.rm = T)))
    })
  ) %>%
  unnest(bin_max, bin_sd) %>%
  mutate(bin_max = ifelse(is.na(bin_sd), NA, bin_max)) %>%
  ggplot() +
  geom_point(aes(x = location, y = Species, color = bin_max, size = bin_sd)) +
  scale_color_gradientn(name = r, colors = c("purple", "orange")) +
  facet_grid(phytogroups~., scales = "free_y", space = "free_y")




testing <- left_join(metadata_filter, asv_table_filter) %>%
  filter(location == "Petermann Island") %>%
  ggplot() +
  geom_density(aes(x = Surface_Salinity_, weights = log10(rare_reads))) +
  geom_point(aes(x = Surface_Salinity_, y = log10(rare_reads)))+
  facet_wrap(~Species)


#### by location: bins, max in bins, loess to max, max of loess, sd of loess
species.bin.max <- left_join(metadata_filter, asv_table_filter)  %>%
  #filter(location %in% stations) %>%
  mutate(sp_rel = sum(reads)/7000,
         bins =  as.numeric(as.character(cut(!!sym(r),
                                             breaks = seq(min(!!sym(r)),
                                                          max(!!sym(r)),
                                                          length.out = 15),
                                             labels = seq(min(!!sym(r)),
                                                          max(!!sym(r)),
                                                          length.out = 15)[-15],
                                             include.lowest = T)))) %>%
  group_by(phytogroups, Species, location, bins) %>%
  reframe(max_in_bin = max(sp_rel)) %>%
  group_by(phytogroups, Species, location, bins) %>%
  select(all_of(c("bins", "max_in_bin"))) %>%
  nest(data = c(bins, max_in_bin)) %>%
  mutate(
    model = map(data, ~ {
      if(length(unique(.$bins[.$max_in_bin >0])) < 5 | nrow(.) == 0){
        NULL
      }else{
        .x %>%
          loess(max_in_bin ~ bins, ., span = 0.75, degree = 2)
      }}),
    fit = map(model, ~{
      if(is.null(.)){
        vector(length = length(1:140))
      }else{
        predict(., 1:140)
      }
    }%>% as_tibble() %>% mutate(bins = 1:140, value = ifelse(value<0,0,value))),
    bin_max = map(fit, ~{.$bins[which.max(.$value)[1]]}),
    bin_sd = map(fit, ~{
      sqrt(
        sum(.$value*(.$bins - .$bins[which.max(.$value)[1]])^2, na.rm = T)/
          (((length(.$value>0)-1)/length(.$value))*sum(.$value, na.rm = T)))
    })
  ) %>%
  unnest(bin_max, bin_sd)


species.loc.bin.df <- left_join(metadata_filter, asv_table_filter) %>%
  group_by(across(all_of(
    c("location", "phytogroups", "Species", env_vars)))) %>%
  reframe(reads = sum(rare_reads, na.rm = T),
          prop_reads = reads/7000) %>%
  mutate(bins =  cut(!!sym(r),
                     breaks = seq(min(!!sym(r)),
                                  max(!!sym(r)),
                                  length.out = 15),
                     labels = seq(min(!!sym(r)),
                                  max(!!sym(r)),
                                  length.out = 15)[-15],
                     include.lowest = T)) %>%
  group_by(location, phytogroups, Species, bins) %>%
  reframe(max_prop = max(prop_reads)) %>%
  group_by(location) %>%
  mutate(n = n()) %>%
  filter(n > 50) %>% pull(location) %>% unique()

full_bins <- as.numeric(as.character(levels(species.loc.bin.df$bins)))

species.loc.bin.df %>%
  mutate(bins = as.numeric(as.character(bins))) %>%
  group_by(phytogroups, Species, location) %>%
  nest(data = c(bins, max_prop)) %>%
  mutate(
    model = map(data, ~ {
      if(length(unique(.$bins[.$max_prop >0])) < 5 | nrow(.) == 0){
        NULL
      }else{
        .x %>%
          loess(max_prop ~ bins, ., span = 0.75, degree = 2)
      }
    }),
    fit = map(model, ~{
      if(is.null(.)){
        vector(length = length(full_bins))
      }else{
        predict(., full_bins)
      }
    }%>%
      as_tibble() %>%
      mutate(bins = full_bins,
             value = ifelse(value<0,0,value))),
    bin_max = map(fit, ~{.$bins[which.max(.$value)[1]]}),
    bin_sd = map(fit, ~{
      sqrt(
        sum(.$value*(.$bins - .$bins[
          which.max(.$value)[1]])^2, na.rm = T)/
          (((length(.$value>0)-1)/length(.$value))*sum(.$value, na.rm = T)))
    })) %>% 
  unnest(bin_max, bin_sd) %>%
  mutate(bin_group = case_when(bin_max <= 32 ~ "Low",
                               bin_max <= 33.5 &
                                 bin_max > 32 ~ "Mid",
                               TRUE ~ "High")) %>%
  ggplot() +
  geom_point(aes(x = location, y = Species, color = bin_group))


## For each top 15 species, plot the maximum proportion of reads within a window of 10 days
env_names <- c("Surface_Temperature_", "Surface_Salinity_", "days_since")

lapply(env_names, function(r){
  plot_a <- left_join(metadata_filter, asv_table_filter) %>%
    group_by(across(all_of(c("sample", "Genus", env_vars)))) %>%
    reframe(reads = sum(rare_reads, na.rm = T),
            prop_reads = reads/7000) %>%
    #filter(Species %in% unique(top_15$Species)) %>%
    mutate(bins =  cut(!!sym(r),
                       breaks = seq(min(!!sym(r)),
                                    max(!!sym(r)),
                                    length.out = 15),
                       labels = seq(min(!!sym(r)),
                                    max(!!sym(r)),
                                    length.out = 15)[-15],
                       include.lowest = T)) %>%
    group_by(Genus, bins) %>%
    reframe(max_prop = max(prop_reads)) %>%
    mutate(bins = as.numeric(as.character(bins))) %>%
    ggplot(data = ., aes(x = bins, y = max_prop)) +
    geom_point() +
    geom_smooth(, method = "loess") +
    theme_bw() +
    my_theme +
    theme(strip.text = element_text(size = 8)) +
    facet_wrap(~Genus) +
    labs(y = "Max proportion of reads", x = capitalize(gsub("_", " ", r)))
  
  plot_name <- paste0("figures/genus_max-loess_over_", r, ".png")
  ggsave(plot = plot_a,
         filename = plot_name,
         width = 8, height = 8, dpi = 300)
})

left_join(metadata_filter, asv_table_filter) %>%
  group_by(across(all_of(c("sample", "Genus", env_vars)))) %>%
  reframe(reads = sum(rare_reads, na.rm = T),
          prop_reads = reads/7000) %>%
  group_by(across(all_of(c("Genus", env_vars)))) %>%
  reframe(max_prop = max(prop_reads)) %>%
  mutate(sal_bins = case_when(Surface_Salinity_ <= 32 ~ "Low",
                              Surface_Salinity_ <= 33.5 &
                                Surface_Salinity_ > 32 ~ "Mid",
                              TRUE ~ "High"),
         temp_bins = case_when(Surface_Temperature_ <= 0 ~ "Low",
                               Surface_Temperature_ <= 2 &
                                 Surface_Temperature_ > 0 ~ "Mid",
                               TRUE ~ "High")) %>%
  ggplot(data = .) +
  geom_point(aes(x = days_since, y = log10(max_prop), color = sal_bins)) +
  facet_wrap(~Genus)



left_join(metadata_filter, asv_table_filter) %>%
  group_by(across(all_of(c("sample", "Genus", env_vars)))) %>%
  reframe(reads = sum(rare_reads, na.rm = T),
          prop_reads = reads/7000) %>%
  mutate(across(c(env_vars), ~as.numeric(as.character(cut(.x,
                                                          breaks = seq(min(.x),
                                                                       max(.x),
                                                                       length.out = 15),
                                                          labels = seq(min(.x),
                                                                       max(.x),
                                                                       length.out = 15)[-15],
                                                          include.lowest = T))),
                .names = "{.col}bin")) %>% 
  group_by(across(all_of(c("Genus", paste0(env_vars, "bin"))))) %>%
  reframe(max_prop = max(prop_reads)) %>%
  mutate(sal_bins = case_when(Surface_Salinity_bin <= 32 ~ "Low",
                              Surface_Salinity_bin <= 33.5 &
                                Surface_Salinity_bin > 32 ~ "Mid",
                              TRUE ~ "High"),
         temp_bins = case_when(Surface_Temperature_bin <= 0 ~ "Low",
                               Surface_Temperature_bin <= 2 &
                                 Surface_Temperature_bin > 0 ~ "Mid",
                               TRUE ~ "High")) %>%
  ggplot(data = .) +
  geom_point(aes(x = days_sincebin, y = log10(max_prop), color = sal_bins)) +
  facet_wrap(~Genus)



left_join(metadata_filter, asv_table_filter) %>%
  group_by(across(all_of(c("sample", "Species", env_vars)))) %>%
  reframe(reads = sum(rare_reads, na.rm = T),
          prop_reads = reads/7000) %>%
  filter(Species %in% unique(top_15$Species)[4]) %>%
  mutate(across(c(env_vars), ~as.numeric(as.character(cut(.x,
                                                          breaks = seq(min(.x),
                                                                       max(.x),
                                                                       length.out = 15),
                                                          labels = seq(min(.x),
                                                                       max(.x),
                                                                       length.out = 15)[-15],
                                                          include.lowest = T))),
                .names = "{.col}bin")) %>% 
  group_by(across(all_of(c("Species", paste0(env_vars, "bin"))))) %>%
  reframe(max_prop = max(prop_reads)) %>%
  plot_ly(., x= ~days_sincebin,
          y= ~Surface_Temperature_bin, z= ~Surface_Salinity_bin,
          type="scatter3d", mode="markers", color = ~max_prop,
          colors= colorRamp(c("purple","orange")))

max_3d_plot<- left_join(metadata_filter, asv_table_filter) %>%
  group_by(across(all_of(c("sample", "Species", env_vars)))) %>%
  reframe(reads = sum(rare_reads, na.rm = T),
          prop_reads = reads/7000) %>%
  filter(Species %in% unique(top_15$Species)) %>%
  mutate(across(c(env_vars),
                ~as.numeric(as.character(
                  cut(.x,
                      breaks = seq(min(.x),
                                   max(.x),
                                   length.out = 15),
                      labels = seq(min(.x),
                                   max(.x),
                                   length.out = 15)[-15],
                      include.lowest = T))),
                .names = "{.col}bin")) %>% 
  group_by(across(all_of(c("Species", paste0(env_vars, "bin"))))) %>%
  reframe(max_reads = max(reads)) %>%
  mutate(days_since_bin = case_when(days_sincebin <= 50 ~ "Early", 
                                    days_sincebin <= 100 &
                                      days_sincebin > 50 ~ "Mid",
                                    TRUE ~ "Late")) %>%
  ggplot(data = ., aes(x = Surface_Temperature_bin,
                       y = Surface_Salinity_bin,
                       color = days_since_bin,
                       size = max_reads)) +
  geom_point() +
  theme_bw() +
  my_theme +
  theme(strip.text = element_text(size = 8)) +
  facet_wrap(~Species) +
  labs(x = "Temperature", y = "Salinity", color = "Season phase",
       size = "Max reads") +
  scale_color_manual(values = c(
    "Early" = "#A0DDFF", "Mid" = "#EF8354", "Late" = "#624CAB"),
    labels = c("Early (< 50 days)", "Mid (51 - 100 days)",
               "Late (> 100 days)"))


ggsave(plot = max_3d_plot,
       filename = "figures/top_15-maximum_bins_3d.png",
       width = 10, height = 8, dpi = 300)

n_bins <- 15

metadata_filter %>%
  select(sample, 
         Surface_Temperature_, Surface_Salinity_, Potential_Density_) %>%
  group_by(sample) %>%
  mutate_at("Potential_Density_", ~as.numeric(mean(unlist(.x)[1:10], na.rm = T))) %>%
  distinct() %>%
  ggplot() +
  geom_contour(aes(x = Surface_Temperature_,
                   y = Surface_Salinity_,
                   fill = Potential_Density_))

max_3d_plot<- left_join(metadata_filter, asv_table_filter) %>%
  group_by(across(all_of(c("sample", "phytogroups", env_vars)))) %>%
  reframe(reads = sum(rare_reads, na.rm = T),
          prop_reads = reads/7000) %>%
  mutate(across(c(env_vars),
                ~as.numeric(as.character(
                  cut(.x,
                      breaks = seq(min(.x),
                                   max(.x),
                                   length.out = n_bins),
                      labels = seq(min(.x),
                                   max(.x),
                                   length.out = n_bins)[-n_bins],
                      include.lowest = T))),
                .names = "{.col}bin")) %>%
  group_by(across(all_of(c("phytogroups", paste0(env_vars, "bin"))))) %>%
  mutate(max_reads = max(reads)) %>%
  mutate(days_since_bin = case_when(days_sincebin <= 50 ~ "Early", 
                                    days_sincebin <= 100 &
                                      days_sincebin > 50 ~ "Mid",
                                    TRUE ~ "Late")) %>%
  ggplot(data = .) +
  geom_isopycnal(aes(x = Surface_Temperature_,
                     y = Surface_Salinity_)) +
  geom_point(aes(x = Surface_Temperature_bin,
                 y = Surface_Salinity_bin,
                 color = days_since_bin,
                 size = max_reads)) +
  theme_bw() +
  my_theme +
  theme(strip.text = element_text(size = 8)) +
  facet_wrap(~phytogroups) +
  labs(x = "Temperature", y = "Salinity", color = "Season phase",
       size = "Max reads") +
  scale_color_manual(values = c(
    "Early" = "#A0DDFF", "Mid" = "#EF8354", "Late" = "#624CAB"),
    labels = c("Early (< 50 days)", "Mid (51 - 100 days)",
               "Late (> 100 days)"))


ggsave(plot = max_3d_plot,
       filename = "figures/phytogroups-maximum_bins_3d.png",
       width = 10, height = 8, dpi = 300)






#### NMDS samples
bio <- asv_table_filter %>%
  select(Species, rare_reads, sample) %>%
  pivot_wider(.,
              id_cols = "sample", names_from = "Species", values_from = "rare_reads", values_fn = function(r)sum(r,na.rm = T), values_fill = 0) %>%
  column_to_rownames(var = "sample")

correlation_matrix <- cor(t(bio), method = "spearman")

bio_NMS1 <-  metaMDS(bio,
                     k = 2,
                     maxit = 200, 
                     trymax = 100,
                     try = 50,
                     wascores= TRUE,
                     weakties = T,
                     expand = TRUE,
                     previous.best = T,
                     autotransform = FALSE)

bio_spp.fit <- envfit(bio_NMS1, bio, permutations = 999) # this fits species vectors
spp.scrs <- as.data.frame(scores(bio_spp.fit, display = "vectors", arrow.mul=2.5))
spp.scrs <- cbind(spp.scrs, vars = rownames(spp.scrs)) #add species names to dataframe
spp.scrs <- cbind(spp.scrs, pval = bio_spp.fit$vectors$pvals) #add pvalues to dataframe so you can select species which are significant
#spp.scrs<- cbind(spp.scrs, abrev = abbreviate(spp.scrs$Species, minlength = 6)) #abbreviate species names
sig.spp.scrs <- subset(spp.scrs, pval<=0.05) #subset data to show species

env_df <- metadata_filter %>%
  select(sample, all_of(env_vars))

bio.env.fit <- envfit(bio_NMS1, env_df, permutations = 999, na.rm = T) 

env.scrs <- as.data.frame(scores(bio.env.fit, display = "vectors", arrow.mul=2.5)) 
env.scrs <- cbind(env.scrs, vars = rownames(env.scrs)) #add species names to dataframe
env.scrs <- cbind(env.scrs, pval = bio.env.fit$vectors$pvals) #add pvalues to dataframe so you can select species which are significant
#spp.scrs<- cbind(spp.scrs, abrev = abbreviate(spp.scrs$Species, minlength = 6)) #abbreviate species names
sig.env.scrs <- subset(env.scrs, pval<=0.05) #subset data to show species 


#Get the vectors for env.fit
df_envfit<-scores(bio.env.fit,display=c("vectors"))
df_envfit<-df_envfit*vegan:::ordiArrowMul(df_envfit)
df_envfit<-as.data.frame(df_envfit)

env_mds <- as.data.frame(scores(bio_NMS1$points)) %>%
  rownames_to_column("sample") %>%
  left_join(env_df)

bio_mds <- sig.spp.scrs %>%
  rownames_to_column("Species") %>%
  left_join(asv_table_filter) %>%
  distinct(Species, NMDS1, NMDS2, phytogroups)

group_colors <- c("#F06400", "#00F064","#008CF0","gold2","#6400F0","#F0008C","grey10")
names(group_colors) <- c("Cryptophytes", "Diatoms", "Dinoflagellates", "Haptophytes",
                         "MAST", "Rhodophytes", "Greenalgae")

require(ggrepel)
ggplot(bio_mds) +
  geom_point(data = env_mds,
             aes(x = MDS1, y = MDS2, fill = Surface_Salinity_), pch = 21, size = 3) +
  geom_segment(data = bio_mds,
               aes(x = 0, xend=NMDS1, y=0, yend=NMDS2, color = phytogroups),
               arrow = arrow(length = unit(0.25, "cm")), lwd = 0.3) +
  geom_text_repel(data = bio_mds,
                  aes(x=NMDS1, y=NMDS2, label = Species, color = phytogroups),
                  cex = 3, direction = "both") +
  scale_fill_gradientn(name = "Days since", colors = c("pink", "purple")) +
  scale_color_manual(name = "", values = group_colors) +
  theme_bw() +
  my_theme

lapply(env_vars, function(x){
  ggplot(bio_mds) +
    geom_point(data = env_mds,
               aes(x = MDS1, y = MDS2, fill = !!sym(x)), pch = 21, size = 3) +
    geom_segment(data = bio_mds,
                 aes(x = 0, xend=NMDS1, y=0, yend=NMDS2, color = phytogroups),
                 arrow = arrow(length = unit(0.25, "cm")), lwd = 0.3) +
    geom_text_repel(data = bio_mds,
                    aes(x=NMDS1, y=NMDS2,
                        label = Species, color = phytogroups),
                    cex = 3, direction = "both", show.legend = F) +
    scale_fill_gradientn(name = capitalize(gsub("_", " ", x)),
                         colors = c("pink", "purple")) +
    scale_color_manual(name = "", values = group_colors) +
    theme_bw() +
    my_theme
  
  plot_name <- paste0("figures/NMDS_", x, ".png")
  ggsave(plot = last_plot(),
         filename = plot_name,
         width = 8, height = 9, dpi = 300)
})


corr_plot <- rcorr(as.matrix(bio), type = "pearson")
diag(corr_plot$P) <- 0
corrplot(corr_plot$r, type = "upper", tl.col = "black", tl.srt = 45, addCoef.col = "grey30", number.cex = 0.5, tl.cex = 0.5,
         p.mat = corr_plot$P, sig.level = 0.05, insig = "blank")

bio_phyto <- bio %>%
  rownames_to_column(var = "sample") %>%
  left_join(env_df) %>%
  column_to_rownames(var = "sample") %>%
  mutate_all(~as.numeric(.x))

env_names <- colnames(env.fit_df)[-1]

corr_plot <- rcorr(as.matrix(bio_phyto), type = "pearson")
diag(corr_plot$P) <- 0
corrplot(corr_plot$r, type = "upper", tl.col = "black", tl.srt = 45,
         number.cex = 0.5, tl.cex = 0.5,
         p.mat = corr_plot$P, sig.level = 0.05, insig = "blank")

corr_species <- as.data.frame(corr_plot$r) %>%
  rownames_to_column("var2") %>%
  pivot_longer(cols = -var2, names_to = "var1", values_to = "r") %>%
  mutate(P = as.vector(corr_plot$P)) %>%
  filter(var1 != var2) %>%
  filter(var1 %in% env_vars) %>%
  filter(var2 %in% unique(colnames(bio))) %>%
  ggplot() +
  geom_point(aes(x = var1, y = var2, color = r,
                 size = if_else(P < 0.05, P, NA))) +
  scale_color_gradient2(low = "blue", high = "red", limits = c(-0.75,0.75)) +
  scale_x_discrete(labels = ~capitalize(str_replace_all(., "_", " "))) +
  scale_size_continuous(guide = "none", transform = "reverse",
                        range = c(1,5)) +
  theme_bw() +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "", color = "Correlation")

ggsave(plot = corr_species,
       filename = "figures/correlation_species_env.png",
       width = 8, height = 10, dpi = 300)

ggplot() +
  geom_point(data = env_mds,
             aes(x = MDS1, y = MDS2, fill = Surface_Salinity_),
             pch = 21, size = 3) +
  geom_segment(data = bio_mds,
               aes(x = 0, xend=NMDS1, y=0, yend=NMDS2, color = phytogroups),
               arrow = arrow(length = unit(0.25, "cm")), lwd = 0.3) +
  scale_fill_viridis_c() +
  scale_color_manual(name = "", values = group_colors)


ggplot() +
  geom_point(data = env_mds,
             aes(x = days_since, y = Surface_Temperature_, fill = Surface_Salinity_), pch = 21, size = 3) +
  # geom_segment(data = bio_mds,
  #              aes(x = 0, xend=NMDS1, y=0, yend=NMDS2, color = phytogroups),
  #              arrow = arrow(length = unit(0.25, "cm")), lwd = 0.3) +
  scale_fill_viridis_c() +
  scale_color_manual(name = "", values = group_colors)



### How is the environment changing over time and space?

metadata_filter %>%
  unnest() %>%
  ggplot() +
  geom_tile(aes(x = days_since, y = -Pressure_,
                fill = Temperature_),width = 3) +
  facet_wrap(~as.numeric(latitude), ncol = 1) +
  theme_bw() +
  my_theme +
  theme(text = element_text(size = 5))

left_join(diversity_filter,  metadata_filter) %>%
  ggplot() +
  geom_point(aes(x = latitude, y = longitude, color = richness_phytogroups)) +
  facet_wrap(~seasonyear)

### Find clusters
bio <- asv_table_filter %>%
  select(Species, rare_reads, sample) %>%
  pivot_wider(.,
              id_cols = "sample", names_from = "Species", values_from = "rare_reads", values_fn = function(r)sum(r,na.rm = T), values_fill = 0) %>%
  column_to_rownames(var = "sample")

correlation_matrix <- cor(bio, method = "spearman")
hierarchical_result <- hclust(dist(1 - correlation_matrix), method = "ward.D2")
clusters <- cutree(hierarchical_result, k = 3)

bio_NMS1 <-  metaMDS(dist(1 - correlation_matrix),
                     k = 2,
                     maxit = 200, 
                     trymax = 100,
                     try = 50,
                     wascores= TRUE,
                     weakties = T,
                     expand = TRUE,
                     previous.best = T,
                     autotransform = FALSE)

data_scores_cluster <- as.data.frame(scores(bio_NMS1)) %>%
  rownames_to_column("Species") %>%
  mutate(clusters = factor(
    Species, levels = names(clusters), labels = clusters)) %>%
  left_join(asv_table_filter %>% select(Species, phytogroups, sample, reads) %>% distinct(), by = "Species")

env.fit_df <- metadata_filter %>%
  select_if(Negate(is.list)) %>%
  select(sample, Surface_Salinity_, Surface_Temperature_, secchi_depth, latitude, longitude, days_since, month, Meltwater_Fraction_, land_dist_km)


### Correlations
bio_phyto <- asv_table_filter %>%
  select(Species, phytogroups, rare_reads, sample) %>%
  mutate(cluster = factor(Species,
                          levels = names(clusters), labels = clusters),
         ID = paste(cluster, phytogroups, sep = "_")) %>%
  pivot_wider(.,
              id_cols = "sample", names_from = "ID", values_from = "rare_reads", values_fn = function(r)sum(r,na.rm = T), values_fill = 0)

bio_names <- colnames(bio_phyto)[-1]
bio_phyto <- bio_phyto %>%
  left_join(env.fit_df) %>%
  column_to_rownames(var = "sample") %>%
  mutate_all(~as.numeric(.x))

env_names <- colnames(env.fit_df)[-1]

corr_plot <- rcorr(as.matrix(bio_phyto), type = "pearson")
diag(corr_plot$P) <- 0
corrplot(corr_plot$r, type = "upper", tl.col = "black", tl.srt = 45, addCoef.col = "grey30", number.cex = 0.5, tl.cex = 0.5,
         p.mat = corr_plot$P, sig.level = 0.05, insig = "blank")


as.data.frame(corr_plot$r) %>%
  rownames_to_column("var2") %>%
  pivot_longer(cols = -var2, names_to = "var1", values_to = "r") %>%
  mutate(P = as.vector(corr_plot$P)) %>%
  filter(var1 != var2) %>%
  filter(var1 %in% env_names) %>%
  filter(var2 %in% bio_names) %>%
  filter(P < 0.05) %>%
  separate(var2, into = c("cluster", "phytogroup"), sep = "_") %>%
  ggplot() +
  geom_point(aes(x = var1, y = phytogroup, color = r, size = P)) +
  scale_color_gradient2(low = "blue", high = "red", limits = c(-0.75,0.75)) +
  scale_x_discrete(labels = ~capitalize(str_replace_all(., "_", " "))) +
  scale_size_continuous(guide = "none", transform = "reverse",
                        range = c(1,5)) +
  theme_bw() +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "", color = "Correlation") +
  facet_wrap(~cluster)


bio_phyto %>%
  rownames_to_column("sample") %>%
  pivot_longer(cols = -all_of(c("sample", env_names)), names_to = "ID", values_to = "rare_reads") %>%
  separate(ID, into = c("cluster", "phytogroup"), sep = "_") %>%
  ggplot() +
  geom_point(aes(x = days_since, y = rare_reads, color = phytogroup)) +
  facet_grid(phytogroup~cluster)

bio_phyto %>%
  rownames_to_column("sample") %>%
  pivot_longer(cols = -all_of(c("sample", env_names)), names_to = "ID", values_to = "rare_reads") %>%
  separate(ID, into = c("cluster", "phytogroup"), sep = "_") %>%
  filter(cluster == 1) %>%
  lmer(rare_reads~ Surface_Temperature_ + Surface_Salinity_ + Meltwater_Fraction_ + (1|month) + (1|region))



#source("analysis/ggnested_pattern.R")
devtools::install_github("gmteunisse/ggnested")
install.packages("devtools")
install.packages("ggpattern")
require(ggpattern)

asv_table_filter %>%
  select(Species, phytogroups, rare_reads, sample) %>%
  mutate(cluster = factor(Species,
                          levels = names(clusters), labels = clusters)) %>%
  group_by(cluster) %>%
  mutate(clust_sum = sum(rare_reads, na.rm = T)) %>%
  group_by(Species, cluster, phytogroups) %>%
  reframe(sp_rel = sum(rare_reads, na.rm = T)/clust_sum) %>%
  distinct(Species, cluster, phytogroups, .keep_all = T) %>%
  group_by(cluster) %>%
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
    pat = case_when(phytogroups == "Diatoms" & cluster == "2"~"wave",
                    phytogroups == "Diatoms" & cluster == "3"~"stripe",
                    TRUE~"none"),
    pat = if_else(Species == "Porosira_sp.", "circle", pat)) %>%
  ggnested_pattern(data = ., aes_string(main_group = "phytogroups",
                                        sub_group = "Species", 
                                        x = "cluster", y = "sp_rel", pattern = "pat"),
                   main_palette = group_colors) + 
  geom_col_pattern(pattern_fill="white", pattern_color="white", pattern_key_scale_factor=0.25,
                   pattern_size = 0.1) +
  scale_y_continuous(name = "Relative abundance",
                     expand = c(0,0)) +
  scale_x_discrete(name = "Cluster", expand = c(1E-3,1E-3), breaks = c(1,2,3)) +
  theme_nested(theme_linedraw) + 
  theme(
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(size = 12),
    legend.title = element_blank(),
    panel.border = element_rect(color = "black")) +
  guides(fill=guide_legend(ncol = 1))


##### Analysis for Ch. 2 (now 3)

#source()
my_theme = theme_linedraw() + theme(text = element_text(size = 14), strip.background = element_blank(), strip.text = element_text(face = "bold", color = "black"))

# diversity_df, asv_table_rare, metadata_filter
metadata_filter <- metadata_w_samples %>%
  drop_na(secchi_depth) %>%
  rename(sample = Genetics_18sv9_Sample_ID) %>%
  filter(!sample %in% samples_removed)

asv_filter <- asv_table_rare %>%
  drop_na(phytogroups) %>%
  filter(sample %in% unique(metadata_filter$sample)) %>%
  group_by(Feature.ID) %>%
  reframe(total_reads = sum(rare_reads, na.rm = T),
          total_samples = length(unique(sample[rare_reads > 0]))) %>%
  filter(total_reads > 5 & total_samples > 5) %>%
  pull(Feature.ID) %>%
  unique()

asv_table_filter <- asv_table_rare %>%
  filter(rare_reads > 0) %>%
  drop_na(phytogroups) %>%
  filter(sample %in% unique(metadata_filter$sample)) %>%
  filter(Feature.ID %in% asv_filter)

diversity_filter <- diversity_df %>%
  filter(sample %in% unique(metadata_filter$sample))

#### Driving questions: How does salinity and temperature influence the microbial community?

## Pull the 9 most sampled locations and filter the data

stations <- c("Cierva Cove","Cuverville Island", "Danco Island",
              "Neko Harbour", "Orne Harbour", "Paradise Harbour", "Petermann Island", "Wilhelmina Bay")
env_vars <- c("Surface_Temperature_", "Surface_Salinity_", "month", "Meltwater_Fraction_", "days_since", "MLD_20", "secchi_depth")
p_groups <- unique(asv_table_filter$phytogroups)

## Raw values 
lapply(stations, function(s){
  lapply(env_vars, function(e){
    #lapply(p_groups, function(r){
      plot_a <- left_join(metadata_filter, asv_table_filter) %>%
        filter(location == s) %>%
        mutate(value = !!sym(e)) %>%
        ggplot() +
        geom_point(aes(x = value, y = log10(rare_reads/7000))) +
        stat_smooth(aes(x = value, y = log10(rare_reads/7000)),
                    method = "loess") +
        theme_bw() +
        my_theme +
        theme(strip.text = element_text(size = 8)) +
        facet_wrap(~phytogroups) +
        labs(y = "log10(Proportion of reads)",
             x = capitalize(gsub("_", " ", e))) +
        ggtitle(paste0("Station = ", s))
      
      plot_name <- paste0("figures/_loess-log_over_", e, "at-station_", s, ".png")
      
      ggsave(plot = plot_a,
             filename = plot_name,
             width = 14, height = 10, dpi = 300)
    })
  })
#})
## Max bins
lapply(stations, function(s){
  lapply(env_vars, function(e){
    lapply(p_groups, function(r){
      plot_a <- left_join(metadata_filter, asv_table_filter) %>%
        filter(phytogroups == r) %>%
        mutate(bins =  as.numeric(as.character(cut(!!sym(e),
                     breaks = seq(min(!!sym(e)),
                                  max(!!sym(e)),
                                  length.out = 30),
                     labels = seq(min(!!sym(e)),
                                  max(!!sym(e)),
                                  length.out = 30)[-30],
                     include.lowest = T)))) %>%
        group_by(across(all_of(c("phytogroups", "Species", "bins")))) %>%
        reframe(reads = sum(rare_reads, na.rm = T),
                prop_reads = reads/7000,
                max_prop = max(prop_reads)) %>%
        ggplot() +
        geom_point(aes(x = bins, y = log10(max_prop))) +
        stat_smooth(aes(x = bins, y = log10(max_prop)), method = "loess") +
        theme_bw() +
        my_theme +
        theme(strip.text = element_text(size = 8)) +
        facet_wrap(~Species) +
        labs(y = "log10(Max reads)", x = capitalize(gsub("_", " ", e))) +
        ggtitle(paste0("Phytoplankton group: ", r))
      
      plot_name <- paste0("figures/phyto_group_", r,
                          "_loess-log-max_over_", e, ".png")
      
      ggsave(plot = plot_a,
             filename = plot_name,
             width = 14, height = 10, dpi = 300)
    })
  })
})

group_plot <- "All_groups"
pca_meta <- metadata_filter %>%
  mutate(month = as.numeric(month)) %>%
  select(sample, all_of(env_vars)) %>%
  select(-c(Meltwater_Fraction_, month))

pca_asv <- asv_table_filter %>%
  select(sample, phytogroups, rare_reads) %>%
  pivot_wider(id_cols = "sample", names_from = "phytogroups", values_from = "rare_reads", values_fn = function(r){log10(sum(r, na.rm = T))}) 

pca_df <- left_join(pca_meta, pca_asv) %>%
  mutate_at(vars(-sample), as.numeric) %>%
  column_to_rownames("sample")
pca_df <- scale(pca_df)
# require("missMDA")
# require("factoextra")
nb <- estim_ncpPCA(pca_df, ncp.max = 5)
pca_res <- imputePCA(pca_df, ncp = nb$ncp)
res.pca <- prcomp(pca_res$completeObs)
mvdf_pca <- res.pca
table_name <- paste0("figures/pca_summary-", group_plot, ".csv")
write.csv(summary(mvdf_pca)$importance, file = table_name)
plot_a <- fviz_contrib(mvdf_pca, choice = "var", axes = 1:2) +
  ggtitle(paste0("Phytoplankton group: ", group_plot))
plot_name <- paste0("figures/contrib_plot_", group_plot, ".png")
ggsave(filename = plot_name,
       plot = plot_a,
       dpi = 300, 
       width = 5,
       height = 5,
       units = "in")
plot_b <- fviz_pca_var(mvdf_pca, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = T) +
  ggtitle(paste0("Phytoplankton group: ", group_plot))
plot_name <- paste0("figures/pca_biplot_", group_plot, ".png")
ggsave(filename = plot_name,
       plot = plot_b,
       dpi = 300, 
       width = 7,
       height = 7,
       units = "in")


lapply(stations, function(e){
pca_meta <- metadata_filter %>%
  filter(location == e) %>% 
  mutate(month = as.numeric(month)) %>%
  select(sample, all_of(env_vars)) %>%
  drop_na() %>%
  column_to_rownames("sample")

pca_out <- prcomp(pca_meta, scale. = T)
summary(pca_out)

pca_arrows <- as.data.frame(pca_out$rotation) %>%
  rownames_to_column("var")
pca_arrows$x_start <- 0
pca_arrows$y_start <- 0

pca_add <- cbind(pca_meta, pca_out$x) %>%
  rownames_to_column("sample")

left_join(pca_add, asv_table_filter) %>%
ggplot() +
  geom_point(aes(x = PC1, y = PC2))+
  geom_segment(data = pca_arrows,
               aes(x = x_start, y = y_start, xend = PC1*4, yend = PC2*4,
                   color = var),
               arrow = arrow()) +
  ggtitle(paste0("Location = ", e))
})

library(vegan)

cca_test <- cca(pca_asv, pca_meta)
summary(cca_test)

veg_1 = as.data.frame(cca_test$CCA$biplot)
veg_1["env"] = row.names(veg_1)

veg_2 = as.data.frame(cca_test$CCA$v)
veg_2["species"] = row.names(veg_2)

ggplot() +
  geom_point(data = veg_1, aes(x = CCA1, y = CCA2), color = "red") +
  geom_point(data =
               veg_2, aes(x = CCA1, y = CCA2), color = "blue") +
  geom_text_repel(data = veg_1,
                  aes(x = CCA1, y = CCA2, label = veg_1$env),
                  nudge_y = -0.05)

pca_asv <- asv_table_filter %>%
  select(sample, Species, rare_reads) %>%
  pivot_wider(id_cols = "sample", names_from = "Species", values_from = "rare_reads", values_fn = sum, values_fill = 0) %>%
  filter(sample %in% rownames(pca_meta)) %>%
  column_to_rownames("sample")

pca_out_asv <- prcomp(pca_asv, scale. = T)
summary(pca_out_asv)

pca_arrows <- as.data.frame(pca_out$rotation) %>%
  rownames_to_column("var")
pca_arrows$x_start <- 0
pca_arrows$y_start <- 0


pca_add <- cbind(pca_meta, pca_out$x) %>%
  rownames_to_column("sample")

left_join(pca_add, asv_table_filter) %>%
  ggplot() +
  geom_point(aes(x = PC1, y = PC2))+
  geom_segment(data = pca_arrows,
               aes(x = x_start, y = y_start, xend = PC1*4, yend = PC2*4,
                   color = var),
               arrow = arrow())


install.packages("ca")
require(ca)
left_join(metadata_filter, asv_table_filter) %>%
  select(sample, Feature.ID, rare_reads, all_of(env_vars)) %>%
  pivot_wider(id_cols = c(Feature.ID, rare_reads))
  
mytable <- with(mydata, table(metadata_filter,asv_table_filter))


max_e_df_reduce <-Reduce(function(x, y) merge(x, y, by = c("phytogroups","Genus", "Species")), max_e_df)


max_e_df <- lapply(env_vars, function(e){
  do.call("rbind", lapply(p_groups, function(r){
    left_join(metadata_filter, asv_table_filter) %>%
      mutate(prop_reads = rare_reads/7000) %>%
      filter(phytogroups == r) %>%
      mutate(bins =  as.numeric(
        as.character(cut(as.numeric(!!sym(e)),
                 breaks = seq(min(as.numeric(!!sym(e)), na.rm = T),
                              max(as.numeric(!!sym(e)), na.rm = T),
                              length.out = 30),
                 labels = seq(min(as.numeric(!!sym(e)), na.rm = T),
                              max(as.numeric(!!sym(e)), na.rm = T),
                              length.out = 30)[-30],
                 include.lowest = T)))) %>%
      group_by(across(all_of(c(
        "phytogroups", "Genus", "Species", "bins")))) %>%
      reframe(max_prop = max(prop_reads)) %>%
      group_by(phytogroups, Species) %>%
      filter(max_prop == max(max_prop)) %>%
      rename_with(~paste0("max_reads_", e), "max_prop") %>%
      rename_with(~paste0("bins_", e), "bins") %>% View()
      select(-c(prop_reads))
  }))
})


max_e_df <- left_join(metadata_filter, asv_table_filter) %>%
  filter(location %in% stations) %>%
  mutate(prop_reads = rare_reads/7000) %>%
  mutate_at(env_vars, ~as.numeric(
    as.character(cut(as.numeric(.),
                     breaks = seq(min(as.numeric(.), na.rm = T),
                                  max(as.numeric(.), na.rm = T),
                                  length.out = 50),
                     labels = seq(min(as.numeric(.), na.rm = T),
                                  max(as.numeric(.), na.rm = T),
                                  length.out = 50)[-50],
                     include.lowest = T)))) %>%
  group_by(across(all_of(c(
    "phytogroups", "Genus", "Species", "location", env_vars)))) %>%
  reframe(max_prop = max(prop_reads)) %>%
  group_by(phytogroups, Species, location) %>%
  filter(max_prop == max(max_prop))

lapply(env_vars[-4], function(e){
lapply(p_groups, function(r){
plot_a <- max_e_df %>%
  filter(phytogroups == r) %>%
  mutate(value = !!sym(e)) %>%
  ggplot(data = ., aes(x = days_since,
                       y = value,
                       color = Genus)) +
  geom_point() +
  theme_bw() +
  my_theme +
  theme(strip.text = element_text(size = 8)) +
  scale_x_continuous(name = "Days since start of season") +
  scale_y_continuous(name = capitalize(gsub("_", " ", e))) +
  labs(color = "Genus") +
  ggtitle(paste0("Phytoplankton group: ", r))

ggsave(plot = plot_a,
       filename = paste0("figures/max_bin_",e,"over_days-color_",r,".png"),
       width = 10, height = 6, dpi = 300)
})
})

max_e_df %>%
  ggplot(data = ., aes(x = days_since,
                       y = ,
                       color = log10(max_prop))) +
  geom_point() +
  theme_bw() +
  my_theme +
  theme(strip.text = element_text(size = 8)) +
  scale_x_continuous(name = "Days since start of season")
  
left_join(metadata_filter, asv_table_filter) %>%
  group_by(across(all_of(c("sample", "Species", env_vars)))) %>%
  reframe(reads = sum(rare_reads, na.rm = T),
          prop_reads = reads/7000) %>%
  filter(Species %in% unique(top_15$Species)) %>%
  mutate(across(c(env_vars),
                ~as.numeric(as.character(
                  cut(.x,
                      breaks = seq(min(.x),
                                   max(.x),
                                   length.out = 15),
                      labels = seq(min(.x),
                                   max(.x),
                                   length.out = 15)[-15],
                      include.lowest = T))),
                .names = "{.col}bin")) %>% 
  group_by(across(all_of(c("Species", paste0(env_vars, "bin"))))) %>%
  reframe(max_reads = max(reads)) %>%
  mutate(days_since_bin = case_when(days_sincebin <= 50 ~ "Early", 
                                    days_sincebin <= 100 &
                                      days_sincebin > 50 ~ "Mid",
                                    TRUE ~ "Late")) %>%
  ggplot(data = ., aes(x = Surface_Temperature_bin,
                       y = Surface_Salinity_bin,
                       color = days_since_bin,
                       size = max_reads)) +
  geom_point() +
  theme_bw() +
  my_theme +
  theme(strip.text = element_text(size = 8)) +
  facet_wrap(~Species) +
  labs(x = "Temperature", y = "Salinity", color = "Season phase",
       size = "Max reads") +
  scale_color_manual(values = c(
    "Early" = "#A0DDFF", "Mid" = "#EF8354", "Late" = "#624CAB"),
    labels = c("Early (< 50 days)", "Mid (51 - 100 days)",
               "Late (> 100 days)"))
  
  

### Top 15 species
top_15 <- asv_table_filter %>%
  select(Species, rare_reads) %>%
  mutate(total_reads = sum(rare_reads, na.rm = T)) %>%
  group_by(Species) %>%
  reframe(sp_rel = sum(rare_reads, na.rm = T)/total_reads) %>%
  distinct(Species, .keep_all = T) %>%
  slice_max(sp_rel, n = 15)
 
#### Fit a density estimate to each species within each location
r = "Surface_Salinity_"

full_bins <- left_join(metadata_filter, asv_table_filter) %>%
  filter(location %in% stations) %>%
  mutate(bins = cut(!!sym(r),
                     breaks = seq(min(!!sym(r)),
                                  max(!!sym(r)),
                                  by = 0.01),
                     labels = seq(min(!!sym(r)),
                                  max(!!sym(r)),
                                  by = 0.01)[-15],
                     include.lowest = T)) %>%
  pull(bins) %>% levels() %>% as.character() %>% as.numeric()

species.loess.df <- left_join(metadata_filter, asv_table_filter) %>%
  filter(location %in% stations) %>%
  mutate(bins = as.numeric(
    as.character(cut(!!sym(r),
                     breaks = seq(min(!!sym(r)),
                                  max(!!sym(r)),
                                  length.out = 15),
                     labels = seq(min(!!sym(r)),
                                  max(!!sym(r)),
                                  length.out = 15)[-15],
                     include.lowest = T)))) %>%
  group_by(location) %>%
  mutate(loc_sum = sum(rare_reads)) %>%
  group_by(phytogroups, Genus, Species, location, bins) %>%
  reframe(sp_rel = sum(rare_reads)/loc_sum) %>%
  nest(data = c(bins, sp_rel)) %>%
  mutate(
    model = map(data, ~ {
      if(length(unique(.$bins[.$sp_rel >0])) < 5 | nrow(.) == 0){
        NULL
      }else{
        .x %>%
          loess(sp_rel ~ bins, ., span = 0.75, degree = 2)
      }}),
    fit = map(model, ~{
      if(is.null(.)){
        vector(length = length(full_bins))
      }else{
        predict(., full_bins)
      }
    }%>% as_tibble() %>% mutate(bins = full_bins, value = ifelse(value<0,0,value))),
    bin_max = map(fit, ~{.$bins[which.max(.$value)[1]]}),
    bin_sd = map(fit, ~{
      sqrt(
        sum(.$value*(.$bins - .$bins[which.max(.$value)[1]])^2, na.rm = T)/
          (((length(.$value>0)-1)/length(.$value))*sum(.$value, na.rm = T)))
    })
  ) %>%
  unnest(bin_max, bin_sd) %>%
  mutate(bin_max = ifelse(is.na(bin_sd), NA, bin_max)) %>%
  ggplot() +
  geom_point(aes(x = location, y = Species, color = bin_max, size = bin_sd)) +
  scale_color_gradientn(name = r, colors = c("purple", "orange")) +
  facet_grid(phytogroups~., scales = "free_y", space = "free_y")




 testing <- left_join(metadata_filter, asv_table_filter) %>%
    filter(location == "Petermann Island") %>%
    ggplot() +
    geom_density(aes(x = Surface_Salinity_, weights = log10(rare_reads))) +
    geom_point(aes(x = Surface_Salinity_, y = log10(rare_reads)))+
    facet_wrap(~Species)

 
#### by location: bins, max in bins, loess to max, max of loess, sd of loess
species.bin.max <- left_join(metadata_filter, asv_table_filter)  %>%
  #filter(location %in% stations) %>%
  mutate(sp_rel = sum(reads)/7000,
         bins =  as.numeric(as.character(cut(!!sym(r),
                     breaks = seq(min(!!sym(r)),
                                  max(!!sym(r)),
                                  length.out = 15),
                     labels = seq(min(!!sym(r)),
                                  max(!!sym(r)),
                                  length.out = 15)[-15],
                     include.lowest = T)))) %>%
  group_by(phytogroups, Species, location, bins) %>%
  reframe(max_in_bin = max(sp_rel)) %>%
  group_by(phytogroups, Species, location, bins) %>%
  select(all_of(c("bins", "max_in_bin"))) %>%
  nest(data = c(bins, max_in_bin)) %>%
  mutate(
    model = map(data, ~ {
      if(length(unique(.$bins[.$max_in_bin >0])) < 5 | nrow(.) == 0){
        NULL
      }else{
        .x %>%
          loess(max_in_bin ~ bins, ., span = 0.75, degree = 2)
      }}),
    fit = map(model, ~{
      if(is.null(.)){
        vector(length = length(1:140))
      }else{
        predict(., 1:140)
      }
    }%>% as_tibble() %>% mutate(bins = 1:140, value = ifelse(value<0,0,value))),
    bin_max = map(fit, ~{.$bins[which.max(.$value)[1]]}),
    bin_sd = map(fit, ~{
      sqrt(
        sum(.$value*(.$bins - .$bins[which.max(.$value)[1]])^2, na.rm = T)/
          (((length(.$value>0)-1)/length(.$value))*sum(.$value, na.rm = T)))
    })
  ) %>%
  unnest(bin_max, bin_sd)


species.loc.bin.df <- left_join(metadata_filter, asv_table_filter) %>%
  group_by(across(all_of(
    c("location", "phytogroups", "Species", env_vars)))) %>%
  reframe(reads = sum(rare_reads, na.rm = T),
          prop_reads = reads/7000) %>%
  mutate(bins =  cut(!!sym(r),
                     breaks = seq(min(!!sym(r)),
                                  max(!!sym(r)),
                                  length.out = 15),
                     labels = seq(min(!!sym(r)),
                                  max(!!sym(r)),
                                  length.out = 15)[-15],
                     include.lowest = T)) %>%
  group_by(location, phytogroups, Species, bins) %>%
  reframe(max_prop = max(prop_reads)) %>%
  group_by(location) %>%
  mutate(n = n()) %>%
  filter(n > 50) %>% pull(location) %>% unique()

full_bins <- as.numeric(as.character(levels(species.loc.bin.df$bins)))

species.loc.bin.df %>%
  mutate(bins = as.numeric(as.character(bins))) %>%
  group_by(phytogroups, Species, location) %>%
  nest(data = c(bins, max_prop)) %>%
  mutate(
    model = map(data, ~ {
      if(length(unique(.$bins[.$max_prop >0])) < 5 | nrow(.) == 0){
        NULL
      }else{
        .x %>%
          loess(max_prop ~ bins, ., span = 0.75, degree = 2)
      }
      }),
    fit = map(model, ~{
      if(is.null(.)){
        vector(length = length(full_bins))
      }else{
        predict(., full_bins)
      }
    }%>%
      as_tibble() %>%
      mutate(bins = full_bins,
             value = ifelse(value<0,0,value))),
    bin_max = map(fit, ~{.$bins[which.max(.$value)[1]]}),
    bin_sd = map(fit, ~{
      sqrt(
        sum(.$value*(.$bins - .$bins[
          which.max(.$value)[1]])^2, na.rm = T)/
          (((length(.$value>0)-1)/length(.$value))*sum(.$value, na.rm = T)))
    })) %>% 
  unnest(bin_max, bin_sd) %>%
  mutate(bin_group = case_when(bin_max <= 32 ~ "Low",
                              bin_max <= 33.5 &
                                bin_max > 32 ~ "Mid",
                              TRUE ~ "High")) %>%
  ggplot() +
  geom_point(aes(x = location, y = Species, color = bin_group))


## For each top 15 species, plot the maximum proportion of reads within a window of 10 days
env_names <- c("Surface_Temperature_", "Surface_Salinity_", "days_since")

lapply(env_names, function(r){
plot_a <- left_join(metadata_filter, asv_table_filter) %>%
  group_by(across(all_of(c("sample", "Genus", env_vars)))) %>%
  reframe(reads = sum(rare_reads, na.rm = T),
          prop_reads = reads/7000) %>%
  #filter(Species %in% unique(top_15$Species)) %>%
  mutate(bins =  cut(!!sym(r),
                         breaks = seq(min(!!sym(r)),
                                      max(!!sym(r)),
                                      length.out = 15),
                         labels = seq(min(!!sym(r)),
                                      max(!!sym(r)),
                                    length.out = 15)[-15],
                         include.lowest = T)) %>%
  group_by(Genus, bins) %>%
  reframe(max_prop = max(prop_reads)) %>%
    mutate(bins = as.numeric(as.character(bins))) %>%
  ggplot(data = ., aes(x = bins, y = max_prop)) +
  geom_point() +
  geom_smooth(, method = "loess") +
  theme_bw() +
  my_theme +
  theme(strip.text = element_text(size = 8)) +
  facet_wrap(~Genus) +
  labs(y = "Max proportion of reads", x = capitalize(gsub("_", " ", r)))
  
  plot_name <- paste0("figures/genus_max-loess_over_", r, ".png")
  ggsave(plot = plot_a,
         filename = plot_name,
         width = 8, height = 8, dpi = 300)
})

left_join(metadata_filter, asv_table_filter) %>%
  group_by(across(all_of(c("sample", "Genus", env_vars)))) %>%
  reframe(reads = sum(rare_reads, na.rm = T),
          prop_reads = reads/7000) %>%
  group_by(across(all_of(c("Genus", env_vars)))) %>%
  reframe(max_prop = max(prop_reads)) %>%
  mutate(sal_bins = case_when(Surface_Salinity_ <= 32 ~ "Low",
                              Surface_Salinity_ <= 33.5 &
                                Surface_Salinity_ > 32 ~ "Mid",
                              TRUE ~ "High"),
         temp_bins = case_when(Surface_Temperature_ <= 0 ~ "Low",
                               Surface_Temperature_ <= 2 &
                                 Surface_Temperature_ > 0 ~ "Mid",
                               TRUE ~ "High")) %>%
  ggplot(data = .) +
  geom_point(aes(x = days_since, y = log10(max_prop), color = sal_bins)) +
  facet_wrap(~Genus)



  left_join(metadata_filter, asv_table_filter) %>%
    group_by(across(all_of(c("sample", "Genus", env_vars)))) %>%
    reframe(reads = sum(rare_reads, na.rm = T),
            prop_reads = reads/7000) %>%
    mutate(across(c(env_vars), ~as.numeric(as.character(cut(.x,
                                    breaks = seq(min(.x),
                                                 max(.x),
                                                 length.out = 15),
                                    labels = seq(min(.x),
                                                 max(.x),
                                                 length.out = 15)[-15],
                                    include.lowest = T))),
           .names = "{.col}bin")) %>% 
    group_by(across(all_of(c("Genus", paste0(env_vars, "bin"))))) %>%
    reframe(max_prop = max(prop_reads)) %>%
    mutate(sal_bins = case_when(Surface_Salinity_bin <= 32 ~ "Low",
                                Surface_Salinity_bin <= 33.5 &
                                  Surface_Salinity_bin > 32 ~ "Mid",
                                TRUE ~ "High"),
           temp_bins = case_when(Surface_Temperature_bin <= 0 ~ "Low",
                                 Surface_Temperature_bin <= 2 &
                                   Surface_Temperature_bin > 0 ~ "Mid",
                                 TRUE ~ "High")) %>%
  ggplot(data = .) +
  geom_point(aes(x = days_sincebin, y = log10(max_prop), color = sal_bins)) +
  facet_wrap(~Genus)

  
  
   left_join(metadata_filter, asv_table_filter) %>%
    group_by(across(all_of(c("sample", "Species", env_vars)))) %>%
    reframe(reads = sum(rare_reads, na.rm = T),
            prop_reads = reads/7000) %>%
    filter(Species %in% unique(top_15$Species)[4]) %>%
    mutate(across(c(env_vars), ~as.numeric(as.character(cut(.x,
                                    breaks = seq(min(.x),
                                                 max(.x),
                                                 length.out = 15),
                                    labels = seq(min(.x),
                                                 max(.x),
                                                 length.out = 15)[-15],
                                    include.lowest = T))),
           .names = "{.col}bin")) %>% 
    group_by(across(all_of(c("Species", paste0(env_vars, "bin"))))) %>%
    reframe(max_prop = max(prop_reads)) %>%
    plot_ly(., x= ~days_sincebin,
      y= ~Surface_Temperature_bin, z= ~Surface_Salinity_bin,
      type="scatter3d", mode="markers", color = ~max_prop,
      colors= colorRamp(c("purple","orange")))
  
  max_3d_plot<- left_join(metadata_filter, asv_table_filter) %>%
     group_by(across(all_of(c("sample", "Species", env_vars)))) %>%
     reframe(reads = sum(rare_reads, na.rm = T),
             prop_reads = reads/7000) %>%
     filter(Species %in% unique(top_15$Species)) %>%
     mutate(across(c(env_vars),
                 ~as.numeric(as.character(
                   cut(.x,
                        breaks = seq(min(.x),
                                     max(.x),
                                     length.out = 15),
                        labels = seq(min(.x),
                                     max(.x),
                                     length.out = 15)[-15],
                                     include.lowest = T))),
                   .names = "{.col}bin")) %>% 
     group_by(across(all_of(c("Species", paste0(env_vars, "bin"))))) %>%
     reframe(max_reads = max(reads)) %>%
     mutate(days_since_bin = case_when(days_sincebin <= 50 ~ "Early", 
                                      days_sincebin <= 100 &
                                        days_sincebin > 50 ~ "Mid",
                                      TRUE ~ "Late")) %>%
      ggplot(data = ., aes(x = Surface_Temperature_bin,
                         y = Surface_Salinity_bin,
                         color = days_since_bin,
                         size = max_reads)) +
    geom_point() +
    theme_bw() +
    my_theme +
    theme(strip.text = element_text(size = 8)) +
    facet_wrap(~Species) +
    labs(x = "Temperature", y = "Salinity", color = "Season phase",
         size = "Max reads") +
     scale_color_manual(values = c(
       "Early" = "#A0DDFF", "Mid" = "#EF8354", "Late" = "#624CAB"),
       labels = c("Early (< 50 days)", "Mid (51 - 100 days)",
                  "Late (> 100 days)"))
  

  ggsave(plot = max_3d_plot,
         filename = "figures/top_15-maximum_bins_3d.png",
         width = 10, height = 8, dpi = 300)

n_bins <- 15

metadata_filter %>%
  select(sample, 
         Surface_Temperature_, Surface_Salinity_, Potential_Density_) %>%
  group_by(sample) %>%
  mutate_at("Potential_Density_", ~as.numeric(mean(unlist(.x)[1:10], na.rm = T))) %>%
  distinct() %>%
  ggplot() +
  geom_contour(aes(x = Surface_Temperature_,
                    y = Surface_Salinity_,
                   fill = Potential_Density_))

  max_3d_plot<- left_join(metadata_filter, asv_table_filter) %>%
    group_by(across(all_of(c("sample", "phytogroups", env_vars)))) %>%
    reframe(reads = sum(rare_reads, na.rm = T),
            prop_reads = reads/7000) %>%
    mutate(across(c(env_vars),
                  ~as.numeric(as.character(
                    cut(.x,
                        breaks = seq(min(.x),
                                     max(.x),
                                     length.out = n_bins),
                        labels = seq(min(.x),
                                     max(.x),
                                     length.out = n_bins)[-n_bins],
                        include.lowest = T))),
                  .names = "{.col}bin")) %>%
    group_by(across(all_of(c("phytogroups", paste0(env_vars, "bin"))))) %>%
    mutate(max_reads = max(reads)) %>%
    mutate(days_since_bin = case_when(days_sincebin <= 50 ~ "Early", 
                                      days_sincebin <= 100 &
                                        days_sincebin > 50 ~ "Mid",
                                      TRUE ~ "Late")) %>%
    ggplot(data = .) +
    geom_isopycnal(aes(x = Surface_Temperature_,
                       y = Surface_Salinity_)) +
    geom_point(aes(x = Surface_Temperature_bin,
                   y = Surface_Salinity_bin,
                   color = days_since_bin,
                   size = max_reads)) +
    theme_bw() +
    my_theme +
    theme(strip.text = element_text(size = 8)) +
    facet_wrap(~phytogroups) +
    labs(x = "Temperature", y = "Salinity", color = "Season phase",
         size = "Max reads") +
    scale_color_manual(values = c(
      "Early" = "#A0DDFF", "Mid" = "#EF8354", "Late" = "#624CAB"),
      labels = c("Early (< 50 days)", "Mid (51 - 100 days)",
                 "Late (> 100 days)"))
  
  
  ggsave(plot = max_3d_plot,
         filename = "figures/phytogroups-maximum_bins_3d.png",
         width = 10, height = 8, dpi = 300)
  
  
  
  
  
  
#### NMDS samples
bio <- asv_table_filter %>%
  select(Species, rare_reads, sample) %>%
  pivot_wider(.,
              id_cols = "sample", names_from = "Species", values_from = "rare_reads", values_fn = function(r)sum(r,na.rm = T), values_fill = 0) %>%
  column_to_rownames(var = "sample")

correlation_matrix <- cor(t(bio), method = "spearman")

bio_NMS1 <-  metaMDS(bio,
                     k = 2,
                     maxit = 200, 
                     trymax = 100,
                     try = 50,
                     wascores= TRUE,
                     weakties = T,
                     expand = TRUE,
                     previous.best = T,
                     autotransform = FALSE)

bio_spp.fit <- envfit(bio_NMS1, bio, permutations = 999) # this fits species vectors
spp.scrs <- as.data.frame(scores(bio_spp.fit, display = "vectors", arrow.mul=2.5))
spp.scrs <- cbind(spp.scrs, vars = rownames(spp.scrs)) #add species names to dataframe
spp.scrs <- cbind(spp.scrs, pval = bio_spp.fit$vectors$pvals) #add pvalues to dataframe so you can select species which are significant
#spp.scrs<- cbind(spp.scrs, abrev = abbreviate(spp.scrs$Species, minlength = 6)) #abbreviate species names
sig.spp.scrs <- subset(spp.scrs, pval<=0.05) #subset data to show species

env_df <- metadata_filter %>%
  select(sample, all_of(env_vars))

bio.env.fit <- envfit(bio_NMS1, env_df, permutations = 999, na.rm = T) 

env.scrs <- as.data.frame(scores(bio.env.fit, display = "vectors", arrow.mul=2.5)) 
env.scrs <- cbind(env.scrs, vars = rownames(env.scrs)) #add species names to dataframe
env.scrs <- cbind(env.scrs, pval = bio.env.fit$vectors$pvals) #add pvalues to dataframe so you can select species which are significant
#spp.scrs<- cbind(spp.scrs, abrev = abbreviate(spp.scrs$Species, minlength = 6)) #abbreviate species names
sig.env.scrs <- subset(env.scrs, pval<=0.05) #subset data to show species 


#Get the vectors for env.fit
df_envfit<-scores(bio.env.fit,display=c("vectors"))
df_envfit<-df_envfit*vegan:::ordiArrowMul(df_envfit)
df_envfit<-as.data.frame(df_envfit)

env_mds <- as.data.frame(scores(bio_NMS1$points)) %>%
  rownames_to_column("sample") %>%
  left_join(env_df)

bio_mds <- sig.spp.scrs %>%
  rownames_to_column("Species") %>%
  left_join(asv_table_filter) %>%
  distinct(Species, NMDS1, NMDS2, phytogroups)

group_colors <- c("#F06400", "#00F064","#008CF0","gold2","#6400F0","#F0008C","grey10")
names(group_colors) <- c("Cryptophytes", "Diatoms", "Dinoflagellates", "Haptophytes",
                         "MAST", "Rhodophytes", "Greenalgae")

require(ggrepel)
ggplot(bio_mds) +
  geom_point(data = env_mds,
             aes(x = MDS1, y = MDS2, fill = Surface_Salinity_), pch = 21, size = 3) +
  geom_segment(data = bio_mds,
               aes(x = 0, xend=NMDS1, y=0, yend=NMDS2, color = phytogroups),
               arrow = arrow(length = unit(0.25, "cm")), lwd = 0.3) +
  geom_text_repel(data = bio_mds,
                  aes(x=NMDS1, y=NMDS2, label = Species, color = phytogroups),
                  cex = 3, direction = "both") +
  scale_fill_gradientn(name = "Days since", colors = c("pink", "purple")) +
  scale_color_manual(name = "", values = group_colors) +
  theme_bw() +
  my_theme

lapply(env_vars, function(x){
  ggplot(bio_mds) +
    geom_point(data = env_mds,
               aes(x = MDS1, y = MDS2, fill = !!sym(x)), pch = 21, size = 3) +
    geom_segment(data = bio_mds,
                 aes(x = 0, xend=NMDS1, y=0, yend=NMDS2, color = phytogroups),
                 arrow = arrow(length = unit(0.25, "cm")), lwd = 0.3) +
    geom_text_repel(data = bio_mds,
                    aes(x=NMDS1, y=NMDS2,
                        label = Species, color = phytogroups),
                    cex = 3, direction = "both", show.legend = F) +
    scale_fill_gradientn(name = capitalize(gsub("_", " ", x)),
                         colors = c("pink", "purple")) +
    scale_color_manual(name = "", values = group_colors) +
    theme_bw() +
    my_theme
  
  plot_name <- paste0("figures/NMDS_", x, ".png")
  ggsave(plot = last_plot(),
         filename = plot_name,
         width = 8, height = 9, dpi = 300)
})


corr_plot <- rcorr(as.matrix(bio), type = "pearson")
diag(corr_plot$P) <- 0
corrplot(corr_plot$r, type = "upper", tl.col = "black", tl.srt = 45, addCoef.col = "grey30", number.cex = 0.5, tl.cex = 0.5,
         p.mat = corr_plot$P, sig.level = 0.05, insig = "blank")

bio_phyto <- bio %>%
  rownames_to_column(var = "sample") %>%
  left_join(env_df) %>%
  column_to_rownames(var = "sample") %>%
  mutate_all(~as.numeric(.x))

env_names <- colnames(env.fit_df)[-1]

corr_plot <- rcorr(as.matrix(bio_phyto), type = "pearson")
diag(corr_plot$P) <- 0
corrplot(corr_plot$r, type = "upper", tl.col = "black", tl.srt = 45,
         number.cex = 0.5, tl.cex = 0.5,
         p.mat = corr_plot$P, sig.level = 0.05, insig = "blank")

corr_species <- as.data.frame(corr_plot$r) %>%
  rownames_to_column("var2") %>%
  pivot_longer(cols = -var2, names_to = "var1", values_to = "r") %>%
  mutate(P = as.vector(corr_plot$P)) %>%
  filter(var1 != var2) %>%
  filter(var1 %in% env_vars) %>%
  filter(var2 %in% unique(colnames(bio))) %>%
  ggplot() +
  geom_point(aes(x = var1, y = var2, color = r,
                 size = if_else(P < 0.05, P, NA))) +
  scale_color_gradient2(low = "blue", high = "red", limits = c(-0.75,0.75)) +
  scale_x_discrete(labels = ~capitalize(str_replace_all(., "_", " "))) +
  scale_size_continuous(guide = "none", transform = "reverse",
                        range = c(1,5)) +
  theme_bw() +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "", color = "Correlation")

ggsave(plot = corr_species,
       filename = "figures/correlation_species_env.png",
       width = 8, height = 10, dpi = 300)

ggplot() +
  geom_point(data = env_mds,
             aes(x = MDS1, y = MDS2, fill = Surface_Salinity_),
             pch = 21, size = 3) +
  geom_segment(data = bio_mds,
               aes(x = 0, xend=NMDS1, y=0, yend=NMDS2, color = phytogroups),
               arrow = arrow(length = unit(0.25, "cm")), lwd = 0.3) +
  scale_fill_viridis_c() +
  scale_color_manual(name = "", values = group_colors)


ggplot() +
  geom_point(data = env_mds,
             aes(x = days_since, y = Surface_Temperature_, fill = Surface_Salinity_), pch = 21, size = 3) +
  # geom_segment(data = bio_mds,
  #              aes(x = 0, xend=NMDS1, y=0, yend=NMDS2, color = phytogroups),
  #              arrow = arrow(length = unit(0.25, "cm")), lwd = 0.3) +
  scale_fill_viridis_c() +
  scale_color_manual(name = "", values = group_colors)



### How is the environment changing over time and space?

metadata_filter %>%
  unnest() %>%
ggplot() +
  geom_tile(aes(x = days_since, y = -Pressure_,
                fill = Temperature_),width = 3) +
  facet_wrap(~as.numeric(latitude), ncol = 1) +
  theme_bw() +
  my_theme +
  theme(text = element_text(size = 5))

left_join(diversity_filter,  metadata_filter) %>%
ggplot() +
  geom_point(aes(x = latitude, y = longitude, color = richness_phytogroups)) +
  facet_wrap(~seasonyear)

### Find clusters
bio <- asv_table_filter %>%
  select(Species, rare_reads, sample) %>%
  pivot_wider(.,
              id_cols = "sample", names_from = "Species", values_from = "rare_reads", values_fn = function(r)sum(r,na.rm = T), values_fill = 0) %>%
  column_to_rownames(var = "sample")

correlation_matrix <- cor(bio, method = "spearman")
hierarchical_result <- hclust(dist(1 - correlation_matrix), method = "ward.D2")
clusters <- cutree(hierarchical_result, k = 3)

bio_NMS1 <-  metaMDS(dist(1 - correlation_matrix),
                     k = 2,
                     maxit = 200, 
                     trymax = 100,
                     try = 50,
                     wascores= TRUE,
                     weakties = T,
                     expand = TRUE,
                     previous.best = T,
                     autotransform = FALSE)

data_scores_cluster <- as.data.frame(scores(bio_NMS1)) %>%
  rownames_to_column("Species") %>%
  mutate(clusters = factor(
    Species, levels = names(clusters), labels = clusters)) %>%
  left_join(asv_table_filter %>% select(Species, phytogroups, sample, reads) %>% distinct(), by = "Species")

env.fit_df <- metadata_filter %>%
  select_if(Negate(is.list)) %>%
  select(sample, Surface_Salinity_, Surface_Temperature_, secchi_depth, latitude, longitude, days_since, month, Meltwater_Fraction_, land_dist_km)


### Correlations
bio_phyto <- asv_table_filter %>%
  select(Species, phytogroups, rare_reads, sample) %>%
  mutate(cluster = factor(Species,
                          levels = names(clusters), labels = clusters),
         ID = paste(cluster, phytogroups, sep = "_")) %>%
  pivot_wider(.,
              id_cols = "sample", names_from = "ID", values_from = "rare_reads", values_fn = function(r)sum(r,na.rm = T), values_fill = 0)

bio_names <- colnames(bio_phyto)[-1]
bio_phyto <- bio_phyto %>%
  left_join(env.fit_df) %>%
  column_to_rownames(var = "sample") %>%
  mutate_all(~as.numeric(.x))

env_names <- colnames(env.fit_df)[-1]

corr_plot <- rcorr(as.matrix(bio_phyto), type = "pearson")
diag(corr_plot$P) <- 0
corrplot(corr_plot$r, type = "upper", tl.col = "black", tl.srt = 45, addCoef.col = "grey30", number.cex = 0.5, tl.cex = 0.5,
         p.mat = corr_plot$P, sig.level = 0.05, insig = "blank")


as.data.frame(corr_plot$r) %>%
  rownames_to_column("var2") %>%
  pivot_longer(cols = -var2, names_to = "var1", values_to = "r") %>%
  mutate(P = as.vector(corr_plot$P)) %>%
  filter(var1 != var2) %>%
  filter(var1 %in% env_names) %>%
  filter(var2 %in% bio_names) %>%
  filter(P < 0.05) %>%
  separate(var2, into = c("cluster", "phytogroup"), sep = "_") %>%
  ggplot() +
  geom_point(aes(x = var1, y = phytogroup, color = r, size = P)) +
  scale_color_gradient2(low = "blue", high = "red", limits = c(-0.75,0.75)) +
  scale_x_discrete(labels = ~capitalize(str_replace_all(., "_", " "))) +
  scale_size_continuous(guide = "none", transform = "reverse",
                        range = c(1,5)) +
  theme_bw() +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "", color = "Correlation") +
  facet_wrap(~cluster)


bio_phyto %>%
  rownames_to_column("sample") %>%
  pivot_longer(cols = -all_of(c("sample", env_names)), names_to = "ID", values_to = "rare_reads") %>%
  separate(ID, into = c("cluster", "phytogroup"), sep = "_") %>%
  ggplot() +
  geom_point(aes(x = days_since, y = rare_reads, color = phytogroup)) +
  facet_grid(phytogroup~cluster)
  
bio_phyto %>%
  rownames_to_column("sample") %>%
  pivot_longer(cols = -all_of(c("sample", env_names)), names_to = "ID", values_to = "rare_reads") %>%
  separate(ID, into = c("cluster", "phytogroup"), sep = "_") %>%
  filter(cluster == 1) %>%
  lmer(rare_reads~ Surface_Temperature_ + Surface_Salinity_ + Meltwater_Fraction_ + (1|month) + (1|region))



#source("analysis/ggnested_pattern.R")
devtools::install_github("gmteunisse/ggnested")
install.packages("devtools")
install.packages("ggpattern")
require(ggpattern)

asv_table_filter %>%
  select(Species, phytogroups, rare_reads, sample) %>%
  mutate(cluster = factor(Species,
                          levels = names(clusters), labels = clusters)) %>%
  group_by(cluster) %>%
  mutate(clust_sum = sum(rare_reads, na.rm = T)) %>%
  group_by(Species, cluster, phytogroups) %>%
  reframe(sp_rel = sum(rare_reads, na.rm = T)/clust_sum) %>%
  distinct(Species, cluster, phytogroups, .keep_all = T) %>%
  group_by(cluster) %>%
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
    pat = case_when(phytogroups == "Diatoms" & cluster == "2"~"wave",
                    phytogroups == "Diatoms" & cluster == "3"~"stripe",
                    TRUE~"none"),
    pat = if_else(Species == "Porosira_sp.", "circle", pat)) %>%
  ggnested_pattern(data = ., aes_string(main_group = "phytogroups",
                                        sub_group = "Species", 
                                        x = "cluster", y = "sp_rel", pattern = "pat"),
                   main_palette = group_colors) + 
  geom_col_pattern(pattern_fill="white", pattern_color="white", pattern_key_scale_factor=0.25,
                   pattern_size = 0.1) +
  scale_y_continuous(name = "Relative abundance",
                     expand = c(0,0)) +
  scale_x_discrete(name = "Cluster", expand = c(1E-3,1E-3), breaks = c(1,2,3)) +
  theme_nested(theme_linedraw) + 
  theme(
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(size = 12),
    legend.title = element_blank(),
    panel.border = element_rect(color = "black")) +
  guides(fill=guide_legend(ncol = 1))



#### Binning tS-days
# mutate_at(env_vars, ~as.numeric(
#   as.character(cut(as.numeric(.),
#                    breaks = seq(min(as.numeric(.), na.rm = T),
#                                 max(as.numeric(.), na.rm = T),
#                                 length.out = n_bins),
#                    labels = seq(min(as.numeric(.), na.rm = T),
#                                 max(as.numeric(.), na.rm = T),
#                                 length.out = n_bins)[-n_bins],
#                    include.lowest = T)))) %>%
#   group_by(across(all_of(c(
#     "phytogroups", "Genus", "Species", "seasonyear", "location", env_vars)))) %>%
#   reframe(max_prop = max(prop_reads)) %>%
#   group_by(phytogroups, Genus, Species, seasonyear, location) %>%
#   filter(max_prop == max(max_prop))


#### Set up phyloseq
###Phyloseq this
asv_wide <- asv_table_filter %>%
  dplyr::select(Feature.ID, sample, rare_reads) %>%
  pivot_wider(values_from = rare_reads,
              id_cols = Feature.ID, names_from = sample)
asv_wide_df <- as.data.frame(asv_wide[,-1]) #gets rid of Feature.ID column
asv_wide_df[is.na(asv_wide_df)] <- 0 #if there are NA reads, make 0

#Make otu_table
rownames(asv_wide_df) <- asv_wide$Feature.ID #make FeatureID the rownames
in_biom <- otu_table(asv_wide_df, taxa_are_rows = T)
taxa_names(in_biom) <- asv_wide$Feature.ID

#Make sample_data
in_biom_metadata <- sample_data(metadata_filter)
sample_names(in_biom_metadata) <- metadata_filter$sample

#Make tax_table
taxatable <- asv_table_filter %>% 
  select(Feature.ID, Kingdom, Supergroup, Division, Class, Order, Family, Genus, Species, phytogroups) %>%
  distinct() %>%
  as.data.frame()
in_biom_tax <- tax_table(taxatable %>% select(-Feature.ID) %>% as.matrix()) 
taxa_names(in_biom_tax) <- taxatable$Feature.ID

#Make phyloseq object
physeq <- phyloseq(in_biom,
                   in_biom_metadata,
                   in_biom_tax)

# Get the top taxa
top_level <- "phytogroups"
nested_level <- "Species"
sample_order <- NULL
top_asv <- top_taxa(physeq, n_taxa = 10)

# Create names for NA taxa
ps_tmp <- top_asv$ps_obj %>%
  name_na_taxa()

# Add labels to taxa with the same names
ps_tmp <- ps_tmp %>%
  label_duplicate_taxa(tax_level = nested_level)

# Generate a palette basedon the phyloseq object
pal <- taxon_colours(ps_tmp,
                     tax_level = top_level)

# Convert physeq to df
psdf <- psmelt(ps_tmp)
nested_top_taxa(physeq,
                top_tax_level = "Class",
                nested_tax_level = "Species",
                n_top_taxa = 3, 
                n_nested_taxa = 3)

in_biom_tax
















