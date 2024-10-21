#### Metadata Chlorophyll

## Load chl_a data 
find_files <- list.files("data/Chla_Rick/Chla_TimeSeries/", pattern = "*.csv")
chl_load <- do.call("rbind", lapply(find_files,function(r){
  read_csv(paste0("data/Chla_Rick/Chla_TimeSeries/", r))[-1,]
}))
chl_load$time <- as_date(chl_load$time)
chl_load[,2:ncol(chl_load)] <- lapply(chl_load[,2:ncol(chl_load)], as.numeric)

## Pull unique ID codes, Dates, lat/lon, and secchi depths
subset_meta <- metadata %>%
  mutate(Date = as_date(date_local)) %>%
  group_by(UNIQUE_ID_CODE, Date, latitude, longitude) %>%
  group_keys()

## Add the metadata to the chla file by date
meta_w_chl <- left_join(subset_meta, chl_load, by = join_by(Date == time))

## The resolution of the latitude and longitude in the metadata is higher than in the chla
## file, so we need to find the closest latitude and longitude in the metadata to the chla file
## Steps: Calculate the aboslute difference between the lats and lons in each data set, group by
## the unique ID code, time, and secchi depth, arrange by the differences, and filter by the first
## row of each group
meta_match_chla <- meta_w_chl %>%
  mutate(diff_lat = abs(latitude.x - latitude.y),
         diff_long = abs(longitude.x - longitude.y)) %>%
  group_by(UNIQUE_ID_CODE, Date, latitude.x, longitude.x) %>%
  arrange(diff_lat, diff_long) %>%
  filter(row_number()==1) %>%
  select(UNIQUE_ID_CODE, Date, latitude.x, longitude.x, chlor_a) %>%
  rename(latitude = latitude.x, longitude = longitude.x)

## Add the chla data to the metadata
metadata_w_chl <- left_join(metadata, meta_match_chla,
                            by = c("UNIQUE_ID_CODE", "date_local" = "Date", "latitude", "longitude"))

# Use a linear model to predict chlorophull from secchi depth
chl_fit <- lm(log10(chlor_a) ~ log10(secchi_depth), data = metadata_w_chl %>% filter(seasonyear == "2018-2019"))

secchi_2017 <- metadata_w_chl %>%
  filter(year == 2017) %>%
  pull(secchi_depth)

metadata_w_chl$chlor_a[metadata_w_chl$year == 2017] <- 10^predict(chl_fit, data.frame(secchi_depth = secchi_2017))

library("ggpmisc")

lm_eqn <- function(m){
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

label_df <- data.frame(seasonyear = "2018-2019", x = 15, y = 0, label = lm_eqn(chl_fit))

chla_regress <- metadata_w_chl  %>%
  filter(seasonyear%in%c("2017-2018", "2018-2019","2019-2020","2021-2022")) %>%
  ggplot(aes(x = log10(secchi_depth),
             y = log10(chlor_a))) +
  geom_point(aes(color = ifelse(year == "2017", "red", "black")),
             show.legend = FALSE) +
  scale_color_manual(name = "",
                     values = c("red" = "red", "black" = "black")) +
  geom_text(data = label_df, x = 15, y = 0, label = lm_eqn(chl_fit), parse = TRUE) +
  geom_smooth(data = metadata_w_chl %>% filter(seasonyear == "2018-2019"),
              aes(x = log10(secchi_depth),
                  y = log10(chlor_a)),
              method = "lm") +
  facet_wrap(~seasonyear, ncol = 2) +
  scale_x_continuous(name = "Obs. Secchi Depth (m)") +
  scale_y_continuous(name = "Chl-a (log10 mg/m^3)") +
  theme_bw() +
  my_theme

ggsave(
  filename =  "data/Chla_Rick/chla_regressions.jpg",
  chla_regress,
  width = 7,
  height = 5,
  units = "in",
  dpi = 300
)


### Chla-part 2 : 293.9/(Zsd^2.345)
metadata_w_chl$chlor_a_2 <- 293.9/(metadata_w_chl$secchi_depth^2.345)


## Add the CTD data to the metadata
metadata_full <- left_join(metadata_w_chl,
                           ctd_w_id %>% select(-c(latitude,longitude)),
                           by = "UNIQUE_ID_CODE")


metadata_full %>%
  filter(seasonyear == "2021-2022") %>%
  group_by(UNIQUE_ID_CODE) %>%
  mutate(chlor_rbr = mean(unlist(Chlorophyll_a_)[1:10], na.rm = T)) %>%
  pull(chlor_rbr)


chla_options <- metadata_full  %>%
  filter(seasonyear%in%c("2017-2018", "2018-2019","2019-2020","2021-2022")) %>%
  group_by(UNIQUE_ID_CODE) %>%
  mutate(chlor_rbr = mean(unlist(Chlorophyll_a_)[1:5], na.rm = T)) %>%
  ggplot() +
  geom_point(aes(x = log10(secchi_depth),
                 y = log10(chlor_a),
                 color = ifelse(year == "2017", "Fitted", "Sat."))) +
  geom_point(aes(x = log10(secchi_depth),
                 y = log10(chlor_a_2),
                 color = "Z Eqn.")) +
  geom_point(aes(x = log10(secchi_depth),
                 y = log10(chlor_rbr),
                 color = "RBR (1-5m)")) +
  scale_color_manual(name = "",
                     values = c(
                       "Fitted" = "blue",
                       "Sat." = "black",
                       "Z Eqn." = "red",
                       "RBR (1-5m)" = "green")) +
  facet_wrap(~seasonyear) +
  scale_x_continuous(name = "Obs. Secchi Depth (log10 m)") +
  scale_y_continuous(name = "Chl-a (log10 mg/m^3)") +
  theme_bw() +
  my_theme 

ggsave(
  filename =  "data/Chla_Rick/chla_all_options.jpg",
  chla_options,
  width = 7,
  height = 5,
  units = "in",
  dpi = 300
)