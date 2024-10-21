
library("readxl")
library("plyr")
library("dplyr")
library("readr")
library("stringr")
library("powerjoin")
library("lubridate")
library("ggplot2")
library("tidyr")

my_theme = theme_linedraw() + theme(text = element_text(size = 14), strip.background = element_blank(), strip.text = element_text(face = "bold", color = "black"))

# intro <- paste("On this date: ", mini_master$date_format, " the coordinates or the times \n
#                    of the CTD/RBR files and the master sheet don't agree")
# ctdInfo <- paste("The CTD/RBR File_ID: ", sub_ctd$File_ID)
# metaInfo <- paste("The Mastersheet  ", sub_ctd$File_ID)
# 
# user <- dlg_input("Who are you?", Sys.info()["user"])$res

files_home <- c("~/Documents/Fjord_Phyto/Chapter_2/Data/")
#### Lat lon switched in rbr_2022 - potentially
ctd_load <- read_csv(paste0(files_home, "All_final_df.csv"))

rbr_2022 <- read_xlsx(paste0(files_home, "Nik-All-CTD-Merged.xlsx"))
missing_lats <- read_csv("~/Desktop/missing_lat.csv") %>%
  mutate(File_ID = gsub("\\\\", "-", File_ID))
master_sheet <- read_xlsx(paste0(files_home, "mastersheet_06042024.xlsx"), sheet = 1)
master_sheet <- master_sheet[-which(apply(master_sheet, 2, function(r)sum(is.na(r))) == nrow(master_sheet))]
master_sheet$day <- as.numeric(master_sheet$day)
master_sheet$month <- as.numeric(master_sheet$month)
master_sheet$year <- as.numeric(master_sheet$year)
master_sheet$posix <- as.POSIXlt(paste(master_sheet$year, master_sheet$month, master_sheet$day, sep = "-"), format = "%Y-%m-%d")
master_sheet$time_convert <- master_sheet$posix + (as.numeric(master_sheet$time_local)*24*3600)
master_sheet$time_local <- strftime(master_sheet$time_convert, "%H:%M:%S")
master_sheet$latitude <- as.numeric(master_sheet$latitude)
master_sheet$longitude <- as.numeric(master_sheet$longitude)

ch1_meta <- read_csv("~/Downloads/metadatafix_fixed.csv")
ch1_meta$time_local <- as.character(ch1_meta$time_local)
ch1_meta$Site_Name[ch1_meta$Site_Name == "Enterprise Island" & !is.na(ch1_meta$Site_Name)] <- 
  "Wilhelmina Bay"
ch1_meta$Operator[is.na(ch1_meta$Operator)] <- "-999"
### Match sample_ids to mastersheet
match_master_samples <- left_join(ch1_meta, master_sheet, by = c("year", "month", "day", "samples" = "Genetics_18sv9_Sample ID", "Operator" = "ship operator"))

#write_csv(match_master_samples, file = "~/Downloads/master_samples_match.csv")

colnames(ctd_load) <- gsub(" ", "_", str_trim(gsub("\\(.*","",colnames(ctd_load))))
colnames(ctd_load)[9:12] <- c("Time", "Latitude", "Longitude", "File_ID")
ctd_load <- ctd_load[,-c(13:15)]
ctd_load$Time <- as.character(ctd_load$Time)
colnames(rbr_2022) <- gsub(" ", "_", colnames(rbr_2022))
colnames(rbr_2022)[14] <- "File_ID"
rbr_2022$Density <- rbr_2022$Density_anomaly + 1000 ## Convert anomaly to measured

rbr_ctd_all <- rbind.fill(ctd_load, rbr_2022) %>%
  filter(Depth > 0) %>%
  #mutate(Date = format(Time, "%Y-%m-%d")) %>%
  mutate(File_ID = gsub("\\\\", "-", File_ID)) %>%
  left_join(., missing_lats, by = "File_ID") %>%
  mutate(Longitude = if_else(is.na(Longitude.x), Longitude.y, Longitude.x),
         Latitude = if_else(is.na(Latitude.x), Latitude.y, Latitude.x)) %>%
  select(-contains(c(".x", ".y"))) %>%
  filter(File_ID != "208040_20211216_1133_examplecrapdates_gps")

### Match master with sample ids to CTD/RBR files
### Turn lat lon into signif
# rbr_ctd_all$Longitude <- signif(rbr_ctd_all$Longitude, 3)
# rbr_ctd_all$Latitude <- signif(rbr_ctd_all$Latitude, 3)
# match_master_samples$Long <- signif(match_master_samples$Long, 3)
# match_master_samples$Lat <- signif(match_master_samples$Lat, 3)

list_cols <- colnames(rbr_ctd_all)[!colnames(rbr_ctd_all) %in% c("Time", "File_ID", "cruise_season", "Longitude", "Latitude")]
unique_ctd_id <- rbr_ctd_all %>%
  group_by(Time, File_ID, cruise_season, Longitude, Latitude) %>%
  reframe(across(list_cols,~list(.x))) %>%
  mutate(day = day(Time),
         month = month(Time),
         year = year(Time),
         hour = hour(Time))

## Match based off day, month, year, UTC hour, lat, long
match_master_samples$hour <- hour(match_master_samples$time_utc)


for (i in seq_along(unique(match_master_samples$sample_id))){
  mini_master <- match_master_samples %>%
    filter(sample_id == unique(match_master_samples$sample_id)[i])
  
  if (all(is.na(mini_master$Lat))){
    mini_master$date_format <- as_datetime(paste(mini_master$Date, mini_master$time_utc),
                                                                   format = "%m/%d/%Y %H:%M:%S")
    na_ctd <- sub_ctd
    na_ctd[,!colnames(na_ctd) %in% c("day","month","year")] <- NA
    na_ctd$date_format <- NA
    combinedOutput <- left_join(mini_master, na_ctd, by = c("day", "month", "year"))
  }else{
  ## Day, month, year should match 
  sub_ctd <- unique_ctd_id %>%
    filter(day == mini_master$day &
             month == mini_master$month &
             year == mini_master$year)
  
  if (nrow(sub_ctd) == 0){
    mini_master$date_format <- as_datetime(paste(mini_master$Date, mini_master$time_utc),
                                           format = "%m/%d/%Y %H:%M:%S")
    na_ctd <- sub_ctd
    na_ctd[,!colnames(na_ctd) %in% c("day","month","year")] <- NA
    na_ctd$date_format <- NA
    combinedOutput <- left_join(mini_master, na_ctd, by = c("day", "month", "year"))
  }else{
  # Are the latitude and longitude within 0.15 degrees?
  lat_test <- abs(sub_ctd$Latitude - mini_master$Lat) < 0.15
  lon_test <- abs(sub_ctd$Longitude - mini_master$Long) < 0.15
  
  if (!all(lat_test, lon_test) & any(lat_test, lon_test)){
    lat_lon_id <- which(lat_test == TRUE & lon_test == TRUE)
    sub_ctd <- sub_ctd[lat_lon_id,]
    lat_test <- lat_test[lat_lon_id]
    lon_test <- lon_test[lat_lon_id]
  }
  # Are the times within one hour?
  mini_master$date_format <- as_datetime(paste(mini_master$Date, mini_master$time_utc),
                                     format = "%m/%d/%Y %H:%M:%S")
  sub_ctd$date_format <- as_datetime(sub_ctd$Time,
                                      format = "%Y-%m-%d %H:%M:%S")
  
  timeDiff <- abs(difftime(mini_master$date_format,sub_ctd$date_format, units = "hours")) < 1 
  
  if (!all(timeDiff) & any(timeDiff)){
    sub_ctd <- sub_ctd[which(timeDiff == TRUE),]
    timeDiff <- timeDiff[which(timeDiff == TRUE)]
  }
  
  if (all(lat_test, lon_test, timeDiff) & !is.na(all(lat_test, lon_test, timeDiff))){
  combinedOutput <- left_join(mini_master, sub_ctd, by = c("day", "month", "year"))
  }else{
  na_ctd <- sub_ctd
  na_ctd[,!colnames(na_ctd) %in% c("day","month","year")] <- NA
  combinedOutput <- left_join(mini_master, na_ctd, by = c("day", "month", "year"))
  }
  }
  }
  if (i == 1){
    combined_df <- combinedOutput
  } else{
    combined_df <- rbind(combined_df, combinedOutput)
  }
}

#combined_df %>% filter(sample_id == "ManifestSample_004") %>% View()
#df_new['Meltwater Fraction (%)'] = (-0.021406* df_new.mean_sal + 0.740392) * 100
# Meltwater calc
surface_calcs <- combined_df %>% 
  tidyr::drop_na(File_ID) %>%
  group_by(File_ID) %>%
  mutate(surface_salinity = mean(unlist(Salinity)[unlist(Depth) <= 1], na.rm = T),
          surface_temp = mean(unlist(Temperature)[unlist(Depth) <= 1], na.rm = T),
          surface_density = mean(unlist(Density)[unlist(Depth) <= 1], na.rm = T),
          meltwater = (-0.021406 * surface_salinity + 0.740392) * 100,
         days_since = yday(date_format.x) - yday("2000-11-01"),
         days_since = ifelse(days_since < 0, days_since + 365, days_since),
         period = factor(case_when(days_since < 50 ~ "Early",
                                   days_since >= 51 & days_since < 100 ~ "Mid",
                                   days_since >= 100 ~ "Late"), levels = c("Early", "Mid", "Late")))

### Load in the chl_a files and match to the time and 
find_files <- list.files("~/Documents/Fjord_Phyto/Chapter_2/Data/chl_a/", pattern = "*.csv")
chl_load <- do.call("rbind", lapply(find_files,function(r){
  read_csv(paste0("~/Documents/Fjord_Phyto/Chapter_2/Data/chl_a/", r))[-1,]
}))
chl_load$time <- as.Date(chl_load$time, format = "%Y-%m-%dT%H:%M:%SZ")
find_chla_match <- function(row){
  chl_load %>%
    filter(time == row$Date &
             latitude == latitude[which.min(abs(as.numeric(latitude) - row$Lat))] &
             longitude == longitude[which.min(abs(as.numeric(longitude) - row$Long))]) %>%
    pull(chlor_a)
}

subset_to_chla <- surface_calcs %>%
  mutate(Date = as.Date(Date, format = "%m/%d/%Y")) %>%
  group_by(File_ID, Date, Lat, Long) %>%
  group_keys() %>%
  left_join(., chl_load, by = join_by(Date == time)) %>%
  mutate(diff_lat = abs(Lat - as.numeric(latitude)),
         diff_long = abs(Long - as.numeric(longitude))) %>%
  group_by(File_ID, Date, Lat, Long) %>%
  arrange(diff_lat, diff_long) %>%
  slice(1) %>%
  select(File_ID, Date, Lat, Long, chlor_a)
  
surface_with_chla <- surface_calcs %>%
  mutate(Date = as.Date(Date, format = "%m/%d/%Y")) %>%
  left_join(., subset_to_chla, by = c("File_ID", "Date", "Lat", "Long"))

surface_with_chla %>%
  ggplot() +
  geom_point(aes(x = Lat, y = Long, color = as.numeric(chlor_a))) +
  scale_color_viridis_c() +
  facet_grid(~year)




 
surface_with_chla %>% 
  group_by(File_ID) %>%
  mutate(lat = mean(Latitude, na.rm = T),
         lon = mean(Longitude, na.rm = T),
          site = paste0(lat, " ", lon)) %>%
  unnest(cols = c(Depth, Temperature, Salinity, Density)) %>%
  group_by(File_ID) %>%
  pivot_longer(cols = c(Temperature, Salinity, Density),names_to = "var", values_to = "vals") %>%
  ggplot() + 
  geom_point(aes(x = vals, y = -Depth, color = period)) +
  facet_grid(~var, scales = "free_x")

{prelim_chla_ts <- surface_with_chla %>% 
  group_by(File_ID) %>%
  mutate(lat = mean(Latitude, na.rm = T),
         lon = mean(Longitude, na.rm = T),
         site = paste0(lat, " ", lon)) %>%
  unnest(cols = c(Depth, Temperature, Salinity, Density)) %>%
  filter(Salinity>30) %>%
  ggplot() +
  geom_point(aes(x = Salinity, y = Temperature, color = as.numeric(chlor_a))) +
  scale_color_viridis_c(name = "Surface\nChl-a") +
  facet_grid(~period) +
  my_theme

ggsave(
  filename = "~/Documents/Fjord_Phyto/Chapter_2/Figures/chlor_ts_prelim.png",
  prelim_chla_ts,
  height = 3,
  width = 7,
  units = "in"
)}
#A) Shallow mld, low salinity -> low density at surface, increasde with depth
#B) Normal well mixed, 
rbr_ctd_all %>% 
  group_by(Latitude, Longitude, Time) %>%
  arrange(Depth) %>%
  mutate(diff = Temperature - Temperature[which.min(Depth)]) %>%
  reframe(mld = Depth[which.max(diff)])

plot(do.call("cbind", approx(subset_test$Density, y = -subset_test$Depth, method = "linear", n = 1000, ties = mean,  na.rm = TRUE)))

mini_profile <- do.call("cbind", approx(subset_test$Density, y = -subset_test$Depth, method = "linear", n = 1000, ties = mean,  na.rm = TRUE))

mini_profile[,1] - mini_profile[which.max(mini_profile[,2]),1]

plot(mini_profile[,1] - mini_profile[which.max(mini_profile[,2]),1], mini_profile[,2])

rbind.fill(ctd_load, rbr_2022) %>%
  filter(Depth > 0) %>%
  mutate(Date = format(Time, "%Y-%m-%d")) %>%
  mutate(File_ID = gsub("\\\\", "-", File_ID)) %>%
  group_by(File_ID) %>%
  reframe(n = n()) %>%
  filter(n == 1) %>% pull(File_ID)


