# Quick metadata processing script (Convert from Excel to .Rdata data frame)
# Takes in the metadata file and processes it to output "metadata" data frame

# Load Libraries ----
packages <- c("tidyverse", "dplyr", "tidyr", "readxl", "geosphere")


funlist <-  lapply(packages, function(x) {
  if (x %in% rownames(installed.packages())) {
    require(x, character.only = TRUE)
  } else {
    install.packages(x, character.only = TRUE); require(x, character.only = TRUE)
  }
})
todays_date <- format(Sys.Date(), "%m%d%Y")
options(max.print = 100)
# Set Working Directory ----
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/FjordPhyto/")

# Load and process metadata -----
metadata_raw <- read_excel(
    "data/metadata/FjordPhyto-MASTER-accessed_06122024.xlsx",
    sheet = 1, na = c(" ", "NA", "#N/A", "N/A", "-999")) %>%
    as_tibble()

# Change column names to match with R conventions
colnames(metadata_raw) <- gsub(" ", "_", colnames(metadata_raw))

# Add in columns for Antarctic season ID, season year, region names, and drop NA
metadata <- metadata_raw %>%
    mutate(
        season = case_when(
            month %in% c("11", "12", "01") ~ "summer",
            month %in% c("02", "03", "04") ~ "fall",
            month %in% c("05", "06", "07") ~ "winter",
            month %in% c("08", "09", "10") ~ "spring"),
        seasonyear = case_when(
            season == "summer" & month == "01" ~ paste0(
                as.numeric(year) - 1, "-", as.numeric(year)),
            season == "summer" & month %in% c("11", "12") ~ paste0(
                as.numeric(year), "-", as.numeric(year) + 1),
            season == "fall" ~ paste0(
                as.numeric(year) - 1, "-", as.numeric(year)),
            season == "winter" ~ paste0(
                as.numeric(year) - 1, "-", as.numeric(year)),
            season == "spring" ~ paste0(
                as.numeric(year), "-", as.numeric(year) + 1)),
          days_since = yday(date_local) - yday("2000-11-01"),
          days_since = ifelse(days_since < 0, days_since + 365, days_since)) %>%
    drop_na(longitude, latitude)

metadata$region <-"shetlands"
metadata$region[
  which(metadata$latitude < -63 & metadata$latitude > -66)] <-"middle"
metadata$region[
  which(metadata$latitude < -66 & metadata$latitude > -73)] <-"southern"
metadata$region[
  which(metadata$longitude < -50 & metadata$longitude > -58)] <-"northern"

dist_vec <- matrix(ncol = 2, nrow = length(unique(metadata$location)))
for (r in 1:length(unique(metadata$location))) {
  lon1 <- unique(metadata$longitude[metadata$location == "Stonington Island"])
  lat1 <- unique(metadata$latitude[metadata$location == "Stonington Island"])
  lon2 <- mean(unique(metadata$longitude[metadata$location == unique(metadata$location)[r]]))
  lat2 <- mean(unique(metadata$latitude[metadata$location == unique(metadata$location)[r]]))
  dist_vec[r,2] <- distm(c(lon1, lat1), c(lon2, lat2), fun = distHaversine)
  dist_vec[r,1] <- unique(metadata$location)[r]
}

dist_df <- arrange(as.data.frame(dist_vec), V2) %>%
  mutate(site_num = 1:nrow(dist_vec))
colnames(dist_df) <- c("location", "distance_to_stonington", "site_num")

metadata$site_ids <- factor(metadata$location, levels = dist_df$location, labels = dist_df$site_num)

# metadata$land_dist_km <- dist2land(metadata[,c("latitude", "longitude")], bind = F)

# Metadata Output ----
save(metadata, file = paste0("data/metadata/fjord_phyto_processed-", todays_date, ".Rdata"))


## Load CTD Data
ctd_load <- read_csv("data/ctd_rbr/ctdfjorder_data_20240728212609_nik.csv")
colnames(ctd_load) <- gsub(" ", "_", colnames(ctd_load))
colnames(ctd_load) <- gsub("\\(.*?\\)", "", colnames(ctd_load))
colnames(ctd_load) <- gsub("-", "", colnames(ctd_load))

### Figure out which variables are raw values (i.e. measured with depth and time) versus computed values (i.e. one value for each CTD run)
unique_keys <- ctd_load %>% 
  mutate(date = as_date(timestamp)) |> 
  filter(filename == unique(filename)[1]) %>%
  summarise(across(everything(), ~length(unique(.x)))) %>% 
  select(where(~.x == 1)) %>% names()

### Turn the full data frame into one where each row is a unique CTD run with the computed values and turn the depth columns into list columns 
ctd_rows <- ctd_load |> 
  mutate(date = as_date(timestamp)) |> 
  group_by_at(unique_keys) %>%
  reframe(across(everything(), ~list(.x))) %>%
  rename_with(~paste0(.x, "_ctd"), 
              all_of(names(ctd_load)[names(ctd_load) %in% names(metadata)])) %>%
  drop_na(Surface_Salinity_) %>%
  filter(Profile_ID == 0) %>%
  distinct(Unique_ID, .keep_all = TRUE)

# Now we know there are 293 unique CTD files
metadata_w_ctd <- left_join(ctd_rows, metadata,
                            by = c("Unique_ID" = "UNIQUE_ID_CODE")) %>%
  filter(!is.na(time_local))

metadata_w_samples <- left_join(ctd_rows, metadata,
                                by = c("Unique_ID" = "UNIQUE_ID_CODE")) %>%
  filter(!is.na(Genetics_18sv9_Sample_ID) &
         longitude != 0 & !is.na(time_local)) 

### There are 150 unique samples with CTD casts, secchi depths, and associated metadata
# allison_output <- asv_load %>%
#   filter(sample %in% paste0("ManifestSample_", sprintf('%0.3d', 93:102))) %>%
# left_join(., metadata_w_samples, by = c("sample" = "Genetics_18sv9_Sample_ID")) %>%
#   mutate_if(is.list, unlist)
# 
# write.csv(allison_output, "data/arctic_output.csv", row.names = FALSE)
