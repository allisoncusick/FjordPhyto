# Quick metadata processing script (Convert from Excel to .Rdata data frame)
# Takes in the metadata file and processes it to output "metadata" data frame

# Load Libraries ----
packages <- c("tidyverse", "dplyr", "tidyr", "readxl")


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
setwd("~/Documents/GitHub/FjordPhyto/")

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
        region = case_when(
            latitude < -63 & latitude > -66 ~ "middle",
            latitude < -66 & latitude > -73 ~ "southern",
            longitude < -50 & longitude > -58 ~ "northern",
            TRUE ~ "shetlands")) %>%
    drop_na(longitude, latitude) %>%
    filter(location != "")

# Metadata Output ----
save(metadata, file = paste0("data/metadata/fjord_phyto_processed-", todays_date, ".Rdata"))
## Load CTD Data
ctd_load <- read_csv("data/ctd_rbr/outputclean.csv")
colnames(ctd_load) <- gsub(" ", "_", colnames(ctd_load))
colnames(ctd_load) <- gsub("\\(.*?\\)", "", colnames(ctd_load))

## Make a simple data frame with all the CTD files and their corresponding ID
ctd_to_ID <- metadata %>%
  select(UNIQUE_ID_CODE, CTD_cast_file_name) %>%
  separate(CTD_cast_file_name, into = c("CTD1", "CTD2", "CTD3"), sep = ",") %>%
  #mutate_at(vars(CTD1, CTD2, CTD3), sub, pattern = "CC", replacement = "") %>%
  mutate_at(
    vars(CTD1, CTD2, CTD3), sub, pattern = "\\..*", replacement = "") %>%
  pivot_longer(cols = c(CTD1, CTD2, CTD3),
  names_to = "CTD", values_to = "CTD_file_name") %>%
  select(-CTD) %>%
  drop_na(CTD_file_name)


#write_csv(ctd_to_ID, "data/ctd_rbr/CTD_ID-list.csv")
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
  reframe(across(everything(), ~list(.x)))
# Now we know there are 294 unique CTD files

# Add in the UNIQUE_ID_CODE using the simple data frame we made above
ctd_w_id <- left_join(ctd_to_ID, ctd_rows, by = c("CTD_file_name" = "filename")) %>%
filter(!is.na(UNIQUE_ID_CODE) & !is.na(latitude)) %>%
distinct(UNIQUE_ID_CODE, .keep_all = T)
# There were 27 UNIQUE ID CODES that have a CTD file associated with them

## Load chl_a data 
find_files <- list.files("data/Chla_Rick/Chla_TimeSeries/", pattern = "*.csv")
chl_load <- do.call("rbind", lapply(find_files,function(r){
  read_csv(paste0("data/Chla_Rick/Chla_TimeSeries/", r))[-1,]
}))
chl_load$time <- as_date(chl_load$time)

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
  mutate(diff_lat = abs(latitude.x - as.numeric(latitude.y)),
         diff_long = abs(longitude.x - as.numeric(longitude.y))) %>%
  group_by(UNIQUE_ID_CODE, Date, latitude.x, longitude.x) %>%
  arrange(diff_lat, diff_long) %>%
  filter(row_number()==1) %>%
  select(UNIQUE_ID_CODE, Date, latitude.x, longitude.x, chlor_a) %>%
  rename(latitude = latitude.x, longitude = longitude.x)

## Add the chla data to the metadata
metadata_w_chl <- left_join(metadata, meta_match_chla,
  by = c("UNIQUE_ID_CODE", "date_local" = "Date", "latitude", "longitude"))

# Use a linear model to predict chlorophull from secchi depth
chl_fit <- lm(as.numeric(chlor_a) ~ as.numeric(secchi_depth), data = metadata_w_chl %>% filter(year == 2018))
secchi_2017 <- metadata_w_chl %>%
filter(year == 2017) %>%
pull(secchi_depth)

metadata_w_chl$chlor_a[metadata_w_chl$year == 2017] <- predict(chl_fit, data.frame(secchi_depth = secchi_2017))


## Add the CTD data to the metadata
metadata_full <- left_join(metadata_w_chl,
  ctd_w_id %>% select(-c(latitude,longitude)),
  by = "UNIQUE_ID_CODE")

metadata_full %>%
select(UNIQUE_ID_CODE, date_local, secchi_depth, Genetics_18sv9_Sample_ID, filename.y) %>%
View()

metadata_w_samples <- metadata_full %>%
filter(!is.na(Genetics_18sv9_Sample_ID)) %>%
filter(!is.na(CTD_file_name))



