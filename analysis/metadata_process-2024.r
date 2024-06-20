# Quick metadata processing script (Convert from Excel to .Rdata data frame)
# Takes in the metadata file and processes it to output "metadata" data frame

# Load Libraries ----
packages <- c("tidyverse", "dplyr", "tidyr", "readxl")


funlist <-  lapply(packages, function(x) {
  if (x %in% rownames(installed.packages())) {
    require(x, character.only = TRUE)
  } else {
    install.packages(x, character.only = TRUE); require(x, character.only = TUE)
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

ctd_to_ID <- metadata %>%
  select(UNIQUE_ID_CODE, CTD_cast_file_name) %>%
  separate(CTD_cast_file_name, into = c("CTD1", "CTD2", "CTD3"), sep = ",") %>%
  mutate_at(vars(CTD1, CTD2, CTD3), sub, pattern = "CC", replacement = "") %>%
  mutate_at(
    vars(CTD1, CTD2, CTD3), sub, pattern = "\\..*", replacement = "") %>%
  pivot_longer(cols = c(CTD1, CTD2, CTD3),
  names_to = "CTD", values_to = "CTD_file_name") %>%
  select(-CTD) %>%
  drop_na(CTD_file_name)

#write_csv(ctd_to_ID, "data/ctd_rbr/CTD_ID-list.csv")

add_id_df <- ctd_load %>%
  mutate(date = as_date(timestamp)) %>%
  group_by(date, filename, latitude, longitude) %>%
  reframe(across(everything(), ~map_if(.x, length(unique(.x)) > 1, ~list(.x)))) %>%
  left_join(., ctd_to_ID, by = c("filename" = "CTD_file_name")) %>%
  left_join(metadata, ., by = c("UNIQUE_ID_CODE"))

## Load chl_a data 
find_files <- list.files("data/Chla_Rick/Chla_TimeSeries/", pattern = "*.csv")
chl_load <- do.call("rbind", lapply(find_files,function(r){
  read_csv(paste0("data/Chla_Rick/Chla_TimeSeries/", r))[-1,]
}))
chl_load$time <- as_date(chl_load$time)

subset_to_chla <- metadata %>%
  mutate(Date = as_date(date_local)) %>%
  group_by(UNIQUE_ID_CODE, Date, latitude, longitude, secchi_depth) %>%
  group_keys() %>%
  left_join(., chl_load, by = join_by(Date == time)) %>%
  mutate(diff_lat = abs(latitude.x - as.numeric(latitude.y)),
         diff_long = abs(longitude.x - as.numeric(longitude.y))) %>%
  group_by(UNIQUE_ID_CODE, Date, latitude.x, longitude.x, secchi_depth) %>%
  arrange(diff_lat, diff_long) %>%
  slice(1) %>%
  select(UNIQUE_ID_CODE, Date, latitude.x, longitude.x, chlor_a, secchi_depth)

# subset_to_chla %>%
# rename(latitude = latitude.x, longitude = longitude.x) %>%
# write_csv(., "data/Chla_Rick/chla_ID-list.csv")

add_chla_ctd <- add_id_df %>%
  left_join(., subset_to_chla,
  by = c("UNIQUE_ID_CODE" = "UNIQUE_ID_CODE")) %>%
  select(-matches("\\.y$")) %>%
  rename(latitude = latitude.x.x, longitude = longitude.x.x) %>%
  group_by(UNIQUE_ID_CODE)
  
unique(unlist(add_chla_ctd[7, "MLD 1"]))


