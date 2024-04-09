#### Script for loading in and cleaning up data 

#### Load in the data -----
#see the README file for how this was re-created
load("~/Documents/GitHub/FjordPhyto/data/WAP_set_meta_rerun.Rdata")
fix_names <- read_csv(paste0(data_home, "site_names_edit_manual.csv"), skip = 0)

#### Metadata ----
#### Change mistake in data
metadata$year[metadata$year == "2020" & metadata$month == "12" & metadata$day == "19"] <- "2019"
metadata$Date <-  as.Date(paste(metadata$day,metadata$month,metadata$year,sep = '/'),format = '%d/%m/%Y')

#make SEASONS brackets 
metadata$month <- factor(as.character(metadata$month), levels = c("11","12","1","2","3","4"))
metadata$season <- "2017-2018"
metadata$season[which(mdy("06-01-2019") > metadata$Date & metadata$Date > mdy("09-01-2018"))] <- "2018-2019"
metadata$season[which(mdy("06-01-2020") > metadata$Date & metadata$Date > mdy("09-01-2019"))] <- "2019-2020"
metadata$season[which(mdy("06-01-2021") > metadata$Date & metadata$Date > mdy("09-01-2020"))] <- "2020-2021" #word of note on 2020-2021 code, there should not be any seasons since it was PANDEMIC year no travel
metadata$season[which(mdy("06-01-2022") > metadata$Date & metadata$Date > mdy("09-01-2021"))] <- "2021-2022"

#make REGIONS brackets "northern" (anything north of -63) "middle" (anything between -63 and -65) "southern" (anything south of -65)
metadata$region <-"shetlands" #assign everything something, start with northern then parse/ID other labels
#metadata$region[which( latitude less than -63 == northern, latitude 63 - 65 = middle, latitude >65 southern)]
metadata$region[which(metadata$Lat < -63 & metadata$Lat > -66)] <-"middle"
metadata$region[which(metadata$Lat < -66 & metadata$Lat > -73)] <-"southern"
metadata$region[which(metadata$Long < -50 & metadata$Long > -58)] <-"northern"

metadata <-  metadata %>%
  mutate(Site_Name_2 = ifelse(
    grepl("Transect", Site_Name) & grepl("Neko", Site_Name),
    "Neko Harbour",
    Site_Name),
    Site_Name_2 = ifelse(Site_Name_2 == "Mikkelson Harbour", "Mikkelsen Harbour", Site_Name_2),
    Site_Name_2 = ifelse(
      grepl("Ple", Site_Name_2),
      "Pleneau Bay",
      Site_Name_2)) %>%
  drop_na(Long, Lat) %>%
  filter(Site_Name_2 != "")

metadata_raw <- metadata


### Change site transect site names and spelling errors ---
metadata <- metadata %>%
  filter(Site_Name_2 != "Whalers Bay") %>%
  mutate(Site_Name_2 = ifelse(
    Site_Name_2 == "Andvord Bay Transect Station 1 Useful",
                              "Useful Island", Site_Name_2),
         Site_Name_2 = ifelse(
           Site_Name_2 == "Andvord Bay Transect Station 2 Errera",
                              "Errera Channel", Site_Name_2),
         Site_Name_2 = ifelse(
           Site_Name_2 == "Andvord Bay Transect Station 3 Middle",
                              "Mid Andvord Bay", Site_Name_2),
         Site_Name_2 = ifelse(
           Site_Name_2 == "Andvord Bay Transect Station 5 Interior",
                              "Bagshawe Glacier", Site_Name_2),
         Site_Name_2 = ifelse(
           Site_Name_2 == "Horshoe Island",
                              "Horseshoe Island", Site_Name_2),
         Site_Name_2 = factor(Site_Name_2, levels = fix_names$Rename),
         site_ids = factor(Site_Name_2, levels = fix_names$Rename, labels = fix_names$site_num)) %>%
  group_by(Site_Name_2) %>%
  mutate(lon = mean(Long),
         lat = mean(Lat),
         region = factor(region, levels = c("southern", "middle", "northern", "shetlands")))


#### Taxonomy ----
##### Split taxa into known taxonomic groups ----
taxa_table_split <- separate(taxa_table, pr2_Taxon, sep = ";",
                             into = c("Kingdom",
                                      "Supergroup",
                                      "Division",
                                      "Class",
                                      "Order",
                                      "Family",
                                      "Genus",
                                      "Species"))

## Change the naming conventions to removes NAs in taxonomy
## (i.e., so each Feature.ID has a string for each taxonomic level)
## The number of Xs represent the number of NA levels between the species level and the top known level.
## For each Feature.ID, figure out where the NAs are, determine the number of Xs depending on
## the number of NAs. Add in the "_sp." if the species name is unknown
for (i in 1:length(unique(taxa_table_split$Feature.ID))){
  taxa_vec <- subset(taxa_table_split, Feature.ID == unique(taxa_table_split$Feature.ID)[i])
  taxa_na <- which(is.na(taxa_vec))
  
  if (length(taxa_na) > 0){
    number_of_X <- 1:length(taxa_na)
    first_taxa <- colnames(taxa_table_split)[taxa_na[1]-1]
    if (grepl("_X", taxa_vec[first_taxa])){
      number_of_X <- number_of_X + 1
    }
    last_taxa <- colnames(taxa_table_split)[taxa_na[length(taxa_na)]]
    X_vector <- unlist(lapply(number_of_X, function(r){paste0(rep("X",r), collapse = "")}))
    new_list <- paste(gsub("_X", "", taxa_vec[first_taxa]), X_vector, sep = "_")
    names(new_list) <- colnames(taxa_table_split)[taxa_na]
    
    taxa_vec[names(new_list)] <- new_list
    taxa_table_split[
      taxa_table_split$Feature.ID == unique(taxa_table_split$Feature.ID)[i],
      names(new_list)] <- new_list
    
    if (last_taxa == "Species"){
      taxa_table_split[
        taxa_table_split$Feature.ID == unique(taxa_table_split$Feature.ID)[i],
        "Species"]  <- paste0(taxa_table_split[
          taxa_table_split$Feature.ID == unique(taxa_table_split$Feature.ID)[i],
          "Genus"], "_sp.")
    } 
    
  }
  
}

##### Define phytoplankton groups ----
## Based on Class and Division 
taxa_table_split <- taxa_table_split %>%
  mutate(
    phytogroups = case_when(
      Class == "Bacillariophyta" ~ "Diatoms",
      Division == "Dinoflagellata" ~ "Dinoflagellates",
      Division %in% c("Chlorophyta", "Prasinodermophyta") ~ "Greenalgae", 
      Division == "Haptophyta" ~ "Haptophytes",
      Division == "Rhodophyta" ~ "Rhodophytes",
      Division == "Cryptophyta"~ "Cryptophytes",
      grepl("MAST-", Class) ~ "MAST"))

## Pivot the asv table and fill in NAs with 0s
## This ensures that every ASV has a value in every sample even if that value
## is 0
asv_table_wide <- pivot_wider(asv_table,
                              id_cols = "Feature.ID",
                              names_from = "samples",
                              values_from = "reads",
                              values_fill = 0)

## Pivot back to long format as a preference
asv_table_fix <- pivot_longer(
  asv_table_wide,
  cols = colnames(asv_table_wide)[-1],
  names_to = "samples",
  values_to = "reads") 


## ASVs that were present in other samples and run with those samples but not present in WAP
not_present <- asv_table_fix %>%
  group_by(Feature.ID) %>%
  reframe(total_samples = length(unique(samples[reads > 0]))) %>%
  filter(total_samples == 0) %>%
  pull(Feature.ID) %>%
  unique()

asv_table_raw <- asv_table_fix %>%
  filter(!Feature.ID %in% not_present)  


#### Rarefy data set ------
asv_table_wide_2 <- pivot_wider(
  asv_table_raw, id_cols = "samples", names_from = "Feature.ID", values_from = "reads")

raremax <- 7000 #can change this cut off to something else if justified 
## Keep samples (rows) where the sum of all the reads is greater than raremax
samples_keep <- asv_table_wide_2[rowSums(asv_table_wide_2[,-1]) >= raremax,1]$samples
## Remove samples (rows) where the sum of all the reads is less than raremax
samples_removed <- asv_table_wide_2[rowSums(asv_table_wide_2[,-1]) < raremax,1]$samples
metadata <- metadata %>%
  filter(samples %in% samples_keep)

## Rarefy using vegan package
Srare <- rrarefy(asv_table_wide_2[rowSums(asv_table_wide_2[,-1]) >= raremax,-1], raremax)
df_rare <- Srare %>%
  as.data.frame() %>%
  mutate(samples = samples_keep) %>%
  pivot_longer(cols = !samples, names_to = "Feature.ID", values_to = "reads")
colnames(asv_table_raw)[3] <- "raw_reads"

asv_table_rare <- left_join(asv_table_raw, df_rare, by = c("Feature.ID", "samples"))

#### File outputs ----
## This creates: metadata and asv_table_rare





