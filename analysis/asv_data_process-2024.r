# Code to process ASV tables

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

# Load in data ----
# 2017 - 2020 have one raw file,
# 2021 - 2022 have two raw files

files_with_tax <- list.files("data/asv_data/", pattern = "asv_count_tax.*\\.csv$", full.names = TRUE)

files_wo_tax <- list.files("data/asv_data/", pattern = "asv_table.*\\.csv$", full.names = TRUE)
tax_files_add <- list.files("data/asv_data/", pattern = "pr2_taxonomy.*\\.csv$", full.names = TRUE)
tax_files_match <- tax_files_add[
    substr(tax_files_add, nchar(tax_files_add) - 11, nchar(tax_files_add)) %in%
    substr(files_wo_tax, nchar(files_wo_tax) - 11, nchar(files_wo_tax))]

## With taxonomy inluded ---
lapply(files_with_tax, function(r){
    input <- read.csv(r)
    colnames(input)[c(1, (ncol(input) - 1):ncol(input))] <- c("Feature.ID", "pr2_Taxon", "pr2_Confidence")
    input$file <- r
    output <- input %>%
    pivot_longer(cols = -c("Feature.ID", "pr2_Taxon", "pr2_Confidence", "file"), names_to = "sample", values_to = "reads")
    return(output)
}) %>% bind_rows() -> asv_with_tax

## Without taxonomy included ---
lapply(seq_along(files_wo_tax), function(r){
    input <- read.csv(files_wo_tax[r])
    colnames(input)[1] <- c("Feature.ID")
    input$file <- files_wo_tax[r]

    taxa_input <- read.csv(tax_files_match[r])
    colnames(taxa_input)[2:3] <- c("pr2_Taxon", "pr2_Confidence")
    input_merge <- left_join(input, taxa_input, by = "Feature.ID")
    output <- input_merge %>%
        pivot_longer(cols = -c("Feature.ID", "pr2_Taxon", "pr2_Confidence", "file"), names_to = "sample", values_to = "reads")
    return(output)
}) %>% bind_rows() -> asv_wo_tax


full_asv <- bind_rows(asv_with_tax, asv_wo_tax) %>% # bind the two together
    mutate(sample = str_replace(sample, ".*\\.", "")) %>% # remove the double "sample" in the sample name
    mutate(sample = str_replace(sample, "^([^_]*_[^_]*)_.*$", "\\1")) %>% # remove the strings after the second underscore in the sample names
    filter(pr2_Confidence > 0.97) # filter only samples with > 97% confidence in taxa

# ASV Table Output ----
save(full_asv, file = paste0(
    "data/asv_data/fjord_phyto_ASV-Table_processed-", todays_date, ".Rdata"))


##### Split taxa into known taxonomic groups ----
taxa_table_split <- full_asv %>%
    select(Feature.ID, pr2_Taxon) %>%
    distinct() %>%
    separate(., pr2_Taxon, sep = ";",
            into = c("Kingdom",
                    "Supergroup",
                    "Division",
                    "Class",
                    "Order",
                    "Family",
                    "Genus",
                    "Species")) %>%
    as.data.frame()

## Change the naming conventions to removes NAs in taxonomy
## (i.e., so each Feature.ID has a string for each taxonomic level)
## The number of Xs represent the number of NA levels between the species level and the top known level.
## For each Feature.ID, figure out where the NAs are, determine the number of Xs depending on
## the number of NAs. Add in the "_sp." if the species name is unknown
for (i in seq_along(unique(taxa_table_split$Feature.ID))){
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


#### Rarefy data set ------
asv_table_wide_2 <- full_asv %>%
select(Feature.ID, sample, reads) %>%
  dplyr::group_by(sample, Feature.ID) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L) 
pivot_wider(., id_cols = "sample", names_from = "Feature.ID", values_from = "reads")

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
