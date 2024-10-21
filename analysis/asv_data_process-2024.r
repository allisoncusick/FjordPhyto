# Code to process ASV tables

# Load Libraries ----
packages <- c("tidyverse", "dplyr", "tidyr", "readxl", "vegan")


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
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/FjordPhyto/")

# Load in data ----

raw_asv_files <- list.files("data/asv_data/raw_files", full.names = TRUE)
## Function that reads in files listed above and if its a tsv file reads in as such and if its a csv reads in as such

asv_load <- lapply(raw_asv_files, function(r){
  if (grepl(".tsv", r)){
    input <- read_tsv(r)
  }else{
    input <- read.csv(r)
  }
  colnames(input)[c(1, (ncol(input) - 1):ncol(input))] <- c("Feature.ID", "pr2_Taxon", "pr2_Confidence")
  input$file <- r
  output <- input %>%
    pivot_longer(cols = -c("Feature.ID", "pr2_Taxon", "pr2_Confidence", "file"), names_to = "sample", values_to = "reads")
  
  return(output)
}) %>% bind_rows() %>% # remove the double "sample" in the sample name
    mutate(sample = str_replace(sample, "^([^_]*_[^_]*)_.*$", "\\1")) %>%  # remove the strings after the second underscore in the sample names
  mutate(sample = ifelse(grepl("*2022_08_30_*",file) &
                           sample %in% paste0(
                             "ManifestSample_",
                             sprintf('%0.3d', 1:20)),
                         paste0("MdLP_", sample),
                         sample))

# ASV Table Output ----
save(asv_load, file = paste0(
    "data/asv_data/fjord_phyto_ASV-Table_processed-", todays_date, ".Rdata"))


##### Split taxa into known taxonomic groups ----
taxa_table_split <- asv_load %>%
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

## Add in the phytoplankton groups to the ASV table
full_asv <- left_join(asv_load, taxa_table_split, by = "Feature.ID") %>%
  filter(pr2_Confidence > 0.97)

#### Rarefy data set ------
asv_table_wide_2 <- full_asv %>%
  select(Feature.ID, sample, reads) %>%
  pivot_wider(., id_cols = "sample", names_from = "Feature.ID", values_from = "reads") %>%
  replace(is.na(.), 0)

raremax <- 7000 #can change this cut off to something else if justified

## Keep samples (rows) where the sum of all the reads is greater than raremax
samples_keep <- asv_table_wide_2[rowSums(asv_table_wide_2[,-1]) >= raremax,1]$sample
## Remove samples (rows) where the sum of all the reads is less than raremax
samples_removed <- asv_table_wide_2[rowSums(asv_table_wide_2[,-1],) < raremax,1]$sample

## Rarefy using vegan package
Srare <- vegan::rrarefy(asv_table_wide_2[rowSums(asv_table_wide_2[,-1]) >= raremax,-1], raremax)
df_rare <- Srare %>%
  as.data.frame() %>%
  mutate(sample = samples_keep) %>%
  pivot_longer(cols = !sample,
               names_to = "Feature.ID", values_to = "rare_reads")%>%
  replace(is.na(.), 0)

### Remove Feature.IDs with no rarefied reads
present_filter <- df_rare %>%
  group_by(Feature.ID) %>%
  reframe(total_reads = sum(rare_reads),
          total_samples = length(unique(sample[rare_reads>0]))) %>% 
  filter(total_reads > 5 & total_samples > 10) %>%
  pull(Feature.ID) %>%
  unique()

### Join the full asv table with the rarefied data and remove the 
### Feature.IDs with no rarefied reads and samples with less than raremax reads
asv_table_rare <- left_join(df_rare, full_asv,
                            by = c("Feature.ID", "sample")) %>%
  filter(Feature.ID %in% present_filter & !sample %in% samples_removed)

### Run simple diversity metrics for each sample
diversity_group <- c("all", "phytogroups", 
                     unique(taxa_table_split$phytogroups))
diversity_group <- diversity_group[!is.na(diversity_group)]

for (i in 1:length(diversity_group)){
  
  if (diversity_group[i] == "all"){
    taxa_pull <- taxa_table_split
  }else if (diversity_group[i] == "phytogroups"){
    taxa_pull <- taxa_table_split %>% filter(!is.na(phytogroups))
    }else{
    taxa_pull <- taxa_table_split %>% filter(phytogroups == diversity_group[i])
  }
  
  piv_all <- asv_table_rare %>%
    filter(Species %in% taxa_pull$Species) %>%
    dplyr::select(c(Species, sample, rare_reads)) %>%
    pivot_wider(names_from = "Species",
                values_from = "rare_reads",
                values_fn = sum,
                values_fill = 0)
  
  #vegan::specnumber - # species with non-zero reads
  richness <- specnumber(piv_all[,-1]) 
  #vegan::diversity - Shannon diversity index
  shannon <- diversity(piv_all[,-1], MARGIN = 1, index = "shannon")
  #vegan::evenness - H/log(S)
  evenness <- shannon/log(richness)
  #vegan::diversity - Simpson diversity index
  simpson <- diversity(piv_all[,-1], MARGIN = 1, index = "simpson")
  # Relative abundance of each taxa in each sample
  rel_abun <- rowSums(piv_all[,-1])/raremax 
  # total reads
  total_reads <- rowSums(piv_all[,-1])
  
  ### Create column names for each diversity metric that change relative to the
  ### group
  shannon_name <- paste0("shannon_", diversity_group[i])
  even_name <- paste0("evenness_", diversity_group[i])
  simpson_name <- paste0("simpson_", diversity_group[i])
  richness_name <- paste0("richness_", diversity_group[i])
  rel_abun_name <- paste0("rel_abun_", diversity_group[i])
  total_reads_name <- paste0("total_reads_", diversity_group[i])
  
  ### Create a data frame with the sample names and the diversity metrics
  df_output <- data.frame(sample = piv_all[,1],
                          shannon = shannon,
                          evenness = evenness,
                          simpson = simpson,
                          richness = richness,
                          rel_abun = rel_abun,
                          total_reads = total_reads)
  ### Change the column names to the names created above
  colnames(df_output)[2:ncol(df_output)] <- str_replace(
    c(shannon_name, even_name, simpson_name,
      richness_name, rel_abun_name, total_reads_name), " ", "_")
  
  if (i == 1){
    diversity_df <- df_output 
  }else{
    diversity_df <- left_join(diversity_df, df_output, by = "sample")
  }
}





