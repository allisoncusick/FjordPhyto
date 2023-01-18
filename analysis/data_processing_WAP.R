library(tidyverse)
library(lubridate)

# data processing file

# all metadata 2017-2020

early_meta <- read.csv("data/Metadata_2017_2020.csv")

early_meta <- early_meta %>% filter(SampleName_18sv9_fastq != "")

meta_20_22 <- read.csv("data/Metadata_2022.csv")

meta_20_22 <- meta_20_22 %>% filter(project_name == "WAP")

early_meta$Date <- dmy(early_meta$Date)

meta_20_22$date <- c(dmy(meta_20_22$date[1:39]),mdy(meta_20_22$date[40:62]))

# LOAD IN DATA

# 17-18

asv_table_17_18 <- read.csv("data/asv_count_tax1718.csv")

asv_table_17_18 <- asv_table_17_18 %>% filter(pr2_Confidence > 0.97)

taxa_17_18 <- asv_table_17_18[,c(1,89,90)]
asv_table_17_18 <- asv_table_17_18[,-c(89:90)]

asv_table_17_18 <- asv_table_17_18 %>% pivot_longer(-Feature.ID, names_to = "samples", values_to = "reads")

# 18-19

asv_table_18_19 <- read.csv("data/asv_count_tax1819.csv")

asv_table_18_19 <- asv_table_18_19 %>% filter(pr2_Confidence > 0.97)

taxa_18_19 <- asv_table_18_19[,c(1,82,83)]
asv_table_18_19 <- asv_table_18_19[,-c(82:83)]

asv_table_18_19 <- asv_table_18_19 %>% pivot_longer(-Feature.ID, names_to = "samples", values_to = "reads")

# 19-20

asv_table_19_20 <- read.csv("data/asv_count_tax1920.csv")

asv_table_19_20 <- asv_table_19_20 %>% filter(Confidence > 0.97)

taxa_19_20 <- asv_table_19_20[,c(1,35,36)]
asv_table_19_20 <- asv_table_19_20[,-c(35:36)]

asv_table_19_20 <- asv_table_19_20 %>% pivot_longer(-Feature.ID, names_to = "samples", values_to = "reads")

# 19-20/21-22

asv_table_20_22 <- read.csv("data/asv_table_2022_seq.csv")
colnames(asv_table_20_22)[1] <- "Feature.ID"

taxa_20_22 <- read.csv("data/pr2_taxonomy_2022_seq.csv")

asv_table_20_22 <- asv_table_20_22 %>% pivot_longer(-Feature.ID, names_to = "samples", values_to = "reads")

taxa_20_22 <- taxa_20_22 %>% filter(Confidence > 0.97)

asv_table_20_22 <- asv_table_20_22 %>% filter(Feature.ID %in% taxa_20_22$Feature.ID)

# filter this table to just WAP samples

asv_table_20_22$samples <- gsub("^([^_]*_[^_]*)_.*$", "\\1", asv_table_20_22$samples)

asv_table_20_22 <- asv_table_20_22 %>% filter(samples %in% meta_20_22$sample.id)

taxa_20_22 <- taxa_20_22 %>% filter(Feature.ID %in% asv_table_20_22$Feature.ID)

# Combine datasets

asv_table_17_18$samples <- gsub("sample\\.","",asv_table_17_18$samples) 
asv_table_17_18$samples <- gsub("^([^_]*_[^_]*)_.*$", "\\1", asv_table_17_18$samples)

asv_table_18_19$samples <- gsub("sample\\.","",asv_table_18_19$samples) 
asv_table_18_19$samples <- gsub("^([^_]*_[^_]*)_.*$", "\\1", asv_table_18_19$samples)

asv_table_19_20$samples <- gsub("sample\\.","",asv_table_19_20$samples) 
asv_table_19_20$samples <- gsub("^([^_]*_[^_]*)_.*$", "\\1", asv_table_19_20$samples)

asv_table <- bind_rows(asv_table_17_18,asv_table_18_19, asv_table_19_20, asv_table_20_22)

# change taxa column names

colnames(taxa_19_20) <- colnames(taxa_17_18)
colnames(taxa_20_22) <- colnames(taxa_17_18)

taxa_table <- bind_rows(taxa_17_18, taxa_18_19, taxa_19_20, taxa_20_22)

taxa_table <- taxa_table %>% distinct_all()

# metadata 

early_meta <- early_meta[,c(1,7:12)]
colnames(early_meta) <- c("sample_id","Date","Time","Site_Name","Lat","Long","Operator")

meta_20_22 <- meta_20_22[,c(1,7:12)]
colnames(meta_20_22) <- c("sample_id","Date","Time","Site_Name","Lat","Long","Operator")

metadata <- bind_rows(early_meta, meta_20_22)

metadata$Site_Name <- str_to_title(metadata$Site_Name)

metadata$Site_Name[which(metadata$Site_Name == "Danco")] <- "Danco Island"
metadata$Site_Name[which(metadata$Site_Name == "Cuverville")] <- "Cuverville Island"
metadata$Site_Name[which(metadata$Site_Name == "Paradise Harbour")] <- "Paradise Bay"
metadata$Site_Name[which(metadata$Site_Name == "Paradise Bay/Skornthrop Cove")] <- "Paradise Bay"
metadata$Site_Name[which(metadata$Site_Name == "Base Brown")] <- "Paradise Bay"
metadata$Site_Name[which(metadata$Site_Name == "Orne")] <- "Orne Harbour"
metadata$Site_Name[which(metadata$Site_Name == "Peterman Island")] <- "Petermann Island"
metadata$Site_Name[which(metadata$Site_Name == "Plenneau")] <- "Plenneau Bay"
metadata$Site_Name[which(metadata$Site_Name == "Pleaneau")] <- "Plenneau Bay"
metadata$Site_Name[which(metadata$Site_Name == "Halfmoon")] <- "Halfmoon Island"
metadata$Site_Name[which(metadata$Site_Name == "Half-Moon Island")] <- "Halfmoon Island"
metadata$Site_Name[which(metadata$Site_Name == "Kinnes Cove/Madder Cliffs")] <- "Kinnes Cove"
metadata$Site_Name[which(metadata$Site_Name == "Mikkelson Harbour")] <- "Mikkelsen Harbour"

sites <- metadata %>% group_by(Site_Name) %>%
  summarise(Lat = mean(Lat, na.rm = TRUE),
            Long = mean(Long, na.rm = TRUE))

metadata$Lat <- sites$Lat[match(metadata$Site_Name,sites$Site_Name)]
metadata$Long <- sites$Long[match(metadata$Site_Name,sites$Site_Name)]

save(asv_table,taxa_table,metadata, file = "data/WAP_clean.Rdata")




