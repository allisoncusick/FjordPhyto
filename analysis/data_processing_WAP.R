library(tidyverse)
library(lubridate)

###if you wanna clear and start over rm(list=ls())

# data processing file

# all metadata 2017-2020 samples 1-200 in a submission during 2020

meta_2020_1_200 <- read.csv("data/metadata/2020_pool_18S.csv")

#meta_2020_1_200 <- meta_1_200 %>% filter(sample.id != "")
#is that code above needed? 

#all metadata mixed 2019 2021 samples 1 - 190 in a submission during 2022

meta_2022_1_190 <- read.csv("data/metadata/2022_08_30_pool1_18S.csv")

meta_2022_1_190 <- meta_2022_1_190 %>% filter(project_name == "WAP")

#no need to fix dates anymore
#early_meta$Date <- dmy(early_meta$Date)

#meta_20_22$date <- dmy(meta_20_22$date)

meta_2022_extra <- read.csv("data/metadata/2022_10_18_pool1_18S.csv")

meta_2022_extra <- meta_2022_extra %>% filter(project_name == "WAP")
#a dino sample pops up ... filter out? 

#no need to change dates, did it in original file
#meta_19_22_extra$date <-dmy(meta_19_22_extra$date)


# LOAD IN DATA

# 17-18

asv_table_17_18 <- read.csv("data/asv_count_tax1718.csv")

asv_table_17_18 <- asv_table_17_18 %>% filter(pr2_Confidence > 0.97)

colnames(asv_table_17_18)
taxa_17_18 <- asv_table_17_18[,c(1,89,90)]
asv_table_17_18 <- asv_table_17_18[,-c(89:90)]

asv_table_17_18 <- asv_table_17_18 %>% pivot_longer(-Feature.ID, names_to = "samples", values_to = "reads")

# 18-19

asv_table_18_19 <- read.csv("data/asv_count_tax1819.csv")

asv_table_18_19 <- asv_table_18_19 %>% filter(pr2_Confidence > 0.97)

colnames(asv_table_18_19)
taxa_18_19 <- asv_table_18_19[,c(1,82,83)]
asv_table_18_19 <- asv_table_18_19[,-c(82:83)]

asv_table_18_19 <- asv_table_18_19 %>% pivot_longer(-Feature.ID, names_to = "samples", values_to = "reads")

# 19-20

asv_table_19_20 <- read.csv("data/asv_count_tax1920.csv")

asv_table_19_20 <- asv_table_19_20 %>% filter(Confidence > 0.97)

colnames(asv_table_19_20)
taxa_19_20 <- asv_table_19_20[,c(1,35,36)]
asv_table_19_20 <- asv_table_19_20[,-c(35:36)]

asv_table_19_20 <- asv_table_19_20 %>% pivot_longer(-Feature.ID, names_to = "samples", values_to = "reads")

# 19-20/21-22 from the 08_30_22_p1_18S run

asv_table_20_22 <- read.csv("data/asv_table_2022_seq.csv")
colnames(asv_table_20_22)[1] <- "Feature.ID"
colnames(asv_table_20_22)

taxa_20_22 <- read.csv("data/pr2_taxonomy_2022_seq.csv")

taxa_20_22 <- taxa_20_22 %>% filter(Confidence > 0.97)

asv_table_20_22 <- asv_table_20_22 %>% pivot_longer(-Feature.ID, names_to = "samples", values_to = "reads")

asv_table_20_22 <- asv_table_20_22 %>% filter(Feature.ID %in% taxa_20_22$Feature.ID)

# filter this table to just WAP samples

asv_table_20_22$samples <- gsub("^([^_]*_[^_]*)_.*$", "\\1", asv_table_20_22$samples)

asv_table_20_22 <- asv_table_20_22 %>% filter(samples %in% meta_2022_1_190$sample.id)

taxa_20_22 <- taxa_20_22 %>% filter(Feature.ID %in% asv_table_20_22$Feature.ID)

# add 19-20/21-22 from the 10_18_22_p1_18S run, ADD TO THE OTHER TABLE ABOVE

asv_table_20_22_extra <- read.csv("data/asv_count_tax_10_18_22_p1_18S.csv")
#colnames(asv_table_20_22_extra)[1] <- "Feature.ID"
asv_table_20_22_extra <-asv_table_20_22_extra[,-c(23,24)]
colnames(asv_table_20_22_extra)
taxa_20_22_extra <- read.csv("data/pr2_taxonomy_2022_10_18_p1_18S.csv")

asv_table_20_22_extra <- asv_table_20_22_extra %>% pivot_longer(-Feature.ID, names_to = "samples", values_to = "reads")

taxa_20_22_extra <- taxa_20_22_extra %>% filter(Confidence > 0.97)
colnames(asv_table_20_22_extra)

asv_table_20_22_extra <- asv_table_20_22_extra %>% filter(Feature.ID %in% taxa_20_22_extra$Feature.ID)


# Combine datasets and fix sample labels

asv_table_17_18$samples <- gsub("sample\\.","",asv_table_17_18$samples) 
asv_table_17_18$samples <- gsub("^([^_]*_[^_]*)_.*$", "\\1", asv_table_17_18$samples)

asv_table_18_19$samples <- gsub("sample\\.","",asv_table_18_19$samples) 
asv_table_18_19$samples <- gsub("^([^_]*_[^_]*)_.*$", "\\1", asv_table_18_19$samples)

asv_table_19_20$samples <- gsub("sample\\.","",asv_table_19_20$samples) 
asv_table_19_20$samples <- gsub("^([^_]*_[^_]*)_.*$", "\\1", asv_table_19_20$samples)

asv_table_20_22_extra$samples <- gsub("sample\\.","",asv_table_20_22_extra$samples) 
asv_table_20_22_extra$samples <- gsub("^([^_]*_[^_]*)_.*$", "\\1", asv_table_20_22_extra$samples)

#combine asv_table_20_22 and asv_table_20_22_extra
#then bind_rows including the update

asv_table <- bind_rows(asv_table_17_18,asv_table_18_19, asv_table_19_20, asv_table_20_22, asv_table_20_22_extra)
#does that look funny? 

# change taxa column names

colnames(taxa_19_20) <- colnames(taxa_17_18)
colnames(taxa_20_22) <- colnames(taxa_17_18)
colnames(taxa_20_22_extra) <- colnames(taxa_17_18)
taxa_table <- bind_rows(taxa_17_18, taxa_18_19, taxa_19_20, taxa_20_22, taxa_20_22_extra)

#filtered anything that is a repeat
taxa_table <- taxa_table %>% distinct_all()

# metadata pull out desired columns and relabel them what i want
colnames(meta_2020_1_200)
meta_2020_1_200 <- meta_2020_1_200[,c(1,7:12)]

colnames(meta_2020_1_200) <- c("sample_id","Date","Time","Site_Name","Lat","Long","Operator")

colnames(meta_2022_1_190)
meta_2022_1_190 <- meta_2022_1_190[,c(1,7:12)]
colnames(meta_2022_1_190) <- c("sample_id","Date","Time","Site_Name","Lat","Long","Operator")

colnames(meta_2022_extra)
meta_2022_extra <- meta_2022_extra[,c(1, 7:12)]
colnames(meta_2022_extra) <- c("sample_id","Date","Time","Site_Name","Lat","Long","Operator")

metadata <- bind_rows(meta_2020_1_200, meta_2022_1_190, meta_2022_extra)
#why is line 201 - 426 missing SO MUCH DATA?? 

#check the site names as they're all OVER THE PLACE, label as desired

metadata$Site_Name <- str_to_title(metadata$Site_Name)
metadata$Site_Name
#seems ok except the lines 201 - 426 missing 
#i think the code below isn't needed any longer cuz I changed them all in the original file ... 

#metadata$Site_Name[which(metadata$Site_Name == "Danco")] <- "Danco Island"
#metadata$Site_Name[which(metadata$Site_Name == "Cuverville")] <- "Cuverville Island"
#metadata$Site_Name[which(metadata$Site_Name == "Paradise Harbour")] <- "Paradise Harbour"
#metadata$Site_Name[which(metadata$Site_Name == "Paradise Bay")] <- "Paradise Harbour"
#metadata$Site_Name[which(metadata$Site_Name == "Paradise Bay/Skornthrop Cove")] <- "Paradise Harbour"
#metadata$Site_Name[which(metadata$Site_Name == "Base Brown")] <- "Paradise Harbour"
#metadata$Site_Name[which(metadata$Site_Name == "Orne")] <- "Orne Harbour"
#metadata$Site_Name[which(metadata$Site_Name == "Peterman Island")] <- "Petermann Island"
#metadata$Site_Name[which(metadata$Site_Name == "Plenneau")] <- "Plenneau Bay"
#metadata$Site_Name[which(metadata$Site_Name == "Pleaneau")] <- "Plenneau Bay"
#metadata$Site_Name[which(metadata$Site_Name == "Halfmoon")] <- "Halfmoon Island"
#metadata$Site_Name[which(metadata$Site_Name == "Half-Moon Island")] <- "Halfmoon Island"
#metadata$Site_Name[which(metadata$Site_Name == "Kinnes Cove/Madder Cliffs")] <- "Kinnes Cove"
#metadata$Site_Name[which(metadata$Site_Name == "Mikkelson Harbour")] <- "Mikkelsen Harbour"

#metadata$Site_Name[which(metadata$Site_Name == "Andvord Bay Station 1 Gerlache Useful ")] <- "Andvord Bay Transect Station 1 Useful"
#metadata$Site_Name[which(metadata$Site_Name == "Andvord Bay Transect 1 -Useful")] <- "Andvord Bay Transect Station 1 Useful"

#metadata$Site_Name[which(metadata$Site_Name == "Andvord Bay Transect 2 -Errera" )] <- "Andvord Bay Transect Station 2 Errera"
#metadata$Site_Name[which(metadata$Site_Name == "Andvord Bay Station 2 Errera ")] <- "Andvord Bay Transect Station 2 Errera"

#metadata$Site_Name[which(metadata$Site_Name == "Andvord Bay Transect 3 -Mid"  )] <- "Andvord Bay Transect Station 3 Middle"
#metadata$Site_Name[which(metadata$Site_Name == "Andvord Bay Station 3 Mid " )] <- "Andvord Bay Transect Station 3 Middle"

#metadata$Site_Name[which(metadata$Site_Name == "Anvord Bay Transect 4-Neko Harbour")] <- "Anvord Bay Transect Station 4 Neko Harbour"
#metadata$Site_Name[which(metadata$Site_Name == "Andvord Bay Station 5 Neko ")] <- "Anvord Bay Transect Station 4 Neko Harbour"

#metadata$Site_Name[which(metadata$Site_Name == "Andvord Bay Transect 5-Int")] <- "Andvord Bay Transect Station 5 Interior"
#metadata$Site_Name[which(metadata$Site_Name == "Andvord Bay Station 4 Int ")] <- "Andvord Bay Transect Station 5 Interior"

metadata$Site_Name[which(metadata$Site_Name == "")] <- NA

sites <- metadata %>% group_by(Site_Name) %>%
  summarise(Lat = mean(Lat, na.rm = TRUE),
            Long = mean(Long, na.rm = TRUE))

#need to fix the sites NA GPS data, I added lat/long to original csv file loaded .. still NA .. ignore those
#i went into the original files and changed the info so its all correct.

metadata$Lat <- sites$Lat[match(metadata$Site_Name,sites$Site_Name)]
metadata$Long <- sites$Long[match(metadata$Site_Name,sites$Site_Name)]


#save my lat long problem then save this file this will be a new file and use code load("data/WAP_clean.Rdata")
save(asv_table,taxa_table,metadata, file = "data/WAP_clean.Rdata")




