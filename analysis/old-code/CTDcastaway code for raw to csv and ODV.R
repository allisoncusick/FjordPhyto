####CODE IDEAL script for processing raw CastAwayCTD files, applying instrument drift correction, merging into one df

library(dplyr)
library(stringr)
library(tidyverse)
library(fs)

# 1. have the folder where all raw csv files will go from instrument. 
#need to read multiple files CastAwayCTD .csv files into R
setwd("~/Documents/Fjord_Phyto/FjordPhyto/data/castawayCTD /down/down_raw/")


#what is the point of the list.files below? to list all files in folder with extension CSV
files<-list.files(full.names = TRUE, include.dirs = TRUE)
head(files)

file_names<-list.files(full.names=FALSE, include.dirs = FALSE)
head(file_names)

#2. create a loop that will look inside each file at line 5 %sample type and if it says "invalid" or "point cast" do not use, 
#if it says "Cast" use that file
for(i in 1:length(file_names)){
  
  meta<-read.csv(files[i],sep=",",nrow=27)
  ctd<-read.csv(files[i],sep=",",skip = 28,header = TRUE)
  ctd$datetime_utc<-as.character(meta[3,2])
  ctd$file_id<-as.character(meta[2,2])
  ctd$sample_type<-as.character(meta[5,2])
  ctd$latitude<-as.character(meta[9,2])
  ctd$longitude<-as.character(meta[10,2])
  setwd("~/Documents/Fjord_Phyto/FjordPhyto/data/castawayCTD /down/down_raw/")
  write.csv(x = ctd, row.names=FALSE, file = paste("mod",file_names[i],sep = "_"))
} 

#4. If line 10 "%latitude"  and line 11 "% longitude" are empty, need to go to mastersheet to hand enter GPS where these were cast were taken (recorded into mastersheet from paper data sheet)
#5. then make new columns and fill with: create a loop that will change all files and add file_id, datetime_utc, latitude, longitude to columns
for(i in 1:length(file_names)){
  setwd("~/Documents/Fjord_Phyto/FjordPhyto/data/castawayCTD /down/down_raw/")
  meta<-read.csv(files[i],sep=",",nrow=27)
  ctd<-read.csv(files[i],sep=",",skip = 28,header = TRUE)
  ctd$datetime_utc<-as.character(meta[3,2])
  ctd$file_id<-as.character(meta[2,2])
  ctd$sample_type<-as.character(meta[5,2])
  ctd$latitude<-as.character(meta[9,2])
  ctd$longitude<-as.character(meta[10,2])
  setwd("~/Documents/Fjord_Phyto/FjordPhyto/data/castawayCTD /down/down_raw/")
  write.csv(x = ctd, row.names=FALSE, file = paste("mod",file_names[i],sep = "_"))
} 

#combining all the csv files into one file
ctd<-lapply(Sys.glob("./mod for ODV/mod*.csv"),read.csv)
#this combines them all into one file
ODVmasterCTDdown<-do.call(rbind,ctd)
write.csv(ODVmasterCTDdown, file = "mod for ODV/ODVmasterCTDdown.csv")

#run through the code, apply the instrument drift correction equation by pulling in the metadata thing
# produce corrected files then merge all into 1 giant df, not all instruments have corrections ... mark with an asterix? 
