#Read in clean 18S taxa asv metadata files from 'data_processing_WAP.R' 
load("data/WAP_clean.Rdata")
#THE HOLY GRAIL
#NEVER HAVE TO PROCESS DATA AGAIN!!!!

#Bring in other data sets

#1. READ IN CTD FILE
CTD <-read.csv("data/All_final_CTD.csv")
head(CTD)

#2. READ IN MASTERSHEET (this is downloaded from https://docs.google.com/spreadsheets/d/1UfYEAQ_BIWsEqeNQmZGTfzNZz1alBYffue3LVZ9NBeE/edit#gid=0
MASTERSHEET <- read.csv("data/FjordPhyto MASTER SHEET_pulled5april2023.csv")
#see what columns are in MASTERSHEET
colnames(MASTERSHEET)
#keep only columns I want: "Season.Year" "Season.Phase" "YYYY" "MM" "DD"  "DD_MON_YY"  "TIME.LOCAL"                     
                         #  "LOC.ID"  "Site.Nickname" "Site.Name_OK" "Operator.ID" "Operator"              
                         #   "Latitude.DMS" "Longitude.DMS" "Latitude.DD" "Longitude.DD"       
                         #   "Ds_m"  "Zeu_calc_m" 
                       
# subset df
MASTERSHEETsub <- subset(MASTERSHEET, select = c(2:19))




