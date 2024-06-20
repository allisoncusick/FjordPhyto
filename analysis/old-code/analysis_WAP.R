library(patchwork)
library(lubridate)
library(ggplot2)
#Read in clean 18S taxa asv metadata files from 'data_processing_WAP.R' 
load("data/WAP_clean.Rdata")
#asv_table
#metadata
#taxa_table
#THE HOLY GRAIL
#NEVER HAVE TO PROCESS DATA AGAIN!!!!

#####cHASE:  WHY DOES METADATA HAVE SO MANY ROWS OF NA?!?!?! some weird thing with importing metdada and merging in data_processing_WAP.R 
#this got fixed with my new cleaned set in data_processing_WAP.R line metadata<-read.csv("data/metadata/metadatafix_fixed.csv")

#see the Phytoplankton Diversity.R script cuz i morphed the overall final merged df into a monster

#Bring in other data sets
#noticed that the metadata, ASV, and taxa table aren't aligned in the column for 'samples' 
#change 'asv_table' column for ''samples' to 'sample_id'
#this code above no longer needed cuz i went to the 'data_processing_WAP.R' script and uploaded a new metadata file that had that added
#colnames(asv_table_new)<-c('Feature.ID', 'sample_id', 'reads')
#colnames(asv_table_new)
#unique(metadata$sample_id)

#1. READ IN CTD FILE ~~~~~~~~~~~~~~~~~~~

CTD <-read.csv("data/All_final_CTD.csv")
library(dplyr)



#remove columns that i dont want  from the original CTD file (keep depth, time, lat/long, salinity, temp, density, filename)
#I DONT ACTUALLY THINK THATS A GOOD IDEA, What if i wanna  calculate salinity from conductivity later? or look at pressure ? 
CTDsubset = subset(CTD, select = c(Depth..Meter., Temperature..Celsius., Salinity..Practical.Salinity.Scale., datetime_utc, latitude, longitude, Density..Kilograms.per.Cubic.Meter., file_id, cruise_season))

###chase code summarise vs mutate, (like the 'distinct' file), easier/faster
ctdmeans<-CTDsubset %>% 
  group_by(file_id, datetime_utc) %>%
  summarise(mean_salinity = mean(Salinity..Practical.Salinity.Scale.[Depth..Meter.[0-5] <= 5]),
            mean_temperature = mean(Temperature..Celsius.[Depth..Meter.[0-5] <= 5]),
            mean_density = mean(Density..Kilograms.per.Cubic.Meter.[Depth..Meter.[0-5] <= 5]))
########

#filter by file id, date/time, find surface 0-5 meter depth mean and median salinity
ctdmeansal<-CTDsubset %>% 
  group_by(file_id, datetime_utc) %>%
  mutate(mean_salinity = mean(Salinity..Practical.Salinity.Scale.[Depth..Meter.[0-5] <= 5]))

ctdmediansal <-ctdmeansal %>%
  mutate (median_salinity = median(Salinity..Practical.Salinity.Scale.[Depth..Meter.[0-5]<=5]))

#filter by file id, date/time, find surface 0-5 meter depth mean and median temperature

ctdsaltemp<-ctdmediansal %>% 
  group_by(file_id, datetime_utc) %>%
  mutate(mean_temperature = mean(Temperature..Celsius.[Depth..Meter.[0-5] <= 5]))

ctdmediansaltemp<-ctdsaltemp %>%
  mutate(median_temperature = median(Temperature..Celsius.[Depth..Meter.[0-5]<=5]))

colnames(ctdmediansaltemp) 

#filter by file id, date/time, find surface 0-5 meter depth mean and median density

ctdprefinal<-ctdmediansaltemp %>% 
  group_by(file_id, datetime_utc) %>%
  mutate(mean_density = mean(Density..Kilograms.per.Cubic.Meter.[Depth..Meter.[0-5] <= 5]))

ctdfinal<-ctdprefinal %>%
  mutate(median_density = median(Density..Kilograms.per.Cubic.Meter.[Depth..Meter.[0-5]<=5]))

colnames(ctdfinal) #ctdfinal shows full profile cast info w mean and median sal, temp, density over depth

#write ctd final into sandbox to play around with stitching in Tableau

write.csv(ctdfinal, file = "data/sandbox/ctdfull.csv")
#THIS GETS SQUASHED TO ONLY SHOW ONE 'surface' value (so no profiles depth info, just surface)
#that file gets stitched to dna data


#do i wanna do this for calculate BV? MLD? 
#load in OCE or TEOS package https://www.teos-10.org/software.htm
#https://teos-10.github.io/GSW-R/articles/gsw.html
install.packages('gsw')
library(gsw)
#Brunt-Vaisala Frequency [cycl/h] and potential temp  [degC]
#i went into ODV to derive this variable and then tried to export, not easy to export
#the BV in ODV is negative, weird. 
#Try in the GSW package, calculate psu from conductivity and compare to CTDs calculated psu

SA <- gsw_SP_from_C(C = CTD$Conductivity..MicroSiemens.per.Centimeter., t = CTD$Temperature..Celsius., p =CTD$Pressure..Decibar.)
plot(SA, CTD$Salinity..Practical.Salinity.Scale.)
#why is there 2 lines

#split the datetime_utc into separate columns
ctdfinal_datesplit<- ctdfinal %>% separate(col=datetime_utc, into = c("month", "day", "year"), sep = "/")
ctdfinal <-ctdfinal_datesplit %>% separate(col=year, into = c("year", "time_utc"), sep = " ")
#ctdfix$asDate %>% mutate(Date = as.Date(Date, format= "%dd/%mm/%YYYY")) 


#####CTDFINAL is the full CTD dataset with columns of interest (omitted conductivity, pressure,etc). 
#####NOW only choose individual files (no profile depth info) so  Ican do the next calculations

#show only unique row entries for surface, so no depth profile info, just flattened 2D info filename, mean and median temp salinity density, gps, so no profile depth info this is 2D info
colnames(ctdfinal)

#i broke the command below, try to fix it somehow. 
ctdsurface <- ctdfinal %>% distinct(file_id, month, day, year, time_utc, latitude, longitude, cruise_season, mean_salinity, median_salinity, mean_temperature, median_temperature, mean_density, median_density) #do i need to include other columns? 
#shows 283 "casts" so some were invalid then, cuz originally there were 273. 
#ctddistinct shows df with one entry surface means salinity and temperature and density per cast file
#so write this ctdsurface df out and merge with the metadata/DNA sample id i made in excel
write.csv(ctdsurface, file = "data/sandbox/ctdsurface.csv")
#I opened the ctdsurface file and I copy/pasted the median temp, sal, dens into the "ctd_import_mergetaxa csv and renamed it ctd_import_mergetaxa_2


##I haven't yet omitted any NAs 
##i need to filter out GPS that are in south america/falklands? 
#FIND THE GPS RANGE FOR WAP ONLY see emilys code FHL
#or go in by hand to excel and create working datasheets (UGH)
#arg the time is UTC and i want it to -3 for GMT-3
#create a column for time_local

##~~~~~~~try plotting stuff

#plot stuff with the surface depth information using ctddistinct
plot(x = ctdfinal$mean_salinity, y = ctdfinal$mean_temperature)
plot(x = ctdfinal$median_salinity, y = ctdfinal$median_temperature)
plot(x = ctdfinal$month, y = ctdfinal$median_salinity)
plot(factor(ctddistinct$month), ctddistinct$mean_salinity, xlab="month", ylab="mean salinity") 
plot(factor(ctddistinct$month), ctddistinct$mean_temperature, xlab="month", ylab="mean temperature") 
plot(factor(ctddistinct$year), ctddistinct$mean_salinity, xlab="year", ylab="mean salinity") 
plot(factor(ctddistinct$year), ctddistinct$mean_temperature, xlab="year", ylab="mean temperature") 
plot(ctdfinal$Density..Kilograms.per.Cubic.Meter.[1:68]-ctdfinal$Density..Kilograms.per.Cubic.Meter.[2:69])
plot(ctdfinal$Temperature..Celsius.[1:68]-ctdfinal$Temperature..Celsius.[2:69])
plot(ctdfinal$Salinity..Practical.Salinity.Scale.[1:68]-ctdfinal$Salinity..Practical.Salinity.Scale.[2:69], ylim = c(-.1, .1))
plot(ctdfinal$Temperature..Celsius.[1:68]-ctdfinal$Temperature..Celsius.[2:69], ylim = c(-.1, .1))
plot(ctdfinal$Density..Kilograms.per.Cubic.Meter.[1:68]-ctdfinal$Density..Kilograms.per.Cubic.Meter.[2:69], ylim = c(-.1, .1))
plot(ctdfinal$Temperature..Celsius.[1:68]-ctdfinal$Temperature..Celsius.[2:69], ctdfinal$Salinity..Practical.Salinity.Scale.[1:68]-ctdfinal$Salinity..Practical.Salinity.Scale.[2:69], xlim = c(-.1, .1), ylim = c(-.1, .1))


#plot some full CTD cast profiles vs depth using ctdfinal
plot(x=ctdfinal$Density..Kilograms.per.Cubic.Meter., y=ctdfinal$Depth..Meter., ylim=rev(range(ctdfinal$Depth..Meter.)))
plot(x=ctdfinal$Temperature..Celsius., y=ctdfinal$Depth..Meter., ylim=rev(range(ctdfinal$Depth..Meter.)))
plot(x=ctdfinal$Salinity..Practical.Salinity.Scale., y=ctdfinal$Depth..Meter., ylim=rev(range(ctdfinal$Depth..Meter.)))
plot(x = ctdfinal$mean_salinity, y = ctdfinal$mean_temperature)
plot(x = ctdfinal$month, y = ctdfinal$mean_salinity)
plot(factor(ctdfinal$month), ctdfinal$mean_salinity, xlab="month", ylab="mean salinity") 
plot(factor(ctdfinal$month), ctdfinal$mean_temperature, xlab="month", ylab="mean temperature") 
plot(factor(ctdfinal$cruise_season), ctdfinal$mean_salinity, xlab="season", ylab="mean salinity") 
plot(factor(ctdfinal$cruise_season), ctdfinal$mean_temperature, xlab="season", ylab="mean temperature") 

plot(x=ctdfinal$Temperature..Celsius., y=ctdfinal$Depth..Meter., ylim=rev(range(ctdfinal$Depth..Meter.)))
plot(x=ctdfinal$Salinity..Practical.Salinity.Scale., y=ctdfinal$Depth..Meter., ylim=rev(range(ctdfinal$Depth..Meter.)))
plot(x = ctdfinal$mean_salinity, y = ctdfinal$mean_temperature)
plot(x = ctdfinal$month, y = ctdfinal$mean_salinity)
plot(factor(ctdfinal$month), ctdfinal$mean_salinity, xlab="month", ylab="mean salinity") 
plot(factor(ctdfinal$month), ctdfinal$mean_temperature, xlab="month", ylab="mean temperature") 
plot(factor(ctdfinal$cruise_season), ctdfinal$mean_salinity, xlab="season", ylab="mean salinity") 
plot(factor(ctdfinal$cruise_season), ctdfinal$mean_temperature, xlab="season", ylab="mean temperature") 



#~~~~~~~~~~~~~~~
#write files for the distinct subset and work in excel
#write.csv(ctddistinct, file = "data/ctddistinct.csv")
#colnames (ctddistinct)
#manually open that csv file in excel and associate file id with taxa manifest sample id
###PROBLEMS: lots of NAs in GPS, need to manually fix. then some GPS are NA and dates dont have any matches in MASTERSHEET
#what to do with those? keep? delete? useless...for this goal, I removed casts that were not associated with any sample id
#this also removed the falkland island/drake GPS locations.



#re-import the completed hand manipulated corrected ctd 'master' file
CTDtomergewtaxa<-read.csv("data/ctd_import_mergetaxa.csv")
#to check, throw this file below into map.R - seems ok? 
#added Zeu in metadata file and add the BV and MLD calcs here above? and jacks meltwater? 

#now how to attach this to the taxa and metadata info? 
#re-import the completed hand manipulated corrected ctd 'master' file
CTDtomergewtaxa2<-read.csv("data/ctd_import_mergetaxa_2.csv")
#Associate each ctd file name with a sample/location from metadata (is it better to make each CTD cast a separate file, not in one giant one?)
#add a column thats the linker - the unique.id.code from MASTERSHEET ? do this by hand for 262 casts associate to 'science boats'? 
#for now just link wiht the sampleid. any sample that does not have a ctd cast, is tossed out in the CTDtomergewtaxa file


#DO THE MERGE CTD WITH DNA 
x<-data.frame(CTDtomergewtaxa) 
y<-data.frame(merged_final_phyto_group_ids)
z<-merge(x,y,by="samples")
z
#z is massive .... is this right?! 
#to be safe make a 'sandbox' folder with the 'final files' i am playing with, ease of loading later
#write files 
write.csv(z, file = "data/sandbox/CTD_and_taxa_Merged.csv")
#load this file in if i mess up later .... 

#repeat above code but with new ctd_import_mergetaxa_2 that has median info stitched in
x2<-data.frame(CTDtomergewtaxa2) 
y2<-data.frame(merged_final_phyto_group_ids)
z2<-merge(x2,y2,by="samples")
z2
#z is massive .... is this right?! 
#to be safe make a 'sandbox' folder with the 'final files' i am playing with, ease of loading later
#write files 
write.csv(z2, file = "data/sandbox/CTD_and_taxa_Merged_2.csv")
#ok but this z2 file doesn't have prop_reads yet, so gotta go add that ... and then later add jacks meltwater
#jump down to code for z_prop_2



#plot some stuff for CTD distinct (those only associated with taxa samples) 
plot(x=z$mean_sal, y=z$phyto_groups)
plot(x=ctdfinal$Temperature..Celsius., y=ctdfinal$Depth..Meter., ylim=rev(range(ctdfinal$Depth..Meter.)))
plot(x=ctdfinal$Salinity..Practical.Salinity.Scale., y=ctdfinal$Depth..Meter., ylim=rev(range(ctdfinal$Depth..Meter.)))
plot(x = ctdfinal$mean_salinity, y = ctdfinal$mean_temperature)
plot(x = ctdfinal$month, y = ctdfinal$mean_salinity)
plot(factor(ctdfinal$month), ctdfinal$mean_salinity, xlab="month", ylab="mean salinity") 
plot(factor(ctdfinal$month), ctdfinal$mean_temperature, xlab="month", ylab="mean temperature") 
plot(factor(ctdfinal$cruise_season), ctdfinal$mean_salinity, xlab="season", ylab="mean salinity") 
plot(factor(ctdfinal$cruise_season), ctdfinal$mean_temperature, xlab="season", ylab="mean temperature") 



#take the z dataframe, remove phyto_groups with NA, look at each sample and sum total reads, then find proportion of phytogroups reads per total
#word of note: go through phyto-groups and make sure I assigned things correctly
z_prop_2<- z2 %>%
  filter(!is.na(phyto_groups))%>%
  group_by(samples) %>%
  mutate(totalreads = sum(reads), prop_reads = reads/totalreads) 

# make seasons

z_prop_2$Date <- mdy(z_prop_2$Date)

#make one big label for all then filter out - by seasons 
z_prop_2$season <- "2019-2020"

z_prop_2$season[z_prop_2$Date >= mdy("05/01/2018") & z_prop_2$Date <= mdy("05/01/2019")] <- "2018-2019"
z_prop_2$season[z_prop_2$Date < mdy("05/01/2018")] <- "2017-2018"

#check if all unique no NAs etc 
unique(z_prop_2$median_salinity)
unique(z_prop_2$median_temperature)


#unique(ctdfinal$file_id) #283 file ids
#unique(totalreadspersample$file_id_ctd) #166 file ids


#write files for z_prop_2
write.csv(z_prop_2, file = "data/sandbox/CTD_and_taxa_Merged_z_prop_2.csv")

#now group so I can plot some things
grouped <- z_prop_2 %>% 
           group_by(samples,phyto_groups,season) %>%
         summarise(prop = sum(prop_reads),
          median_temperature = median(median_temperature),
         median_salinity = median(median_salinity))

a <- ggplot(grouped %>% filter(phyto_groups =="diatoms"), aes(x = median_temperature, y = median_salinity, color = prop)) +
  geom_point() + theme_classic() + labs(x = "Temp_median", y = "Sal_median", title = "Diatoms") 
a

b <- ggplot(grouped %>% filter(phyto_groups =="cryptophytes"), aes(x = median_temperature, y = median_salinity, color = prop)) +
  geom_point() + theme_classic() + labs(x = "Temp_median", y = "Sal_median", title = "Cryptophytes")

layout <- "A
           B"

a + b + plot_layout(design = layout)
#take this dataframe and go plot stuff in Phytoplankton Diversity R script but first need to add jacks meltwater

# Basic box plot
p <- ggplot(z2, aes(x=median_temperature, y=phyto_groups)) + 
  geom_boxplot()
p
#ummm that looks like there a very narrow range of temp ... is this right? 

#bring in jacks meltwater
#add a column that calculates jacks meltwater based on the equation 
#Meltwater Fraction (%) = (-0.021406 * df$median_salinity + 0.740392) * 100
#make sure you're using the right file, should have been created above and labeld as
#CTDtaxawmeltwater <-read.csv("data/sandbox/CTD_and_taxa_Merged_z_prop_2.csv")
z_prop_2_mw <- z_prop_2  %>%
  mutate(mw_fraction_perc = (-0.021406 * median_salinity + 0.740392) * 100)
#check it worked  
colnames(z_prop_2_mw) 
unique(z_prop_2_mw$mw_fraction_perc)
#YAY it worked .... 
working_df <- z_prop_2_mw
mw <- ggplot(working_df, aes(x=mw_fraction_perc, y=phyto_groups)) + 
  geom_boxplot()
mw

#MY FINAL DATAFRAME JUNE 1 IS WORKING_DF
#write files out for later use or to share
write.csv(working_df, file = "data/sandbox/working_df_2023JUN01.csv")
#I made new dataframes June 15 and 19 with diversity index columns
#see sandbox/working_df_2023JUN19 

# ^ I need to add MLD and BV (N and N2) here ..... 
#https://cerweb.ifremer.fr/deboyer/mld/Surface_Mixed_Layer_Depth.php
#read in the June1 file
working_df <-read.csv("data/sandbox/working_df_2023JUN01.csv")


#run wilcox test
dino_sal <- working_df %>%
  filter(!is.na(phyto_groups))%>%
  filter(phyto_groups =="dinoflagellates")

d <-dino_sal[,"median_salinity"]

prasino_sal <-working_df %>%
  filter(!is.na(phyto_groups))%>%
  filter(phyto_groups =="prasinophytes")

p <-prasino_sal[,"median_salinity"]
w <- wilcox.test(d, p, alternative = "two.sided")

w 

#playing with some plotting, nothing fancy
dinosal_mw<- working_df %>% filter(!is.na(phyto_groups))%>%
  filter(phyto_groups =="dinoflagellates")

plot(x=dinosal_mw$median_salinity, y=dinosal_mw$prop_reads)


#~~~~~~~~~~~~
#DO NOT HAVE TO DO THE BELOW CUZ I DID IT BY HAND IN EXCEL 
#2. READ IN MASTERSHEET for Zeu (this is downloaded from https://docs.google.com/spreadsheets/d/1UfYEAQ_BIWsEqeNQmZGTfzNZz1alBYffue3LVZ9NBeE/edit#gid=0
#MASTERSHEET <- read.csv("data/FjordPhyto MASTER SHEET_pulled5april2023.csv")
#see what columns are in MASTERSHEET
#colnames(MASTERSHEET)
#keep only columns I want: "Season.Year" "Season.Phase" "YYYY" "MM" "DD"  "DD_MON_YY"  "TIME.LOCAL"                     
                         #  "LOC.ID"  "Site.Nickname" "Site.Name_OK" "Operator.ID" "Operator"              
                         #   "Latitude.DMS" "Longitude.DMS" "Latitude.DD" "Longitude.DD"       
                         #   "Ds_m"  "Zeu_calc_m" 
                       
# subset df
#MASTERSHEETsub <- subset(MASTERSHEET, select = c(2:19))
#Look at the map.R script to make cool interactive bubble plot map and see Zeu with GPS
#but now my CTD file to merge with taxa has the Zeu added there instead .... 
#write a script that cleans up the MASTERSHEET when pulled from offline

