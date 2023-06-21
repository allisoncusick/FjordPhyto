#Phytoplankton Diversity

#load in packages (do i need all these?)
library(tidyverse)
library(lubridate)
library(phyloseq)
library(ape)
library(vegan)
library(dplyr)
library(tidyr)
library(readr)


#load data from 'data_processing_WAP.R'
load("data/WAP_clean.Rdata")
#asv_table
#metadata
#taxa_table

#split taxa_table to individual categories
# A - "Domain" B - "Kingdom"  C- "Phylum"   D- "Class"  E-  "Order"  F = "Family" G- "Genus" H - "Species" 
taxa_table_split <- separate(taxa_table, pr2_Taxon, sep = ";", into = c("Domain","Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

#write files 
#write.csv(taxa_table_split, file = "data/taxa_table_split.csv")
#write.csv(merged_final_phyto_group_ids, file ="data/taxaphyto.csv")

#the point is to make asv so it can plot all the feature IDs as taxon
#I wrote similar script to make diatom plots for SCAR 2022 - use that as a guideline
#pivot asv_table so that each cell is filled with the read counts, each row is a sample, and each column is a feature ID
#do i need to tho if asv_table is in the format i want already? 
#pivot_asv <- asv_tae %>% pivot_wider(names_from = Feature.ID, values_from = reads, values_fill = 0) 
#change column name Feature.ID to hash 
asv_table$hash<-asv_table$Feature.ID 

#must load in metadata to attribute to each sample ID now
#noticed that the metadata, asv, and taxa table aren't named the same in the column for 'samples' 
#change 'metadata' file column for ''sample_id' to 'samples' 

metadata$samples <-metadata$sample_id #well that added a new column with that name ... 

#merge the data with metadata
merged_df <- merge(metadata, asv_table, by="samples")
#merge the metadata + asv, with taxa 
merged_final <-merge(merged_df, taxa_table_split, by="Feature.ID")

#what is the point of this code below here? 
#asv_table$Phylum <- taxa_table_split$Phylum [match(asv_table$hash, taxa_table_split$Feature.ID)]
#asv_table$Class <- taxa_table_split$Class[match(asv_table$hash, taxa_table_split$Feature.ID)]
#asv_table$Order <- taxa_table_split$Order[match(asv_table$hash, taxa_table_split$Feature.ID)]
merged_final$Phylum[is.na(merged_final$Phylum)] <- "Other Eukaryotes"

#see how many "phytoplankton"groups / genera to categorize 
unique(taxa_table_split$Phylum) #37
unique(merged_final$Phylum) #this is the same code for above, just using different data frames
unique(taxa_table_split$Class) #115

#in the taxa_table_split dataframe, categorize groups: diatoms, dinoflagellates, autotrophs, mixotrophs (hamilton)
#classify following adl et al  https://onlinelibrary.wiley.com/doi/full/10.1111/jeu.12691 
#diatoms <- taxa_table_split %>% filter(Class == "Bacillariophyta")
#dinoflagellates <- taxa_table_split %>% filter(Phylum == "Dinoflagellata")
#greenalgae <- taxa_table_split %>% filter(Phylum == "Chlorophyta")
#haptophytes <- taxa_table_split %>% filter(Phylum == "Haptophyta")
#rhodophytes <-taxa_table_split %>% filter(Phylum == "Rhodophyta")
#cryptophytes <- taxa_table_split %>% filter(Phylum == "Cryptophyta") 
#prasinophytes <- taxa_table_split %>% filter(Phylum == "Prasinodermophyta") 
#MAST <-taxa_table_split %>% filter(Class == "MAST") #fix the all MAST-numbers

#in the merged_final dataframe, categorize groups, this will make new dataframes per group
diatoms <- merged_final %>% filter(Class == "Bacillariophyta")
dinoflagellates <- merged_final %>% filter(Phylum == "Dinoflagellata")
greenalgae <- merged_final %>% filter(Phylum == "Chlorophyta")
haptophytes <- merged_final %>% filter(Phylum == "Haptophyta")
rhodophytes <-merged_final %>% filter(Phylum == "Rhodophyta")
cryptophytes <- merged_final %>% filter(Phylum == "Cryptophyta") 
prasinophytes <- merged_final %>% filter(Phylum == "Prasinodermophyta") 
MAST <- merged_final %>% filter(grepl("MAST-",Class))

#or, Alaina showed me that I can take the merged_final dataframe, and make phyto_groups column in one big dataframe

merged_final_phyto_group_ids <- merged_final %>% 
  mutate(
    phyto_groups = case_when(
      Class == "Bacillariophyta" ~ "diatoms",
      Phylum == "Dinoflagellata" & Class != "Syndiniales" ~ "dinoflagellates",
      Phylum == "Chlorophyta" ~ "green algae",
      Phylum == "Haptophyta" ~ "haptophytes",
      Phylum == "Rhodophyta" ~ "rhodophytes",
      Phylum == "Cryptophyta"~ "cryptophytes",
      Phylum == "Prasinodermophyta" ~"prasinophytes",
      grepl("MAST-", Class) ~ "MAST"))


unique(merged_final_phyto_group_ids$phyto_groups) #9 (with NA)
unique(merged_final_phyto_group_ids$Site_Name) #38 with NA
unique(merged_final_phyto_group_ids$sample_id) #271

#write this to sandbox merged final file in case i mess up later
write.csv(merged_final_phyto_group_ids, file = "data/sandbox/merged_final_phyto_group_ids.csv")
#if i get stuck, read in that file 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#now try to plot the big merged df based on groups

#assign colors to 9 

color_phyto_groups<-c("#372564",
                        "#80da3a",
                        "#6c38c4",
                        "#68d76e",
                        "#ca4ecd",
                        "#cad44c",
                        "#656ad5",
                        "#daa83e",
                        "#8b7bc4")
merged_final_phyto_group_ids %>%
  filter(!is.na(phyto_groups)) %>% 
  ggplot( aes(x = Date, y = reads, fill = phyto_groups)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = color_phyto_groups)
#how do i check if this is in date order? or look at each season at a time? 


#now try to plot merged_df  just for each group! see code from the old 2a_Analysis to look at Diatoms.R script
#colors from https://medialab.github.io/iwanthue/

###DIATOMS

unique(diatoms$Genus) #31
unique(diatoms$Species) #44

#FOR MYSELF DIATOM MASTER COLOR SHEET (do i want to make color codes based on "pennate vs centric" ?)
#something is strange, because i see navicula, rhizoselenia, biddulphia, and others in microscope that aren't in this list, but they were in my SCAR 2022 results
diatoms_genus <-as.factor(c("Actinocyclus", "Arcocellulus","Attheya","Amphora","Achnanthes",
                              "Bacteriastrum",
                            "Chaetoceros", "Coscinodiscus", "Contricribra", "Corethron",   
                            "Ditylum", "Eucampia",  "Hemialus", "Leptocylindrus",
                            "Melosira","Minutocellus", "Odontella", 
                                "Paralia", "Pleurosigma", "Polar-centric-Mediophyceae_X",
                               "Porosira", "Proboscia", "Pseudo-nitzschia", "Pleurosigma","Pseudogomphonema",
                              "Radial-centric-basal-Coscinodiscophyceae_X", "Raphid-pennate_X",
                               "Skeletonema", "Shionodiscus", "Stellarima", "Thalassiosira", "NA"))

diatoms_colors <-c("#372564",
                    "#80da3a",
                    "#6c38c4",
                    "#68d76e",
                    "#ca4ecd",
                    "#cad44c",
                    "#656ad5",
                    "#daa83e",
                    "#8b7bc4",
                    "#578d3d",
                    "#d83f8f",
                    "#6dd7a4",
                    "#db412e",
                    "#66d4d5",
                    "#c84960",
                    "#c5d191",
                    "#853576",
                    "#b9d1c4",
                    "#49212e",
                    "#6fa8d9",
                    "#cc713a",
                    "#505f7f",
                    "#8f7b3b",
                    "#d887c0",
                    "#405228",
                    "#cab0cd",
                    "#25342f",
                    "#d99e88",
                    "#53897d",
                    "#793522",
                    "#906767")
#without NAs (caution, date is NOT in correct order, need to order by year then month then day)
diatoms %>%
  filter(!is.na(Genus)) %>% 
  ggplot( aes(x = Date, y = reads, fill = Genus)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = diatoms_colors)


#now I want to plot above but sort x-axis based on site name THEN date (caution, dates not in order)
diatoms %>%
  filter(!is.na(Genus)) %>% 
ggplot(aes(x = Date, y = reads, fill = Genus)) +
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = diatoms_colors) #+ facet_wrap(.~Site_Name, scales="free") 
#if you dont want location, remove the + facet_wrap command

#### DINOFLAGELLATES 

unique(dinoflagellates$Order) #15
unique(dinoflagellates$Species) #78
unique(dinoflagellates$Class) #5 


#plot by class
dino_colors_class <-c("#372564",
                   "#80da3a",
                   "#6c38c4",
                   "#68d76e",
                   "#ca4ecd")
#this has NA in it
ggplot(dinoflagellates, aes(x = Date, y = reads, fill = Class)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = dino_colors_class)

#plot by order
dino_colors_order <-c("#773a30",
                      "#6fd054",
                      "#6a4bcc",
                      "#c7cc57",
                      "#ca4ac5",
                      "#87d1ae",
                      "#d34c3c",
                      "#789fc3",
                      "#c98b39",
                      "#523273",
                      "#57713d",
                      "#c44b7e",
                      "#3a3941",
                      "#bc8aca",
                      "#cc9f8e")

dinoflagellates %>%
  filter(!is.na(Order)) %>% 
  ggplot( aes(x = Date, y = reads, fill = Order)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = dino_colors_order)

#now I want to plot above but sort x-axis based on site name THEN date (caution, dates not in order)
dinoflagellates %>%
  filter(!is.na(Order)) %>% 
  ggplot(aes(x = Date, y = reads, fill = Order)) +
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = dino_colors_order) + facet_wrap(.~Site_Name, scales="free") 
#if you dont want location, remove the + facet_wrap command

###CRYPTOPHYTES
unique(cryptophytes$Order) #1
unique(cryptophytes$Species) #2
unique(cryptophytes$Class) #1 

#plot by species
crypto_colors_species <-c("#773a30",
                      "#6fd054")

cryptophytes %>%
  filter(!is.na(Species)) %>% 
  ggplot( aes(x = Site_Name, y = reads, fill = Species)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = crypto_colors_species)
#by site then date
cryptophytes %>%
  filter(!is.na(Species)) %>% 
  ggplot(aes(x = Date, y = reads, fill = Species)) +
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = crypto_colors_species) + facet_wrap(.~Site_Name, scales="free") 
#if you dont want location, remove the + facet_wrap command

#hmmm cant really see the hemiselmis, by site or date ....

###GREENALGAE
unique(greenalgae$Order) #11
unique(greenalgae$Species) #26
unique(greenalgae$Genus) #26 

green_colors_genus <-c("#6a3bc5",
                       "#7ad846",
                       "#c942cd",
                       "#d2cf3d",
                       "#6370d8",
                       "#63d98c",
                       "#d74887",
                       "#5c903e",
                       "#c573d3",
                       "#cbd986",
                       "#382971",
                       "#d37a39",
                       "#677ab5",
                       "#dc453b",
                       "#61c5b9",
                       "#85377a",
                       "#b39548",
                       "#452236",
                       "#b6caad",
                       "#923638",
                       "#91b5d4",
                       "#694b2e",
                       "#cb92c3",
                       "#3f5b37",
                       "#cb928c",
                       "#405564")

#plot by Genus
greenalgae %>%
  filter(!is.na(Genus)) %>% 
  ggplot( aes(x = Date, y = reads, fill = Genus)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = green_colors_genus)
#by site then date
greenalgae %>%
  filter(!is.na(Genus)) %>% 
  ggplot(aes(x = Date, y = reads, fill = Genus)) +
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = green_colors_genus) + facet_wrap(.~Site_Name, scales="free")

####PRASINO
unique(prasinophytes$Order) #2

prasi_colors_order <-c("#773a30",
                          "#6fd054")

prasinophytes %>%
  filter(!is.na(Order)) %>% 
  ggplot( aes(x = Date, y = reads, fill = Order)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = prasi_colors_order) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#by site then date
prasinophytes %>%
  filter(!is.na(Order)) %>% 
  ggplot(aes(x = Date, y = reads, fill = Order)) +
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = prasi_colors_order) + facet_wrap(.~Site_Name, scales="free")

###HAPTO
unique(haptophytes$Order) #4
unique(haptophytes$Species) #8
unique(haptophytes$Genus) #6

hapto_colors_species <-c("#6a3bc5",
                       "#7ad846",
                       "#c942cd",
                       "#d2cf3d",
                       "#6370d8",
                       "#63d98c",
                       "#d74887", "#b39548")
#by date
haptophytes %>%
  filter(!is.na(Species)) %>% 
  ggplot( aes(x = Date, y = reads, fill = Species)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = hapto_colors_species) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#by site 
haptophytes %>%
  filter(!is.na(Species)) %>% 
  ggplot( aes(x = Site_Name, y = reads, fill = Species)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = hapto_colors_species) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#by site then date
haptophytes %>%
  filter(!is.na(Species)) %>% 
  ggplot(aes(x = Date, y = reads, fill = Species)) +
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = hapto_colors_species) + facet_wrap(.~Site_Name, scales="free")


##~~~~~~~~~~~ Now bring in the cleaned CTD data that only has casts associated with samples
#use merged_final_phyto_group_ids to associate to df 'z' from the analysis_WAP.R file
#need to still fix date and NAs etc

unique(z$phyto_groups) #watch out for the NAs - how do i fix that? 

color_phyto_groups<-c("#372564",
                      "#80da3a",
                      "#6c38c4",
                      "#68d76e",
                      "#ca4ecd",
                      "#cad44c",
                      "#656ad5",
                      "#daa83e",
                      "#8b7bc4")
z %>%
  filter(!is.na(phyto_groups)) %>% 
  ggplot( aes(x = mean_sal, y = reads, fill = phyto_groups)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = color_phyto_groups)

#plot for mean temp  - these bar plots mean nothing
z %>%
  filter(!is.na(phyto_groups)) %>% 
  ggplot( aes(x = mean_temp, y = reads, fill = phyto_groups)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = color_phyto_groups)

#plot for mean density - these bar plots mean nothing
z %>%
  filter(!is.na(phyto_groups)) %>% 
  ggplot( aes(x = mean_dens, y = reads, fill = phyto_groups)) +
  geom_dotplot(binaxis = "", position = "fill") +
  scale_fill_manual(values = color_phyto_groups)

# Basic dot plot taxa over mean salinity
library(ggplot2)
z %>%
  filter(!is.na(phyto_groups)) %>% 
  ggplot(aes(x = mean_sal, y = reads, colour = phyto_groups)) + 
  geom_point ()

# Basic dot plot taxa over mean temperature
library(ggplot2)
z %>%
  filter(!is.na(phyto_groups)) %>% 
  ggplot(aes(x = mean_temp, y = reads, colour = phyto_groups)) + 
  geom_point ()

#crap, thats showing reads, not a fraction of the total. ....

#make z have a column with fraction of total ... 
library(ggplot2)
z_prop %>%
  ggplot(aes(x = mean_sal, y = prop_reads, colour = phyto_groups)) + 
  geom_point () + theme_classic()

ctdfinal %>%
  ggplot(aes(x = median_salinity, y = prop_reads, colour = phyto_groups)) + 
  geom_point () + theme_classic()

# Basic dot plot taxa over mean temperature
library(ggplot2)
z %>%
  filter(!is.na(phyto_groups)) %>% 
  ggplot(aes(x = mean_temp, y = reads, colour = phyto_groups)) + 
  geom_point ()
f

##################################################
##### DIVERSITY ANALYSIS (TAMMY AND NCOG Chase AND BEN C)
#load in the working df and manipulate 
which(colnames(working_df)=="samples")
df_diversity <-working_df[,c(2,27,38)]
summary(df_diversity$phyto_groups)
df_diversity$phyto_groups <-as.factor(df_diversity$phyto_groups)

colnames(sumdf_wider)

# Group by sample and then by phytoplankton group and sum prop_reads
summary_df <-working_df %>% 
  group_by(samples, phyto_groups) %>% 
  summarize(counts = sum(prop_reads))
#if i group by sample, then sum the counts - it should give me the totalreads column
#i'm in long form, so to do diversity analysis i have to do wide form, use pivot

sumdf_wider<- summary_df %>% pivot_wider(names_from = phyto_groups, values_from = counts, values_fill = 0)
#had to add a values_fill 0 cuz prasino was being funky making NAs

#look at it

#use this data to do the shannon diversity analysis on phyto_groups
matrix_data<-as.matrix(sumdf_wider[,-1])
# If doesn't work, remove sample id column, then you can add it back later, thats what the [,-1] means

rownames(matrix_data)<-sumdf_wider$samples
matrix_data_BC<-metaMDS(matrix_data, distance = "bray", na.rm=FALSE)
orditorp(matrix_data_BC, "sites")

#need to get rid of NAs or it will mess things up. prasinophytes have NAs ...
temp2 <- replace(temp2, is.na(temp2), 0)
#try chases code Shannon Diversity
temp$shannon_index <-vegan::diversity(temp, MARGIN = 1, index = "shannon")
#above will not work cuz column 'samples' is not numeric so remove 
temp2$shannon_index <-vegan::diversity(temp2, MARGIN = 1, index = "shannon")



#try chases code evenness
S <- apply(temp2>0,1,sum)
temp2$evenness <-vegan::diversity(temp2, index = "shannon")/log(S)

#richness
richness <- apply(temp2, 1, function(x) length(which(x != 0)))
temp2$richness <- richness
#why are there numbers in richness but not evenness and shannon
#help can I make the NAs say 0 instead

#inv_simp

temp2$inv_simp <- vegan::diversity(temp2, MARGIN = 1, index = "invsimpson")



# chao1
#chao1 <- estimateR(temp)

#asv_copy$chao1 <- chao1[2,]

# gini

#gini <- apply(temp2, 1, gini)

#temp2$gini <- gini

# Dissimilarity
#below doesn't work ? open that file, its gibberish

dissimilar <- vegdist(temp2, method = "bray", binary = FALSE)
dissimilar <- as.matrix(dissimilar)

save(dissimilar, file = 'dissimilar_matrix') #ummm thats a garbage file in my fjordphyto folder


#add the samples column back in to temp2
temp3<-bind_cols(temp[,1], temp2[,2:12])
#yahoo now sample, phyto_groups, indexes are together
#merge into the working_df
working_df$shannon <-temp3$shannon_index[match(working_df$samples, temp3$samples)]
working_df$evenness <-temp3$evenness[match(working_df$samples, temp3$samples)]
working_df$richness <-temp3$richness[match(working_df$samples, temp3$samples)]
working_df$inv_simp <-temp3$inv_simp[match(working_df$samples, temp3$samples)]


#save the working df for chase as a r.file then send 
write.csv(working_df, file = "data/sandbox/working_df_2023JUN15.csv") #for chase
write.csv(working_df, file = "data/sandbox/working_df_2023JUN19.csv") #updated with other indexes
#is it OK that the prasinophytes were written to be 0?

#sumdf_wider is phyto_group level I want to look at all species, so need to redo this with that level of column naming
#how? 
#I would need to ask it to pull and make columns the following, down to the lowest taxonomic resolution, if possible to species: 
#Class == "Bacillariophyta" ~ "diatoms",
#Phylum == "Dinoflagellata" & Class != "Syndiniales" ~ "dinoflagellates",
#Phylum == "Chlorophyta" ~ "green algae",
#Phylum == "Haptophyta" ~ "haptophytes",
#Phylum == "Rhodophyta" ~ "rhodophytes",
#Phylum == "Cryptophyta"~ "cryptophytes",
#Phylum == "Prasinodermophyta" ~"prasinophytes",
#grepl("MAST-", Class) ~ "MAST"

#but if there were NAs down to a lower classification level, I would need to tell it to take the last named level and ignor the NA
#column names: samples, all diatoms, all dino, all green algae, all x y z 
#cells filled with prop_reads of those respective to the total reads per sample (need to recalculate that?)
#so filter dataframe, and only keep Class Bacillariophyto, Phylum Dinoflagellata and Class Syndiniales, and Phylum Chlorophyta, and etc

newdf<-working_df %>% 
  group_by(samples) %>%
