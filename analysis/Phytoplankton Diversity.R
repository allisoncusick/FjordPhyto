#Phytoplankton Diversity

#load data
load("data/WAP_clean.Rdata")

# A - "Domain" B - "Kingdom"  C- "Phylum"   D- "Class"  E-  "Order"  F = "Family" G- "Genus" H - "Species" 

#example
taxa_table_split <- separate(taxa_table, pr2_Taxon, sep = ";", into = c("Domain","Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

#write files 
#write.csv(taxa_table_split, file = "data/taxa_table_split.csv")

#see how many groups, genera there are. 
#categorize groups: diatoms, dinoflagellates, autotrophs, mixotrophs (hamilton)
#classify following adl et al  https://onlinelibrary.wiley.com/doi/full/10.1111/jeu.12691 
unique(taxa_table_split$Phylum)
unique(taxa_table_split$Class)

diatoms <- taxa_table_split %>% filter(Class == "Bacillariophyta")
dinoflagellates <- taxa_table_split %>% filter(Phylum == "Dinoflagellata")
greenalgae <- taxa_table_split %>% filter(Phylum == "Chlorophyta")
haptophytes <- taxa_table_split %>% filter(Phylum == "Haptophyta")
rhodophytes <-taxa_table_split %>% filter(Phylum == "Rhodophyta")
cryptophytes <- taxa_table_split %>% filter(Phylum == "Cryptophyta") 
prasinophytes <- taxa_table_split %>% filter(Phylum == "Prasinodermophyta") 
MAST <-taxa_table_split %>% filter(Class == "MAST-\\.") #fix the all MAST-number


         