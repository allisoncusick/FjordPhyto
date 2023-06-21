#process the manifest files and the returned meltwater dO18 data from Oregon State Labs

manifestGRN <-read.csv("data/meltwater dO18/d18O-bottle 1_52-manifest Greenland.csv")
manifestWAP1 <-read.csv("data/meltwater dO18/d18O-bottle-53_70_manifest-WAP FjordPhyto.csv")
manifestWAP2 <- read.csv("data/meltwater dO18/d18O-bottle-2022 JULY Round2_manifest-WAP FjordPhyto.csv")

#edited files from Oregon State Labs original
meltwaterGRNWAP1 <-read.csv("data/meltwater dO18/meltwater_1_70_220418_Cusick_edited.csv")
#filter above file to split into Greenland 1 - 52 and WAP1 bottles 53-70

meltwaterWAP2 <-read.csv("data/meltwater dO18/meltwater_2022 JULY Round 2_221117_Cusick_edited.csv")

#merge manifest with meltwater results 
#do some fun plots ... 
#use Jacks meltwater algorithm to compare Salinity to deltaO18 ...