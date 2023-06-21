library(tidyverse)
library(ggmap)
library(rgdal)
library(mapproj)

load("data/WAP_clean.Rdata")

#make a map for greenland sampling sites
Greenland <-read.csv("data/FN Greenland 2023 map.csv")
# Create a color palette with handmade bins for zeu.
mybins <- seq(0, 35, by=5)
mypalette <- colorBin( palette="YlOrBr", domain=Greenland$Date, na.color="transparent", bins=mybins)

mytext <- paste(
  "Country ", Greenland$Country, "<br/>", 
  "Site: ", Greenland$Port, "<br/>", 
  "Date ", Greenland$Date, sep="") %>%
  lapply(htmltools::HTML)


# Final Map for greenland sites
m <- leaflet(Greenland) %>% 
  addTiles()  %>% 
  setView( lat=65, lng=-55 , zoom=3) %>%
  addProviderTiles("Esri.WorldImagery") %>%
  addCircleMarkers(~long, ~lat, fillOpacity = 0.7, color="blue", radius=8, stroke=FALSE, label=mytext)
                   
m 



#Plot a basic Map
map <- map_data("world")

ggplot() + 
geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black", linewidth = 0.25) +
  coord_map(projection = "ortho",xlim = c(-75,-55), ylim = c(-72,-60),orientation = c(-100,-80,-12.5)) +
  geom_point(data = metadata, aes(x = Long, y = Lat), size = 4, pch = 4) 

##### MAKE A MAP USING LEAFLET https://r-graph-gallery.com/19-map-leafletr.html  ######
library(leaflet)

# Note: if you do not already installed it, install it with:
# install.packages("leaflet")

# Example 1: NASA
m <- leaflet() %>% 
  addTiles() %>% 
  setView( lng = -65, lat = -65, zoom = 5 ) %>% 
  addProviderTiles("NASAGIBS.ViirsEarthAtNight2012") #here you can change the map types, see below
m

# Example 2: World Imagery
m <- leaflet() %>% 
  addTiles() %>% 
  setView( lng = -65, lat = -65, zoom = 5 ) %>% 
  addProviderTiles("Esri.WorldImagery")
m

#you can try these other tiles and find more here (https://github.com/leaflet-extras/leaflet-providers)
#Nasa: NASAGIBS.ViirsEarthAtNight2012
#Google map: Esri.WorldImagery
#Gray: Esri.WorldGrayCanvas
#Terrain: Esri.WorldTerrain
#Topo Map: Esri.WorldTopoMap
#OceanBase Map: Esri.OceanBasemap 

# save the widget in a html file if needed.
#library(htmlwidgets)
#saveWidget(m, file=paste0( getwd(), "/HtmlWidget/backgroundMapTile.html", width="1000px"))
#########

#now try with my data plotting Zeu from MASTERSHEET
#also try with CTD data below zeu 

library(leaflet)
MASTERSHEET <- read.csv("data/FjordPhyto MASTER SHEET_pulled5april2023.csv")

#relabel GPS Decimal Degrees  

MASTERSHEET$long <-MASTERSHEET$Longitude.DD
MASTERSHEET$lat <- MASTERSHEET$Latitude.DD


#MASTERSHEET <- MASTERSHEET %>% 
 # rename("lat" = "Latitude.DD") 
#MASTERSHEET <- MASTERSHEET %>% 
 # rename("lng" = "Longitude.DD")

colnames(MASTERSHEET)


#FIGURE OUT HOW TO REMOVE Zeu NAs ... 
#ADD THAT CODE HERE .... 

# Create a color palette with handmade bins for zeu.
mybins <- seq(0, 35, by=5)
mypalette <- colorBin( palette="YlOrBr", domain=MASTERSHEET$Zeu_calc_m, na.color="transparent", bins=mybins)

mytext <- paste(
  "Operator: ", MASTERSHEET$Operator, "<br/>", 
  "Site: ", MASTERSHEET$Site.Name_OK, "<br/>", 
  "Zeu: ", MASTERSHEET$Zeu_calc_m, sep="") %>%
  lapply(htmltools::HTML)


# Final Map for Zeu
m <- leaflet(MASTERSHEET) %>% 
  addTiles()  %>% 
  setView( lat=-65, lng=-65 , zoom=6) %>%
  addProviderTiles("Esri.WorldImagery") %>%
  addCircleMarkers(~long, ~lat, 
                   fillColor = ~mypalette(Zeu_calc_m), fillOpacity = 0.7, color="white", radius=8, stroke=FALSE,
                   label = mytext,
                   labelOptions = labelOptions( style = list("font-weight" = "normal", padding = "3px 8px"), textsize = "13px", direction = "auto")
  ) %>%
  addLegend( pal=mypalette, values=~Zeu_calc_m, opacity=0.9, title = "Euphotic Depth", position = "bottomright" )

m 

###### MAP FOR CTD SALINITY AND TEMPERATURE 
CTD <-read.csv("data/All_final_CTD.csv")
library(dplyr)


ctdtest<-CTD %>% 
  group_by(file_id, datetime_utc) %>%
  mutate(mean_salinity = mean(Salinity..Practical.Salinity.Scale.[Depth..Meter.[0-5] <= 5]))

ctdtest<-ctdtest %>% 
  group_by(file_id, datetime_utc) %>%
  mutate(mean_temperature = mean(Temperature..Celsius.[Depth..Meter.[0-5] <= 5]))

colnames(ctdtest)

#keep filename, mean temp, mean salinity, gps
ctdfilt<-ctdtest %>% distinct(file_id, mean_salinity, mean_temperature, latitude, longitude)

ctdfilt$long<-ctdfilt$longitude
ctdfilt$lat <- ctdfilt$latitude

colnames(ctdfilt)
na_ctdfilt <-na.omit(ctdfilt)
#i downloaded this and hand entered missing entries then reimported in the 'analysis_WAP' ....
#see analysis_WAP.R line importing 'CTDtomergewtaxa'

#alaina says another way to code the above is: 
#ctdfilt<-CTD %>%
  # group_by(file_id, datetime_utc) %>%
  # mutate(
   # mean_salinity = mean(Salinity..Practical.Salinity.Scale.[Depth..Meter.[0-5] <= 5]),
  #  mean_temperature = mean(Temperature..Celsius.[Depth..Meter.[0-5] <= 5]),
  #  mean_density = mean(Density..Kilograms.per.Cubic.Meter.[Depth..Meter.[0-5] <= 5]),
  #  long = longitude,
  #  lat = latitude) %>%
  # distinct(file_id, mean_salinity, mean_temperature, mean_density, latitude, longitude)

#colnames(ctdfilt)

#bring in a different file from analysis_WAP called 'CTD to merge w taxa ....
CTDtomergewtaxa$long<-CTDtomergewtaxa$Long_ctd
CTDtomergewtaxa$lat <- CTDtomergewtaxa$Lat_ctd
na_ctd <-na.omit(CTDtomergewtaxa)

# Create a color palette with handmade bins for CTD salinity values.
mybins <- seq(25, 34, by=0.5)
mypalette <- colorBin( palette="YlOrBr", domain=na_ctd$mean_salinity, na.color="transparent", bins=mybins)

mytext <- paste(
  "Site: ", na_ctd$datetime_utc, "<br/>",
  "meanTemperature: ", na_ctd$mean_temp, "<br/>",
  "meanSalinity: ", na_ctd$mean_sal, sep="") %>%
  lapply(htmltools::HTML)


# Final Map for ctd salinity

m <- leaflet(na_ctd) %>% 
  addTiles()  %>% 
  setView( lat=-65, lng=-65 , zoom=6) %>%
  addProviderTiles("Esri.WorldImagery") %>%
  addCircleMarkers(~long, ~lat, 
                   fillColor = ~mypalette(mean_sal), fillOpacity = 0.7, color="white", radius=8, stroke=FALSE,
                   label = mytext,
                   labelOptions = labelOptions( style = list("font-weight" = "normal", padding = "3px 8px"), textsize = "13px", direction = "auto")
  ) %>%
  addLegend( pal=mypalette, values=~mean_sal, opacity=0.9, title = "Mean Salinity", position = "bottomright" )

m 

# Create a color palette with handmade bins for CTD temperature.
range(ctdfilt$mean_temperature)

mybins <- seq(-2, 11, by=1)
mypalette <- colorBin( palette="YlOrBr", domain=ctdfilt$mean_temperature, na.color="transparent", bins=mybins)

mytext <- paste(
  "Site: ", ctdfilt$datetime_utc, "<br/>",
  "meanTemperature: ", ctdfilt$mean_temperature, "<br/>",
  "meanSalinity: ", ctdfilt$mean_salinity, sep="") %>%
  lapply(htmltools::HTML)

m <- leaflet(na_ctdfilt) %>% 
  addTiles()  %>% 
  setView( lat=-65, lng=-65 , zoom=6) %>%
  addProviderTiles("Esri.WorldImagery") %>%
  addCircleMarkers(~long, ~lat, 
                   fillColor = ~mypalette(mean_temperature), fillOpacity = 0.7, color="white", radius=8, stroke=FALSE,
                   label = mytext,
                   labelOptions = labelOptions( style = list("font-weight" = "normal", padding = "3px 8px"), textsize = "13px", direction = "auto")
  ) %>%
  addLegend( pal=mypalette, values=~mean_temperature, opacity=0.9, title = "Mean Temperature", position = "bottomright" )

m 
#the temp range doesn't depict the proper scale...like if change range to -2 to 4 the colors show way more in 4 C range. 
#amazing interactive, love it

#############################

#create a map showing Shannon Diversity (high number, higher diversity)
mapShannon<-read.csv("data/sandbox/working_df_2023JUN15.csv")
mapshannon_noNA <-mapShannon %>% filter(!is.na(shannon))

# Create a color palette with handmade bins for shannon.
mybins <- seq(0, 2, by=0.2)
mypalette <- colorBin( palette="YlOrBr", domain=mapshannon_noNA$shannon, na.color="transparent", bins=mybins)

mytext <- paste(
  "Site: ", mapshannon_noNA$Site_Name, "<br/>", 
  "Shannon: ", mapshannon_noNA$shannon, sep="") %>%
  lapply(htmltools::HTML)

# Final Map for Shannon
m <- leaflet(mapshannon_noNA) %>% 
  addTiles()  %>% 
  setView( lat=-65, lng=-65 , zoom=6) %>%
  addProviderTiles("Esri.WorldImagery") %>%
  addCircleMarkers(~Long, ~Lat, 
                   fillColor = ~mypalette(shannon), fillOpacity = 0.7, color="white", radius=8, stroke=FALSE,
                   label = mytext,
                   labelOptions = labelOptions( style = list("font-weight" = "normal", padding = "3px 8px"), textsize = "13px", direction = "auto")
  ) %>%
  addLegend( pal=mypalette, values=~shannon, opacity=0.9, title = "Shannon Diversity", position = "bottomright" )

m 
