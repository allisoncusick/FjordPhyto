library(tidyverse)
library(ggmap)
library(rgdal)
library(mapproj)

load("data/WAP_clean.Rdata")

#Plot a basic Map
map <- map_data("world")

ggplot() + 
geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black", linewidth = 0.25) +
  coord_map(projection = "ortho",xlim = c(-75,-55), ylim = c(-72,-60),orientation = c(-100,-80,-12.5)) +
  geom_point(data = metadata, aes(x = Long, y = Lat), size = 4, pch = 4) 

##### MAKE A MAP USING LEAFLET ######
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
  addProviderTiles("Esri.OceanBasemap")
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


library(leaflet)
MASTERSHEET <- read.csv("data/FjordPhyto MASTER SHEET_pulled5april2023.csv")
#relabel latdd long dd 

MASTERSHEET <- MASTERSHEET %>% 
  rename("lat" = "Latitude.DD") 
MASTERSHEET <- MASTERSHEET %>% 
  rename("long" = "Longitude.DD")

colnames(MASTERSHEET)

# Create a color palette with handmade bins.
mybins <- seq(0, 40, by=10)
mypalette <- colorBin( palette="YlOrBr", domain=MASTERSHEET$Zeu_calc_m, na.color="transparent", bins=mybins)

mytext <- paste(
  "Operator: ", MASTERSHEET$Operator, "<br/>", 
  "Site: ", MASTERSHEET$Site.Name_OK, "<br/>", 
  "Zeu: ", MASTERSHEET$Zeu_calc_m, sep="") %>%
  lapply(htmltools::HTML)


# Final Map
m <- leaflet(MASTERSHEET) %>% 
  addTiles()  %>% 
  setView( lat=-65, lng=-65 , zoom=4) %>%
  addProviderTiles("Esri.WorldImagery") %>%
  addCircleMarkers(~long, ~lat, 
                   fillColor = ~mypalette(Zeu_calc_m), fillOpacity = 0.7, color="white", radius=8, stroke=FALSE,
                   label = mytext,
                   labelOptions = labelOptions( style = list("font-weight" = "normal", padding = "3px 8px"), textsize = "13px", direction = "auto")
  ) %>%
  addLegend( pal=mypalette, values=~Zeu_calc_m, opacity=0.9, title = "Euphotic Depth", position = "bottomright" )

m 


