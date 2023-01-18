library(tidyverse)
library(ggmap)
library(rgdal)
library(mapproj)

load("data/WAP_clean.Rdata")

map <- map_data("world")

ggplot() + 
geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black", linewidth = 0.25) +
  coord_map(projection = "ortho",xlim = c(-80,-55), ylim = c(-72,-60),orientation = c(-100,-80,-12.5)) +
  geom_point(data = metadata, aes(x = Long, y = Lat), size = 4, pch = 4) 


