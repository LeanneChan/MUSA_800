##### standarlize #####
setwd("C:/Users/HP/Desktop/fire")

library(tidyverse)
library(sf)
library(QuantPsyc)
library(RSocrata)
library(viridis)
library(caret)
library(spatstat)
library(spdep)
library(FNN)
library(grid)
library(gridExtra)
library(knitr)
library(kableExtra)
library(tidycensus)

# create a mapTheme() function that is used to standardized the map outputs 
mapTheme <- function(base_size = 12) {
  theme(
    text = element_text( color = "black"),
    plot.title = element_text(size = 14,colour = "black"),
    plot.subtitle=element_text(face="italic"),
    plot.caption=element_text(hjust=0),
    axis.ticks = element_blank(),
    panel.background = element_blank(),axis.title = element_blank(),
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=2)
  )
}

##### create fishnet #####
# read the boundary data
philly <- 
  st_read("http://data.phl.opendata.arcgis.com/datasets/405ec3da942d4e20869d4e1449a2be48_0.geojson")%>%
  st_transform(crs=2272)
# create the 500ft grid cell fishnet
fishnet <- 
  st_make_grid(philly, cellsize = 500) %>%
  st_sf()
# clip the fishnet by the boundary
fishnet <- 
  fishnet[philly,] %>%
  mutate(uniqueID = rownames(.)) %>%
  dplyr::select(uniqueID)
# plot the fishnet
ggplot() +
  geom_sf(data=philly,fill=NA) +
  geom_sf(data=fishnet,fill=NA) +
  labs(title = "Fishnet in Philly") +
  mapTheme()

##### load DV data#####
# load the fire data
fire <- read.csv("2015 to 2019 fires for u of pa report run 1620.csv")%>%
        subset(.,longitude>0|latitude>0)%>%
        st_as_sf(.,coords = c("longitude", "latitude"), crs = 4326, agr = "constant") %>%
        st_transform(st_crs(fishnet))

firecount <- st_intersection(fire, philly)

ggplot() + 
  geom_sf(data =firecount, colour="#006699", size=0.7) +
  geom_sf(data = philly,fill=NA) +
  labs(title= "firecount , Philly") +
  mapTheme()

##### load Predictors data#####
