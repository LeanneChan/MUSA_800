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
library(lubridate)
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

##### load the boundary data#####
philly <- 
  st_read("http://data.phl.opendata.arcgis.com/datasets/405ec3da942d4e20869d4e1449a2be48_0.geojson")%>%
  st_transform(crs=2272)

# load the neighborhood data
neighbor <- 
  st_read("data/Neighborhoods_Philadelphia/Neighborhoods_Philadelphia.shp")%>%
  st_transform(crs=2272)
# plot the neighbourhood
ggplot() +
  geom_sf(data=philly,fill=NA) +
  geom_sf(data=neighbor,fill=NA) +
  labs(title = "neighborhoods in Philly") +
  mapTheme()

# load the census tracts data
tracts <- 
  st_read("http://data.phl.opendata.arcgis.com/datasets/ccdc0d022f4949a7b11333bd37231aef_0.geojson")%>%
  st_transform(crs=2272)
# plot the census tracts
ggplot() +
  geom_sf(data=philly,fill=NA) +
  geom_sf(data=tracts,fill=NA) +
  labs(title = "tracts in Philly") +
  mapTheme()

##### load original data#####
# load the fire data
fire <- read.csv("origin_data/2015 to 2019 fires for u of pa report run 1620.csv")%>%
  subset(.,longitude>0|latitude>0)%>%
  st_as_sf(.,coords = c("longitude", "latitude"), crs = 4326, agr = "constant") %>%
  st_transform(st_crs(fishnet))%>%
  st_intersection(., philly)

ggplot() + 
  geom_sf(data =fire, colour="#006699", size=0.7) +
  geom_sf(data = philly,fill=NA) +
  labs(title= "firecount , Philly") +
  mapTheme()

# load the inspection data
inspect <- 
  st_read("origin_data/PDF_Hydrant_Inspections_2019_for_MUSA/PDF_Hydrant_Inspections_2019.shp") %>%
  st_transform(crs=2272) 

ggplot() + 
  geom_sf(data =inspect, colour="#006699", size=0.7) +
  geom_sf(data = philly,fill=NA) +
  labs(title= "fireinspect , Philly") +
  mapTheme()

# load the hydrant data
hydrants <- 
  st_read("origin_data/Hydrants_2019/Hydrants_2019.shp") %>%
  st_transform(crs=2272) 

ggplot() + 
  geom_sf(data =hydrants, colour="#006699", size=0.7) +
  geom_sf(data = philly,fill=NA) +
  labs(title= "hydrant , Philly") +
  mapTheme()

# load the engine locals data
engine_new <- 
  st_read("origin_data/EngineLocals/EngineLocals_New.shp") %>%
  st_transform(crs=2272) 
ggplot() + 
  geom_sf(data =engine_new, colour="#006699", size=0.7) +
  geom_sf(data = philly,fill=NA) +
  labs(title= "engine_new , Philly") +
  mapTheme()

engine_old <- 
  st_read("origin_data/EngineLocals/EngineLocals_Old.shp") %>%
  st_transform(crs=2272) 
ggplot() + 
  geom_sf(data =engine_old, colour="#006699", size=0.7) +
  geom_sf(data = philly,fill=NA) +
  labs(title= "engine_old , Philly") +
  mapTheme()



##### explore fire data by time (to be solved) #####
fire <- 
  fire %>%
  mutate(intervaday=floor_date(ymd_hms(alm_date),unit='day'),
         year=year(interval60),
         month=month(interval60),
         dotw=wday(interval60,label=TRUE),
         weekend = ifelse(dotw %in% c("Sun", "Sat"), "Weekend", "Weekday")) 

##### explore leading cause of fires across philly #####
reorder_size <- function(x) {
  factor(x, levels = names(sort(table(x),decreasing=T)))
}
ggplot(fire, aes(reorder_size(inci_type))) + 
  geom_bar() +  
  theme(axis.text.x = element_text(angle = 60))

#visualize type 113
fire113 <- 
  fire %>%
  filter(., inci_type == 113)
ggplot() + 
  geom_sf(data =fire113, colour="#006699", size=0.7) +
  geom_sf(data = philly,fill=NA) +
  labs(title= "fire113count , Philly") +
  mapTheme()

#visualize type 1110
fire1110 <- 
  fire %>%
  filter(., inci_type == 1110)
ggplot() + 
  geom_sf(data =fire1110, colour="#006699", size=0.7) +
  geom_sf(data = philly,fill=NA) +
  labs(title= "fire1110count , Philly") +
  mapTheme()

#visualize type 150
fire150 <- 
  fire %>%
  filter(., inci_type == 150)
ggplot() + 
  geom_sf(data =fire150, colour="#006699", size=0.7) +
  geom_sf(data = philly,fill=NA) +
  labs(title= "fire150count , Philly") +
  mapTheme()





##### visualize the fire & fire hydrants data by neighborhoods#####
# count the fire by neighborhood
fire_neighbor <- 
  fire %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., neighbor, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count))
# plot the fire neighborhood map
ggplot() +
  geom_sf(data =fire_neighbor, aes(fill = count)) +
  scale_fill_viridis(option = "A") +
  labs(title = "Count of fire for the neighborhood") +
  mapTheme()

# count the hydrants by neighborhood
hydrants_neighbor <- 
  hydrants %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., neighbor, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count))
# plot the hydrants neighborhood map
ggplot() +
  geom_sf(data =hydrants_neighbor, aes(fill = count)) +
  scale_fill_viridis(option = "A") +
  labs(title = "Count of hydrants for the neighborhood") +
  mapTheme()

##### visualize the fire data by census tracts#####
# count the fire by tracts
fire_tracts <- 
  fire %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., tracts, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count))
# plot the fire tracts map
ggplot() +
  geom_sf(data =fire_tracts, aes(fill = count)) +
  scale_fill_viridis(option = "A") +
  labs(title = "Count of fire for the census tracts") +
  mapTheme()

# count the hydrants by tracts
hydrants_tracts <- 
  hydrants %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., tracts, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count))
# plot the hydrants tracts map
ggplot() +
  geom_sf(data =hydrants_tracts, aes(fill = count)) +
  scale_fill_viridis(option = "A") +
  labs(title = "Count of fire for the census tracts") +
  mapTheme()

##### visualize the fire data by 800*800 fishnet#####
# Fire hydrants shall be provided for detached one- and two-family dwellings in accordance with both of the following: 
# (1) The maximum distance to a fire hydrant from the closest point on the building shall not exceed 600 ft (122183 m).
# (2) The maximum distance between fire hydrants shall not exceed 800 ft (244 m).
# create the 800ft grid cell fishnet
fishnet <- 
  st_make_grid(philly, cellsize = 800) %>%
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

# join fire data to the fishnet
fire_fishnet <- 
  fire %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., fishnet, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count),
         uniqueID = rownames(.),
         cvID = sample(round(nrow(fishnet) / 24), size=nrow(fishnet), replace = TRUE))
# map the fire fishnet map
ggplot() +
  geom_sf(data = fire_fishnet, aes(fill = count)) +
  scale_fill_viridis(option = "A") +
  labs(title = "Count of fire for the fishnet") +
  mapTheme()

# join hydrant data to the fishnet
hydrants_fishnet <- 
  hydrants %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., fishnet, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count),
         uniqueID = rownames(.),
         cvID = sample(round(nrow(fishnet) / 24), size=nrow(fishnet), replace = TRUE))
# map the hydrants fishnet map
ggplot() +
  geom_sf(data = hydrants_fishnet, aes(fill = count)) +
  scale_fill_viridis(option = "A") +
  labs(title = "Count of hydrants for the fishnet") +
  mapTheme()









##### load Predictors data#####
