
########## PART 1 ##########
##### standarlize #####
setwd("C:/Users/HP/Desktop/MUSA800_practium")

library(tigris)
options(tigris_use_cache = TRUE) 
library(tidycensus)
library(viridis)
library(riem)
library(gridExtra)
library(knitr)
library(kableExtra)
library(RSocrata)
library(gifski)
library(png)
library(ggsn)
library(gganimate)
library(tidyverse)
library(sf)
library(spdep)
library(caret)
library(ckanr)
library(FNN)
library(grid)
library(gridExtra)
library(knitr)
library(kableExtra)
library(car) #test colinarity
library(ggplot2)
library(dplyr)
library(gmodels)
library(MASS)
library(e1071)
library(Hmisc)
library(Boruta)
library(corrplot)
library(viridis)
library(stargazer)
library(tigris)
library(GoodmanKruskal)
library(spatstat)
library(riem)
library(readr)
library(readxl)
library(QuantPsyc)
library(lubridate)
library(jsonlite)

# create a mapTheme() and plotTheme()function that is used to standardized the map outputs 
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
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
}

plotTheme <- function(base_size = 12) {
  theme(
    text = element_text( color = "black"),
    plot.title = element_text(size = 14,colour = "black"),
    plot.subtitle = element_text(face="italic"),
    plot.caption = element_text(hjust=0),
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_line("grey80", size = 0.1),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=2),
    strip.background = element_rect(fill = "grey80", color = "white"),
    strip.text = element_text(size=12),
    axis.title = element_text(size=12),
    axis.text = element_text(size=10),
    plot.background = element_blank(),
    legend.background = element_blank(),
    legend.title = element_text(colour = "black", face = "italic"),
    legend.text = element_text(colour = "black", face = "italic"),
    strip.text.x = element_text(size = 14)
  )
}

#Nearest neighbor (NND) function
nn_function <- function(measureFrom,measureTo,k) {
  measureFrom_Matrix <-
    as.matrix(measureFrom)
  measureTo_Matrix <-
    as.matrix(measureTo)
  nn <-   
    get.knnx(measureTo, measureFrom, k)$nn.dist
  output <-
    as.data.frame(nn) %>%
    rownames_to_column(var = "thisPoint") %>%
    gather(points, point_distance, V1:ncol(.)) %>%
    arrange(as.numeric(thisPoint)) %>%
    group_by(thisPoint) %>%
    dplyr::summarize(pointDistance = mean(point_distance)) %>%
    arrange(as.numeric(thisPoint)) %>% 
    dplyr::select(-thisPoint) %>%
    pull()
  
  return(output)  
}

##### DATA (external sources)#####
#load the boundary data
philly<- st_read("http://data.phl.opendata.arcgis.com/datasets/405ec3da942d4e20869d4e1449a2be48_0.geojson") %>%
  st_set_crs(4326) %>%
  na.omit() %>%
  st_transform(2272)

# load the neighborhood data
neighbor <- 
  st_read("C:/Users/HP/Desktop/MUSA800_practium/data/Neighborhoods_Philadelphia/Neighborhoods_Philadelphia.shp")%>%
  na.omit() %>%
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
  st_set_crs(4326) %>%
  na.omit() %>%
  st_transform(crs=2272)

##### DATA (original)#####
# load the fire data
fire_excel<-read_excel('C:/Users/HP/Desktop/MUSA800_practium/origin_data/2015 to 2019 fires for u of pa report run 1620.xlsx')
fire.sf <- 
  fire_excel%>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, agr = "constant") %>%
  st_set_crs(4326) %>%
  st_transform(2272)
fire.sf <- 
  st_intersection(fire.sf, philly)

#Plot all fire points (2015-2019)
#ggplot() +
#geom_sf(data = philly, fill = "grey40") +
#geom_sf(data = fire.sf , colour = "#FA7800", size = .75) +
#labs(title = "Philly fires") +
#mapTheme()

# load hydrant inspection data
inspect <- 
  st_read("C:/Users/HP/Desktop/MUSA800_practium/origin_data/PDF_Hydrant_Inspections_2019_for_MUSA/PDF_Hydrant_Inspections_2019.shp") %>%
  st_transform(crs=2272) 

# load the hydrant data
hydrants <- 
  st_read("C:/Users/HP/Desktop/MUSA800_practium/origin_data/Hydrants_2019/Hydrants_2019.shp") %>%
  st_transform(crs=2272) 

# load the engine locals data
engine_new <- 
  st_read("C:/Users/HP/Desktop/MUSA800_practium/origin_data/EngineLocals/EngineLocals_New.shp") %>%
  st_transform(crs=2272) 
engine_old <- 
  st_read("C:/Users/HP/Desktop/MUSA800_practium/origin_data/EngineLocals/EngineLocals_Old.shp") %>%
  st_transform(crs=2272) 










########## PART 2 ##########
##### /new/ visualize the fire & fire hydrants data by neighborhoods#####
# count the fire by neighborhood
fire_neighbor <- 
  fire.sf %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., neighbor, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count))

# count the hydrants by neighborhood
hydrants_neighbor <- 
  hydrants %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., neighbor, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count))

# plot the the fire neighborhood map & hydrants neighborhood map
grid.arrange(ggplot() +
               geom_sf(data =fire_neighbor, aes(fill = count)) +
               scale_fill_viridis(option = "A") +
               labs(title = "Count of fire for the neighborhood") +
               mapTheme(), 
             ggplot() +
               geom_sf(data =hydrants_neighbor, aes(fill = count)) +
               scale_fill_viridis(option = "A") +
               labs(title = "Count of hydrants for the neighborhood") +
               mapTheme(),ncol=2)

##### /new/ visualize the fire & fire hydrants data by census tracts#####
# count the fire by tracts
fire_tracts <- 
  fire.sf %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., tracts, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count))

# count the hydrants by tracts
hydrants_tracts <- 
  hydrants %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., tracts, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count))

# plot the fire tracts map &  the hydrants tracts map
grid.arrange(ggplot() +
               geom_sf(data =fire_tracts, aes(fill = count)) +
               scale_fill_viridis(option = "A") +
               labs(title = "Count of fire for the census tracts") +
               mapTheme(),
             ggplot() +
               geom_sf(data =hydrants_tracts, aes(fill = count)) +
               scale_fill_viridis(option = "A") +
               labs(title = "Count of hydrants for the census tracts") +
               mapTheme(),ncol=2)

##### /new/ visualize the fire & fire hydrants data by 800*800 fishnet#####
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
  geom_sf(data=fishnet,fill=NA) +
  labs(title = "Philadelphia fishnet") +
  mapTheme()

# join fire data to the fishnet
fire_fishnet <- 
  fire.sf %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., fishnet, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count),
         uniqueID = rownames(.),
         cvID = sample(round(nrow(fishnet) / 24), size=nrow(fishnet), replace = TRUE))

# join hydrant data to the fishnet
hydrants_fishnet <- 
  hydrants %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., fishnet, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count),
         uniqueID = rownames(.),
         cvID = sample(round(nrow(fishnet) / 24), size=nrow(fishnet), replace = TRUE))

#Map the fire fishnet map and hydrants fishnet map
grid.arrange(ggplot() +
               geom_sf(data = fire_fishnet, aes(fill = count)) +
               scale_fill_viridis(option = "A") +
               labs(title = "Count of fire incidents for the fishnet 2015-2019") +
               mapTheme(),ggplot() +
               geom_sf(data = hydrants_fishnet, aes(fill = count)) +
               scale_fill_viridis(option = "A") +
               labs(title = "Count of hydrants for the fishnet") +
               mapTheme(),
             ncol=2)






########## PART 3 ##########
##### /new/ explore leading cause of fires across philly #####
reorder_size <- function(x) {
  factor(x, levels = names(sort(table(x),decreasing=T)))
}
ggplot(fire.sf, aes(reorder_size(inci_type))) + 
  geom_bar() +  
  theme(axis.text.x = element_text(angle = 60))

# fire113
fire113 <- 
  fire.sf %>%
  filter(., inci_type == 113)

fire113_fishnet <- 
  fire113 %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., fishnet, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count),
         uniqueID = rownames(.),
         cvID = sample(round(nrow(fishnet) / 24), size=nrow(fishnet), replace = TRUE))

fire113_ppp <- as.ppp(st_coordinates(fire113), W = st_bbox(fire113_fishnet))
fire113_KD <- spatstat::density.ppp(fire113_ppp, 1000)

# fire1110
fire1110 <- 
  fire.sf %>%
  filter(., inci_type == 1110)
fire1110_fishnet <- 
  fire1110 %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., fishnet, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count),
         uniqueID = rownames(.),
         cvID = sample(round(nrow(fishnet) / 24), size=nrow(fishnet), replace = TRUE))

fire1110_ppp <- as.ppp(st_coordinates(fire1110), W = st_bbox(fire1110_fishnet))
fire1110_KD <- spatstat::density.ppp(fire1110_ppp, 1000)


#fire 150
fire150 <- 
  fire.sf %>%
  filter(., inci_type == 150)

fire150_fishnet <- 
  fire150 %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., fishnet, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count),
         uniqueID = rownames(.),
         cvID = sample(round(nrow(fishnet) / 24), size=nrow(fishnet), replace = TRUE))

fire150_ppp <- as.ppp(st_coordinates(fire150), W = st_bbox(fire150_fishnet))
fire150_KD <- spatstat::density.ppp(fire150_ppp, 1000)


grid.arrange(ggplot() + 
               geom_sf(data =fire113, colour="#006699", size=0.7) +
               geom_sf(data = philly,fill=NA) +
               labs(title= "fire150count , Philly") +
               mapTheme(),
             ggplot() + 
               geom_sf(data =fire1110, colour="#006699", size=0.7) +
               geom_sf(data = philly,fill=NA) +
               labs(title= "fire150count , Philly") +
               mapTheme(),
             ggplot() + 
               geom_sf(data =fire150, colour="#006699", size=0.7) +
               geom_sf(data = philly,fill=NA) +
               labs(title= "fire150count , Philly") +
               mapTheme(),nrow=1)

grid.arrange(as.data.frame(fire113_KD) %>%
               st_as_sf(coords = c("x", "y"), crs = st_crs(hydrants_fishnet)) %>%
               aggregate(.,fire113_fishnet, mean) %>%
               ggplot() +
               geom_sf(aes(fill=value),colour=NA) +
               ##geom_sf(data = sample_n(fire.sf, 1500), size = .5) +
               labs(title = "Kernel density of fire 113") +
               scale_fill_viridis() +
               mapTheme(),
             as.data.frame(fire1110_KD) %>%
               st_as_sf(coords = c("x", "y"), crs = st_crs(hydrants_fishnet)) %>%
               aggregate(.,fire1110_fishnet, mean) %>%
               ggplot() +
               geom_sf(aes(fill=value),colour=NA) +
               ##geom_sf(data = sample_n(fire.sf, 1500), size = .5) +
               labs(title = "Kernel density of fire 1110") +
               scale_fill_viridis() +
               mapTheme() ,
             as.data.frame(fire150_KD) %>%
               st_as_sf(coords = c("x", "y"), crs = st_crs(hydrants_fishnet)) %>%
               aggregate(.,fire150_fishnet, mean) %>%
               ggplot() +
               geom_sf(aes(fill=value),colour=NA) +
               labs(title = "Kernel density of fire 150") +
               scale_fill_viridis() +
               mapTheme(), nrow=1)

##### /new/ explore fire count & fire kernel density by year and month #####
#0.Standarize date
fire.sf <- 
  fire.sf %>%
  mutate(date=dmy(alm_date),
         Year=year(date),
         Month=month(date))

#1. Histogram for total count 2015-2019
ggplot(fire_fishnet, aes(count)) + 
  geom_histogram(binwidth = 1) +
  labs(title = "Distribution of fire incidents by grid cell 2015-2019")

#2. Visualize by year 
#2015
fire2015<-
  fire.sf %>%
  filter(Year == 2015)
# join fire data to the fishnet
fire2015_fishnet <- 
  fire2015 %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., fishnet, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count),
         uniqueID = rownames(.),
         cvID = sample(round(nrow(fishnet) / 24), size=nrow(fishnet), replace = TRUE))

#2016
fire2016<-
  fire.sf %>%
  filter(Year == 2016)
# join fire data to the fishnet
fire2016_fishnet <- 
  fire2016 %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., fishnet, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count),
         uniqueID = rownames(.),
         cvID = sample(round(nrow(fishnet) / 24), size=nrow(fishnet), replace = TRUE))

#2017
fire2017<-
  fire.sf %>%
  filter(Year == 2017)
# join fire data to the fishnet
fire2017_fishnet <- 
  fire2017 %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., fishnet, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count),
         uniqueID = rownames(.),
         cvID = sample(round(nrow(fishnet) / 24), size=nrow(fishnet), replace = TRUE))

#2018
fire2018<-
  fire.sf %>%
  filter(Year == 2018)
# join fire data to the fishnet
fire2018_fishnet <- 
  fire2018 %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., fishnet, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count),
         uniqueID = rownames(.),
         cvID = sample(round(nrow(fishnet) / 24), size=nrow(fishnet), replace = TRUE))

#2019
fire2019<-
  fire.sf %>%
  filter(Year == 2019)
# join fire data to the fishnet
fire2019_fishnet <- 
  fire2019 %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., fishnet, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count),
         uniqueID = rownames(.),
         cvID = sample(round(nrow(fishnet) / 24), size=nrow(fishnet), replace = TRUE))

#Map the fire fishnet map (2015-2019)
grid.arrange(
  ggplot() +
    geom_sf(data = fire2015_fishnet, aes(fill = count),colour=NA) +
    scale_fill_viridis(option = "A") +
    labs(title = "Count of fire incidents for the fishnet 2015") +
    mapTheme(),
  
  ggplot() +
    geom_sf(data = fire2016_fishnet, aes(fill = count),colour=NA) +
    scale_fill_viridis(option = "A") +
    labs(title = "Count of fire incidents for the fishnet 2016") +
    mapTheme(),
  
  ggplot() +
    geom_sf(data = fire2017_fishnet, aes(fill = count),colour=NA) +
    scale_fill_viridis(option = "A") +
    labs(title = "Count of fire incidents for the fishnet 2017") +
    mapTheme(),
  
  ggplot() +
    geom_sf(data = fire2018_fishnet, aes(fill = count),colour=NA) +
    scale_fill_viridis(option = "A") +
    labs(title = "Count of fire incidents for the fishnet 2018") +
    mapTheme(),
  
  ggplot() +
    geom_sf(data = fire2019_fishnet, aes(fill = count),colour=NA) +
    scale_fill_viridis(option = "A") +
    labs(title = "Count of fire incidents for the fishnet 2019") +
    mapTheme(),
  ncol=3)
``

#Total fire incidents 2015-2019
fire_ppp <- as.ppp(st_coordinates(fire.sf), W = st_bbox(fire_fishnet))
fire_KD <- spatstat::density.ppp(fire_ppp, 1000)
as.data.frame(fire_KD) %>%
  st_as_sf(coords = c("x", "y"), crs = st_crs(fire_fishnet)) %>%
  aggregate(.,fire_fishnet, mean) %>%
  ggplot() +
  geom_sf(aes(fill=value),colour=NA) +
  labs(title = "Kernel density of total fire incidents 2015-2019") +
  scale_fill_viridis(option='inferno') +
  mapTheme()

#By each year
fire_ppp_2015 <- as.ppp(st_coordinates(fire2015), W = st_bbox(fire2015_fishnet))
fire_KD_2015 <- spatstat::density.ppp(fire_ppp_2015, 1000)

fire_ppp_2016 <- as.ppp(st_coordinates(fire2016), W = st_bbox(fire2016_fishnet))
fire_KD_2016 <- spatstat::density.ppp(fire_ppp_2016, 1000)

fire_ppp_2017 <- as.ppp(st_coordinates(fire2017), W = st_bbox(fire2017_fishnet))
fire_KD_2017<- spatstat::density.ppp(fire_ppp_2017, 1000)

fire_ppp_2018 <- as.ppp(st_coordinates(fire2018), W = st_bbox(fire2018_fishnet))
fire_KD_2018 <- spatstat::density.ppp(fire_ppp_2018, 1000)

fire_ppp_2019 <- as.ppp(st_coordinates(fire2019), W = st_bbox(fire2019_fishnet))
fire_KD_2019 <- spatstat::density.ppp(fire_ppp_2019, 1000)

grid.arrange( #Plot
  as.data.frame(fire_KD_2015) %>% #2015
    st_as_sf(coords = c("x", "y"), crs = st_crs(fire2015_fishnet)) %>%
    aggregate(.,fire2015_fishnet, mean) %>%
    ggplot() +
    geom_sf(aes(fill=value),colour=NA) +
    labs(title = "Kernel density of fire incidents 2015") +
    scale_fill_viridis(option='inferno') +
    mapTheme(), 
  
  as.data.frame(fire_KD_2016) %>% #2016
    st_as_sf(coords = c("x", "y"), crs = st_crs(fire2016_fishnet)) %>%
    aggregate(.,fire2016_fishnet, mean) %>%
    ggplot() +
    geom_sf(aes(fill=value),colour=NA) +
    labs(title = "Kernel density of fire incidents 2016") +
    scale_fill_viridis(option='inferno') +
    mapTheme(), 
  
  as.data.frame(fire_KD_2017) %>% #2017
    st_as_sf(coords = c("x", "y"), crs = st_crs(fire2017_fishnet)) %>%
    aggregate(.,fire2017_fishnet, mean) %>%
    ggplot() +
    geom_sf(aes(fill=value),colour=NA) +
    labs(title = "Kernel density of fire incidents 2017") +
    scale_fill_viridis(option='inferno') +
    mapTheme(), 
  
  as.data.frame(fire_KD_2018) %>% #2018
    st_as_sf(coords = c("x", "y"), crs = st_crs(fire2018_fishnet)) %>%
    aggregate(.,fire2018_fishnet, mean) %>%
    ggplot() +
    geom_sf(aes(fill=value),colour=NA) +
    labs(title = "Kernel density of fire incidents 2018") +
    scale_fill_viridis(option='inferno') +
    mapTheme(), 
  
  as.data.frame(fire_KD_2019) %>% #2019
    st_as_sf(coords = c("x", "y"), crs = st_crs(fire2019_fishnet)) %>%
    aggregate(.,fire2019_fishnet, mean) %>%
    ggplot() +
    geom_sf(aes(fill=value),colour=NA) +
    labs(title = "Kernel density of fire incidents 2019") +
    scale_fill_viridis(option='inferno') +
    mapTheme(), 
  ncol=3)

#3. Visualize by month
#Jan
fireJan<-
  fire.sf %>%
  filter(Month == 1)
fireJan_fishnet <- 
  fireJan %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., fishnet, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count),
         uniqueID = rownames(.),
         cvID = sample(round(nrow(fishnet) / 24), size=nrow(fishnet), replace = TRUE))
fire_ppp_Jan <- as.ppp(st_coordinates(fireJan), W = st_bbox(fireJan_fishnet))
fire_KD_Jan <- spatstat::density.ppp(fire_ppp_Jan, 1000)

#Feb
fireFeb<-
  fire.sf %>%
  filter(Month == 2)

fireFeb_fishnet <- 
  fireFeb %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., fishnet, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count),
         uniqueID = rownames(.),
         cvID = sample(round(nrow(fishnet) / 24), size=nrow(fishnet), replace = TRUE))
fire_ppp_Feb <- as.ppp(st_coordinates(fireFeb), W = st_bbox(fireFeb_fishnet))
fire_KD_Feb <- spatstat::density.ppp(fire_ppp_Feb, 1000)

#Mar
fireMar<-
  fire.sf %>%
  filter(Month == 3)

fireMar_fishnet <- 
  fireMar %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., fishnet, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count),
         uniqueID = rownames(.),
         cvID = sample(round(nrow(fishnet) / 24), size=nrow(fishnet), replace = TRUE))
fire_ppp_Mar <- as.ppp(st_coordinates(fireMar), W = st_bbox(fireMar_fishnet))
fire_KD_Mar <- spatstat::density.ppp(fire_ppp_Mar, 1000)

#Apr
fireApr<-
  fire.sf %>%
  filter(Month == 4)

fireApr_fishnet <- 
  fireApr %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., fishnet, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count),
         uniqueID = rownames(.),
         cvID = sample(round(nrow(fishnet) / 24), size=nrow(fishnet), replace = TRUE))
fire_ppp_Apr <- as.ppp(st_coordinates(fireApr), W = st_bbox(fireApr_fishnet))
fire_KD_Apr <- spatstat::density.ppp(fire_ppp_Apr, 1000)

#May
fireMay<-
  fire.sf %>%
  filter(Month == 5)

fireMay_fishnet <- 
  fireMay %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., fishnet, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count),
         uniqueID = rownames(.),
         cvID = sample(round(nrow(fishnet) / 24), size=nrow(fishnet), replace = TRUE))
fire_ppp_May <- as.ppp(st_coordinates(fireMay), W = st_bbox(fireMay_fishnet))
fire_KD_May <- spatstat::density.ppp(fire_ppp_May, 1000)

#Jun
fireJun<-
  fire.sf %>%
  filter(Month == 6)

fireJun_fishnet <- 
  fireJun %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., fishnet, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count),
         uniqueID = rownames(.),
         cvID = sample(round(nrow(fishnet) / 24), size=nrow(fishnet), replace = TRUE))
fire_ppp_Jun <- as.ppp(st_coordinates(fireJun), W = st_bbox(fireJun_fishnet))
fire_KD_Jun <- spatstat::density.ppp(fire_ppp_Jun, 1000)

#Jul
fireJul<-
  fire.sf %>%
  filter(Month == 7)

fireJul_fishnet <- 
  fireJul %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., fishnet, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count),
         uniqueID = rownames(.),
         cvID = sample(round(nrow(fishnet) / 24), size=nrow(fishnet), replace = TRUE))
fire_ppp_Jul <- as.ppp(st_coordinates(fireJul), W = st_bbox(fireJul_fishnet))
fire_KD_Jul <- spatstat::density.ppp(fire_ppp_Jul, 1000)

#Aug
fireAug<-
  fire.sf %>%
  filter(Month == 8)

fireAug_fishnet <- 
  fireAug %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., fishnet, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count),
         uniqueID = rownames(.),
         cvID = sample(round(nrow(fishnet) / 24), size=nrow(fishnet), replace = TRUE))
fire_ppp_Aug <- as.ppp(st_coordinates(fireAug), W = st_bbox(fireAug_fishnet))
fire_KD_Aug <- spatstat::density.ppp(fire_ppp_Aug, 1000)

#Sep
fireSep<-
  fire.sf %>%
  filter(Month == 9)

fireSep_fishnet <- 
  fireSep %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., fishnet, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count),
         uniqueID = rownames(.),
         cvID = sample(round(nrow(fishnet) / 24), size=nrow(fishnet), replace = TRUE))
fire_ppp_Sep <- as.ppp(st_coordinates(fireSep), W = st_bbox(fireSep_fishnet))
fire_KD_Sep <- spatstat::density.ppp(fire_ppp_Sep, 1000)

#Oct
fireOct<-
  fire.sf %>%
  filter(Month == 10)

fireOct_fishnet <- 
  fireOct %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., fishnet, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count),
         uniqueID = rownames(.),
         cvID = sample(round(nrow(fishnet) / 24), size=nrow(fishnet), replace = TRUE))
fire_ppp_Oct <- as.ppp(st_coordinates(fireOct), W = st_bbox(fireOct_fishnet))
fire_KD_Oct <- spatstat::density.ppp(fire_ppp_Oct, 1000)

#Nov
fireNov<-
  fire.sf %>%
  filter(Month == 11)

fireNov_fishnet <- 
  fireNov %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., fishnet, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count),
         uniqueID = rownames(.),
         cvID = sample(round(nrow(fishnet) / 24), size=nrow(fishnet), replace = TRUE))
fire_ppp_Nov <- as.ppp(st_coordinates(fireNov), W = st_bbox(fireNov_fishnet))
fire_KD_Nov <- spatstat::density.ppp(fire_ppp_Nov, 1000)

#Dec
fireDec<-
  fire.sf %>%
  filter(Month == 12)

fireDec_fishnet <- 
  fireDec %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., fishnet, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count),
         uniqueID = rownames(.),
         cvID = sample(round(nrow(fishnet) / 24), size=nrow(fishnet), replace = TRUE))
fire_ppp_Dec <- as.ppp(st_coordinates(fireDec), W = st_bbox(fireDec_fishnet))
fire_KD_Dec <- spatstat::density.ppp(fire_ppp_Dec, 1000)

grid.arrange(
  as.data.frame(fire_KD_Jan) %>% #1
    st_as_sf(coords = c("x", "y"), crs = st_crs(fireJan_fishnet)) %>%
    aggregate(.,fireJan_fishnet, mean) %>%
    ggplot() +
    geom_sf(aes(fill=value),colour=NA) +
    labs(title = "Kernel density of fire incidents Jan") +
    scale_fill_viridis(option='inferno') +
    mapTheme(),
  as.data.frame(fire_KD_Feb) %>% #2
    st_as_sf(coords = c("x", "y"), crs = st_crs(fireJan_fishnet)) %>%
    aggregate(.,fireJan_fishnet, mean) %>%
    ggplot() +
    geom_sf(aes(fill=value),colour=NA) +
    labs(title = "Kernel density of fire incidents Feb") +
    scale_fill_viridis(option='inferno') +
    mapTheme(),
  as.data.frame(fire_KD_Mar) %>% #3
    st_as_sf(coords = c("x", "y"), crs = st_crs(fireJan_fishnet)) %>%
    aggregate(.,fireJan_fishnet, mean) %>%
    ggplot() +
    geom_sf(aes(fill=value),colour=NA) +
    labs(title = "Kernel density of fire incidents Mar") +
    scale_fill_viridis(option='inferno') +
    mapTheme(),
  as.data.frame(fire_KD_Apr) %>% #4
    st_as_sf(coords = c("x", "y"), crs = st_crs(fireJan_fishnet)) %>%
    aggregate(.,fireJan_fishnet, mean) %>%
    ggplot() +
    geom_sf(aes(fill=value),colour=NA) +
    labs(title = "Kernel density of fire incidents Apr") +
    scale_fill_viridis(option='inferno') +
    mapTheme(),
  as.data.frame(fire_KD_May) %>% #5
    st_as_sf(coords = c("x", "y"), crs = st_crs(fireJan_fishnet)) %>%
    aggregate(.,fireJan_fishnet, mean) %>%
    ggplot() +
    geom_sf(aes(fill=value),colour=NA) +
    labs(title = "Kernel density of fire incidents May") +
    scale_fill_viridis(option='inferno') +
    mapTheme(),
  as.data.frame(fire_KD_Jun) %>% #6
    st_as_sf(coords = c("x", "y"), crs = st_crs(fireJan_fishnet)) %>%
    aggregate(.,fireJan_fishnet, mean) %>%
    ggplot() +
    geom_sf(aes(fill=value),colour=NA) +
    labs(title = "Kernel density of fire incidents Jun") +
    scale_fill_viridis(option='inferno') +
    mapTheme(),
  as.data.frame(fire_KD_Jul) %>% #7
    st_as_sf(coords = c("x", "y"), crs = st_crs(fireJan_fishnet)) %>%
    aggregate(.,fireJan_fishnet, mean) %>%
    ggplot() +
    geom_sf(aes(fill=value),colour=NA) +
    labs(title = "Kernel density of fire incidents Jul") +
    scale_fill_viridis(option='inferno') +
    mapTheme(),
  as.data.frame(fire_KD_Aug) %>% #8
    st_as_sf(coords = c("x", "y"), crs = st_crs(fireJan_fishnet)) %>%
    aggregate(.,fireJan_fishnet, mean) %>%
    ggplot() +
    geom_sf(aes(fill=value),colour=NA) +
    labs(title = "Kernel density of fire incidents Aug") +
    scale_fill_viridis(option='inferno') +
    mapTheme(),
  as.data.frame(fire_KD_Sep) %>% #9
    st_as_sf(coords = c("x", "y"), crs = st_crs(fireJan_fishnet)) %>%
    aggregate(.,fireJan_fishnet, mean) %>%
    ggplot() +
    geom_sf(aes(fill=value),colour=NA) +
    labs(title = "Kernel density of fire incidents Sep") +
    scale_fill_viridis(option='inferno') +
    mapTheme(),
  as.data.frame(fire_KD_Oct) %>% #10
    st_as_sf(coords = c("x", "y"), crs = st_crs(fireJan_fishnet)) %>%
    aggregate(.,fireJan_fishnet, mean) %>%
    ggplot() +
    geom_sf(aes(fill=value),colour=NA) +
    labs(title = "Kernel density of fire incidents Oct") +
    scale_fill_viridis(option='inferno') +
    mapTheme(),
  as.data.frame(fire_KD_Nov) %>% #11
    st_as_sf(coords = c("x", "y"), crs = st_crs(fireJan_fishnet)) %>%
    aggregate(.,fireJan_fishnet, mean) %>%
    ggplot() +
    geom_sf(aes(fill=value),colour=NA) +
    labs(title = "Kernel density of fire incidents Nov") +
    scale_fill_viridis(option='inferno') +
    mapTheme(),
  as.data.frame(fire_KD_Dec) %>% #12
    st_as_sf(coords = c("x", "y"), crs = st_crs(fireJan_fishnet)) %>%
    aggregate(.,fireJan_fishnet, mean) %>%
    ggplot() +
    geom_sf(aes(fill=value),colour=NA) +
    labs(title = "Kernel density of fire incidents Dec") +
    scale_fill_viridis(option='inferno') +
    mapTheme(),
  ncol=3
)


##### Local Moran's I statistics of fire count #####
fire_fishnet.nb <- poly2nb(fire_fishnet, queen=TRUE)
fire_fishnet.weights <- nb2listw(fire_fishnet.nb, style="W", zero.policy=TRUE)

fire_fishnet.localMorans <- 
  cbind(
    as.data.frame(localmoran(fire_fishnet$count, fire_fishnet.weights)),
    as.data.frame(fire_fishnet, NULL)) %>% 
  st_sf() %>%
  dplyr::select(Fire_TotalCount = count, 
                Local_Morans_I = Ii, 
                P_Value = `Pr(z > 0)`) %>%
  mutate(Significant_Hotspots = ifelse(P_Value <= 0.05, 1, 0)) %>%
  gather(Variable, Value, -geometry)

vars <- unique(fire_fishnet.localMorans$Variable)
varList <- list()

for(i in vars){
  varList[[i]] <- 
    ggplot() +
    geom_sf(data = filter(fire_fishnet.localMorans, Variable == i), aes(fill = Value), colour=NA) +
    scale_fill_viridis(name="") +
    labs(title=i) +
    mapTheme()}
#Plot
do.call(grid.arrange,c(varList, ncol = 2, top = "Local Morans I statistics, Fire incidents in Philadelphia 2015-2019"))

##### Hotspots of fire count #####
fire_fishnet <-
  fire_fishnet %>% 
  mutate(fire.isSig = ifelse(localmoran(fire_fishnet$count, 
                                        fire_fishnet.weights)[,5] <= 0.0000001, "Statistically significant", "Not significant")) 
#%>%
#mutate(fire.isSig.dist = nn_function(st_coordinates(st_centroid(fire_fishnet)),
#st_coordinates(st_centroid(
#filter(fire_fishnet, fire.isSig == 1))), 1 ))

#Plot
ggplot() + 
  geom_sf(data = fire_fishnet, aes(fill = fire.isSig),colour=NA) +
  labs(title = "Significantly clustered fire incidents (local hotspots)") +
  scale_fill_manual(values = c("grey", "red"), name = "Values")+
  mapTheme()
##### kernel density of hydrants #####
#Kernel density of hydrants
hydrants_ppp <- as.ppp(st_coordinates(hydrants), W = st_bbox(hydrants_fishnet))
hydrants_KD <- spatstat::density.ppp(hydrants_ppp, 1000)
as.data.frame(hydrants_KD) %>%
  st_as_sf(coords = c("x", "y"), crs = st_crs(hydrants_fishnet)) %>%
  aggregate(.,hydrants_fishnet, mean) %>%
  ggplot() +
  geom_sf(aes(fill=value),colour=NA) +
  ##geom_sf(data = sample_n(fire.sf, 1500), size = .5) +
  labs(title = "Kernel density of fire hydrants") +
  scale_fill_viridis() +
  mapTheme()

##### kernel density map of NND to hydrants #####
fire_fishnet$hydrants.nn =
  nn_function(st_coordinates(st_centroid(fire_fishnet)), st_coordinates(hydrants), 10)

ggplot() +
  geom_sf(data=fire_fishnet, aes(fill=hydrants.nn),colour=NA) +
  scale_fill_viridis() +
  labs(title="Nearest Neighbor Distance to Hydrants (10)") +
  mapTheme()





########## PART 4 ##########
##### /new/  explore the relationship between year/month and fire count #####
ggplot(fire.sf, aes(Year)) + 
  geom_bar() +
  theme(axis.text.x = element_text(angle = 60))

ggplot(fire.sf, aes(Month)) + 
  geom_bar() +
  theme(axis.text.x = element_text(angle = 60))
##### /new/ explore the relationship between zipcode and fire count #####
ggplot(fire.sf, aes(reorder_size(zip))) + 
  geom_bar() +
  theme(axis.text.x = element_text(angle = 60))

##### /new/ explore the relationship between zoning type and fire count#####
#load the zoning type data
landuse <- st_read("http://data-phl.opendata.arcgis.com/datasets/e433504739bd41049de5d8f4a22d34ba_0.geojson") %>%
  st_set_crs(4326) %>%
  st_transform(2272)
#join the landuse data to the fire data
fire.sf.zone <- 
  fire.sf %>%
  st_join(.,landuse,join=st_intersects, left= TRUE) %>%
  dplyr::select(inci_no,C_DIG1DESC)%>%
  na.omit()
# visualize the relationship between zoning type and fire count
ggplot(fire.sf.zone, aes(reorder_size(C_DIG1DESC)))+  geom_bar() +  
  theme(axis.text.x = element_text(angle = 60))

##### /new/  explore the relationship between hydrant density and fire density#####
fire_hydrant_fishnet <- merge(x=as.data.frame(fire_fishnet),y=as.data.frame(hydrants_fishnet),by="uniqueID",all=TRUE)
ggplot(data=fire_hydrant_fishnet,aes(x=count.y,y=count.x)) +
  geom_point(colour ="#006699",size=0.8) +
  geom_smooth(method = "lm", se=F, colour = "#336699") +
  labs(title = "Fire counts as a function of Hydrants counts") +
  plotTheme()

