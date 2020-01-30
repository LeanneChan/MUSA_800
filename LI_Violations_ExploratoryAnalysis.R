library(lubridate)
library(viridis)
library(gridExtra)
library(tidyverse)
library(sf)
library(jsonlite)
library(readxl)
library(esri2sf) # could not install 
install.packages("remotes")
remotes::install_github("CityOfPhiladelphia/rphl") # could not be installed 

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

# load the fire data
fire_excel<-read_excel('2015 to 2019 fires for u of pa report run 1620.xlsx')

fire.sf <- 
  fire_excel%>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, agr = "constant") %>%
  st_set_crs(4326) %>%
  st_transform(2272)

fire.sf <- 
  st_intersection(fire.sf, philly)

fire_fishnet <- 
  fire.sf %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., fishnet, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count),
         uniqueID = rownames(.),
         cvID = sample(round(nrow(fishnet) / 24), size=nrow(fishnet), replace = TRUE))
#Map the fire fishnet map (2015-2019)
fire<-ggplot() +
  geom_sf(data = fire_fishnet, aes(fill = count), color=NA) +
  scale_fill_viridis(option = "A", direction=-1) +
  labs(title = "Count of fire incidents for the fishnet 2015-2019") +
  mapTheme()


# city limits
philly <- 
  st_read("http://data.phl.opendata.arcgis.com/datasets/405ec3da942d4e20869d4e1449a2be48_0.geojson")%>%
  st_transform(crs=2272)

fishnet <- 
  st_make_grid(philly, cellsize = 800) %>%
  st_sf()
# clip the fishnet by the boundary
fishnet <- 
  fishnet[philly,] %>%
  mutate(uniqueID = rownames(.)) %>%
  dplyr::select(uniqueID)

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

# License and Inspection Violations Data 

unsafe_structureJSON<- fromJSON("https://phl.carto.com/api/v2/sql?q=SELECT%20*%20FROM%20li_violations%20WHERE%20(violationtype%20=%20%27PM15-108.1%27)")
unsafeDF <- 
  as.data.frame(unsafe_structureJSON$rows) %>%
  select(censustract,violationdate,violationdescription, status, casestatus,casepriority, geocode_x,geocode_y) %>%
  na.omit(., cols=c("geocode_x", "geocode_y"))%>%
  st_as_sf(coords = c("geocode_x", "geocode_y"), crs = 2272) %>%
  mutate(
    X = st_coordinates(.)[,1],
    Y = st_coordinates(.)[,2])

ggplot() + 
  geom_sf(data=unsafeDF) +
  geom_sf(data=tracts) +
  geom_point(data=unsafeDF,
             aes(x=X,y=Y), size=0.2, color='red') +
  labs(title="Unsafe Buildings 2007-2020",
       subtitle=paste("Source: L&I Violations", sep="; ")) +
  mapTheme()

# FISNET PLOT 

unsafeStructure_fishnet <- 
  unsafeDF %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., fishnet, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count),
         uniqueID = rownames(.),
         cvID = sample(round(nrow(fishnet) / 24), size=nrow(fishnet), replace = TRUE))

unsafe<- ggplot() +
  geom_sf(data = unsafeStructure_fishnet, aes(fill = count),colour=NA) +
  scale_fill_viridis(direction=-1,option="A") +
  labs(title="Count of Unsafe Buildings ALL code=PM15-108.1 ",
       subtitle=paste("Source: L&I Violations", sep="; ")) +
  mapTheme()

# from past 365 days 
unsafe_structureJSON_1Yr<- fromJSON("https://phl.carto.com/api/v2/sql?q=SELECT%20*%20FROM%20li_violations%20WHERE%20(violationtype%20=%20%27PM15-108.1%27)%20AND%20(violationdate%20%3E=%20current_date%20-%20365)")  
unsafeDF_1Yr <- 
  as.data.frame(unsafe_structureJSON_1Yr$rows) %>%
  select(censustract,violationdate,violationdescription, status, casestatus,casepriority, geocode_x,geocode_y) %>%
  na.omit(., cols=c("geocode_x", "geocode_y"))%>%
  st_as_sf(coords = c("geocode_x", "geocode_y"), crs = 2272) %>%
  mutate(
    X = st_coordinates(.)[,1],
    Y = st_coordinates(.)[,2])

ggplot() +
  geom_sf(data=unsafeDF_1Yr) +
  geom_sf(data=tracts) +
  geom_point(data=unsafeDF_1Yr,
             aes(x=X,y=Y), size=0.2, color='red') +
  labs(title="Unsafe Structure Past 365 Days, code=PM15-108.1",
       subtitle=paste("Source: L&I Violations", sep="; ")) +
  mapTheme()

unsafeStructure1YR_fishnet <- 
  unsafeDF_1Yr %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., fishnet, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count),
         uniqueID = rownames(.),
         cvID = sample(round(nrow(fishnet) / 24), size=nrow(fishnet), replace = TRUE))
ggplot() +
  geom_sf(data = unsafeStructure1YR_fishnet, aes(fill = count),colour=NA) +
  scale_fill_viridis(option="A",direction=-1) +
  labs(title = "Count of unsafe structures Past Year") +
  mapTheme()

# lack smoke alarm 
lack_smokeAlarm<- fromJSON("https://phl.carto.com/api/v2/sql?q=SELECT%20*%20FROM%20li_violations%20WHERE%20(violationtype%20=%20%27FC-907.3/20%27)")
lack_smokeAlarm_DF<-as.data.frame(lack_smokeAlarm$rows) %>%
  select(censustract,violationdate,violationdescription, status, casestatus,casepriority, geocode_x,geocode_y) %>%
  na.omit(., cols=c("geocode_x", "geocode_y"))%>%
  st_as_sf(coords = c("geocode_x", "geocode_y"), crs = 2272) %>%
  mutate(
    X = st_coordinates(.)[,1],
    Y = st_coordinates(.)[,2])

ggplot() +
  geom_sf(data=lack_smokeAlarm_DF) +
  geom_sf(data=tracts) +
  geom_point(data=lack_smokeAlarm_DF,
             aes(x=X,y=Y), size=0.2, color='red') +
  labs(title="Lack Smoke Alarm ALL, code=FC-907.3/20",
       subtitle=paste("Source: L&I Violations", sep="; ")) +
  mapTheme()

LackSmokeAlarm_fishnet <- 
  lack_smokeAlarm_DF %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., fishnet, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count),
         uniqueID = rownames(.),
         cvID = sample(round(nrow(fishnet) / 24), size=nrow(fishnet), replace = TRUE))
smokeAlm<- ggplot() +
  geom_sf(data = LackSmokeAlarm_fishnet, aes(fill = count),colour=NA) +
  scale_fill_viridis(option="A",direction=-1) +
  labs(title="Count of Lack Smoke Alarm ALL, code=FC-907.3/20",
       subtitle=paste("Source: L&I Violations", sep="; "))+
  mapTheme()


# Fire extinguisher not certifies 
fire_ext_notcertified<- fromJSON("https://phl.carto.com/api/v2/sql?q=SELECT%20*%20FROM%20li_violations%20WHERE%20(violationtype%20=%20%27FC-906.2/2%27)")
fire_ext_notcertified_DF<-as.data.frame(fire_ext_notcertified$rows) %>%
  select(censustract,violationdate,violationdescription, status, casestatus,casepriority, geocode_x,geocode_y) %>%
  na.omit(., cols=c("geocode_x", "geocode_y"))%>%
  st_as_sf(coords = c("geocode_x", "geocode_y"), crs = 2272) %>%
  mutate(
    X = st_coordinates(.)[,1],
    Y = st_coordinates(.)[,2])

ggplot() +
  geom_sf(data=fire_ext_notcertified_DF) +
  geom_sf(data=tracts) +
  geom_point(data=fire_ext_notcertified_DF,
             aes(x=X,y=Y), size=0.1, color='red') +
  labs(title="Fire Extinguisher Needs Inspection (ALL), code=FC-906.2/2",
       subtitle=paste("Source: L&I Violations", sep="; ")) +
  mapTheme()

fireExtNotCertified_fishnet <- 
  fire_ext_notcertified_DF %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., fishnet, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count),
         uniqueID = rownames(.),
         cvID = sample(round(nrow(fishnet) / 24), size=nrow(fishnet), replace = TRUE))
fireExt<-ggplot() +
  geom_sf(data = fireExtNotCertified_fishnet, aes(fill = count),colour=NA) +
  scale_fill_viridis(option="A", direction=-1) +
  labs(title="Count of Fire Extinguisher Lacking Inspection (ALL), code=FC-906.2/2",
       subtitle=paste("Source: L&I Violations", sep="; "))+
  mapTheme()




grid.arrange(fire,unsafe, smokeAlm, fireExt, nrow = 2, ncol = 2)




