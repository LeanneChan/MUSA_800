
########## Preperation ##########
##### standarlize function #####
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
library(readr)
library(ggplot2)
library(gganimate)
library(scales)
library(dplyr)
library(gifski)
library(png)
library(ggsn)
library(ggpmisc) # to get R squared

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
    #panel.border = element_rect(colour = "black", fill=NA, size=1)
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

##### count/kernel density/knn functions #####
# count function 
count_net <- function(data,fishnetGrid, name){
  count<- data %>% #get count by fishnet
    dplyr::select() %>% 
    mutate(count = 1) %>% 
    aggregate(., fishnetGrid, sum) %>%
    mutate(count = ifelse(is.na(count), 0, count),uniqueID = rownames(.))%>%
    plyr::rename(.,c('count' = name))
  return (count)
}

# kernel density function 
kernelDensity <- function(data, fishnetGrid, name){
  count<- data %>% #get count by fishnet
    dplyr::select() %>% 
    mutate(count = 1) %>% 
    aggregate(., fishnetGrid, sum) %>%
    mutate(count = ifelse(is.na(count), 0, count),
           uniqueID = rownames(.),
           cvID = sample(round(nrow(fishnetGrid) / 24), size=nrow(fishnetGrid), replace = TRUE))
  
  point_ppp <- as.ppp(st_coordinates(data), W = st_bbox(count))
  KD <- spatstat::density.ppp(point_ppp, 1000)
  
  KD_dataframe<-as.data.frame(KD) %>% 
    st_as_sf(coords = c("x", "y"), crs = st_crs(fishnetGrid)) %>%
    aggregate(.,fishnetGrid, mean)%>%plyr::rename(.,c('value' = name))
  
  return(KD_dataframe)
}

# knn function 
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

##### Creat 800*800 fishnet #####
#load the philly boundary data
philly<- st_read("http://data.phl.opendata.arcgis.com/datasets/405ec3da942d4e20869d4e1449a2be48_0.geojson") %>%
  st_set_crs(4326) %>%
  na.omit() %>%
  st_transform(2272)

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

##### Load 2018 fire data to the fishnet#####
# load the fire data
fire_excel<-read_excel('C:/Users/HP/Desktop/MUSA800_practium/origin_data/2015 to 2019 fires for u of pa report run 1620.xlsx')
fire.sf <- 
  fire_excel%>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, agr = "constant") %>%
  st_set_crs(4326) %>%
  st_transform(2272)
fire.sf <- 
  st_intersection(fire.sf, philly)
fire.sf <- 
  fire.sf %>%
  mutate(date=dmy(alm_date),
         Year=year(date),
         Month=month(date))
fire_18 <- fire.sf[fire.sf$Year==2018,]

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


# join fire data to the fishnet
fire_fishnet <- 
  fire_18 %>% 
  dplyr::select() %>% 
  mutate(count = 1) %>% 
  aggregate(., fishnet, sum) %>%
  mutate(count = ifelse(is.na(count), 0, count),
         uniqueID = rownames(.),
         cvID = sample(round(nrow(fishnet) / 24), size=nrow(fishnet), replace = TRUE))




########## Feature engineering ##########
##### Fire is.Sig & Fire distance #####
# define the significant fire spot (set the significant level as 0.05)
fire_fishnet.nb <- poly2nb(fire_fishnet, queen=TRUE)
fire_fishnet.weights <- nb2listw(fire_fishnet.nb, style="W", zero.policy=TRUE)
# add isSig and isSig.dist variables
fire_fishnet <-
  fire_fishnet %>% 
  mutate(fire.isSig = ifelse(localmoran(fire_fishnet$count, 
                                        fire_fishnet.weights)[,5] <= 0.05, 1, 0)) %>%
  mutate(fire.isSig.dist = nn_function(st_coordinates(st_centroid(fire_fishnet)),
                                       st_coordinates(st_centroid(
                                         filter(fire_fishnet, fire.isSig == 1))), 1 ))

##### Hydrant (count/kernel density---KD is better) #####
# add hydrant count
var_fishnet <- 
  count_net(hydrants,fishnet,"hydrant_count") %>%
  as.data.frame(.)%>%
  dplyr::select(hydrant_count,uniqueID)%>%
  merge(as.data.frame(fire_fishnet),. , by="uniqueID") %>%
  st_sf()
# add hydrant kernel density
var_fishnet <- 
  kernelDensity(hydrants, fishnet, 'hydrant_KD') %>%
  st_join(.,var_fishnet,join=st_equals)

##### 2018 L& I violations (count/kernel density/knn---knn is better) #####
# load L&I violations data, following carto API. 
readData<-function(link){
  raw<-fromJSON(link)
  formated<- as.data.frame(raw$rows) %>%
    dplyr::select(censustract,violationdate,violationdescription, status, casestatus, casepriority, geocode_x,geocode_y) %>%
    mutate('Year' = year(ymd_hms(violationdate)))%>%
    na.omit(., cols=c("geocode_x", "geocode_y"))%>%
    st_as_sf(coords = c("geocode_x", "geocode_y"), crs = 2272) %>%
    mutate(
      X = st_coordinates(.)[,1],
      Y = st_coordinates(.)[,2])
  return(formated)
}
# unsafe structures
unsafe_structures<-readData("https://phl.carto.com/api/v2/sql?q=SELECT%20*%20FROM%20li_violations%20WHERE%20(violationtype%20=%20%27PM15-108.1%27)%20AND%20(violationdate%20%3E=%20%272015-01-01%27)")
unsafe_18 <- unsafe_structures[unsafe_structures$Year==2018,]
# Imminently Dangerous
ID_structures<-readData("https://phl.carto.com/api/v2/sql?q=SELECT%20*%20FROM%20li_violations%20WHERE%20(violationtype%20=%20%27PM15-110.1%27)%20AND%20(violationdate%20%3E=%20%272015-01-01%27)")
ID_18 <- ID_structures[ID_structures$Year==2018,]
# Lacking Proper Fire Extinguisher for Commercial Cooking
CE_structures<-readData("https://phl.carto.com/api/v2/sql?q=SELECT%20*%20FROM%20li_violations%20WHERE%20(violationtype%20=%20%27FC-904.11/5%27)%20AND%20(violationdate%20%3E=%20%272015-01-01%27)")
CE_18 <- CE_structures[CE_structures$Year==2018,]
# Lacking smoke Alarm
ALM_structures<-readData("https://phl.carto.com/api/v2/sql?q=SELECT%20*%20FROM%20li_violations%20WHERE%20(violationtype%20=%20%27FC-907.3/20%27)%20AND%20(violationdate%20%3E=%20%272015-01-01%27)")
ALM_18 <- ALM_structures[ALM_structures$Year==2018,]

# # add unsafe_18 count
# var_fishnet <- 
#   count_net(unsafe_18,fishnet,"unsafe_18_count") %>%
#   as.data.frame(.)%>%
#   dplyr::select(unsafe_18_count,uniqueID)%>%
#   merge(as.data.frame(var_fishnet),. , by="uniqueID") %>%
#   st_sf()
# # add unsafe_18 kernel density
# var_fishnet <- 
#   kernelDensity(unsafe_18, fishnet, 'unsafe_18_KD') %>%
#   st_join(.,var_fishnet,join=st_equals)
# 
# # add ID_18 count
# var_fishnet <- 
#   count_net(ID_18,fishnet,"ID_18_count") %>%
#   as.data.frame(.)%>%
#   dplyr::select(ID_18_count,uniqueID)%>%
#   merge(as.data.frame(var_fishnet),. , by="uniqueID") %>%
#   st_sf()
# # add ID_18 kernel density
# var_fishnet <- 
#   kernelDensity(ID_18, fishnet, 'ID_18_KD') %>%
#   st_join(.,var_fishnet,join=st_equals)

# # add CE_18 count
# var_fishnet <-
#   count_net(CE_18,fishnet,"CE_18_count") %>%
#   as.data.frame(.)%>%
#   dplyr::select(CE_18_count,uniqueID)%>%
#   merge(as.data.frame(var_fishnet),. , by="uniqueID") %>%
#   st_sf()
# # add CE_18 kernel density
# var_fishnet <-
#   kernelDensity(CE_18, fishnet, 'CE_18_KD') %>%
#   st_join(.,var_fishnet,join=st_equals)

# # add ALM_18 count
# var_fishnet <- 
#   count_net(ALM_18,fishnet,"ALM_18_count") %>%
#   as.data.frame(.)%>%
#   dplyr::select(ALM_18_count,uniqueID)%>%
#   merge(as.data.frame(var_fishnet),. , by="uniqueID") %>%
#   st_sf()
# # add ALM_18 kernel density
# var_fishnet <- 
#   kernelDensity(ALM_18, fishnet, 'ALM_18_KD') %>%
#   st_join(.,var_fishnet,join=st_equals)


# knn method
var_fishnet$unsafe_18.nn3 <-nn_function(st_coordinates(st_centroid(var_fishnet)),st_coordinates(st_centroid(unsafe_18)),3) 
var_fishnet$ID_18.nn3 <-nn_function(st_coordinates(st_centroid(var_fishnet)),st_coordinates(st_centroid(ID_18)),3) 
var_fishnet$ALM_18.nn3 <-nn_function(st_coordinates(st_centroid(var_fishnet)),st_coordinates(st_centroid(ALM_18)),3) 
var_fishnet$CE_18.nn1 <-nn_function(st_coordinates(st_centroid(var_fishnet)),st_coordinates(st_centroid(CE_18)),1) 


##### 2018 311 requests(count/kernel density/knn---knn is better) #####
# load 311 data 
req <- read.csv("C:/Users/HP/Desktop/public_cases_fc.csv") 
req <- 
  req %>%
  mutate(Year=year(requested_datetime)) %>%
  dplyr::select(lon,lat,service_name,Year)%>%
  na.omit()%>%
  st_as_sf(.,coords=c("lon","lat"),crs=2272)

table(req$service_name)

hydrant_Down18 <- req %>%
  filter(service_name=="Hydrant Knocked Down (No Water)" & Year==2018)   
hydrant_req <- req %>%
  filter(service_name=="Hydrant Request" & Year==2018)  
smoke_dector18 <- req%>%
  filter(service_name=="Smoke Detector" & Year==2018)    
complainfire18 <- req %>%
  filter(.,service_name=="Complaints against Fire or EMS" & Year==2018)
FireRes_or_Com18 <- req %>%
  filter(service_name=="Fire Residential or Commercial" & Year==2018) 

# knn method
var_fishnet$hydrant_Down18.nn3 <-nn_function(st_coordinates(st_centroid(var_fishnet)),st_coordinates(st_centroid(hydrant_Down18)),3) 
var_fishnet$hydrant_req.nn3 <-nn_function(st_coordinates(st_centroid(var_fishnet)),st_coordinates(st_centroid(hydrant_req)),3) 
var_fishnet$smoke_dector18.nn3 <-nn_function(st_coordinates(st_centroid(var_fishnet)),st_coordinates(st_centroid(smoke_dector18)),3) 
var_fishnet$complainfire18.nn3 <-nn_function(st_coordinates(st_centroid(var_fishnet)),st_coordinates(st_centroid(complainfire18)),3) 
var_fishnet$FireRes_or_Com18.nn3 <-nn_function(st_coordinates(st_centroid(var_fishnet)),st_coordinates(st_centroid(FireRes_or_Com18)),3) 

##### 2018 census tracts ####
# load the 2018 census data
census <- get_acs(geography = "tract",
                  variables=c("B01003_001", "B19013_001", 
                              "B02001_002", "B01002_001","B06009_005"),
                  key="d72b594e4d0f9b9b34217cdea8a4bcbc60354e21",
                  state=42,
                  year=2018,
                  geometry=TRUE,
                  county=101,
                  output="wide") %>%
  rename(Total_Pop =  B01003_001E,
         Med_Inc = B19013_001E,
         Med_Age = B01002_001E,
         White_Pop = B02001_002E,
         Bachelor_Pop=B06009_005E) %>%
  dplyr::select(Total_Pop, Med_Inc, White_Pop,Bachelor_Pop, Med_Age,GEOID, geometry) %>%
  mutate(Percent_White = White_Pop / Total_Pop,Percent_Bachelor=Bachelor_Pop/ Total_Pop)%>%
  dplyr::select(- White_Pop,-Bachelor_Pop)

var_fishnet <- 
  st_join(var_fishnet,census%>% st_transform(2272),st_intersects,left=TRUE) 

var_fishnet[is.na(var_fishnet)] <- 0

##### correlation test#####
# correlation Matirix between predictors
Cor <- 
  var_fishnet %>% 
  st_set_geometry(NULL) %>%
  dplyr::select(-uniqueID, -cvID,-GEOID,-hydrant_count,-count) 

col <- colorRampPalette(c("#6699CC","white","#FFCCCC"))

corrplot(cor(Cor), method = "color", col = col(20),
         type = "upper", order = "alphabet", number.cex = .7,
         addCoef.col = "white", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         diag = TRUE)

# correlation between predictors and DV
correlation.long <-
  st_set_geometry(var_fishnet, NULL) %>%
  dplyr::select(-uniqueID, -cvID,-GEOID,-hydrant_count) %>%
  gather(Variable, Value, -count)

correlation.cor <-
  correlation.long %>%
  group_by(Variable) %>%
  summarise(correlation = cor(Value, count, use = "complete.obs"))

ggplot(correlation.long, aes(Value, count)) +
  geom_point(size = 0.2,colour = "#6699CC") +
  geom_text(data = correlation.cor, aes(label = paste("r =", round(correlation, 2))),
            x=-Inf, y=Inf, vjust = 1, hjust = -.1) +
  geom_smooth(method = "lm", se = FALSE, colour = "#333366") +
  facet_wrap(~Variable, ncol = 4, scales = "free") +
  labs(title = "Fire count as a function of potential factors")



########## Modeling & Validation ##########
##### cross validation by LOGO and K-fold #####
# load the neighborhood data
neighbor <- 
  st_read("C:/Users/HP/Desktop/MUSA800_practium/data/Neighborhoods_Philadelphia/Neighborhoods_Philadelphia.shp")%>%
  na.omit() %>%
  st_transform(crs=2272)

# create the final fishnet
final_fishnet <-
  st_centroid(var_fishnet) %>%
  st_join(., dplyr::select(neighbor, NAME)) %>%
  st_set_geometry(NULL) %>%
  left_join(dplyr::select(var_fishnet, geometry, uniqueID)) %>%
  st_sf() %>%
  na.omit()

# cross validation
str(final_fishnet)
reg.ss.vars <- c('fire.isSig','fire.isSig.dist','hydrant_KD',
                 'unsafe_18.nn3','ID_18.nn3','ALM_18.nn3','CE_18.nn1',
                 "Total_Pop","Med_Inc","Med_Age","Percent_White","Percent_Bachelor",
                 "hydrant_Down18.nn3","hydrant_req.nn3","smoke_dector18.nn3","complainfire18.nn3","FireRes_or_Com18.nn3")

crossValidate <- function(dataset, id, dependentVariable, indVariables) {
  
  allPredictions <- data.frame()
  cvID_list <- unique(dataset[[id]])
  
  for (i in cvID_list) {
    
    thisFold <- i
    cat("This hold out fold is", thisFold, "\n")
    
    fold.train <- filter(dataset, dataset[[id]] != thisFold) %>% as.data.frame() %>% 
      dplyr::select(id, geometry, indVariables, dependentVariable)
    fold.test  <- filter(dataset, dataset[[id]] == thisFold) %>% as.data.frame() %>% 
      dplyr::select(id, geometry, indVariables, dependentVariable)
    
    regression <-
      glm(count ~ ., family = "poisson", 
          data = fold.train %>% 
            dplyr::select(-geometry, -id))
    
    thisPrediction <- 
      mutate(fold.test, Prediction = predict(regression, fold.test, type = "response"))
    
    allPredictions <-
      rbind(allPredictions, thisPrediction)
    
  }
  return(st_sf(allPredictions))
}

# k fold validation
reg.ss.cv <- crossValidate(
  dataset = final_fishnet,
  id = "cvID",
  dependentVariable = "count",
  indVariables = reg.ss.vars) %>%
  dplyr::select(cvID = cvID, count, Prediction, geometry)

# spatial cross validation
reg.ss.spatialCV <- crossValidate(
  dataset = final_fishnet,
  id = "NAME",
  dependentVariable = "count",
  indVariables = reg.ss.vars) %>%
  dplyr::select(cvID = NAME, count, Prediction, geometry)

##### visualize the result #####
# calculate the error  
reg.summary <- 
  rbind(
    mutate(reg.ss.cv,        Error = Prediction - count,
           Regression = "Random k-fold CV: Spatial Structure"),
    
    mutate(reg.ss.spatialCV, Error = Prediction - count,
           Regression = "Spatial LOGO-CV: Spatial Structure")) %>%
  st_sf()

# visulize the predicted and observed fire count
st_set_geometry(reg.summary, NULL) %>%
  group_by(Regression) %>%
  mutate(fire_Decile = ntile(count, 10)) %>%
  group_by(Regression, fire_Decile) %>%
  summarise(meanObserved = mean(count, na.rm=T),
            meanPrediction = mean(Prediction, na.rm=T)) %>%
  gather(Variable, Value, -Regression, -fire_Decile) %>%          
  ggplot(aes(fire_Decile, Value, shape = Variable)) +
  geom_point(size = 2) + geom_path(aes(group =fire_Decile), colour = "black") +
  scale_shape_manual(values = c(2, 17)) +
  facet_wrap(~Regression) + xlim(0,10) +
  labs(title = "Predicted and observed fire count")

# check the mae and the sd_mae
st_set_geometry(reg.summary, NULL) %>%
  group_by(Regression) %>% 
  summarise(MAE = round(mean(abs(Prediction - count), na.rm = T),2),
            SD_MAE = round(sd(abs(Prediction - count), na.rm = T),2)) %>% 
  kable(caption = "MAE by regression") %>%
  kable_styling("striped", full_width = F) %>%
  row_spec(2, color = "black", background = "white")

mean(var_fishnet$count)
