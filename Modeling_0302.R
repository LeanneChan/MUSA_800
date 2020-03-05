########## Preperation ##########
##### standarlize function #####
setwd("C:/Users/HP/Desktop/MUSA800_practium")

colors=c("#61c0bf", "#bbded6", "#fae3d9","#ffb6b9","#fc7978")

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
    panel.grid.major = element_blank(),
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

##### Load fire data to the fishnet#####
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
  fire.sf %>% 
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

##### 2015-2019 L& I violations (count/kernel density/knn---knn is better) #####
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
unsafe <- unsafe_structures[unsafe_structures$Year>=2015 && unsafe_structures$Year<=2019,]
# Imminently Dangerous
ID_structures<-readData("https://phl.carto.com/api/v2/sql?q=SELECT%20*%20FROM%20li_violations%20WHERE%20(violationtype%20=%20%27PM15-110.1%27)%20AND%20(violationdate%20%3E=%20%272015-01-01%27)")
ID <- ID_structures[ID_structures$Year>=2015 && ID_structures$Year<=2019]
# Lacking Proper Fire Extinguisher for Commercial Cooking
CE_structures<-readData("https://phl.carto.com/api/v2/sql?q=SELECT%20*%20FROM%20li_violations%20WHERE%20(violationtype%20=%20%27FC-904.11/5%27)%20AND%20(violationdate%20%3E=%20%272015-01-01%27)")
CE <- CE_structures[CE_structures$Year>=2015 && CE_structures$Year<=2019,]
# Lacking smoke Alarm
ALM_structures<-readData("https://phl.carto.com/api/v2/sql?q=SELECT%20*%20FROM%20li_violations%20WHERE%20(violationtype%20=%20%27FC-907.3/20%27)%20AND%20(violationdate%20%3E=%20%272015-01-01%27)")
ALM <- ALM_structures[ALM_structures$Year>=2015 && ALM_structures$Year<=2019,]

# knn method
var_fishnet$unsafe.nn3 <-nn_function(st_coordinates(st_centroid(var_fishnet)),st_coordinates(st_centroid(unsafe)),3) 
var_fishnet$ID.nn3 <-nn_function(st_coordinates(st_centroid(var_fishnet)),st_coordinates(st_centroid(ID)),3) 
var_fishnet$ALM.nn3 <-nn_function(st_coordinates(st_centroid(var_fishnet)),st_coordinates(st_centroid(ALM)),3) 
var_fishnet$CE.nn3 <-nn_function(st_coordinates(st_centroid(var_fishnet)),st_coordinates(st_centroid(CE)),3) 

# # add unsafe count
# var_fishnet <- 
#   count_net(unsafe,fishnet,"unsafe_18_count") %>%
#   as.data.frame(.)%>%
#   dplyr::select(unsafe_count,uniqueID)%>%
#   merge(as.data.frame(var_fishnet),. , by="uniqueID") %>%
#   st_sf()
# # add unsafe kernel density
# var_fishnet <- 
#   kernelDensity(unsafe, fishnet, 'unsafe_18_KD') %>%
#   st_join(.,var_fishnet,join=st_equals)
# 
# # add ID_18 count
# var_fishnet <- 
#   count_net(ID,fishnet,"ID_count") %>%
#   as.data.frame(.)%>%
#   dplyr::select(ID_count,uniqueID)%>%
#   merge(as.data.frame(var_fishnet),. , by="uniqueID") %>%
#   st_sf()
# # add ID_18 kernel density
# var_fishnet <- 
#   kernelDensity(ID, fishnet, 'ID_KD') %>%
#   st_join(.,var_fishnet,join=st_equals)

# # add CE count
# var_fishnet <-
#   count_net(CE,fishnet,"CE_count") %>%
#   as.data.frame(.)%>%
#   dplyr::select(CE_count,uniqueID)%>%
#   merge(as.data.frame(var_fishnet),. , by="uniqueID") %>%
#   st_sf()
# # add CE kernel density
# var_fishnet <-
#   kernelDensity(CE, fishnet, 'CE_KD') %>%
#   st_join(.,var_fishnet,join=st_equals)

# # add ALM count
# var_fishnet <- 
#   count_net(ALM,fishnet,"ALM_count") %>%
#   as.data.frame(.)%>%
#   dplyr::select(ALM_count,uniqueID)%>%
#   merge(as.data.frame(var_fishnet),. , by="uniqueID") %>%
#   st_sf()
# # add ALM kernel density
# var_fishnet <- 
#   kernelDensity(ALM, fishnet, 'ALM_KD') %>%
#   st_join(.,var_fishnet,join=st_equals)


##### 2015-2019 311 requests(count/kernel density/knn---knn is better) #####
# load 311 data 
req <- read.csv("C:/Users/HP/Desktop/MUSA800_practium/public_cases_fc.csv") 
req <- 
  req %>%
  mutate(Year=year(requested_datetime)) %>%
  dplyr::select(lon,lat,service_name,Year)%>%
  na.omit()%>%
  st_as_sf(.,coords=c("lon","lat"), crs=4326)%>%
  st_transform(2272)
table(req$service_name)

hydrant_Down<- req %>%
  filter(.,service_name=="Hydrant Knocked Down (No Water)" & Year>=2015 & Year<=2019)   
hydrant_req<- req %>%
  filter(.,service_name=="Hydrant Request" & Year>=2015 & Year<=2019)  
smoke_dector<- req%>%
  filter(.,service_name=="Smoke Detector" & Year>=2015 & Year<=2019)    
complainfire<- req %>%
  filter(.,service_name=="Complaints against Fire or EMS" & Year>=2015 & Year<=2019)
FireRes_or_Com<- req %>%
  filter(.,service_name=="Fire Residential or Commercial" & Year>=2015 & Year<=2019) 

# knn method
var_fishnet$hydrant_Down.nn3 <-nn_function(st_coordinates(st_centroid(var_fishnet)),st_coordinates(st_centroid(hydrant_Down)),3) 
var_fishnet$hydrant_req.nn3 <-nn_function(st_coordinates(st_centroid(var_fishnet)),st_coordinates(st_centroid(hydrant_req)),3) 
var_fishnet$smoke_dector.nn3 <-nn_function(st_coordinates(st_centroid(var_fishnet)),st_coordinates(st_centroid(smoke_dector)),3) 
var_fishnet$complainfire.nn3 <-nn_function(st_coordinates(st_centroid(var_fishnet)),st_coordinates(st_centroid(complainfire)),3) 
var_fishnet$FireRes_or_Com.nn3 <-nn_function(st_coordinates(st_centroid(var_fishnet)),st_coordinates(st_centroid(FireRes_or_Com)),3) 

##### census tracts ####
# load the census data
census <- get_acs(geography = "tract",
                  variables=c("B01003_001", "B19013_001", 
                              "B02001_002", "B01002_001","B06009_005"),
                  key="d72b594e4d0f9b9b34217cdea8a4bcbc60354e21",
                  state=42,
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

col <- colorRampPalette(c("#61c0bf", "#bbded6", "#f8f8f8","#ffb6b9","#fc7978"))
 
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
  geom_point(size = 0.2,colour = "#5eb7b7") +
  geom_text(data = correlation.cor, aes(label = paste("r =", round(correlation, 2))),
            x=-Inf, y=Inf, vjust = 1, hjust = -.1) +
  geom_smooth(method = "lm", se = FALSE, colour = "#fc7978") +
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
                 'ALM.nn3','CE.nn3',
                 "Total_Pop","Med_Inc","Med_Age","Percent_White","Percent_Bachelor",
                 "hydrant_Down.nn3","hydrant_req.nn3","complainfire.nn3","FireRes_or_Com.nn3")
#'unsafe.nn3',"smoke_dector.nn3",'ID.nn3'
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

# map the error 
colors1=c("#5eb7b7", "#96d1c7", "#f8f8f8","#ffafb0","#fc7978")
ggplot() +
  geom_sf(aes(fill = Error),data=reg.summary) +
  facet_wrap(~Regression) +
  #scale_fill_viridis(begin =1, end =0.2,limits=c(-70,70),option='inferno') +
  scale_fill_gradientn(colors = colors1,limits=c(-30,30)) +
  labs(title = "Fire Count Errors by Regression") +
  mapTheme()

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
  summarise(MAPE = round(mean(abs(Prediction - count), na.rm = T),2),
            SD_MAPE = round(sd(abs(Prediction - count), na.rm = T),2),
            MEAN_firecount=5.98) %>% 
  kable(caption = "MAE by regression") %>%
  kable_styling("striped", full_width = F) %>%
  row_spec(2, color = "black", background = "white")

mean(var_fishnet$count)

##### Compare with the traditional kernel density #####
# Compute kernel density
fire_ppp <- as.ppp(st_coordinates(fire.sf), W = st_bbox(final_fishnet))
fire_KD <- spatstat::density.ppp(fire_ppp, 1000)
# Convert kernel density to grid cells taking the mean
fire_KDE_sf <- as.data.frame(fire_KD) %>%
  st_as_sf(coords = c("x", "y"), crs = st_crs(final_fishnet)) %>%
  aggregate(., final_fishnet, mean) %>%
  #Mutate the Risk_Category field as defined below.
  mutate(label = "Kernel Density",
         Risk_Category = ntile(value, 100),
         Risk_Category = case_when(
           Risk_Category >= 90 ~ "90% to 100%",
           Risk_Category >= 70 & Risk_Category <= 89 ~ "70% to 89%",
           Risk_Category >= 50 & Risk_Category <= 69 ~ "50% to 69%",
           Risk_Category >= 30 & Risk_Category <= 49 ~ "30% to 49%",
           Risk_Category >= 1 & Risk_Category <= 29 ~ "1% to 29%")) %>%
  # Bind to a layer where test set crime counts are spatially joined to the fisnnet.
  bind_cols(
    aggregate(
      dplyr::select(fire.sf) %>% mutate(Count = 1), ., length) %>%
      mutate(Count = replace_na(Count, 0))) %>%
  #Select the fields we need
  dplyr::select(label, Risk_Category, Count)
head(fire_KDE_sf)

## risk predictors
#Spatial LOGO-CV: Spatial Structure
fire_risk_sf <-
  filter(reg.summary, Regression == "Random k-fold CV: Spatial Structure") %>%
  mutate(label = "Risk Predictions",
         Risk_Category = ntile(Prediction, 100),
         Risk_Category = case_when(
           Risk_Category >= 90 ~ "90% to 100%",
           Risk_Category >= 70 & Risk_Category <= 89 ~ "70% to 89%",
           Risk_Category >= 50 & Risk_Category <= 69 ~ "50% to 69%",
           Risk_Category >= 30 & Risk_Category <= 49 ~ "30% to 49%",
           Risk_Category >= 1 & Risk_Category <= 29 ~ "1% to 29%")) %>%
  bind_cols(
    aggregate(
      dplyr::select(fire.sf) %>% mutate(Count = 1), ., length) %>%
      mutate(Count = replace_na(Count, 0))) %>%
  dplyr::select(label,Risk_Category,Count)

## The map comparing kernel density to risk predictions
rbind(fire_KDE_sf, fire_risk_sf) %>%
  gather(Variable, Value, -label, -Risk_Category, -geometry) %>%
  ggplot() +
  geom_sf(aes(fill = Risk_Category), colour = NA) +
  geom_sf(data = sample_n(fire.sf, 1500), size = .1, colour = "black") +
  facet_wrap(~label, ) +
  #scale_fill_viridis(discrete = TRUE,option='inferno')+
  scale_fill_manual(values = colors) +
  labs(title="Comparison of Kernel Density and Risk Predictions",
       subtitle="Relative to test set points (in black)") +
  mapTheme()

## bar plot
rbind(fire_KDE_sf, fire_risk_sf) %>%
  st_set_geometry(NULL) %>%
  gather(Variable, Value, -label, -Risk_Category) %>%
  group_by(label, Risk_Category) %>%
  summarise(count = sum(Value)) %>%
  ungroup() %>%
  group_by(label) %>%
  mutate(Rate_of_test_set_fires = count / sum(count)) %>%
  ggplot(aes(Risk_Category,Rate_of_test_set_fires)) +
  geom_bar(aes(fill=label),position="dodge", stat="identity") +
  labs(title="Comparison of Kernel Density and Risk Predictions",
       subtitle="rate of test set fires") +
  scale_fill_manual(values = c("#61c0bf", "#fc7978"))
