## -------------------------------------------------
##           Join French detections - traps           
## ------------------------------------------------- 

rm(list = ls())

library(tidyverse)
library(sf)
library(rgdal)
library(mapview)
library(lubridate)

## ---- Load function ----

dist_nearest <- function(point, piege){
  point <- point %>%
    mutate(trap = st_nearest_feature(point,piege))
  
  A <- rep(NA, nrow(point))
  
  for (i in 1 : nrow(point)){
    A[i] <- st_distance(point[i,], piege[point$trap[i],])}
  
  point <- point %>%
    mutate(dist = A)
  
  return(point)
}

## ---- Load background map ----

map <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/clip_pyros2.shp") %>%
  st_transform(map, crs = 32631) # WGS84 31N

## ---- Load traps ----

## -- A. Transects --

itineraires <- st_read("D:/MargSalas/Oso/Datos/Effort_raw/France/Itineraires_2020/Itineraire.2020.1.shp")%>%
  distinct(.keep_all = TRUE) %>%
  rename("NOM"="nom") %>%
  dplyr::select(NOM) %>%
  mutate(id = row_number())

## -- B. Hair traps --

# We distinguish hair trap that are on an itinerary (because they are more frequently visited) 
# from hair trap outside an itinerary. 
# Hair trap on an itinerary = at least 500m from the closest itinerary.

## DON'T HAVE THEM YET. DO IT USE ANY

## -- C. Camera traps --

# Make one camera traps file by year, i.e. remove camera traps that aren't close to a hair trap. 
# Only camera traps can't allow the identification of individuals so we remove them from the study. 

piege_photos_2020 <- st_read("D:/MargSalas/Oso/Datos/Effort_raw/France/Cameras_2020/AppareilsAutos_2020.shp") %>%
  rename("NOM" = "Nom") %>%
  st_drop_geometry()
