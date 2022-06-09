## -------------------------------------------------
##           Join spanish detections - traps           
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

## Trampas Catalonia

trap_onlycat <- readxl::read_xlsx("D:/MargSalas/Oso/Datos/Effort_raw/Spain/Revisions trampes 2020 V2.xlsx", sheet = 1) %>% 
  janitor::clean_names() %>%
  filter(!is.na(coord_x)) %>%
  select(codi_tr,tipus_tr,coord_x,coord_y) %>%
  mutate(site = "hair") %>% # poils seul
  st_as_sf(coords = c("coord_x","coord_y"), 
           crs = CRS("+proj=utm +zone=31 +datum=WGS84"))

  trap_onlycat$site <- ifelse(trap_onlycat$tipus_tr == "Mixte", "both","hair")

## Trampas Aran
#  De momento cojo las de Aran 2021, porque en teoría es un sistema
## nuevo de cuadrículas en el que las trampas se repiten y son las mismas para 2020 y 2021

trap_aran <- readxl::read_xlsx("D:/MargSalas/Oso/Datos/Effort_raw/Spain/Trampes_2021_SFF_CAR_CGA.xlsx") %>% 
  janitor::clean_names() %>%
  filter(comarca == "Vall d'Aran") %>%
  select(codi_tr,tipus_tr,coord_x,coord_y) %>%
  #mutate(trap_id = paste(row_number(), "c")) %>% 
  mutate(site = "hair") %>% # poils seul
  st_as_sf(coords = c("coord_x","coord_y"), 
           crs = CRS("+proj=utm +zone=31 +datum=WGS84")) 
#mutate(trap = row_number())

trap_aran$site <- ifelse(trap_aran$tipus_tr == "Mixte", "both","hair")

##ATENCIÓN!! AQUÍ FALTA LA INFORMACIÓN DE NAVARRA
# No la añado porque faltan info en los datos, pero hay que añadirlo
# en la nueva base de datos Seguiment 96-21


# Join Catalonia and Aran
trap_cat <- trap_onlycat %>%
  rbind(trap_aran) %>%
  mutate(trap_id = paste(row_number(), "c")) %>%
  mutate(trap = row_number())

## Creamos un object para cada tipo de trampa

trap_cat_pels <- trap_cat %>%
  filter(site == "hair") 

trap_cat_foto <- trap_cat %>%
  filter(site == "both") 

st_bbox(trap_cat_pels)
#plot(st_geometry(map),xlim = st_bbox(trap_cat_pels)[c(1,3)], ylim = st_bbox(trap_cat_pels)[c(2, 4)])
plot(st_geometry(map))
plot(st_geometry(trap_cat_pels), add = TRUE)

#mapview(trap_cat_pels) + mapview(trap_cat_foto, col.regions = "red")

## ---- Load detections ----

setwd("D:/MargSalas/Oso/Datos/Tablas_finales")
os <- read.csv("Seguiment_Ossos_Pirineus_1996_2020_taula_final.csv", header = TRUE, row.names = NULL) %>%
  filter(Year == 2020 & Probable_Individual != "Indetermined" & Country == "Spain") %>%
  select(-X.1) %>%
  mutate(date = as_date(os$Date_register, format = "%d/%m/%Y"),
         month = month(date))
os$month[157] <- 4 # Correct mannually because date was not exact and we only had month
os$month[158] <- 6
os <- os %>%
  filter(month < 12, month > 4)  # 7 months form mai to november 

## Keep systematic data
dat_cat_syst <- os[which(os$Method %in% c("Sampling_station", "Transect") &
                 os$Obs_type %in% c("Photo","Photo/Video", "Hair")), ]

pts_foto_2020 <- dat_cat_syst %>% 
  st_as_sf(coords = c("x_long","y_lat"), 
           crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) %>%
  st_transform(crs = 32631) %>%
  filter(Obs_type %in% c("Photo","Photo/Video"))

pts_pels_2020 <- dat_cat_syst %>% 
  st_as_sf(coords = c("x_long","y_lat"), 
           crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) %>%
  st_transform(crs = 32631) %>%
  filter(Obs_type %in% c("Hair"))


## ---- Association detection - nearest trap ----

pts_foto_2020 <- dist_nearest(pts_foto_2020, trap_cat_foto) %>%
  #st_drop_geometry() %>%
  left_join(trap_cat_foto %>% st_drop_geometry(), by = "trap") %>%
  dplyr::select(Year, Date_register, Probable_Individual, Sex, Method, Obs_type, trap_id, dist)

pts_pels_2020 <- dist_nearest(pts_pels_2020, trap_cat) %>% ## Here I join to trap_cat that contains hair and mixed (both)
  #st_drop_geometry() %>%
  left_join(trap_cat %>% st_drop_geometry(), by = "trap") %>%
  dplyr::select(Year, Date_register, Probable_Individual, Sex, Method, Obs_type, trap_id, dist)

mapview(trap_cat_pels) + mapview(trap_cat_foto, col.regions = "green") + 
  mapview(pts_foto_2020, col.regions = "darkgreen", cex = 2) + 
  mapview(pts_pels_2020, col.regions = "magenta", cex = 2)

## With a threshold of 500 m there is only one that is saved

# Cases that could be worth checking with santi
# There is an observation but there is no trap associated, it is VERY far from any trap and the coordinates are good
## Esmolet 15/06, Pepito 12/08, New 18-03 05/07 (sist auto?there is even a camera?)
# Should we include a trap in this ones??



