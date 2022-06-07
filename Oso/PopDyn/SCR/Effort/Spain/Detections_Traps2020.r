## -------------------------------------------------
##           Join spanish detections - traps           
## ------------------------------------------------- 

rm(list = ls())

library(tidyverse)
library(sf)
library(rgdal)

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

## Estas son SÓLO las trampas de Catalonia, falta ARAN

trap_cat <- readxl::read_xlsx("D:/MargSalas/Oso/Datos/Effort_raw/Revisions trampes 2020 V2.xlsx", sheet = 1) %>% 
  janitor::clean_names() %>%
  filter(!is.na(coord_x)) %>%
  select(codi_tr,tipus_tr,coord_x,coord_y) %>%
  mutate(trap_id = paste(row_number(), "c")) %>% 
  mutate(site = "hair") %>% # poils seul
  st_as_sf(coords = c("coord_x","coord_y"), 
           crs = CRS("+proj=utm +zone=31 +datum=WGS84"))  %>%
  mutate(trap = row_number())

trap_cat$site <- ifelse(trap_cat$tipus_tr == "Mixte", "both","hair")

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
os <- read.csv("Seguiment_Ossos_Pirineus_1996_2020_taula_final.csv", header = TRUE, row.names = NULL)
os <- os[which(os$Year == 2020 & os$Probable_Individual != "Indetermined" & os$Country == "Spain"),-1]

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
