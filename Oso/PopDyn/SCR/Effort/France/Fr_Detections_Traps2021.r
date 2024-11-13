## -------------------------------------------------
##           Join French detections - traps
##                        2021 
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

itineraires_2021 <- st_read("D:/MargSalas/Oso/Datos/Effort_raw/France/Itineraires_2021/Itineraire2021.shp")%>%
  distinct(.keep_all = TRUE) %>%
  st_transform(crs = CRS("+proj=utm +zone=31 +datum=WGS84")) %>%
  rename("NOM"="nom") %>%
  dplyr::select(NOM) 

Cagateille <- st_read("D:/MargSalas/Oso/Datos/Effort_raw/France/Itineraires_2021/Iti_Shp_a_ajouter/Iti_Shp_a_ajouter/Cagateille_2020_2021.shp") %>%
  mutate(NOM = "Cagateille") %>%
  distinct(.keep_all = TRUE) %>%
  st_zm(drop = TRUE, what = "ZM") %>% # To remove the z dimension, so that we can join it with the rest of transects
  st_transform(crs = CRS("+proj=utm +zone=31 +datum=WGS84"))%>%
  dplyr::select(NOM) 

Luzenac <- st_read("D:/MargSalas/Oso/Datos/Effort_raw/France/Itineraires_2021/Iti_Shp_a_ajouter/Iti_Shp_a_ajouter/Luzenac_2021.shp")%>%
  mutate(NOM = "Luzenac") %>%
  distinct(.keep_all = TRUE) %>%
  st_zm(drop = TRUE, what = "ZM") %>%
  st_transform(crs = CRS("+proj=utm +zone=31 +datum=WGS84"))%>%
  dplyr::select(NOM)  

Savignac <- st_read("D:/MargSalas/Oso/Datos/Effort_raw/France/Itineraires_2021/Iti_Shp_a_ajouter/Iti_Shp_a_ajouter/Savignac_2021.shp")%>%
  mutate(NOM = "Savignac") %>%
  distinct(.keep_all = TRUE) %>%
  st_zm(drop = TRUE, what = "ZM") %>%
  st_transform(crs = CRS("+proj=utm +zone=31 +datum=WGS84"))%>%
  dplyr::select(NOM)  

Soussoueou <- st_read("D:/MargSalas/Oso/Datos/Effort_raw/France/Itineraires_2021/Iti_Shp_a_ajouter/Iti_Shp_a_ajouter/Soussoueou_2020_2021.shp")%>%
  mutate(NOM = "Soussoueou") %>%
  distinct(.keep_all = TRUE) %>%
  st_zm(drop = TRUE, what = "ZM") %>%
  st_transform(crs = CRS("+proj=utm +zone=31 +datum=WGS84"))%>%
  dplyr::select(NOM)  

itineraires_2021 <- rbind(itineraires_2021, Cagateille, Luzenac, Savignac, Soussoueou) %>%
  mutate(id = row_number())

mapview(itineraires_2021)

## -- B. Hair traps --

# We load all hair traps from 2022 (which contain the hairtraps from 2020 and 2021)
piege_poils <- st_read("D:/MargSalas/Oso/Datos/Effort_raw/France/2022/Appat_smola_2022.shp") %>%
  distinct(.keep_all = TRUE) %>%
  st_transform(crs = CRS("+proj=utm +zone=31 +datum=WGS84"))

# We keep the ones that fall within the itineraries of 2021
# Hair trap on an itinerary = at least 500m from the closest itinerary.
# With this threshold we remove the ones that belong to transects that weren't done this year

seuil <- 500

piege_poils_it_2021 <- dist_nearest(piege_poils,itineraires_2021) %>%
  filter(dist < seuil) %>%
  mutate("itineraire" = NA) # Variable trap = transecto al que pertenece la trampa de pelo

#mapview(itineraires_2021) + mapview(piege_poils, col.regions = "magenta", cex = 2) + mapview(piege_poils_it_2021, col.regions = "darkgreen", cex = 2)
# There are two transects without hair traps, but is an error (those transects weren't active in 2021). We just ignore it

## -- C. Camera traps --

piege_photos_2021 <- st_read("D:/MargSalas/Oso/Datos/Effort_raw/France/Cameras_2021/App.2021.shp") %>%
  rename("NOM" = "Nom") %>%
  st_transform(crs = CRS("+proj=utm +zone=31 +datum=WGS84"))

## ---- Create SITES ----

# 1)) De las trampas de pelo que están en los transectos 
# miramos cuales tienen cámara en proximidad y asociamos a cada cámara a menos 
# de 500m de una trampa de pelo, la trampa de pelo más cercana
## ** Importante! Las coordenadas de las trampas de pelo servirán como puntos GPS de nuestros "sites".

piege_photos_2021 <- dist_nearest(piege_photos_2021, piege_poils_it_2021) # Variable trap = trampa de pelo más cercana a la camara (número de fila en piege_poils_it_2020)

mapview(piege_poils_it_2021, col.regions = "magenta", cex = 3) + mapview(piege_photos_2021, cex = 3)

# On prendre les pieges photos qui sont à moins de 500 m d'un piège à poils qui se situe sur un itinéraire. 
piege_photos_it_2021 <- piege_photos_2021 %>%
  filter(dist < seuil) # piège photo sur un itineraire

piege_poils_it_2021 <- piege_poils_it_2021 %>%
  mutate("site" = "hair") # Site possède un piège à poils seul

mapview(piege_poils_it_2021, col.regions = "magenta", cex = 3) + mapview(piege_photos_it_2021, cex = 3)

# 2))
# Si les pièges à poils des itinéraires sont associés à au moins un piège photos alors on modifie la covariable site pour les rattacher au piège à poil le plus proche. 
# Pour chaque piège photo sur un itinéraire on l'associé au piège à poil le plus proche, ce qui forme un site de type 3 (piège à poil + piège photo) 

for (i in 1 : length(piege_photos_it_2021$trap)){ 
  piege_poils_it_2021$site[piege_photos_it_2021$trap[i]] <- "both"
}

# On vérifie que les différents tableaux ont les mêmes variables
piege_poils_it_2021 <- piege_poils_it_2021 %>%
  dplyr::select(Nom_Iti, site)

# 3)) Como piege_poils_it_2021 incluye los sitios con y sin camara y no hay más que añadir, estos son todos los sites

sites_2021 <- piege_poils_it_2021

# On ajoute la variable effort qui indique si le site se situe ou non sur un itinéraire
sites_2021$effort <- "it" # Como no se añaden las cámaras que están lejos de un transecto, todo es dentro de transecto

# On les range par type de site et on leur attribue un numéro
sites_2021 <- sites_2021 %>%
  arrange(site) %>%
  mutate(trap_id = paste(row_number(),"f"))

## ---- Detections ----

