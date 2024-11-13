## -------------------------------------------------
##           Join French detections - traps
##                        2020 
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

itineraires_2020 <- st_read("D:/MargSalas/Oso/Datos/Effort_raw/France/Itineraires_2020/Itineraire.2020.1.shp")%>%
  distinct(.keep_all = TRUE) %>%
  st_transform(crs = CRS("+proj=utm +zone=31 +datum=WGS84")) %>%
  rename("NOM"="nom") %>%
  dplyr::select(NOM) %>%
  mutate(id = row_number())

## -- B. Hair traps --

# We load all hair traps from 2022 (which contain the hairtraps from 2020 and 2021)
piege_poils <- st_read("D:/MargSalas/Oso/Datos/Effort_raw/France/2022/Appat_smola_2022.shp") %>%
  distinct(.keep_all = TRUE) %>%
  st_transform(crs = CRS("+proj=utm +zone=31 +datum=WGS84"))

# We keep the ones that fall within the itineraries of 2020
# Hair trap on an itinerary = at least 500m from the closest itinerary.
# With this threshold we remove the ones that belong to transects that weren't done this year

seuil <- 500

piege_poils_it_2020 <- dist_nearest(piege_poils,itineraires_2020) %>%
  filter(dist < seuil) %>%
  mutate("itineraire" = NA) # Variable trap = transecto al que pertenece la trampa de pelo

#mapview(itineraires_2020) + mapview(piege_poils, col.regions = "magenta", cex = 2) + mapview(piege_poils_it_2020, col.regions = "darkgreen", cex = 2)

## -- C. Camera traps --

piege_photos_2020 <- st_read("D:/MargSalas/Oso/Datos/Effort_raw/France/Cameras_2020/AppareilsAutos_2020.shp") %>%
  rename("NOM" = "Nom") %>%
  st_transform(crs = CRS("+proj=utm +zone=31 +datum=WGS84"))

## ---- Create SITES ----

# 1)) De las trampas de pelo que están en los transectos 
# miramos cuales tienen cámara en proximidad y asociamos a cada cámara a menos 
# de 500m de una trampa de pelo, la trampa de pelo más cercana
## ** Importante! Las coordenadas de las trampas de pelo servirán como puntos GPS de nuestros "sites".

piege_photos_2020 <- dist_nearest(piege_photos_2020, piege_poils_it_2020) # Variable trap = trampa de pelo más cercana a la camara (número de fila en piege_poils_it_2020)

#mapview(piege_poils_it_2020, col.regions = "magenta", cex = 3) + mapview(piege_photos_2020, cex = 3)

# On prendre les pieges photos qui sont à moins de 500 m d'un piège à poils qui se situe sur un itinéraire. 
piege_photos_it_2020 <- piege_photos_2020 %>%
  filter(dist < seuil) # piège photo sur un itineraire

piege_poils_it_2020 <- piege_poils_it_2020 %>%
  mutate("site" = "hair") # Site possède un piège à poils seul

#mapview(piege_poils_it_2020, col.regions = "magenta", cex = 3) + mapview(piege_photos_it_2020, cex = 3)

# 2))
# Si les pièges à poils des itinéraires sont associés à au moins un piège photos alors on modifie la covariable site pour les rattacher au piège à poil le plus proche. 
# Pour chaque piège photo sur un itinéraire on l'associé au piège à poil le plus proche, ce qui forme un site de type 3 (piège à poil + piège photo) 

for (i in 1 : length(piege_photos_it_2020$trap)){ 
  piege_poils_it_2020$site[piege_photos_it_2020$trap[i]] <- "both"
}

# On vérifie que les différents tableaux ont les mêmes variables
piege_poils_it_2020 <- piege_poils_it_2020 %>%
  dplyr::select(Nom_Iti, site)

# 3)) Como piege_poils_it_2020 incluye los sitios con y sin camara y no hay más que añadir, estos son todos los sites

sites_2020 <- piege_poils_it_2020

# On ajoute la variable effort qui indique si le site se situe ou non sur un itinéraire
sites_2020$effort <- "it" # Como no se añaden las cámaras que están lejos de un transecto, todo es dentro de transecto

# On les range par type de site et on leur attribue un numéro
sites_2020 <- sites_2020 %>%
  arrange(site) %>%
  mutate(trap_id = paste(row_number(),"f"))

## ---- Detections ----

# LOAD dataset that has been contrasted with dataset of Maelis (Prospind_2017-2019_Maelisv2.xlsx)
# ONLY for 2017-2019 (so this year 2020 is not contrasted with the french data)
# Whether the hair sample comes from a trap (Poils appât) or not (poils spontanées) is only identified
# for 2017-2019, so I will need to use all hair samples and assume they come from traps.

setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
os <- openxlsx::read.xlsx('Seguiment_Ossos_Pirineus_1996_2021_taula_final_2_cubLocations.xlsx') %>% 
  filter(Year %in% c(2020) & Confirmed_Individual != "Indetermined" & Country == "France") # Removed indetermined

os <- os %>% mutate(date = as_date(os$Date_register, format = "%d/%m/%Y"),
                    month = month(date)) %>%
  filter(month < 12, month > 4) %>% # 7 months form may to november
  filter(!is.na(x_long)) # Remove observations without coordinates


## Keep systematic data
dat2 <- os[which(os$Method %in% c("Sampling_station", "Transect") &
                   os$Obs_type %in% c("Photo","Photo/Video", "Hair", "Video")), ]

# We distinguish the 3 methods 
dat2 <- dat2 %>%
  mutate(method = "itineraire") # Itinerary 

dat2[which(dat2$Obs_type == "Photo" | dat2$Obs_type == "Photo/Video" | dat2$Obs_type == "Video"),"method"] <- "piege photos" # Camera trap
dat2[which(dat2$Obs_type == "Hair"),"method"] <- "piege poils" # Hair trap
