## -------------------------------------------------
##      Join detections - traps France and Spain    
##                      2020-2021
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

## ---- FRANCE ----
## ---- 1. Load and sort out transects and traps ----

#### a) Itineraries ####
# Make one file by year with only itineraries that have been made. 

# 2020

itineraires_2020 <- st_read("D:/MargSalas/Oso/Datos/Effort_raw/France/Itineraires_2020/Itineraire.2020.1.shp")%>%
  distinct(.keep_all = TRUE) %>%
  st_transform(crs = CRS("+proj=utm +zone=31 +datum=WGS84")) %>%
  rename("NOM"="nom") %>%
  dplyr::select(NOM) %>%
  mutate(id = row_number())

# 2021

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

#### b) Hair traps ####
# Make one hair traps file by year. 

# We load all hair traps from 2022 (which contain the hairtraps from 2020 and 2021)
piege_poils <- st_read("D:/MargSalas/Oso/Datos/Effort_raw/France/2022/Appat_smola_2022.shp") %>%
  distinct(.keep_all = TRUE) %>%
  st_transform(crs = CRS("+proj=utm +zone=31 +datum=WGS84"))

# We keep the ones that fall within the itineraries of 2020 / 2021
# Hair trap on an itinerary = at least 500m from the closest itinerary.
# With this threshold we remove the ones that belong to transects that weren't done this year

seuil <- 500

piege_poils_it_2020 <- dist_nearest(piege_poils,itineraires_2020) %>%
  filter(dist < seuil) %>%
  mutate("itineraire" = NA) # Variable trap = transecto al que pertenece la trampa de pelo

#mapview(itineraires_2020) + mapview(piege_poils, col.regions = "magenta", cex = 2) + mapview(piege_poils_it_2020, col.regions = "darkgreen", cex = 2)

piege_poils_it_2021 <- dist_nearest(piege_poils,itineraires_2021) %>%
  filter(dist < seuil) %>%
  mutate("itineraire" = NA) # Variable trap = transecto al que pertenece la trampa de pelo

#mapview(itineraires_2021) + mapview(piege_poils, col.regions = "magenta", cex = 2) + mapview(piege_poils_it_2021, col.regions = "darkgreen", cex = 2)
# There are two transects without hair traps, but is an error (those transects weren't active in 2021). We just ignore it

#### c) Camera traps ####
# Make one camera traps file by year

piege_photos_2020 <- st_read("D:/MargSalas/Oso/Datos/Effort_raw/France/Cameras_2020/AppareilsAutos_2020.shp") %>%
  rename("NOM" = "Nom") %>%
  st_transform(crs = CRS("+proj=utm +zone=31 +datum=WGS84"))

piege_photos_2021 <- st_read("D:/MargSalas/Oso/Datos/Effort_raw/France/Cameras_2021/App.2021.shp") %>%
  rename("NOM" = "Nom") %>%
  st_transform(crs = CRS("+proj=utm +zone=31 +datum=WGS84"))

## ---- 2. Create SITES  ----

# 1)) De las trampas de pelo que están en los transectos 
# miramos cuales tienen cámara en proximidad y asociamos a cada cámara a menos 
# de 500m de una trampa de pelo, la trampa de pelo más cercana
## ** Importante! Las coordenadas de las trampas de pelo servirán como puntos GPS de nuestros "sites".

# 2020

piege_photos_2020 <- dist_nearest(piege_photos_2020, piege_poils_it_2020) # Variable trap = trampa de pelo más cercana a la camara (número de fila en piege_poils_it_2020)

#mapview(piege_poils_it_2020, col.regions = "magenta", cex = 3) + mapview(piege_photos_2020, cex = 3)

# On prendre les pieges photos qui sont à moins de 500 m d'un piège à poils qui se situe sur un itinéraire. 
piege_photos_it_2020 <- piege_photos_2020 %>%
  filter(dist < seuil) # piège photo sur un itineraire

piege_poils_it_2020 <- piege_poils_it_2020 %>%
  mutate("site" = "hair") # Site possède un piège à poils seul

#mapview(piege_poils_it_2020, col.regions = "magenta", cex = 3) + mapview(piege_photos_it_2020, cex = 3)

# 2021 

piege_photos_2021 <- dist_nearest(piege_photos_2021, piege_poils_it_2021) # Variable trap = trampa de pelo más cercana a la camara (número de fila en piege_poils_it_2020)

#mapview(piege_poils_it_2021, col.regions = "magenta", cex = 3) + mapview(piege_photos_2021, cex = 3)

# On prendre les pieges photos qui sont à moins de 500 m d'un piège à poils qui se situe sur un itinéraire. 
piege_photos_it_2021 <- piege_photos_2021 %>%
  filter(dist < seuil) # piège photo sur un itineraire

piege_poils_it_2021 <- piege_poils_it_2021 %>%
  mutate("site" = "hair") # Site possède un piège à poils seul

#mapview(piege_poils_it_2021, col.regions = "magenta", cex = 3) + mapview(piege_photos_it_2021, cex = 3)

# 2))
# Si les pièges à poils des itinéraires sont associés à au moins un piège photos alors on modifie la covariable site pour les rattacher au piège à poil le plus proche. 
# Pour chaque piège photo sur un itinéraire on l'associé au piège à poil le plus proche, ce qui forme un site de type 3 (piège à poil + piège photo) 

# 2020

for (i in 1 : length(piege_photos_it_2020$trap)){ 
  piege_poils_it_2020$site[piege_photos_it_2020$trap[i]] <- "both"
}

# On vérifie que les différents tableaux ont les mêmes variables
piege_poils_it_2020 <- piege_poils_it_2020 %>%
  dplyr::select(Nom_Iti, site)

# 2021

for (i in 1 : length(piege_photos_it_2021$trap)){ 
  piege_poils_it_2021$site[piege_photos_it_2021$trap[i]] <- "both"
}

# On vérifie que les différents tableaux ont les mêmes variables
piege_poils_it_2021 <- piege_poils_it_2021 %>%
  dplyr::select(Nom_Iti, site)

# 3)) Como piege_poils_it_202X incluye los sitios con y sin camara y no hay más que añadir, estos son todos los sites

sites_2020 <- piege_poils_it_2020
sites_2021 <- piege_poils_it_2021

# On ajoute la variable effort qui indique si le site se situe ou non sur un itinéraire
sites_2020$effort <- "it" # Como no se añaden las cámaras que están lejos de un transecto, todo es dentro de transecto
sites_2021$effort <- "it" 

# On les range par type de site et on leur attribue un numéro
sites_2020 <- sites_2020 %>%
  arrange(site) %>%
  mutate(trap_id = paste(row_number(),"f"))

sites_2021 <- sites_2021 %>%
  arrange(site) %>%
  mutate(trap_id = paste(row_number(),"f"))

# Transform to make crs fit with detections

sites_2020 <- sites_2020 %>%
  st_transform(crs = 4326)

sites_2021 <- sites_2021 %>%
  st_transform(crs = 4326)

## ---- 3. Detections ----

# LOAD dataset that has been contrasted with dataset of Maelis (Prospind_2017-2019_Maelisv2.xlsx)
# ONLY for 2017-2019 (so this years 2020 and 2021 is not contrasted with the french data)
# Whether the hair sample comes from a trap (Poils appât) or not (poils spontanées) is only identified
# for 2017-2019, so I will need to use all hair samples and assume they come from traps.

setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
os <- openxlsx::read.xlsx('Seguiment_Ossos_Pirineus_1996_2021_taula_final_2_cubLocations.xlsx') %>% 
  filter(Year %in% c(2020,2021) & Confirmed_Individual != "Indetermined" & Country == "France") # Removed indetermined

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

#### 3.1. Separate the detections of each type  ####

# Détections issues de pièges photos (photographie automatiques, vidéo automatique)
pts_photo <- dat2 %>%
  filter(method =="piege photos") %>%
  st_as_sf(coords = c("x_long", "y_lat"),
           crs = 4326)

pts_photo_2020 <- pts_photo %>%
  filter(Year == 2020) 

pts_photo_2021 <- pts_photo %>%
  filter(Year == 2021) 

# Détections issues de pièges à poils 
pts_poils <- dat2 %>%
  filter(method =="piege poils") %>%
  st_as_sf(coords = c("x_long", "y_lat"),
           crs = 4326)

pts_poils_2020 <-  pts_poils %>%
  filter(Year == 2020)

pts_poils_2021 <-  pts_poils %>%
  filter(Year == 2021)

#### 3.2. We match detection indices with traps  ####

# Association entre la détection et le piège à poil le plus proche en tenant compte des pièges actifs selon les années

# Photo
pts_photo_2020 <- dist_nearest(pts_photo_2020, 
                               sites_2020 %>% 
                                 filter(site != "hair") %>%
                                 mutate("trap" = row_number())) %>%
  #mapview(pts_photo_2020, cex = 2) + mapview(sites_2020, col.regions = "red", cex = 2)
  ## ASP: Distance from detection (photo) to the closest camera trap, and asigns row number of
  # the sites file (which is the same as the id of the camera)
  st_drop_geometry() %>%
  left_join(sites_2020 %>% 
              filter(site != "hair") %>%
              mutate("trap" = row_number()) %>% 
              st_drop_geometry(), by = "trap") %>%
  ## ASP: Joins the detections to the camera trap id (whch is the row number)
  dplyr::select(Year, month, Confirmed_Individual, Sex, Obs_type, trap_id, dist)


pts_photo_2021 <- dist_nearest(pts_photo_2021, 
                               sites_2021 %>% 
                                 filter(site != "hair") %>%
                                 mutate("trap" = row_number())) %>%
  #mapview(pts_photo_2021, cex = 2) + mapview(sites_2021, col.regions = "red", cex = 2)
  ## ASP: Distance from detection (photo) to the closest camera trap, and asigns row number of
  # the sites file (which is the same as the id of the camera)
  st_drop_geometry() %>%
  left_join(sites_2021 %>% 
              filter(site != "hair") %>%
              mutate("trap" = row_number()) %>% 
              st_drop_geometry(), by = "trap") %>%
  ## ASP: Joins the detections to the camera trap id (whch is the row number)
  dplyr::select(Year, month, Confirmed_Individual, Sex, Obs_type, trap_id, dist)

#Poils 

pts_poils_2020 <- dist_nearest(pts_poils_2020, 
                               sites_2020 %>%
                                 mutate("trap" = row_number())) %>%
  #mapview(pts_poils_2020, cex = 2) + mapview(sites_2020, col.regions = "red", cex = 2)
  st_drop_geometry() %>%
  left_join(sites_2020 %>%
              mutate("trap" = row_number()) %>% 
              st_drop_geometry(), by = "trap") %>%
  dplyr::select(Year, month, Confirmed_Individual, Sex, Obs_type, trap_id, dist)

pts_poils_2021 <- dist_nearest(pts_poils_2021, 
                               sites_2021 %>%
                                 mutate("trap" = row_number())) %>%
  #mapview(pts_poils_2021, cex = 2) + mapview(sites_2021, col.regions = "red", cex = 2)
  st_drop_geometry() %>%
  left_join(sites_2021 %>%
              mutate("trap" = row_number()) %>% 
              st_drop_geometry(), by = "trap") %>%
  dplyr::select(Year, month, Confirmed_Individual, Sex, Obs_type, trap_id, dist)


#### 3.3. We remove detection that are to far away a trap, because they must be errors ####

pts_2020 <- rbind(pts_photo_2020,pts_poils_2020) %>%
  filter(dist < seuil) 
pts_2021 <- rbind(pts_photo_2021,pts_poils_2021) %>%
  filter(dist < seuil)


## ---- SPAIN ----
## ---- 1. Load and sort out traps ----
#### a) Hair traps ####