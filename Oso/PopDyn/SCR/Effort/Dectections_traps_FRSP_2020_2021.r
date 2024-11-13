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

pts_2020_FR <- rbind(pts_photo_2020,pts_poils_2020) %>%
  filter(dist < seuil) 
pts_2021_FR <- rbind(pts_photo_2021,pts_poils_2021) %>%
  filter(dist < seuil)


## ---- SPAIN ----
## ---- 1. Load and sort out traps ----
#### a) 2020 ####

## Trampas Catalonia

trap_onlycat_2020 <- readxl::read_xlsx("D:/MargSalas/Oso/Datos/Effort_raw/Spain/Revisions trampes 2020 V2.xlsx", sheet = 1) %>% 
  janitor::clean_names() %>%
  filter(!is.na(coord_x)) %>%
  select(codi_tr,tipus_tr,coord_x,coord_y) %>%
  mutate(site = "hair") %>% # poils seul
  st_as_sf(coords = c("coord_x","coord_y"), 
           crs = CRS("+proj=utm +zone=31 +datum=WGS84"))

trap_onlycat_2020$site <- ifelse(trap_onlycat_2020$tipus_tr == "Mixte", "both","hair")

## Trampas Aran
#  De momento cojo las de Aran 2021, porque en teoría es un sistema
## nuevo de cuadrículas en el que las trampas se repiten y son las mismas para 2020 y 2021

trap_aran_2020 <- readxl::read_xlsx("D:/MargSalas/Oso/Datos/Effort_raw/Spain/Trampes_2021_SFF_CAR_CGA.xlsx") %>% 
  janitor::clean_names() %>%
  filter(comarca == "Vall d'Aran") %>%
  select(codi_tr,tipus_tr,coord_x,coord_y) %>%
  #mutate(trap_id = paste(row_number(), "c")) %>% 
  mutate(site = "hair") %>% # poils seul
  st_as_sf(coords = c("coord_x","coord_y"), 
           crs = CRS("+proj=utm +zone=31 +datum=WGS84")) 
#mutate(trap = row_number())

trap_aran_2020$site <- ifelse(trap_aran_2020$tipus_tr == "Mixte", "both","hair")

## Trampas Navarra

trap_nav_2020 <- readxl::read_xlsx("D:/MargSalas/Oso/Datos/Effort_raw/Spain/Copia de INDICIOS OSO PARDO-CAMARAS FOTOTRAMPEO 2021-22 B-GMA D-3 RONCAL-SALAZAR.xlsx", sheet = 1) %>% 
  janitor::clean_names() %>%
  rename(codi_tr = toponimia) %>%
  select(codi_tr, x_utm, y_utm) %>%
  mutate(tipus_tr = "Mixte") %>% # Lo pongo como Mixto, porque aunque no haya pelo estamos seguras de los individuos (Claverina)
  mutate(site = "both") %>%
  st_as_sf(coords = c("x_utm","y_utm"), 
           crs = CRS("+proj=utm +zone=30 +datum=WGS84") ) %>%
  st_transform(CRS("+proj=utm +zone=31 +datum=WGS84"))

#mapview(trap_nav_2020)

## Trampas añadidas de 2021
# Hay observaciones en 2020 de estaciones de muestreo que no están asociadas a una trampa en 2020, pero coinciden exactamente con una trampa en 2021.
# Añado estas trampas de 2021 porque lo más probable esque yo no las tenga pero ya existiesen.

trap_add2021 <- readxl::read_xlsx("D:/MargSalas/Oso/Datos/Effort_raw/Spain/Trampes_2021_SFF_CAR_CGA.xlsx") %>% 
  janitor::clean_names() %>%
  select(codi_tr,tipus_tr,coord_x,coord_y) %>%
  filter(codi_tr %in% c("D510502", "C730142", "E610189", "D710212")) %>%
  mutate(site = "hair") %>% # poils seul
  st_as_sf(coords = c("coord_x","coord_y"), 
           crs = CRS("+proj=utm +zone=31 +datum=WGS84")) 

trap_add2021$site <- ifelse(trap_add2021$tipus_tr == "Mixte", "both","hair")

## Trampa no registrada con observaciones múltiples en 2020 y 2021

# Calcular centroide entre detecciones de ambos años = trampa
setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
os_add_trap <- openxlsx::read.xlsx('Seguiment_Ossos_Pirineus_1996_2021_taula_final_2.xlsx') %>%
  filter(Year %in% c(2020,2021) & Confirmed_Individual != "Indetermined" & Country == "Spain") %>%
  filter(Confirmed_Individual %in% c("Bonabe", "New 20-14") & Date_register %in% c("17/08/2020", "06/09/2021") & Obs_type %in% c("Hair", "Video"))

mx <- mean(os_add_trap[,colnames(os_add_trap) == "x_long"])
my <- mean(os_add_trap[,colnames(os_add_trap) == "y_lat"])

trap_add <- data.frame(codi_tr = "new", tipus_tr = "Mixte", site = "both", x = mx, y = my) %>%
  st_as_sf(coords = c("x","y"), 
           crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) %>%
  st_transform(crs = 32631) 

# Save to use in 2021
setwd("D:/MargSalas/Oso/Datos/Effort_raw/Spain")
#save(trap_add, file = "trap_add.RData")

# Join Catalonia, Aran and Navarra, new traps 

trap_cat_2020 <- trap_onlycat_2020 %>%
  rbind(trap_aran_2020) %>%
  rbind(trap_nav_2020) %>%
  rbind(trap_add2021) %>%
  rbind(trap_add) %>%
  arrange(by = site) %>% # Very important for order later
  mutate(trap_id = paste(row_number(), "c"))  %>% ## This is the ID of all the traps together
  mutate(trap = row_number())  ## This is different to Maelis, who adds it in the next step**
# ** I do it like this because I will join the hair detections to all the traps (both and hair),
# and I need that the trap number is already in

## Creamos un object para cada tipo de trampa

trap_cat_foto_2020 <- trap_cat_2020 %>% 
  filter(site == "both") # VERY IMPORTANT that "both" is located first in trap_cat, so that trap number is the same as row number. Necesary for dist_nearest

trap_cat_pels_2020 <- trap_cat_2020 %>% # This file is not used to join with detections (I join with all hair and both)
  filter(site == "hair")  # ONLY for ploting (trap is not the row number)

#mapview(trap_cat_pels) + mapview(trap_cat_foto, col.regions = "red")

#### a) 2021 ####

## Trampas Catalonia y Aran

trap_cat_aran_2021 <- readxl::read_xlsx("D:/MargSalas/Oso/Datos/Effort_raw/Spain/Trampes_2021_SFF_CAR_CGA.xlsx") %>% 
  janitor::clean_names() %>%
  select(codi_tr,tipus_tr,coord_x,coord_y) %>%
  #mutate(trap_id = paste(row_number(), "c")) %>% 
  mutate(site = "hair") %>% # poils seul
  st_as_sf(coords = c("coord_x","coord_y"), 
           crs = CRS("+proj=utm +zone=31 +datum=WGS84")) 
#mutate(trap = row_number())

trap_cat_aran_2021$site <- ifelse(trap_cat_aran_2021$tipus_tr == "Mixte", "both","hair")

## Trampas Navarra

trap_nav_2021 <- readxl::read_xlsx("D:/MargSalas/Oso/Datos/Effort_raw/Spain/Copia de INDICIOS OSO PARDO-CAMARAS FOTOTRAMPEO 2021-22 B-GMA D-3 RONCAL-SALAZAR.xlsx", sheet = 2) %>% 
  janitor::clean_names() %>%
  rename(codi_tr = toponimia) %>%
  select(codi_tr, x_utm, y_utm) %>%
  mutate(tipus_tr = "Mixte") %>% # Lo pongo como Mixto, porque aunque no haya pelo estamos seguras de los individuos (Claverina)
  mutate(site = "both") %>%
  st_as_sf(coords = c("x_utm","y_utm"), 
           crs = CRS("+proj=utm +zone=30 +datum=WGS84") ) %>%
  st_transform(CRS("+proj=utm +zone=31 +datum=WGS84"))

#mapview(trap_nav_2021)

## Trampa no registrada con observaciones múltiples en 2020 y 2021 (average location detections)
setwd("D:/MargSalas/Oso/Datos/Effort_raw/Spain")
load("trap_add.RData")


# Join Catalonia, Aran and Navarra

trap_cat_2021 <- trap_cat_aran_2021 %>%
  rbind(trap_nav_2021) %>%
  rbind(trap_add) %>%
  arrange(by = site) %>% # Very important for order later
  mutate(trap_id = paste(row_number(), "c"))  %>% ## This is the ID of all the traps together
  mutate(trap = row_number())  ## This is different to Maelis, who adds it in the next step**
# ** I do it like this because I will join the hair detections to all the traps (both and hair),
# and I need that the trap number is already in

## Creamos un object para cada tipo de trampa

trap_cat_foto_2021 <- trap_cat_2021 %>% 
  filter(site == "both") # VERY IMPORTANT that "both" is located first in trap_cat, so that trap number is the same as row number. Necesary for dist_nearest

trap_cat_pels_2021 <- trap_cat_2021 %>% # This file is not used to join with detections (I join with all hair and both)
  filter(site == "hair")  # ONLY for ploting (trap is not the row number)

#mapview(trap_cat_pels) + mapview(trap_cat_foto, col.regions = "red")

## ---- 2. Detections ----

setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
os <- openxlsx::read.xlsx('Seguiment_Ossos_Pirineus_1996_2021_taula_final_2_cubLocations.xlsx') %>%
  filter(Year %in%  c(2020,2021) & Confirmed_Individual != "Indetermined" & Country == "Spain") 
os <- os %>%
  mutate(date = as_date(os$Date_register, format = "%d/%m/%Y"),
         month = month(date))

os$month[302] <- 6 # Correct mannually because date was not exact and we only had month
os$month[303] <- 4
os$month[452] <- 9
os$month[453] <- 9

os <- os %>%
  filter(month < 12, month > 4)  # 7 months form may to november 


## Keep systematic data
dat_cat_syst <- os[which(os$Method %in% c("Sampling_station", "Transect") &
                           os$Obs_type %in% c("Photo","Photo/Video", "Hair", "Video")), ]

#### 2.1. Separate the detections of each type  ####

## 2020

dat_cat_syst_2020 <- dat_cat_syst %>% 
  filter(Year == 2020)

pts_foto_2020 <- dat_cat_syst_2020 %>% 
  st_as_sf(coords = c("x_long","y_lat"), 
           crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) %>%
  st_transform(crs = 32631) %>%
  filter(Obs_type %in% c("Photo","Photo/Video", "Video"))

pts_pels_2020 <- dat_cat_syst_2020 %>% 
  st_as_sf(coords = c("x_long","y_lat"), 
           crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) %>%
  st_transform(crs = 32631) %>%
  filter(Obs_type %in% c("Hair"))

## 2021

dat_cat_syst_2021 <- dat_cat_syst %>% 
  filter(Year == 2021)

pts_foto_2021 <- dat_cat_syst_2021 %>% 
  st_as_sf(coords = c("x_long","y_lat"), 
           crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) %>%
  st_transform(crs = 32631) %>%
  filter(Obs_type %in% c("Photo","Photo/Video", "Video"))

pts_pels_2021 <- dat_cat_syst_2021 %>% 
  st_as_sf(coords = c("x_long","y_lat"), 
           crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) %>%
  st_transform(crs = 32631) %>%
  filter(Obs_type %in% c("Hair"))

#### 2.2. We match detection indices with traps  ####
# Association detection - nearest trap

##2020

pts_foto_2020 <- dist_nearest(pts_foto_2020, trap_cat_foto_2020) %>%
  #st_drop_geometry() %>%
  left_join(trap_cat_foto_2020 %>% st_drop_geometry(), by = "trap") %>%
  dplyr::select(Confirmed_Individual, month, Sex, trap_id, trap, dist)

pts_pels_2020 <- dist_nearest(pts_pels_2020, trap_cat_2020) %>% ## Here I join to trap_cat that contains hair and mixed (both, this is a bit different than maelis)
  #st_drop_geometry() %>%
  left_join(trap_cat_2020 %>% st_drop_geometry(), by = "trap") %>%
  dplyr::select(Confirmed_Individual, month, Sex, trap_id, trap, dist)

mapview(trap_cat_pels_2020) + mapview(trap_cat_foto_2020, col.regions = "green") + 
  mapview(pts_foto_2020, col.regions = "darkgreen", cex = 2) + 
  mapview(pts_pels_2020, col.regions = "magenta", cex = 2)

##2021

pts_foto_2021 <- dist_nearest(pts_foto_2021, trap_cat_foto_2021) %>%
  #st_drop_geometry() %>%
  left_join(trap_cat_foto_2021 %>% st_drop_geometry(), by = "trap") %>%
  dplyr::select(Confirmed_Individual, month, Sex, trap_id, trap, dist)

pts_pels_2021 <- dist_nearest(pts_pels_2021, trap_cat_2021) %>% ## Here I join to trap_cat that contains hair and mixed (both, this is a bit different than maelis)
  #st_drop_geometry() %>%
  left_join(trap_cat_2021 %>% st_drop_geometry(), by = "trap") %>%
  dplyr::select(Confirmed_Individual, month, Sex, trap_id, trap, dist)

mapview(trap_cat_pels_2021) + mapview(trap_cat_foto_2021, col.regions = "green") + 
  mapview(pts_foto_2021, col.regions = "darkgreen", cex = 2) + 
  mapview(pts_pels_2021, col.regions = "magenta", cex = 2)

# Here trap_id is = trap. I don't know why in MK is different, but this is the only way that works for me.
# Because every year will have a different trap set,

## With a threshold of 500 m there are few saved (IN 2020)

# Cases that could be worth checking with santi
# There is an observation but there is no trap associated, it is VERY far from any trap and the coordinates are good
## Esmolet 15/06, Pepito 12/08, New 18-03 05/07 (sist auto?there is even a camera?)
# Should we include a trap in this ones??

#### 2.3. Join, format, and remove detections further than threshold distance ####

seuil <- 500 #♣ Threshold distance, arbitrary (same as MK)

## 2020

# Format traps
trap_2020 <- trap_cat_2020 %>%
  rename("NOM" = "codi_tr") %>%
  mutate(pays = "Espagne") %>%
  mutate(suivi = "systematic") %>%
  select(NOM, site, trap_id, pays, trap, suivi, geometry)

# Combine and format detections
pts_2020_SP <- rbind(pts_foto_2020,pts_pels_2020) %>%
  rename("id" = "Confirmed_Individual") %>%
  rename("sex" = "Sex") %>%
  mutate(suivi = "systematic") %>%
  filter(dist < seuil) %>%
  select(-dist)

mapview(trap_cat) + mapview(pts_2020, col.regions = "green", cex = 2)

## 2021

# Format traps
trap_2021 <- trap_cat_2021 %>%
  rename("NOM" = "codi_tr") %>%
  mutate(pays = "Espagne") %>%
  mutate(suivi = "systematic") %>%
  select(NOM, site, trap_id, pays, trap, suivi, geometry)

# Combine and format detections
pts_2021_SP <- rbind(pts_foto_2021,pts_pels_2021) %>%
  rename("id" = "Confirmed_Individual") %>%
  rename("sex" = "Sex") %>%
  mutate(suivi = "systematic") %>%
  filter(dist < seuil) %>%
  select(-dist)

## ---- Combine France and Catalunya ----
# 1. Sites

sites_2020_t <- sites_2020 %>%
  rename("NOM" = "Nom_Iti") %>%
  mutate(pays = "France") %>%
  st_transform(map, crs = 32631) %>%
  rbind(trap_cat_2020 %>%
          mutate(pays = "Espagne") %>%
          rename("NOM" = "codi_tr")%>%
          mutate(effort = "it") %>%
          dplyr::select(NOM, site, effort, geometry, trap_id, pays)) %>%
  mutate(trap = row_number()) %>%
  #st_transform(map, crs = 32631) %>% # Transform a UTM to run in SCR
  mutate(X = unlist(map(geometry,1)),
         Y = unlist(map(geometry,2))) %>%
  dplyr::select(trap,trap_id,X,Y,site, pays) %>%
  st_drop_geometry()

sites_2021_t <- sites_2021 %>%
  rename("NOM" = "Nom_Iti") %>%
  mutate(pays = "France") %>%
  st_transform(map, crs = 32631) %>%
  rbind(trap_cat_2021 %>%
          mutate(pays = "Espagne") %>%
          rename("NOM" = "codi_tr")%>%
          mutate(effort = "it") %>%
          dplyr::select(NOM, site, effort, geometry, trap_id, pays)) %>%
  mutate(trap = row_number()) %>%
  #st_transform(map, crs = 32631) %>% # Transform a UTM to run in SCR
  mutate(X = unlist(map(geometry,1)),
         Y = unlist(map(geometry,2))) %>%
  dplyr::select(trap,trap_id,X,Y,site, pays) %>%
  st_drop_geometry()

# 2. Detections

pts_2020_t <- pts_2020_FR %>%
  mutate(id  = as.character(Confirmed_Individual)) %>%
  mutate(month = month %>%
           as.character() %>%
           as.numeric()) %>%
  mutate(sex  = as.character(Sex)) %>%
  mutate(trap_id  = as.character(trap_id)) %>%
  dplyr::select(id, month, sex, trap_id) %>%
  rbind(pts_2020_SP %>%
          st_drop_geometry() %>%
          mutate(id  = as.character(id)) %>%
          mutate(month = month %>%
                   as.character() %>%
                   as.numeric()) %>%
          mutate(sex  = as.character(sex)) %>%
          mutate(trap_id  = as.character(trap_id)) %>%
          dplyr::select(id, month, sex, trap_id)) %>%
  left_join(sites_2020_t, by = "trap_id") %>%
  dplyr::select(id, month, sex, trap, trap_id) %>%
  mutate(sex = as.factor(sex)) 

pts_2021_t <- pts_2021_FR %>%
  mutate(id  = as.character(Confirmed_Individual)) %>%
  mutate(month = month %>%
           as.character() %>%
           as.numeric()) %>%
  mutate(sex  = as.character(Sex)) %>%
  mutate(trap_id  = as.character(trap_id)) %>%
  dplyr::select(id, month, sex, trap_id) %>%
  rbind(pts_2021_SP %>%
          st_drop_geometry() %>%
          mutate(id  = as.character(id)) %>%
          mutate(month = month %>%
                   as.character() %>%
                   as.numeric()) %>%
          mutate(sex  = as.character(sex)) %>%
          mutate(trap_id  = as.character(trap_id)) %>%
          dplyr::select(id, month, sex, trap_id)) %>%
  left_join(sites_2021_t, by = "trap_id") %>%
  dplyr::select(id, month, sex, trap, trap_id) %>%
  mutate(sex = as.factor(sex)) 

## ---- Create EDF and TDF ----

# I DO IT IN THE COMMON SCRIPT :) Detections_Traps_2020_2021_ALL

