## -------------------------------------------------
##      Join detections - traps France and Spain    
##                      2017-2021
## ------------------------------------------------- 

# Divided in three parts:
# 1. Period 2017-2019
# 2. Period 2020-2021
# 3. Join all data in edf and tdf for SCR analysis

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


## -------------------------------------------------
##                      2017-2019
## ------------------------------------------------- 

## ---- FRANCE ----
## ---- 1. Load and sort out transects and traps ----
#### a) Itineraries ####
# Make one file by year with only itineraries that have been made. 

itineraires <- st_read("D:/MargSalas/Oso/Datos/Effort_raw/France/Raw/Raw/shp/Traps2017_2019/Itineraires_2017_2020.shp")%>%
  distinct(.keep_all = TRUE) %>%
  rename("NOM"="nom") %>%
  dplyr::select(NOM) %>%
  mutate(id = row_number())

# On enlève à la main les doublons qu'il reste
# l'itinéraire sans nom est le même que Rouze, Ste Colombe Escouloubre1 correspond à Ste Colombe, Ste Colombe Escouloubre2 correspond à Escouloubre. 
itineraires <- itineraires %>%
  filter(NOM != "Ste Colombe Escouloubre1" & NOM != "Ste Colombe Escouloubre2" & NOM != "Cazaux de Larboust") %>%
  # Il reste à enlever ceux qui ont le même nom mais qui diffère légèrement, soit : Arrioutort2014 : id = 44,    Er_Gazies2014 : id = 41, Larry2014 : id = 38, Sost  : id = 63 et Tuquet2014 : id = 48
  filter(id != 44 & id != 41 & id != 38 & id != 48)

# Les itinéraires qui n'ont pas été parcouru en 2017 : 
# "Bern", "Saoubette, "Orlu", "Formiguere", "Mijanes1", "Majanes2", "Rouze ou id = 64", "Escouloubre", "St Collombe", "Cazaux-Larboust"
itineraires_2017 <- itineraires %>%
  filter(NOM !="Bern" & NOM !="Saoubette" & NOM != "Escouloubre" & NOM != "Formiguere" & NOM != "Orlu" & NOM != "Mijanes1" & NOM != "Majanes2" & NOM != "Cazaux-Larboust" & NOM != "St Collombe" & NOM != "Rouze")  

# Les itinéraires de 2018
# St Collombe, Esouloubre et Formiguères actif à partir de juin
itineraires_2018 <- itineraires %>%
  filter(NOM !="Saoubette")

# Les itinéraires qui n'ont pas été parcouru en 2019 : 
# Prunadière, Burat et Bezin Garraux 
itineraires_2019 <- itineraires %>%
  filter(id != 53) %>% # retire Prunadière : Problème d'accent
  filter(NOM != "Burat" & NOM != "Bezin Garraux 2014")

#### b) Hair traps ####
# Make one hair traps file by year. 

##### b.1) Load hair traps from three french regions and join them
piege_poils1719 <- rgdal::readOGR("D:/MargSalas/Oso/Datos/Effort_raw/France/Raw/Raw/shp/Traps2017_2019/PiegesPoils_2017_2019.shp", use_iconv=TRUE, encoding = "UTF-8")

# On formate les pièges à poils en spatial
piege_poils1719 <- piege_poils1719 %>%
  st_as_sf(coords = coords) %>%
  dplyr::select(Nom) %>%
  rename("NOM"="Nom") %>%
  st_transform(crs = 27572)

# Importer les pièges à poils de Arrioutort 
arrioutort <- read.table("D:/MargSalas/Oso/Datos/Effort_raw/France/Raw/Raw/AppatsSmola_PO_2019Bis.csv", dec=",", header=T, sep=";")  %>%
  janitor::clean_names() %>% 
  rename("NOM" = "field_1",
         "X" = "x", 
         "Y" = "y") %>%
  dplyr::select(NOM, X, Y)  %>%
  filter (str_detect(NOM, "arrioutort")) %>%
  st_as_sf(coords = c("X", "Y"),
           crs = 27572)

# Importer les pièges à poils de Haute Arriège
haute_arriege <- read.table("D:/MargSalas/Oso/Datos/Effort_raw/France/Raw/Raw/AppareilsAutos2019Haute_Ariege_Maelis.csv", dec=",", header=T, sep=";") %>%
  janitor::clean_names() %>% 
  rename("X" = "x_lb_ii_etendu",
         "Y" = "y_lb_ii_etendu",
         "NOM" = "lieu") %>%
  dplyr::select(NOM, X, Y) %>%
  st_as_sf(coords = c("X", "Y"),
           crs = 27572)

# Fusionner tous les fichiers pièges 
piege_poils1719 <- rbind(piege_poils1719,haute_arriege, arrioutort) %>%
  distinct(.keep_all = TRUE) 

##### b.2) Distinguish hair trap that are on an itinerary ( <500 m, because they are more frequently visited) and outside
# And join it with itinerary id

seuil <- 500

piege_poils1719_it <- dist_nearest(piege_poils1719,itineraires) %>%
  filter(dist < seuil) %>%
  mutate("itineraire" = NA) 

# Pour chaque piège a poils on remet le bon numéro correspondant à l'itinéraire associé
for (i in 1 : length(piege_poils1719_it$trap)) {
  piege_poils1719_it$itineraire[i] <- itineraires$id[piege_poils1719_it$trap[i]]
  piege_poils1719_it$NOM[i] <- itineraires$NOM[piege_poils1719_it$trap[i]]
}

# On enlève la variable trap pour ne pas se mélanger par la suite 
piege_poils1719_it <- piege_poils1719_it  %>%
  dplyr::select(geometry,NOM,itineraire)

##### b.3) One hair trap file by year (remove traps that are on an itinerary that haven't been made)

## 2017 
piege_poils_it_2017 <- piege_poils1719_it %>% 
  filter(NOM !="Bern" & NOM !="Saoubette" & NOM != "Escouloubre" & NOM != "Formiguere" & NOM != "Orlu" & NOM != "Mijanes1" & NOM != "Majanes2" & NOM != "Cazaux-Larboust" & NOM != "St Collombe" & NOM != "Rouze") 

## 2018 
piege_poils_it_2018 <- piege_poils1719_it %>% 
  filter(NOM !="Saoubette")

## 2019 
piege_poils_it_2019 <- piege_poils1719_it %>%
  filter(itineraire != 53) %>% # retire Prunadière : Problème d'accent
  filter(NOM != "Burat" & NOM != "Bezin Garraux 2014")


#### c) Camera traps ####
# Make one camera traps file by year
# For each year, by knowledge --> Remove camera traps that aren't close to a hair trap. 
# Only camera traps can't allow the identification of individuals so we remove them from the study. 

# Modification des coordonnées de Bibet :  X = 501083 et Y = 1755754 
# sinon : X = 501292.999994 et Y= 1754615.999986
bibet <- data.frame(NOM = "Bibet", 
                    X = 501083,
                    Y = 1755754) %>%
  st_as_sf(coords = c("X", "Y"),
           crs = 27572)

piege_photos_2017 <- st_read("D:/MargSalas/Oso/Datos/Effort_raw/France/Raw/Raw/shp/Camera_trap_france/AppareilsAuto2017.shp") %>%
  rename("NOM" = "Nom") %>%
  st_drop_geometry() # On arrondit les coordonnnées pour avoir la même précision entre les années et pouvoir joindre le tableau avec le format des appareils (photo, vidéo)

# on évite de perdre de l'info 
piege_photos_2017$X[49] <- 433673 
piege_photos_2017$Y[49] <- 1751301
piege_photos_2017$X[50] <- 362206
piege_photos_2017$Y[50] <- 1763137
piege_photos_2017$X <- round(piege_photos_2017$X,0) 
piege_photos_2017$Y <- round(piege_photos_2017$Y,0) 

# On repasse en spatial
piege_photos_2017 <- piege_photos_2017 %>%
  dplyr::select(NOM, X, Y) %>%
  st_as_sf(coords = c("X", "Y"),
           crs = 27572) %>%
  filter(NOM != "Bibet" | is.na(NOM)) %>% # on enlève Bibet car on va modifier ses coordonnées
  filter(NOM != "Galedre" & NOM != "St Jean" & NOM != "St Mamet" & NOM != "Paloumère" & NOM != "Saubé"  & NOM != "Soulas" | is.na(NOM)) # on enlève les pièges qui n'ont pas de pièges à poils à proximité, ils ont pour objectif de faire de la détection mais ne permettent pas l'identification des individus. 
piege_photos_2017 <-piege_photos_2017[-43,] # Il reste à enlever Coueq

piege_photos_2018 <- st_read("D:/MargSalas/Oso/Datos/Effort_raw/France/Raw/Raw/shp/Camera_trap_france/AppareilsAutos2018.shp")
piege_photos_2018 <- st_set_crs(piege_photos_2018,27572) %>%
  filter(NOM != "Bibet"|is.na(NOM)) %>%
  filter(NOM != "Rib\u0082rot" & NOM != "St Jean" & NOM != "Sarrouges" & NOM != "Paloum\u008are" & NOM != "Coueq"  | is.na(NOM)) %>% # Piège photo de détection
  dplyr::select(NOM)

piege_photos_2019 <- st_read("D:/MargSalas/Oso/Datos/Effort_raw/France/Raw/Raw/shp/Camera_trap_france/AppareilsAutos2019.shp")
piege_photos_2019 <- st_set_crs(piege_photos_2019,27572) %>%
  dplyr::select(NOM) %>%
  filter(NOM != "Bibet"|is.na(NOM)) %>%
  filter(NOM != "Galedre" & NOM != "Rib\u0082rot" & NOM != "St Jean" & NOM != "Sarrouges"  & NOM != "Coueq"& NOM != "Comus"& NOM != "Merial"  | is.na(NOM)) # Pieges photos de détection

# Ajout d'un appareil photo vallée du lys 
vallee_lys <- data.frame(NOM = "Vallee du lys", 
                         X = 456512,
                         Y = 1749685) %>%
  st_as_sf(coords = c("X", "Y"),
           crs = 27572)

# Ajout de l'appareil Bibet et vallee du lys. Puis on ajoute un numéro a chaque piège photo
piege_photos_2017 <- rbind(piege_photos_2017,bibet)
piege_photos_2018 <- rbind(piege_photos_2018,bibet, vallee_lys)
piege_photos_2019 <- rbind(piege_photos_2019, bibet, vallee_lys)



## ---- 2. Create SITES  ----

# Site = zone where different kind of traps (photo, video, hair) can capture an individual in a radius of 500m
# We consider that detections among those trap can't be independent because they are really close to each other.
# Then we consider them together and define a special detection probability for each of these sites. 

#We can differentiates these sites by a categorical covariate which take into account the types of traps present at a site :
#  - camera = only camera trap (NOT USED)
#  - hair = only hair trap (have to be on an itinerary)
#  - both = camera trap/ hair trap

# Para cada año, 3 acciones (más detalle en script Tidy Maelis):

# 1)) Seleccionamos las trampas de pelo que están en los transectos. De estas, 
# miramos cuales tienen cámara en proximidad y asociamos a cada cámara a menos 
# de 500m de una trampa de pelo, la trampa de pelo más cercana
## ** Importante! Las coordenadas de las trampas de pelo servirán como puntos GPS de nuestros "sites".

# 2)) En segundo lugar, nos interesan las cámaras que no están
# en transectos (es decir, cualquier cámara a más de 500 m de una cámara ubicada en un transecto). 
# Para estas cámaras miramos si hay una trampa de pelo a menos de 500 m
# En caso afirmativo, modificamos la covariable "site" (3); de lo contrario, no se hace nada (1).

# 3)) Queda por formar la matriz de "sites". Los "sites" corresponden a:
#   - Trampas de pelo en transectos (asociadas o no a un cámara) 
#   - Cámaras que no están en transectos (asociadas o no a una trampa de pelo). En este caso,
# todas las trampas de pelo que no están en transecto o que no están asociados con una cámara se eliminan de la lista ya que no son visitados regularmente.

#### 2.1. 2017 ####

# 1))
# pour tout appareil photo on cherche le piège à poils sur un itinéraire le plus proche. Et on calcule la distance les séparants
piege_photos_2017 <- dist_nearest(piege_photos_2017, piege_poils_it_2017) 

## ASP: trap es el id (número de fila en piege_poils_it_2017) de la trampa de pelo 
# más cercana a esa cámara

# On sépare les pièges photos qui sont à plus de 500m d'un piège à poils de ceux qui sont à moins de 500 m d'un piège à poils qui se situe sur un itinéraire. 
piege_photos_it_2017 <- piege_photos_2017 %>%
  filter(dist < seuil) # piège photo sur un itineraire

piege_photos_syst_2017 <-  piege_photos_2017 %>%
  filter(dist >= seuil) %>% # piege photo hors d'un itineraire
  mutate("site" = "camera") %>% # piège photo seul
  dplyr::select(NOM, site)

piege_poils_it_2017 <- piege_poils_it_2017 %>%
  mutate("site" = "hair") # Site possède un piège à poils seul

# 2))
# Si les pièges à poils des itinéraires sont associés à au moins un piège photos alors on modifie la covariable site pour les rattacher au piège à poil le plus proche. 
# Pour chaque piège photo sur un itinéraire on l'associé au piège à poil le plus proche, ce qui forme un site de type 3 (piège à poil + piège photo) 
for (i in 1 : length(piege_photos_it_2017$trap)){ 
  piege_poils_it_2017$site[piege_photos_it_2017$trap[i]] <- "both"
}

# Ajout d'un piège à poils combiné à un appareil photo qu'il manque
mines_blende <- data.frame(NOM = "Mines de Blende", 
                           X = 476618,
                           Y = 1765777,
                           site = "both") %>%
  st_as_sf(coords = c("X", "Y"),
           crs = 27572)

# Ajout d'un appareil vidéo à Bouquemont 
bouquemont <- data.frame(NOM = "Bouquemont", 
                         X = 473144,
                         Y = 1762179,
                         site = "both") %>%
  st_as_sf(coords = c("X", "Y"),
           crs = 27572)
# On vérifie que les différents tableaux ont les mêmes variables
piege_poils_it_2017 <- piege_poils_it_2017 %>%
  dplyr::select(NOM, site)

# 3)) On combine tous les sites
sites_2017 <- rbind(piege_poils_it_2017, piege_photos_syst_2017, mines_blende, bouquemont) 

# On ajoute la variable effort qui indique si le site se situe ou non sur un itinéraire
effort <- c(rep("it",length(piege_poils_it_2017$NOM)),
            rep("hors it",length(piege_photos_syst_2017$NOM)),
            "hors it","it")

sites_2017 <- cbind(sites_2017, effort)

# Fos aussi piège à poils -> site == 3 
sites_2017[which(sites_2017$NOM == "Fos"), "site"] <- "both"

# On enlève les sites avec uniquement des cameras car il ne permettent pas l'identification des individus 
sites_2017 <- sites_2017 %>%
  filter(site != "camera")

# On les range par type de site et on leur attribue un numéro
sites_2017 <- sites_2017 %>%
  arrange(site) %>%
  mutate(trap_id = paste(row_number(),"f"))

# Transform to make crs fit with detections
sites_2017 <- sites_2017 %>%
  st_transform(crs = 4326)

#### 2.2. 2018 ####
# 1))
# Identifier les appareils qui sont associé à des pièges à poils qui sont sur des itinéraires  
piege_photos_2018 <- dist_nearest(piege_photos_2018, piege_poils_it_2018) 

# On sépare les pièges photos qui sont à plus de 500m d'une piège à poils de ceux qui sont à moins de 500 m d'un piège à poils qui se situe sur un itinéraire. 
piege_photos_it_2018 <- piege_photos_2018 %>%
  filter(dist < seuil) 

piege_photos_syst_2018 <-  piege_photos_2018 %>%
  filter(dist >= seuil) %>%
  mutate("site" = "camera") %>% # Site possède un piège photo et un piège à poils
  dplyr::select(NOM, site)

piege_poils_it_2018 <- piege_poils_it_2018 %>%
  mutate("site" = "hair") # Site possède au moins un piège à poils seul

# 2)) Si les pièges à poils des itinéraires sont associés à un piège photos alors on modifie la covariable site
for (i in 1 : length(piege_photos_it_2018$trap)){
  piege_poils_it_2018$site[piege_photos_it_2018$trap[i]] <- "both"
}

# On vérifie que les différents tableaux ont les mêmes variables
piege_poils_it_2018 <- piege_poils_it_2018 %>%
  dplyr::select(NOM, site)

# 3)) On combine tous les sites
sites_2018 <- rbind(piege_poils_it_2018, piege_photos_syst_2018, mines_blende, bouquemont) 
effort <- c(rep("it",length(piege_poils_it_2018$NOM)),
            rep("hors it",length(piege_photos_syst_2018$NOM)),
            "hors it","it")
sites_2018 <- cbind(sites_2018, effort)

# Camera trap hors des itinéraires qui possèdent tout de même un piège à poil 
sites_2018[which(sites_2018$NOM == "Fos"), "site"] <- "both"
sites_2018[which(sites_2018$NOM == "Bouquemont"), "site"] <- "both"
sites_2018[which(sites_2018$NOM == "Seridere"), "site"] <- "both"
sites_2018[which(sites_2018$NOM == "Le Mail"), "site"] <- "both"

# On enlève les sites avec uniquement des cameras car il ne permettent pas l'identification des individus 
sites_2018 <- sites_2018 %>%
  filter(site != "camera")

# On les range par type de site et on leur attribue un numéro
sites_2018 <- sites_2018 %>%
  arrange(site) %>%
  mutate(trap_id = paste(row_number(),"f"))

# Transform to make crs fit with detections
sites_2018 <- sites_2018 %>%
  st_transform(crs = 4326)

#### 2.3. 2019 ####
# 1))
# Identifier les appareils qui sont associé à des pièges à poils qui sont sur des itinéraires  
piege_photos_2019 <- dist_nearest(piege_photos_2019, piege_poils_it_2019) 

# On sépare les pièges photos qui sont à plus de 500m d'une piège à poils de ceux qui sont à moins de 500 m d'un piège à poils qui se situe sur un itinéraire. 
piege_photos_it_2019 <- piege_photos_2019 %>%
  filter(dist < seuil) 

piege_photos_syst_2019 <-  piege_photos_2019 %>%
  filter(dist >= seuil) # Cameras far from transects

piege_photos_syst_2019 <- dist_nearest(piege_photos_syst_2019, piege_poils1719) # This is only done this year so a bit weird, maybe check??

piege_photos_syst_both_2019 <-  piege_photos_syst_2019 %>%
  filter(dist < seuil) %>%
  mutate("site" = "both") %>% # Site possède un piège photo et un piège à poils
  dplyr::select(NOM, site)

piege_photos_syst_cam_2019 <-  piege_photos_syst_2019 %>%
  filter(dist >= seuil) %>%
  mutate("site" = "camera") %>% # Site possède un piège photo seulement
  dplyr::select(NOM, site)
piege_photos_syst_cam_2019$NOM[6] <- "Lassas"
piege_photos_syst_2019 <- rbind(piege_photos_syst_both_2019,
                                piege_photos_syst_cam_2019)

piege_poils_it_2019 <- piege_poils_it_2019 %>%
  mutate("site" = "hair") # Site possède au moins un piège à poils seul

# 2)) Si les pièges à poils des itinéraires sont associés à un piège photos alors on modifie la covarible site
for (i in 1 : length(piege_photos_it_2019$trap)){
  piege_poils_it_2019$site[piege_photos_it_2019$trap[i]] <- "both"
}

# On vérifie que les différents tableaux ont les mêmes variables
piege_poils_it_2019 <- piege_poils_it_2019 %>%
  dplyr::select(NOM, site)

# 3)) Join
sites_2019 <- rbind(piege_poils_it_2019, piege_photos_syst_2019, bouquemont)
effort <- c(rep("it",length(piege_poils_it_2019$NOM)),
            rep("hors it",length(piege_photos_syst_2019$NOM)),
            "it")
sites_2019 <- cbind(sites_2019, effort)

# Camera hors it avec un piège photo 
sites_2019[which(sites_2019$NOM == "Fos"), "site"] <- "both"
sites_2019[which(sites_2019$NOM == "Bouquemont"), "site"] <- "both"
sites_2019[which(sites_2019$NOM == "Seridere"), "site"] <- "both"
sites_2019[which(sites_2019$NOM == "Lassas"), "site"] <- "both"

# On enlève les sites avec uniquement des cameras car il ne permettent pas l'identification des individus 
sites_2019 <- sites_2019 %>%
  filter(site != "camera")

# On les range par type de site et on leur attribue un numéro
sites_2019 <- sites_2019 %>%
  arrange(site) %>% 
  mutate(trap_id = paste(row_number(),"f"))

# Transform to make crs fit with detections
sites_2019 <- sites_2019 %>%
  st_transform(crs = 4326)

## ---- 3. Detections ----

# LOAD dataset that has been contrasted with dataset of Maelis (Prospind_2017-2019_Maelisv2.xlsx)
# We have also added if the hair sample comes from a trap (Poils appât) or not (poils spontanées)

setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
os1719_FR <- openxlsx::read.xlsx('Seguiment_Ossos_Pirineus_1996_2021_taula_final.xlsx') %>% 
  filter(Year %in% c(2017,2018,2019) & Confirmed_Individual != "Indetermined" & Country == "France") # Removed indetermined

os1719_FR <- os1719_FR %>% mutate(date = as_date(os1719_FR$Date, format = "%d/%m/%Y"),
                                  month = month(date)) %>%
  filter(month < 12, month > 4) %>% # 7 months form may to november
  filter(!is.na(x_long)) # Remove observations without coordinates


## Keep systematic data
dat2_1719_FR <- os1719_FR[which(os1719_FR$Method %in% c("Sampling_station", "Transect") &
                                  os1719_FR$Obs_type %in% c("Photo","Photo/Video", "Hair", "Video")), ]

# We distinguish the 3 methods 
dat2_1719_FR <- dat2_1719_FR %>%
  mutate(method = "itineraire") # Itinerary 

dat2_1719_FR[which(dat2_1719_FR$Obs_type == "Photo" | dat2_1719_FR$Obs_type == "Photo/Video" | dat2_1719_FR$Obs_type == "Video"),"method"] <- "piege photos" # Camera trap

nrow(dat2_1719_FR[which(dat2_1719_FR$Obs_type == "Hair"),]) - nrow(dat2_1719_FR[which(dat2_1719_FR$Obs_type == "Hair" & dat2_1719_FR$hairtrap == 1),]) # We remove 36 that are hair 
dat2_1719_FR[which(dat2_1719_FR$Obs_type == "Hair" & dat2_1719_FR$hairtrap == 1),"method"] <- "piege poils" # Hair trap

# HERE we can trust that hairtrap = 1 is really a hair trap because it comes from Maelis (only in France 17-19), but only here!
# Then, the information of hairtraps comes from Yon that has assigned it quite randomly (if the hair comes from systematic sampling)

#### 3.1. Separate the detections of each type  ####

# Détections issues de pièges photos (photographie automatiques, vidéo automatique)
pts_photo1719 <- dat2_1719_FR %>%
  filter(method =="piege photos") %>%
  st_as_sf(coords = c("x_long", "y_lat"),
           crs = 4326)

pts_photo_2017 <- pts_photo1719 %>%
  filter(Year == 2017) 

pts_photo_2018 <- pts_photo1719 %>%
  filter(Year == 2018) 

pts_photo_2019 <- pts_photo1719 %>%
  filter(Year == 2019)

# Détections issues de pièges à poils 
pts_poils <- dat2_1719_FR %>%
  filter(method =="piege poils") %>%
  st_as_sf(coords = c("x_long", "y_lat"),
           crs = 4326)

pts_poils_2017 <-  pts_poils %>%
  filter(Year == 2017)

pts_poils_2018 <-  pts_poils %>%
  filter(Year == 2018)

pts_poils_2019 <-  pts_poils %>%
  filter(Year == 2019)

#### 3.2. We match detection indices with traps  ####

# Association entre la détection et le piège à poil le plus proche en tenant compte des pièges actifs selon les années

# Photo
pts_photo_2017 <- dist_nearest(pts_photo_2017, 
                               sites_2017 %>% 
                                 filter(site != "hair") %>%
                                 mutate("trap" = row_number())) %>%
  #mapview(pts_photo_2017, cex = 2) + mapview(sites_2017[c(418,24),], col.regions = "red", cex = 2)
  ## ASP: Distance from detection (photo) to the closest camera trap, and asigns row number of
  # the sites file (which is the same as the id of the camera)
  st_drop_geometry() %>%
  left_join(sites_2017 %>% 
              filter(site != "hair") %>%
              mutate("trap" = row_number()) %>% 
              st_drop_geometry(), by = "trap") %>%
  ## ASP: Joins the detections to the camera trap id (whch is the row number)
  dplyr::select(Year, month, Confirmed_Individual, Sex, Obs_type, trap_id, dist)


pts_photo_2018 <- dist_nearest(pts_photo_2018, 
                               sites_2018 %>% 
                                 filter(site != "hair") %>%
                                 mutate("trap" = row_number())) %>%
  # mapview(pts_photo_2018, cex = 2) + mapview(sites_2018, col.regions = "red", cex = 2)
  ## ASP: Distance from detection (photo) to the closest camera trap, and asigns row number of
  # the sites file (which is the same as the id of the camera)
  st_drop_geometry() %>%
  left_join(sites_2018 %>% 
              filter(site != "hair") %>%
              mutate("trap" = row_number()) %>% 
              st_drop_geometry(), by = "trap") %>%
  ## ASP: Joins the detections to the camera trap id (whch is the row number)
  dplyr::select(Year, month, Confirmed_Individual, Sex, Obs_type, trap_id, dist)

pts_photo_2019 <- dist_nearest(pts_photo_2019, 
                               sites_2019 %>% 
                                 filter(site != "hair") %>%
                                 mutate("trap" = row_number())) %>%
  
  #mapview(pts_photo_2019, cex = 2) + mapview(sites_2019, col.regions = "red", cex = 2)
  
  ## ASP: Distance from detection (photo) to the closest camera trap, and asigns row number of
  # the sites file (which is the same as the id of the camera)
  st_drop_geometry() %>%
  left_join(sites_2019 %>% 
              filter(site != "hair") %>%
              mutate("trap" = row_number()) %>% 
              st_drop_geometry(), by = "trap") %>%
  ## ASP: Joins the detections to the camera trap id (whch is the row number)
  dplyr::select(Year, month, Confirmed_Individual, Sex, Obs_type, trap_id, dist)


#Poils 
pts_poils_2017 <- dist_nearest(pts_poils_2017, 
                               sites_2017 %>%
                                 mutate("trap" = row_number())) %>%
  #mapview(pts_poils_2017, cex = 2) + mapview(sites_2017[c(418,24),], col.regions = "red", cex = 2)
  st_drop_geometry() %>%
  left_join(sites_2017 %>%
              mutate("trap" = row_number()) %>% 
              st_drop_geometry(), by = "trap") %>%
  dplyr::select(Year, month, Confirmed_Individual, Sex, Obs_type, trap_id, dist)

pts_poils_2018 <- dist_nearest(pts_poils_2018, 
                               sites_2018 %>%
                                 mutate("trap" = row_number())) %>%
  st_drop_geometry() %>%
  left_join(sites_2018 %>%
              mutate("trap" = row_number()) %>% 
              st_drop_geometry(), by = "trap") %>%
  dplyr::select(Year, month, Confirmed_Individual, Sex, Obs_type, trap_id, dist)

pts_poils_2019 <- dist_nearest(pts_poils_2019, 
                               sites_2019 %>%
                                 mutate("trap" = row_number())) %>%
  st_drop_geometry() %>%
  left_join(sites_2019 %>%
              mutate("trap" = row_number()) %>% 
              st_drop_geometry(), by = "trap") %>%
  dplyr::select(Year, month, Confirmed_Individual, Sex, Obs_type, trap_id, dist)

#### 3.3. We remove detection that are to far away a trap, because they must be errors ####

pts_2017 <- rbind(pts_photo_2017,pts_poils_2017) %>%
  filter(dist < seuil) 
pts_2018 <- rbind(pts_photo_2018,pts_poils_2018) %>%
  filter(dist < seuil)
pts_2019 <- rbind(pts_photo_2019,pts_poils_2019) %>%
  filter(dist < seuil)

## ---- SPAIN ----
## ---- 1. Load and sort out traps ----
#### a) Hair traps ####

trap_cat1719 <- readxl::read_xlsx("D:/MargSalas/Oso/Datos/Effort_raw/France/Raw/Raw/UTM trampas pelos y camaras.xlsx") %>% 
  janitor::clean_names() %>%
  filter(!is.na(xutm)) %>%
  st_as_sf(coords = c("xutm","yutm"), 
           crs = CRS("+proj=utm +zone=31 +datum=WGS84")) %>% 
  st_transform(crs = 4326) %>%
  rename("NOM" = "itinerari") %>%
  mutate(effort = "it") %>%
  mutate(site = "hair") # poils seul

## ASP: Names all "hair" and then according to foto_video or pels columns assings to camera or both
trap_cat1719[which(!is.na(trap_cat1719$foto_video)),"site"] <- "both" # poils et photo
trap_cat1719[which(is.na(trap_cat1719$pels)),"site"] <- "camera" # photo seul

trap_cat1719 <- trap_cat1719[-which(duplicated(trap_cat1719$geometry)),] # ASP: MK did not do this step, but I do it
# because they are duplicates and it gives problems when assigning detections to traps

# On enlève les pièges photos seul car ils ne permettent pas l'identification des ours 
trap_cat1719 <- trap_cat1719 %>%
  filter(site != "camera")  %>%
  arrange(by = site) %>% # I arrange it so that numbers make more sense (fotos first, then both)
  mutate(trap_id = paste(row_number(), "c"))

# On fait un sous tableau pour les pièges à poils 
trap_cat1719_pels <- trap_cat1719 %>%
  mutate(trap = row_number())

#### b) Camera traps ####

trap_cat1719_foto <- trap_cat1719 %>%
  filter(foto_video != is.na(foto_video)) %>% 
  mutate(trap = row_number()) 

## ---- 2. Detections ----

setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
os1719_SP <- openxlsx::read.xlsx('Seguiment_Ossos_Pirineus_1996_2021_taula_final.xlsx') %>% 
  filter(Year %in% c(2017,2018,2019) & Confirmed_Individual != "Indetermined" & Country == "Spain") # Removed indetermined

os1719_SP <- os1719_SP %>% mutate(date = as_date(os1719_SP$Date, format = "%d/%m/%Y"),
                                  month = month(date)) %>%
  filter(month < 12, month > 4) %>% # 7 months form may to november
  filter(!is.na(x_long)) # Remove observations without coordinates


## Keep systematic data
dat_cat1719_syst <- os1719_SP[which(os1719_SP$Method %in% c("Sampling_station", "Transect") &
                                      os1719_SP$Obs_type %in% c("Photo","Photo/Video", "Hair", "Video")), ]

#### 2.1. Separate the detections of each type  ####

# Camera trap detections
pts_foto <- dat_cat1719_syst %>% 
  st_as_sf(coords = c("x_long","y_lat"), 
           crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) %>%
  filter(Obs_type %in% c("Photo","Photo/Video", "Video"))

pts_foto_2017 <- pts_foto %>%
  filter(Year == 2017) 

pts_foto_2018 <- pts_foto %>%
  filter(Year == 2018) 

pts_foto_2019 <- pts_foto %>%
  filter(Year == 2019) 

# Hair detections
pts_pels <- dat_cat1719_syst %>% 
  st_as_sf(coords = c("x_long","y_lat"), 
           crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) %>%
  filter(Obs_type %in% c("Hair"))

pts_pels_2017 <- pts_pels %>% 
  filter(Year == 2017) 

pts_pels_2018 <- pts_pels %>% 
  filter(Year == 2018)

pts_pels_2019 <- pts_pels %>% 
  filter(Year == 2019)

#### 2.2. We match detection indices with traps  ####

# Association entre la détection et le piège photo le plus proche, 
# les mêmes pièges seraient actifs les 3 années (à vérifier)

# Camera trap detections

pts_foto_2017 <- dist_nearest(pts_foto_2017, trap_cat1719_foto) %>% ## ASP: It assings to the nearest trap. The name of the trap is the variable "trap" (which is created as the row number of trap_cat1719_foto)
  #st_drop_geometry() %>%
  left_join(trap_cat1719_foto %>% st_drop_geometry(), by = "trap") %>% # To get the trapID of the whole trap dataset, it joins it by the variable trap
  dplyr::select(Confirmed_Individual, month, Sex, trap_id, trap, dist)

pts_foto_2018 <- dist_nearest(pts_foto_2018, trap_cat1719_foto) %>% ## ASP: It assings to the nearest trap. The name of the trap is the variable "trap" (which is created as the row number of trap_cat1719_foto)
  #st_drop_geometry() %>%
  left_join(trap_cat1719_foto %>% st_drop_geometry(), by = "trap") %>% # To get the trapID of the whole trap dataset, it joins it by the variable trap
  dplyr::select(Confirmed_Individual, month, Sex, trap_id, trap, dist)

pts_foto_2019 <- dist_nearest(pts_foto_2019, trap_cat1719_foto) %>% ## ASP: It assings to the nearest trap. The name of the trap is the variable "trap" (which is created as the row number of trap_cat1719_foto)
  #st_drop_geometry() %>%
  left_join(trap_cat1719_foto %>% st_drop_geometry(), by = "trap") %>% # To get the trapID of the whole trap dataset, it joins it by the variable trap
  dplyr::select(Confirmed_Individual, month, Sex, trap_id, trap, dist)

# Hair trap detections

pts_pels_2017 <- dist_nearest(pts_pels_2017, trap_cat1719_pels) %>% ## Here I join to trap_cat1719 that contains hair and mixed (both, this is a bit different than maelis)
  #st_drop_geometry() %>%
  left_join(trap_cat1719_pels %>% st_drop_geometry(), by = "trap") %>%
  dplyr::select(Confirmed_Individual, month, Sex, trap_id, trap, dist)

pts_pels_2018 <- dist_nearest(pts_pels_2018, trap_cat1719_pels) %>% ## Here I join to trap_cat1719 that contains hair and mixed (both, this is a bit different than maelis)
  #st_drop_geometry() %>%
  left_join(trap_cat1719_pels %>% st_drop_geometry(), by = "trap") %>%
  dplyr::select(Confirmed_Individual, month, Sex, trap_id, trap, dist)

pts_pels_2019 <- dist_nearest(pts_pels_2019, trap_cat1719_pels) %>% ## Here I join to trap_cat1719 that contains hair and mixed (both, this is a bit different than maelis)
  #st_drop_geometry() %>%
  left_join(trap_cat1719_pels %>% st_drop_geometry(), by = "trap") %>%
  dplyr::select(Confirmed_Individual, month, Sex, trap_id, trap, dist)

#mapview(trap_cat1719_pels) + mapview(trap_cat1719_foto, col.regions = "green") + 
#  mapview(pts_foto_2017, col.regions = "darkgreen", cex = 2) + 
#  mapview(pts_pels_2017, col.regions = "magenta", cex = 2)
# Esto se podría solapar con las trampas de otros años a ver si se puede rescatar alguna donde haya muchas detecciones


#### 2.3. Join and remove detections further than threshold distance ####

pts_cat_2017 <- rbind(pts_foto_2017,pts_pels_2017) %>% 
  filter(dist < seuil) %>%
  mutate(suivi = "systematic")

pts_cat_2018 <- rbind(pts_foto_2018,pts_pels_2018) %>% 
  filter(dist < seuil) %>%
  mutate(suivi = "systematic") 

pts_cat_2019 <- rbind(pts_foto_2019,pts_pels_2019) %>% 
  filter(dist < seuil) %>%
  mutate(suivi = "systematic")

## ---- Combine France and Catalunya ----
# 1. Sites

sites_2017_t <- sites_2017 %>%
  mutate(pays = "France") %>%
  rbind(trap_cat1719 %>%
          mutate(pays = "Espagne") %>%
          dplyr::select(NOM, site, effort, geometry, trap_id, pays)) %>%
  mutate(trap = row_number()) %>%
  st_transform(map, crs = 32631) %>% # Transform a UTM to run in SCR
  mutate(X = unlist(map(geometry,1)),
         Y = unlist(map(geometry,2))) %>%
  dplyr::select(trap,trap_id,X,Y,site, pays) %>%
  st_drop_geometry()

sites_2018_t <- sites_2018 %>%
  mutate(pays = "France") %>%
  rbind(trap_cat1719 %>%
          mutate(pays = "Espagne") %>%
          dplyr::select(NOM, site, effort, geometry, trap_id, pays)) %>%
  mutate(trap = row_number()) %>%
  st_transform(map, crs = 32631) %>% # Transform a UTM to run in SCR
  mutate(X = unlist(map(geometry,1)),
         Y = unlist(map(geometry,2))) %>%
  dplyr::select(trap,trap_id,X,Y,site, pays) %>%
  st_drop_geometry()

sites_2019_t <- sites_2019 %>%
  mutate(pays = "France") %>%
  rbind(trap_cat1719 %>%
          mutate(pays = "Espagne") %>%
          dplyr::select(NOM, site, effort, geometry, trap_id, pays)) %>%
  mutate(trap = row_number()) %>%
  st_transform(map, crs = 32631) %>% # Transform a UTM to run in SCR
  mutate(X = unlist(map(geometry,1)),
         Y = unlist(map(geometry,2))) %>%
  dplyr::select(trap,trap_id,X,Y,site, pays) %>%
  st_drop_geometry()

# 2. Detections

pts_2017_t <- pts_2017 %>%
  mutate(id  = as.character(Confirmed_Individual)) %>%
  mutate(month = month %>%
           as.character() %>%
           as.numeric()) %>%
  mutate(sex  = as.character(Sex)) %>%
  mutate(trap_id  = as.character(trap_id)) %>%
  dplyr::select(id, month, sex, trap_id) %>%
  rbind(pts_cat_2017 %>%
          st_drop_geometry() %>%
          mutate(id  = as.character(Confirmed_Individual)) %>%
          mutate(month = month %>%
                   as.character() %>%
                   as.numeric()) %>%
          mutate(sex  = as.character(Sex)) %>%
          mutate(trap_id  = as.character(trap_id)) %>%
          dplyr::select(id, month, sex, trap_id)) %>%
  left_join(sites_2017_t, by = "trap_id") %>%
  dplyr::select(id, month, sex, trap, trap_id) %>%
  mutate(sex = as.factor(sex)) 


pts_2018_t <- pts_2018 %>%
  mutate(id  = as.character(Confirmed_Individual)) %>%
  mutate(month = month %>%
           as.character() %>%
           as.numeric()) %>%
  mutate(sex  = as.character(Sex)) %>%
  mutate(trap_id  = as.character(trap_id)) %>%
  dplyr::select(id, month, sex, trap_id) %>%
  rbind(pts_cat_2018 %>%
          st_drop_geometry() %>%
          mutate(id  = as.character(Confirmed_Individual)) %>%
          mutate(month = month %>%
                   as.character() %>%
                   as.numeric()) %>%
          mutate(sex  = as.character(Sex)) %>%
          mutate(trap_id  = as.character(trap_id)) %>%
          dplyr::select(id, month, sex, trap_id)) %>%
  left_join(sites_2018_t, by = "trap_id") %>% # ASP: Esto no sería necesario, pero lo hace MK
  dplyr::select(id, month, sex, trap, trap_id) %>%
  mutate(sex = as.factor(sex)) 

pts_2019_t <- pts_2019 %>%
  mutate(id  = as.character(Confirmed_Individual)) %>%
  mutate(month = month %>%
           as.character() %>%
           as.numeric()) %>%
  mutate(sex  = as.character(Sex)) %>%
  mutate(trap_id  = as.character(trap_id)) %>%
  dplyr::select(id, month, sex, trap_id) %>%
  rbind(pts_cat_2019 %>%
          st_drop_geometry() %>%
          mutate(id  = as.character(Confirmed_Individual)) %>%
          mutate(month = month %>%
                   as.character() %>%
                   as.numeric()) %>%
          mutate(sex  = as.character(Sex)) %>%
          mutate(trap_id  = as.character(trap_id)) %>%
          dplyr::select(id, month, sex, trap_id)) %>%
  left_join(sites_2019_t, by = "trap_id") %>%
  dplyr::select(id, month, sex, trap, trap_id) %>%
  mutate(sex = as.factor(sex)) 

## -------------------------------------------------
##                      2020 - 2021
## ------------------------------------------------- 

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
piege_poils2021 <- st_read("D:/MargSalas/Oso/Datos/Effort_raw/France/2022/Appat_smola_2022.shp") %>%
  distinct(.keep_all = TRUE) %>%
  st_transform(crs = CRS("+proj=utm +zone=31 +datum=WGS84"))

# We keep the ones that fall within the itineraries of 2020 / 2021
# Hair trap on an itinerary = at least 500m from the closest itinerary.
# With this threshold we remove the ones that belong to transects that weren't done this year

seuil <- 500

piege_poils_it_2020 <- dist_nearest(piege_poils2021,itineraires_2020) %>%
  filter(dist < seuil) %>%
  mutate("itineraire" = NA) # Variable trap = transecto al que pertenece la trampa de pelo

#mapview(itineraires_2020) + mapview(piege_poils2021, col.regions = "magenta", cex = 2) + mapview(piege_poils_it_2020, col.regions = "darkgreen", cex = 2)

piege_poils_it_2021 <- dist_nearest(piege_poils2021,itineraires_2021) %>%
  filter(dist < seuil) %>%
  mutate("itineraire" = NA) # Variable trap = transecto al que pertenece la trampa de pelo

#mapview(itineraires_2021) + mapview(piege_poils2021, col.regions = "magenta", cex = 2) + mapview(piege_poils_it_2021, col.regions = "darkgreen", cex = 2)
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
os2021_FR<- openxlsx::read.xlsx('Seguiment_Ossos_Pirineus_1996_2021_taula_final.xlsx') %>% 
  filter(Year %in% c(2020,2021) & Confirmed_Individual != "Indetermined" & Country == "France") # Removed indetermined

os2021_FR<- os2021_FR%>% mutate(date = as_date(os2021_FR$Date, format = "%d/%m/%Y"),
                                month = month(date)) %>%
  filter(month < 12, month > 4) %>% # 7 months form may to november
  filter(!is.na(x_long)) # Remove observations without coordinates


## Keep systematic data
dat2_2021_FR <- os2021_FR[which(os2021_FR$Method %in% c("Sampling_station", "Transect") &
                                  os2021_FR$Obs_type %in% c("Photo","Photo/Video", "Hair", "Video")), ]

# We distinguish the 3 methods 
dat2_2021_FR <- dat2_2021_FR %>%
  mutate(method = "itineraire") # Itinerary 

dat2_2021_FR[which(dat2_2021_FR$Obs_type == "Photo" | dat2_2021_FR$Obs_type == "Photo/Video" | dat2_2021_FR$Obs_type == "Video"),"method"] <- "piege photos" # Camera trap

# We don't know if we should take all hair observations, or all hair observations belonging
# to a hair trap in 2020-2021, because we do not know if it's reliable

dat2_2021_FR[which(dat2_2021_FR$Obs_type == "Hair"),"method"] <- "piege poils" # Hair trap

#### 3.1. Separate the detections of each type  ####

# Détections issues de pièges photos (photographie automatiques, vidéo automatique)
pts_photo2021 <- dat2_2021_FR %>%
  filter(method =="piege photos") %>%
  st_as_sf(coords = c("x_long", "y_lat"),
           crs = 4326)

pts_photo_2020 <- pts_photo2021 %>%
  filter(Year == 2020) 

pts_photo_2021 <- pts_photo2021 %>%
  filter(Year == 2021) 

# Détections issues de pièges à poils 
pts_poils <- dat2_2021_FR %>%
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
os_add_trap <- openxlsx::read.xlsx('Seguiment_Ossos_Pirineus_1996_2021_taula_final.xlsx') %>%
  filter(Year %in% c(2020,2021) & Confirmed_Individual != "Indetermined" & Country == "Spain") %>%
  filter(Confirmed_Individual %in% c("Bonabe", "New20-14") & Date %in% c("01/08/2020", "19/08/2021") & Obs_type %in% c("Hair", "Video"))

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

#mapview(trap_cat_pels_2020) + mapview(trap_cat_foto_2020, col.regions = "red")

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

#mapview(trap_cat_pels_2021) + mapview(trap_cat_foto_2021, col.regions = "red")

## ---- 2. Detections ----

setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
os2021_SP <- openxlsx::read.xlsx('Seguiment_Ossos_Pirineus_1996_2021_taula_final.xlsx') %>%
  filter(Year %in%  c(2020,2021) & Confirmed_Individual != "Indetermined" & Country == "Spain") 
os2021_SP <- os2021_SP %>%
  mutate(date = as_date(os2021_SP$Date, format = "%d/%m/%Y"),
         month = month(date))

os2021_SP[which(is.na(os2021_SP$month)),]
os2021_SP$month[350] <- 4 # Correct mannually because date was not exact and we only had month


os2021_SP <- os2021_SP %>%
  filter(month < 12, month > 4)  # 7 months form may to november 


## Keep systematic data
dat_cat2021_syst <- os2021_SP[which(os2021_SP$Method %in% c("Sampling_station", "Transect") &
                                      os2021_SP$Obs_type %in% c("Photo","Photo/Video", "Hair", "Video")), ]

#### 2.1. Separate the detections of each type  ####

## 2020

dat_cat_syst_2020 <- dat_cat2021_syst %>% 
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

dat_cat_syst_2021 <- dat_cat2021_syst %>% 
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

#mapview(trap_cat_pels_2020) + mapview(trap_cat_foto_2020, col.regions = "green") + 
#  mapview(pts_foto_2020, col.regions = "darkgreen", cex = 2) + 
#  mapview(pts_pels_2020, col.regions = "magenta", cex = 2)

##2021

pts_foto_2021 <- dist_nearest(pts_foto_2021, trap_cat_foto_2021) %>%
  #st_drop_geometry() %>%
  left_join(trap_cat_foto_2021 %>% st_drop_geometry(), by = "trap") %>%
  dplyr::select(Confirmed_Individual, month, Sex, trap_id, trap, dist)

pts_pels_2021 <- dist_nearest(pts_pels_2021, trap_cat_2021) %>% ## Here I join to trap_cat that contains hair and mixed (both, this is a bit different than maelis)
  #st_drop_geometry() %>%
  left_join(trap_cat_2021 %>% st_drop_geometry(), by = "trap") %>%
  dplyr::select(Confirmed_Individual, month, Sex, trap_id, trap, dist)

#mapview(trap_cat_pels_2021) + mapview(trap_cat_foto_2021, col.regions = "green") + 
#  mapview(pts_foto_2021, col.regions = "darkgreen", cex = 2) + 
#  mapview(pts_pels_2021, col.regions = "magenta", cex = 2)

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

#mapview(trap_cat) + mapview(pts_2020, col.regions = "green", cex = 2)

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

## -------------------------------------------------
##                      Create EDF and TDF
## ------------------------------------------------- 

# I will need to add the effort before this step, but for now I leave it like this

#EDF
edf <- pts_2017_t %>% 
  mutate(occ = month %>% as.factor() %>% as.numeric) %>% 
  mutate(session = 1) %>%
  mutate(ind = id) %>%
  dplyr::select(session, ind, occ, trap, sex) %>%
  rbind(pts_2018_t %>% 
          mutate(occ = month %>% as.factor() %>% as.numeric) %>% 
          mutate(session = 2) %>%
          mutate(ind = id) %>%
          dplyr::select(session, ind, occ, trap, sex)) %>%
  rbind(pts_2019_t %>% 
          mutate(occ = month %>% as.factor() %>% as.numeric) %>% 
          mutate(session = 3) %>%
          mutate(ind = id) %>%
          dplyr::select(session, ind, occ, trap, sex)) %>%
  rbind(pts_2020_t %>% 
          mutate(occ = month %>% as.factor() %>% as.numeric) %>% 
          mutate(session = 4) %>%
          mutate(ind = id) %>%
          dplyr::select(session, ind, occ, trap, sex)) %>%
  rbind(pts_2021_t %>% 
          mutate(occ = month %>% as.factor() %>% as.numeric) %>% 
          mutate(session = 5) %>%
          mutate(ind = id) %>%
          dplyr::select(session, ind, occ, trap, sex))

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")
save(edf, file = "edf1721_nocubsmother.RData")

#TDF

tdf2017 <- sites_2017_t %>% 
  dplyr::select(trap,X,Y,site, pays)

tdf2018 <- sites_2018_t %>% 
  dplyr::select(trap,X,Y,site, pays)

tdf2019 <- sites_2019_t %>% 
  dplyr::select(trap,X,Y,site, pays)

tdf2020 <- sites_2020_t %>% 
  dplyr::select(trap,X,Y,site, pays)

tdf2021 <- sites_2021_t %>% 
  dplyr::select(trap,X,Y,site, pays)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")

save(tdf2017, file = "tdf2017_nocubsmother.RData")
save(tdf2018, file = "tdf2018_nocubsmother.RData")
save(tdf2019, file = "tdf2019_nocubsmother.RData")
save(tdf2020, file = "tdf2020_nocubsmother.RData")
save(tdf2021, file = "tdf2021_nocubsmother.RData")

## -------------------------------------------------
##                      DATA CHECK
## ------------------------------------------------- 

# Protocol to check data set up:
# - Make sure that all objects in script are UNIQUE (Indexed by year/region) --> OK
# - Take 2 observations from edf each session and track it back

# More thorough protocol (maybe when I have the final checked dataset)
# Check individual cases, to see if several observations are repeated in one place
# over years but do not have a trap near -> maybe create a trap there? Ask?



