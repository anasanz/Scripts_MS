
## -------------------------------------------------
##      Join detections - traps France and Spain    
##                      2017-2019
## ------------------------------------------------- 

## Difference from Maelis: Using Spanish database, as it includes the locations of
## the cubs that were detected with the mother

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
piege_poils <- rgdal::readOGR("D:/MargSalas/Oso/Datos/Effort_raw/France/Raw/Raw/shp/Traps2017_2019/PiegesPoils_2017_2019.shp", use_iconv=TRUE, encoding = "UTF-8")

# On formate les pièges à poils en spatial
piege_poils <- piege_poils %>%
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
piege_poils <- rbind(piege_poils,haute_arriege, arrioutort) %>%
  distinct(.keep_all = TRUE) 

##### b.2) Distinguish hair trap that are on an itinerary ( <500 m, because they are more frequently visited) and outside
# And join it with itinerary id

seuil <- 500

piege_poils_it <- dist_nearest(piege_poils,itineraires) %>%
  filter(dist < seuil) %>%
  mutate("itineraire" = NA) 

# Pour chaque piège a poils on remet le bon numéro correspondant à l'itinéraire associé
for (i in 1 : length(piege_poils_it$trap)) {
  piege_poils_it$itineraire[i] <- itineraires$id[piege_poils_it$trap[i]]
  piege_poils_it$NOM[i] <- itineraires$NOM[piege_poils_it$trap[i]]
}

# On enlève la variable trap pour ne pas se mélanger par la suite 
piege_poils_it <- piege_poils_it  %>%
  dplyr::select(geometry,NOM,itineraire)

##### b.3) One hair trap file by year (remove traps that are on an itinerary that haven't been made)

## 2017 
piege_poils_it_2017 <- piege_poils_it %>% 
  filter(NOM !="Bern" & NOM !="Saoubette" & NOM != "Escouloubre" & NOM != "Formiguere" & NOM != "Orlu" & NOM != "Mijanes1" & NOM != "Majanes2" & NOM != "Cazaux-Larboust" & NOM != "St Collombe" & NOM != "Rouze") 

## 2018 
piege_poils_it_2018 <- piege_poils_it %>% 
  filter(NOM !="Saoubette")

## 2019 
piege_poils_it_2019 <- piege_poils_it %>%
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
  filter(dist >= seuil)

piege_photos_syst_2019 <- dist_nearest(piege_photos_syst_2019, piege_poils)

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
os <- openxlsx::read.xlsx('Seguiment_Ossos_Pirineus_1996_2021_taula_final_2_cubLocations.xlsx') %>% 
  filter(Year %in% c(2017,2018,2019) & Confirmed_Individual != "Indetermined" & Country == "France") # Removed indetermined
  
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

nrow(dat2[which(dat2$Obs_type == "Hair"),]) - nrow(dat2[which(dat2$Obs_type == "Hair" & dat2$Hair_trap == 1),]) # We remove 36 that are hair 
dat2[which(dat2$Obs_type == "Hair" & dat2$Hair_trap == 1),"method"] <- "piege poils" # Hair trap

#### 3.1. Separate the detections of each type  ####

# Détections issues de pièges photos (photographie automatiques, vidéo automatique)
pts_photo <- dat2 %>%
  filter(method =="piege photos") %>%
  st_as_sf(coords = c("x_long", "y_lat"),
           crs = 4326)

pts_photo_2017 <- pts_photo %>%
  filter(Year == 2017) 

pts_photo_2018 <- pts_photo %>%
  filter(Year == 2018) 

pts_photo_2019 <- pts_photo %>%
  filter(Year == 2019)

# Détections issues de pièges à poils 
pts_poils <- dat2 %>%
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
  # mapview(pts_photo_2017, cex = 2) + mapview(sites_2017, col.regions = "red", cex = 2)
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
  #mapview(pts_poils_2017, cex = 2) + mapview(sites_2017, col.regions = "red", cex = 2)
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

trap_cat <- readxl::read_xlsx("D:/MargSalas/Oso/Datos/Effort_raw/France/Raw/Raw/UTM trampas pelos y camaras.xlsx") %>% 
  janitor::clean_names() %>%
  filter(!is.na(xutm)) %>%
  st_as_sf(coords = c("xutm","yutm"), 
           crs = CRS("+proj=utm +zone=31 +datum=WGS84")) %>% 
  st_transform(crs = 4326) %>%
  rename("NOM" = "itinerari") %>%
  mutate(effort = "it") %>%
  mutate(site = "hair") # poils seul

## ASP: Names all "hair" and then according to foto_video or pels columns assings to camera or both
trap_cat[which(!is.na(trap_cat$foto_video)),"site"] <- "both" # poils et photo
trap_cat[which(is.na(trap_cat$pels)),"site"] <- "camera" # photo seul

trap_cat <- trap_cat[-which(duplicated(trap_cat$geometry)),] # ASP: MK did not do this step, but I do it
                  # because they are duplicates and it gives problems when assigning detections to traps

# On enlève les pièges photos seul car ils ne permettent pas l'identification des ours 
trap_cat <- trap_cat %>%
  filter(site != "camera")  %>%
  arrange(by = site) %>% # I arrange it so that numbers make more sense (fotos first, then both)
  mutate(trap_id = paste(row_number(), "c"))

# On fait un sous tableau pour les pièges à poils 
trap_cat_pels <- trap_cat %>%
  mutate(trap = row_number())

#### b) Camera traps ####

trap_cat_foto <- trap_cat %>%
  filter(foto_video != is.na(foto_video)) %>% 
  mutate(trap = row_number()) 

## ---- 2. Detections ----

setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
os <- openxlsx::read.xlsx('Seguiment_Ossos_Pirineus_1996_2021_taula_final_2_cubLocations.xlsx') %>% 
  filter(Year %in% c(2017,2018,2019) & Confirmed_Individual != "Indetermined" & Country == "Spain") # Removed indetermined

os <- os %>% mutate(date = as_date(os$Date_register, format = "%d/%m/%Y"),
                    month = month(date)) %>%
  filter(month < 12, month > 4) %>% # 7 months form may to november
  filter(!is.na(x_long)) # Remove observations without coordinates


## Keep systematic data
dat_cat_syst <- os[which(os$Method %in% c("Sampling_station", "Transect") &
                   os$Obs_type %in% c("Photo","Photo/Video", "Hair", "Video")), ]

#### 2.1. Separate the detections of each type  ####

# Camera trap detections
pts_foto <- dat_cat_syst %>% 
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
pts_pels <- dat_cat_syst %>% 
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

pts_foto_2017 <- dist_nearest(pts_foto_2017, trap_cat_foto) %>% ## ASP: It assings to the nearest trap. The name of the trap is the variable "trap" (which is created as the row number of trap_cat_foto)
  #st_drop_geometry() %>%
  left_join(trap_cat_foto %>% st_drop_geometry(), by = "trap") %>% # To get the trapID of the whole trap dataset, it joins it by the variable trap
  dplyr::select(Confirmed_Individual, month, Sex, trap_id, trap, dist)

pts_foto_2018 <- dist_nearest(pts_foto_2018, trap_cat_foto) %>% ## ASP: It assings to the nearest trap. The name of the trap is the variable "trap" (which is created as the row number of trap_cat_foto)
  #st_drop_geometry() %>%
  left_join(trap_cat_foto %>% st_drop_geometry(), by = "trap") %>% # To get the trapID of the whole trap dataset, it joins it by the variable trap
  dplyr::select(Confirmed_Individual, month, Sex, trap_id, trap, dist)

pts_foto_2019 <- dist_nearest(pts_foto_2019, trap_cat_foto) %>% ## ASP: It assings to the nearest trap. The name of the trap is the variable "trap" (which is created as the row number of trap_cat_foto)
  #st_drop_geometry() %>%
  left_join(trap_cat_foto %>% st_drop_geometry(), by = "trap") %>% # To get the trapID of the whole trap dataset, it joins it by the variable trap
  dplyr::select(Confirmed_Individual, month, Sex, trap_id, trap, dist)

# Hair trap detections

pts_pels_2017 <- dist_nearest(pts_pels_2017, trap_cat_pels) %>% ## Here I join to trap_cat that contains hair and mixed (both, this is a bit different than maelis)
  #st_drop_geometry() %>%
  left_join(trap_cat_pels %>% st_drop_geometry(), by = "trap") %>%
  dplyr::select(Confirmed_Individual, month, Sex, trap_id, trap, dist)

pts_pels_2018 <- dist_nearest(pts_pels_2018, trap_cat_pels) %>% ## Here I join to trap_cat that contains hair and mixed (both, this is a bit different than maelis)
  #st_drop_geometry() %>%
  left_join(trap_cat_pels %>% st_drop_geometry(), by = "trap") %>%
  dplyr::select(Confirmed_Individual, month, Sex, trap_id, trap, dist)

pts_pels_2019 <- dist_nearest(pts_pels_2019, trap_cat_pels) %>% ## Here I join to trap_cat that contains hair and mixed (both, this is a bit different than maelis)
  #st_drop_geometry() %>%
  left_join(trap_cat_pels %>% st_drop_geometry(), by = "trap") %>%
  dplyr::select(Confirmed_Individual, month, Sex, trap_id, trap, dist)

mapview(trap_cat_pels) + mapview(trap_cat_foto, col.regions = "green") + 
  mapview(pts_foto_2017, col.regions = "darkgreen", cex = 2) + 
  mapview(pts_pels_2017, col.regions = "magenta", cex = 2)
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
  rbind(trap_cat %>%
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
  rbind(trap_cat %>%
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
  rbind(trap_cat %>%
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

## ---- Create EDF and TDF ----
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
          dplyr::select(session, ind, occ, trap, sex))

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_final_1719")
save(edf, file = "edf1719.RData")

#TDF: There are a few traps less than MK (5-6 less) because I removed the duplicates

tdf2017 <- sites_2017_t %>% 
  dplyr::select(trap,X,Y,site, pays)

tdf2018 <- sites_2018_t %>% 
  dplyr::select(trap,X,Y,site, pays)

tdf2019 <- sites_2019_t %>% 
  dplyr::select(trap,X,Y,site, pays)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_final_1719")

#save(tdf2017, file = "tdf2017.RData")
#save(tdf2018, file = "tdf2018.RData")
#save(tdf2019, file = "tdf2019.RData")


#### CHECK THAT IS RIGHT --> IT SEEMS OKay ####

# 2017 

# To be able to plot, to convert it in spatial
pts_photo_2017 <- pts_photo %>%
  filter(Year == 2017) 
pts_poils_2017 <-  pts_poils %>%
  filter(Year == 2017)

mapview(trap_cat) + mapview(sites_2017, col.regions = "red") + 
  mapview(pts_cat_2017, col.regions = "green", cex = 2) + # Detections catalonia
  mapview(pts_photo_2017, col.regions = "yellow", cex = 2) + mapview(pts_poils_2017, col.regions = "orange", cex = 2) + # Detections France
  mapview(pts_poils_2017[31,], col.regions = "darkblue", cex = 2)
# Photo Nere France 06/05/2017 (obs 3 yellow) -> Nearest trap 16 f (trap 16)
edf[which(edf$session == 1 & edf$ind == "Nere" & edf$occ == 1),]
sites_2017_t[sites_2017_t$trap == 16,]

#↑ Hair Rodri France 12/11/2017 (obs 43 orange) -> Nearest trap 16 F
edf[which(edf$session == 1 & edf$ind == "Rodri" & edf$occ == 7),]
sites_2017_t[sites_2017_t$trap == 16,]

#♦ Hair Nere France 11/08/2017 (obs 31 orange) -> Nearest trap 161 f (hair)
edf[which(edf$session == 1 & edf$ind == "Nere" & edf$occ == 4),]
sites_2017_t[sites_2017_t$trap == 161,]

#  Hair Hvala (8/2017) (obs 45 green) -> Nearest trap 192 c (192 trap_cat)
edf[which(edf$session == 1 & edf$ind == "Hvala" & edf$occ == 4),]
sites_2017_t[sites_2017_t$trap_id == "192 c",]

# Both Goiat (7/2017) obs 22 green -> Nearest trap 45 c (45 trap_cat)
edf[which(edf$session == 1 & edf$ind == "Goiat" & edf$occ == 3),]
sites_2017_t[sites_2017_t$trap_id == "45 c",]

# 2019 

# To be able to plot, to convert it in spatial
pts_photo_2019 <- pts_photo %>%
  filter(Year == 2019) 
pts_poils_2019 <-  pts_poils %>%
  filter(Year == 2019)

mapview(trap_cat) + mapview(sites_2019, col.regions = "red") + 
  mapview(pts_cat_2019, col.regions = "green", cex = 2) + # Detections catalonia
  mapview(pts_photo_2019, col.regions = "yellow", cex = 2) + mapview(pts_poils_2019, col.regions = "orange", cex = 2) # Detections France

# Photo Canellito France 07/09/2019 (obs 15 yellow) .> Nearest trap 49 f (49)
edf[which(edf$session == 3 & edf$ind == "Canellito" & edf$occ == 5),]
sites_2019_t[sites_2019_t$trap == 49,]

# Hair Cirera cat 7/2019 (obs 9 green) -> Nearest trap 196 c
edf[which(edf$session == 3 & edf$ind == "Cirera" & edf$occ == 3),]
sites_2019_t[sites_2019_t$trap == 712,]

#### COMPARE WITH TDF MK ####

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Effort/France/Data/Systematic")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Effort/France/Data/Systematic")

#---- 1. LOAD THE DETECTION DATA ---- 
edf_new <- edf
tdf2017_new <- tdf2017 %>%
  st_as_sf(coords = c("X", "Y"),
           crs = 32631)
tdf2018_new <- tdf2018 %>%
  st_as_sf(coords = c("X", "Y"),
           crs = 32631) %>%
  st_as_sf(coords = c("X", "Y"),
           crs = 32631)
tdf2019_new <- tdf2019 %>%
  st_as_sf(coords = c("X", "Y"),
           crs = 32631)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Effort/France/Data/Systematic")
load("edf2017_2019_fr.RData")
load("tdf2017_fr.RData")
load("tdf2018_fr.RData")
load("tdf2019_fr.RData")

tdf2017 <- tdf2017 %>%
  st_as_sf(coords = c("X", "Y"),
           crs = 32631)

tdf2018 <- tdf2018 %>%
  st_as_sf(coords = c("X", "Y"),
           crs = 32631)

tdf2018 <- tdf2018 %>%
  st_as_sf(coords = c("X", "Y"),
           crs = 32631)

# It makes sense to have more detections in edf, but in tdf??? I need to check this

mapview(tdf2017_new, col.regions = "red", cex = 2) + mapview(tdf2017, col.regions = "blue", cex = 2)
# In the map they look the same

# I have checked step by step building the script (hers and mine), and the result is exactly
# the same (as is logical). The different number of traps comes from the file ("data.RData"), 
# she has around 30 traps less each year.
# so one day I will need to ask Maelis why and if she did something on that file or deleted traps.
