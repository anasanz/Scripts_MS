
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

## ---- France ----
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




## ---- 3. Detections ----

setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
os <- openxlsx::read.xlsx('Seguiment_Ossos_Pirineus_1996_2021_taula_final_cubLocations.xlsx') %>% 
  filter(Year %in% c(2017,2018,2019) & Confirmed_Individual != "Indetermined" & Country == "France") # Removed indetermined
  
os <- os %>% mutate(date = as_date(os$Date_register, format = "%d/%m/%Y"),
         month = month(date)) %>%
  filter(month < 12, month > 4) %>% # 7 months form may to november
  filter(!is.na(x_long)) # Remove observations without coordinates


## Keep systematic data
dat2 <- os[which(os$Method %in% c("Sampling_station", "Transect") &
                           os$Obs_type %in% c("Photo","Photo/Video", "Hair", "Video")), ]

# CHECK this dataset with dataset Maelis: 
# The biggest difference is that she only keeps the hair samples that come from traps (not spontaneous)
# Todos los que vienen de transecto son espontáneos, y los que vienen de estaciones son trampas de pelo

setwd("D:/MargSalas/Oso/Datos/Data_raw")
df_st <- openxlsx::read.xlsx('Data_Maelis_1719_Syst.xlsx') %>%
  filter(method == "piege poils")
df_trans <- openxlsx::read.xlsx('Data_Maelis_1719_Syst.xlsx') %>%
  filter(method == "itineraire")

# Los que sean pelo en las localicaciones donde Maelis tiene que es una trampa, pongo que vienen de una trampa

filt <- dat2[which(dat2$Obs_type == "Hair" & 
                  dat2$Database %in% c("Seguiment França Ossos 2010-2017.xls", "Seguiment França Ossos 2016-2020.xls")), ]
filt_trans <- filt[which(filt$Method == "Transect"), ]
filt_st <- filt[which(filt$Method == "Sampling_station"), ]

unique(df_st$x) %in% unique(filt_st$X) # Comprobando el primero (excel) me doy cuenta de que no coincide porque
# aunque el tipo de método sea "itineraire", lo ponen como trampa de pelo si la observación viene de una
# trampa de pelo

dat2$hairtrap <- 0
dat2$hairtrap[which(dat2$Obs_type == "Hair" & 
             dat2$Database %in% c("Seguiment França Ossos 2010-2017.xls", "Seguiment França Ossos 2016-2020.xls") &
               dat2$X %in% df_st$x)] <- 1

# 









