## -------------------------------------------------
##           Join spanish detections - traps     
##                      YEAR 2020
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

mapview(trap_nav_2020)

## Trampas 2021
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

trap_cat_foto_2020 <- trap_cat %>% 
  filter(site == "both") # VERY IMPORTANT that "both" is located first in trap_cat, so that trap number is the same as row number. Necesary for dist_nearest

trap_cat_pels_2020 <- trap_cat %>% # This file is not used to join with detections (I join with all hair and both)
  filter(site == "hair")  # ONLY for ploting (trap is not the row number)

#mapview(trap_cat_pels) + mapview(trap_cat_foto, col.regions = "red")

## ---- Load detections ----

# I take only confirmed individuals (we lose 18 detections)

setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
os <- openxlsx::read.xlsx('Seguiment_Ossos_Pirineus_1996_2021_taula_final_cubLocations.xlsx') %>% 
  filter(Year == 2020 & Confirmed_Individual != "Indetermined" & Country == "Spain") # 290 confirmed, 308 probable
os <- os %>%
  mutate(date = as_date(os$Date_register, format = "%d/%m/%Y"),
         month = month(date))
os$month[307] <- 6 # Correct mannually because date was not exact and we only had month
os$month[308] <- 4
os <- os %>%
  filter(month < 12, month > 4)  # 7 months form may to november 


## Keep systematic data
dat_cat_syst <- os[which(os$Method %in% c("Sampling_station", "Transect") &
                 os$Obs_type %in% c("Photo","Photo/Video", "Hair", "Video")), ]

pts_foto_2020 <- dat_cat_syst %>% 
  st_as_sf(coords = c("x_long","y_lat"), 
           crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) %>%
  st_transform(crs = 32631) %>%
  filter(Obs_type %in% c("Photo","Photo/Video", "Video"))

pts_pels_2020 <- dat_cat_syst %>% 
  st_as_sf(coords = c("x_long","y_lat"), 
           crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) %>%
  st_transform(crs = 32631) %>%
  filter(Obs_type %in% c("Hair"))


## ---- Association detection - nearest trap ----

pts_foto_2020 <- dist_nearest(pts_foto_2020, trap_cat_foto) %>%
  #st_drop_geometry() %>%
  left_join(trap_cat_foto %>% st_drop_geometry(), by = "trap") %>%
  dplyr::select(Confirmed_Individual, month, Sex, trap_id, trap, dist)

pts_pels_2020 <- dist_nearest(pts_pels_2020, trap_cat) %>% ## Here I join to trap_cat that contains hair and mixed (both, this is a bit different than maelis)
  #st_drop_geometry() %>%
  left_join(trap_cat %>% st_drop_geometry(), by = "trap") %>%
  dplyr::select(Confirmed_Individual, month, Sex, trap_id, trap, dist)

mapview(trap_cat_pels) + mapview(trap_cat_foto, col.regions = "green") + 
  mapview(pts_foto_2020, col.regions = "darkgreen", cex = 2) + 
  mapview(pts_pels_2020, col.regions = "magenta", cex = 2)

# Here trap_id is = trap. I don't know why in MK is different, but this is the only way that works for me.
# Because every year will have a different trap set, this works for 2020

## With a threshold of 500 m there are few saved

# Cases that could be worth checking with santi
# There is an observation but there is no trap associated, it is VERY far from any trap and the coordinates are good
## Esmolet 15/06, Pepito 12/08, New 18-03 05/07 (sist auto?there is even a camera?)
# Should we include a trap in this ones??

## ---- Join and save data ----
seuil <- 500 #♣ Threshold distance, arbitrary (same as MK)

# Format traps
trap_2020 <- trap_cat %>%
  rename("NOM" = "codi_tr") %>%
  mutate(pays = "Espagne") %>%
  mutate(suivi = "systematic") %>%
  select(NOM, site, trap_id, pays, trap, suivi, geometry)
  
# Combine and format detections
pts_2020 <- rbind(pts_foto_2020,pts_pels_2020) %>%
  rename("id" = "Confirmed_Individual") %>%
  rename("sex" = "Sex") %>%
  mutate(suivi = "systematic") %>%
  filter(dist < seuil) %>%
  select(-dist)

mapview(trap_cat) + mapview(pts_2020, col.regions = "green", cex = 2)


dataSpain20 <- list(det = pts_2020,
                    traps = trap_2020)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Effort/Spain/Data")
save(dataSpain20, file = "dataSpain20.RData")
