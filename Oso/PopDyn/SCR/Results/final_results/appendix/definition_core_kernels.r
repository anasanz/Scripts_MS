
## -------------------------------------------------
##                 CORE AREA
## ------------------------------------------------- 

rm(list = ls())

library(tidyverse)
library(sf)
library(rgdal)
library(mapview)
library(lubridate)
library(adehabitatHR)
library(sp)
library(raster)
library(ggplot2)


setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2021")
os <- read.csv("Data_os_96_21_cubLocations.csv", header = TRUE, row.names = NULL)  %>% 
  filter(Confirmed_Individual != "Indetermined") %>%
  st_as_sf(coords = c("x_long","y_lat"), crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

os <- os[,-which(colnames(os) %in% c("ID_obs"))] # So that the script runs, because this column was added later and it doesn't fit otherwise



map <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/clip_pyros2.shp")


# Core area: Posiciones de hembras reproductoras (eliminar Franska, Palouma, Sarousse)
# -> Tiene poco sentido incluir las localizaciones de hembras reintroducidas no reproductoras? 
#    Sobre todo Franska y Palouma que murieron el mismo año
# -> Quitando estas, todas son hijas de Hvala y Caramelles (hija de Mellba)
# -> Posible delimitación: Distribución de hembras reintroducidas reproductivas,el año de su primera reproducción (asumiendo que están asentadas)
#    * Núcleo central: Distribición Ziva, Mellba y Hvala (tiene sentido porque antes de 2010 no hay seguimiento sistemático)
#    * Núcleo oriental: Distribución Sorita (Claverina no se ha reproducido)

core <- os %>%
  filter(Year %in% c(1997,1998, 2007,2008, 2019, 2020) 
         & Confirmed_Individual %in% c("Mellba", "Ziva", "Hvala", "Sorita") 
         & With_cubs_estimated %in% c("<6month", "1", "2"))


## ---- Calculate kernel utilization distribution and MCP ----

setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2021")
os_spatial <- read.csv("Data_os_96_21_cubLocations.csv", header = TRUE, row.names = NULL)  %>% 
  filter(Confirmed_Individual != "Indetermined") 

os_spatial <- os_spatial[,-which(colnames(os_spatial) %in% c("ID_obs"))]

core_spatial <- os_spatial %>%
  filter(Year %in% c(1997,1998, 2007,2008, 2019, 2020) 
         & Confirmed_Individual %in% c("Mellba", "Ziva", "Hvala", "Sorita") 
         & With_cubs_estimated %in% c("<6month", "1", "2"))

coordinates(core_spatial) <- core_spatial[,c("x_long", "y_lat")]
proj4string(core_spatial) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# 1. Calculate kernel for each individual and join

core_kernel <- core_spatial[,2]

kref <- kernelUD(core_kernel, h="href")
kernel.poly <- getverticeshr(kref, percent = 95) 

core_area_kernelpol <- aggregate(kernel.poly, dissolve = TRUE)

# 2. Calculate kernel for each population nucleus

core_spatial$Nucleus <- ifelse(core_spatial$Confirmed_Individual %in% c("Mellba", "Ziva", "Hvala"), 1, 2)

core_kernel2 <- core_spatial[,24]

kref2 <- kernelUD(core_kernel2, h="href")
kernel.poly2 <- getverticeshr(kref2, percent = 95) 

core_area_kernelpol2 <- aggregate(kernel.poly2, dissolve = TRUE)


# Try also to limit the % to 90 in the occidental core, because it is only one female
# and is over-represented?

core_kernel_n1 <- core_spatial[core_spatial$Nucleus == 1,24]
core_kernel_n2 <- core_spatial[core_spatial$Nucleus == 2,24]

kref <- kernelUD(core_kernel_n1, h="href")
kernel.poly <- getverticeshr(kref, percent = 95) 
core_area_kernelpol_n1 <- aggregate(kernel.poly, dissolve = TRUE)

kref <- kernelUD(core_kernel_n2, h="href")
kernel.poly <- getverticeshr(kref, percent = 80) 
core_area_kernelpol_n2 <- aggregate(kernel.poly, dissolve = TRUE)

########## Plot for methods

sa <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/clip_pyros2_WGS84_31N_all.shp")
eur <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/esp_fr_2.shp") %>%
  st_transform(dpts, crs = crs(sa))
esp <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/ESP_adm/ESP_adm0.shp") %>%
  st_transform(dpts, crs = crs(sa))
fr <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/FRA_adm/FRA_adm0.shp") %>%
  st_transform(dpts, crs = crs(sa))

core <- core %>% st_transform(32631)
core_area_kernelpol_n1 <- spTransform(core_area_kernelpol_n1,CRS("+proj=utm +zone=31 +datum=WGS84 +units=m +no_defs +type=crs"))
core_area_kernelpol_n2 <- spTransform(core_area_kernelpol_n2,CRS("+proj=utm +zone=31 +datum=WGS84 +units=m +no_defs +type=crs"))

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Results_section/Plots")
pdf("SI_core_kernel_95n1_80n2.pdf",7,7)

plot(st_geometry(sa), col = "beige", border = "beige",  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000)) # To add andorra with same color (not present in eur)
plot(st_geometry(eur), col = "beige", border = "beige",  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000), add = TRUE)
plot(st_geometry(esp), col = "beige", border = "white",  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000), add = TRUE)
plot(st_geometry(fr), col = "beige", border = "white",  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000), add = TRUE)
plot(core, pch = 19, col = "brown", add = TRUE)
plot(core_area_kernelpol_n1, pch = 19, col = adjustcolor("yellow", alpha.f = 0.7), add = TRUE)
plot(core_area_kernelpol_n2, pch = 19, col = adjustcolor("yellow", alpha.f = 0.6), add = TRUE)

dev.off()

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Results_section/Plots")
pdf("SI_core_kernel_95n1_80n2_zoom.pdf",5,5)

plot(st_geometry(sa), col = "beige", border = "beige",  xlim = c(st_bbox(sa)[1] + 10000, st_bbox(sa)[3] - 10000),  ylim = c(st_bbox(sa)[2] + 50000, st_bbox(sa)[4] - 50000)) # To add andorra with same color (not present in eur)
plot(st_geometry(eur), col = "beige", border = "beige",  xlim = c(st_bbox(sa)[1]+10000, st_bbox(sa)[3]-10000),  ylim = c(st_bbox(sa)[2] + 50000 , st_bbox(sa)[4] - 50000), add = TRUE)
plot(st_geometry(esp), col = "beige", border = "white",  xlim = c(st_bbox(sa)[1]+10000, st_bbox(sa)[3]-10000),  ylim = c(st_bbox(sa)[2] + 50000 , st_bbox(sa)[4] - 50000), add = TRUE)
plot(st_geometry(fr), col = "beige", border = "white",  xlim = c(st_bbox(sa)[1]+10000, st_bbox(sa)[3]-10000),  ylim = c(st_bbox(sa)[2] + 50000 , st_bbox(sa)[4] - 50000), add = TRUE)
plot(core, pch = 19, cex = 0.2, col = "brown", add = TRUE)
plot(core_area_kernelpol_n1, pch = 19, col = adjustcolor("yellow", alpha.f = 0.7), add = TRUE)
plot(core_area_kernelpol_n2, pch = 19, col = adjustcolor("yellow", alpha.f = 0.6), add = TRUE)

dev.off()
