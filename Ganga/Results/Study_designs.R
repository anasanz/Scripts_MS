

## -------------------------------------------------
##                 Plot study designs
## ------------------------------------------------- 

rm(list=ls())

library(sf)
library(mapview)
library(ggmap)
library(ggplot2)
library(ggmap)
library(gridExtra)
library(rgdal)
library(ggpubr)
library(dplyr)
library(geosphere)
library("viridis")
library(RColorBrewer)

devtools::install_github('oswaldosantos/ggsn')

## ---- Study area delimitation ----
sa <- st_read("D:/MargSalas/Ganga/Data/GIS/Study_area.shp")

## ---- Transects ----

#tr_farm2022 <- st_read("D:/MargSalas/Ganga/Data/GIS/Trans_Ganga_FarmDindis_2022.shp") # These are the analyzed transects in hds model
# Estos no, son en los que la especie ha estado ausente en los ultimos años (categoria pres_abs 0) pero si que los metimos al final

tr_farm <- st_read("D:/MargSalas/Ganga/Data/GIS/Trans_Ganga_FarmDindis.shp") # These are the analyzed transects in hds model
tr_spec <- st_read("D:/MargSalas/Ganga/Data/GIS/Trans_Ganga_500_Specific.shp") # These are the analyzed transects in hds model

tr_spec <- st_transform(tr_spec, crs = st_crs(tr_farm))

bbox <- c(st_bbox(tr_farm)[1], st_bbox(tr_farm)[3], st_bbox(tr_farm)[2], st_bbox(tr_farm)[4])
bbox <- matrix(bbox, ncol = 2, nrow = 2, byrow = TRUE)

m <- mapview(tr_farm, map.types = c("Esri.WorldImagery"), color = "darkblue") +
  mapview(tr_spec, map.types = c("Esri.WorldImagery"), color = "violet")

mapshot(m, file = paste0(getwd(), "/map.png"))

## ---- Fields cmr ----

cr2 <- st_read("D:/MargSalas/Ganga/Data/CMR/effort_22.shp") # These are the analyzed transects in hds model

## ---- Areas density extrapolation ----

zones <- st_read("D:/MargSalas/Ganga/Data/GIS/zonesGanga2021.shp")
zones2 <- st_union(zones)

spa <- st_read("D:/PhD/GIS Ana_14_Mayo/Capes GIS/Special_Protection_Areas_ETRS89.shp")
spa_2 <- st_crop(spa,sa)

# Remove transects outside (not used)
ovr <- c(st_intersects(tr_farm, zones2)) 
subst <- lapply(c(ovr), function(x){x == 1}) == TRUE
ovr[which(is.na(subst))] <- 0

tr_farm2 <- tr_farm[which(unlist(ovr) == 1), ]




## ---- Plot ----

mapview(cr2, map.types = c("Esri.WorldImagery"), col.regions = adjustcolor("#35978f", alpha = 1), col = adjustcolor("black", alpha = 1), lwd = 2, cex = 3) +
  mapview(sa, col.regions = "white", alpha.regions = 0.2, col = adjustcolor("black", alpha = 1), lwd = 2) +
  mapview(tr_farm2, color = "#253494", lwd = 3.5) +
  mapview(tr_spec, color = "#d73027", lwd = 3.5) +
  #mapview(zones2, col.regions = adjustcolor("yellow", alpha = 0.5), col = adjustcolor("yellow", alpha = 1)) +
  mapview(spa_2, col.regions = adjustcolor("khaki", alpha = 0.5), col = adjustcolor("khaki", alpha = 1))
  
plot(1)
plot.new()
legend("center", legend = c("General HDS", "Specific HDS", "CR"), 
       fill = c("#253494", "#d73027", "#35978f"), border = NA, bty = "n", horiz = FALSE, cex = 1.5)

## ---- Plot spain ----

sp <- st_read("D:/MargSalas/Ganga/Data/GIS/ESP_adm1.shp") %>%
  st_transform(crs = st_crs(sa))

setwd("D:/MargSalas/Ganga/Results/Plots/Study_design")
pdf("spain.pdf",7,7)

plot(st_geometry(sp), col = "beige", border = "yellow4")
plot(st_geometry(sa), lwd = 2, add = TRUE)
dev.off()

## ---- Zoom plot study design ----

z <- st_read("D:/MargSalas/Ganga/Data/GIS/Zoom_study_design.shp")

mapview(cr2, map.types = c("Esri.WorldImagery"), col.regions = adjustcolor("#35978f", alpha = 1), col = adjustcolor("black", alpha = 1), lwd = 2, cex = 4) +
  mapview(z, col.regions = "white", alpha.regions = 0, col = adjustcolor("black", alpha = 1), lwd = 2) +
  mapview(tr_farm2, color = "#253494", lwd = 4.5) +
  mapview(tr_spec, color = "#d73027", lwd = 4.5) 

## ---- For GIAE ----

## -------------------------------------------------
##                 Plot study designs
## ------------------------------------------------- 

rm(list=ls())

library(sf)
library(mapview)
library(ggmap)
library(ggplot2)
library(ggmap)
library(gridExtra)
library(rgdal)
library(ggpubr)
library(dplyr)
library(geosphere)
library("viridis")
library(RColorBrewer)

devtools::install_github('oswaldosantos/ggsn')


## ---- Transects ----

tr_farm <- st_read("D:/MargSalas/Ganga/Data/GIS/Trans_Ganga_FarmDindis.shp") # These are the analyzed transects in hds model
tr_spec <- st_read("D:/MargSalas/Ganga/Data/GIS/Trans_Ganga_500_Specific.shp") # These are the analyzed transects in hds model

tr_spec <- st_transform(tr_spec, crs = st_crs(tr_farm))

bbox <- c(st_bbox(tr_farm)[1], st_bbox(tr_farm)[3], st_bbox(tr_farm)[2], st_bbox(tr_farm)[4])
bbox <- matrix(bbox, ncol = 2, nrow = 2, byrow = TRUE)

m <- mapview(tr_farm, map.types = c("Esri.WorldImagery"), color = "darkblue") +
  mapview(tr_spec, map.types = c("Esri.WorldImagery"), color = "violet")

mapshot(m, file = paste0(getwd(), "/map.png"))

## ---- Fields cmr ----

cr <- st_read("D:/MargSalas/Ganga/Data/CMR/XYcoord_g.shp") # These are the analyzed transects in hds model

## ---- Areas density extrapolation ----

zones <- st_read("D:/MargSalas/Ganga/Data/GIS/zonesGanga2021.shp")
zones2 <- st_union(zones)

# Remove transects outside (not used)
ovr <- c(st_intersects(tr_farm, zones2)) 
subst <- lapply(c(ovr), function(x){x == 1}) == TRUE
ovr[which(is.na(subst))] <- 0

tr_farm2 <- tr_farm[which(unlist(ovr) == 1), ]

## ---- Study area delimitation ----
sa <- st_read("D:/MargSalas/Ganga/Data/GIS/Study_area.shp")

## ---- Plot ----

  mapview(sa,  map.types = c("Esri.WorldImagery"), col.regions = "white", alpha.regions = 0.2, col = adjustcolor("black", alpha = 1), lwd = 2) +
  mapview(zones2, col.regions = adjustcolor("yellow", alpha = 0.5), col = adjustcolor("yellow", alpha = 1))
