
## -------------------------------------------------
##      Arrange data DS Farmdindis 2010-2022
## -------------------------------------------------

rm(list=ls())

library(dplyr)
library(stringr)
library(rgdal)
library(sf)

setwd("D:/MargSalas/Ganga/Data/FarmdindisDS")

dat <- read.csv("Transectes_Farmdindis_Ganga_TOTS anys.csv", sep = ";")

colnames(dat)[which(colnames(dat) == "ï..Codi_seca")] <- "Region.Label"
colnames(dat)[which(colnames(dat) == "anys")] <- "Year"
colnames(dat)[which(colnames(dat) == "Hora_inici")] <- "Start_time"
colnames(dat)[which(colnames(dat) == "Hora_final")] <- "End_time"
colnames(dat)[which(colnames(dat) == "Nuvolositat")] <- "Clouds"
colnames(dat)[which(colnames(dat) == "Temperatura")] <- "Temp"
colnames(dat)[which(colnames(dat) == "codiEspecie")] <- "Species"
colnames(dat)[which(colnames(dat) == "Num")] <- "Count" 
dat$Effort <- 500

## ---- Create variable transectID, than matches with the code of the GIS layers (i.e., two digits: 09) ----

#1. Add a 0 before the transect number
for (i in 1:nrow(dat)){ 
  dat$Num_transecte[i] <- paste(0,dat$Num_transecte[i], sep = "")
}

#2. Keep only the last 2 digits (or 3 in the case of the transects that contain 100)

for (i in 1:nrow(dat)) { 
  cent <- substr(dat$Num_transecte[i], 4,4)
  cent <- as.numeric(cent) # NA if it doesnt have 4 digits
  if(is.na(cent)) { # if is NA (has 3 digits)
    dat$Num_transecte[i] <- str_sub(dat$Num_transecte[i], start = -2) # Keep the last 2
  } else { dat$Num_transecte[i] <- str_sub(dat$Num_transecte[i], start = -3)} # Otherwise, keep the last 3
}
# Create variable by pasting it
for (i in 1:nrow(dat)){ 
  dat$transectID[i] <- paste(dat$Region.Label[i],dat$Num_transecte[i], sep = "")
}

## ---- Transect-Year variable ----

for (i in 1:nrow(dat)){ 
  dat$T_Y[i] <- paste(dat$transectID[i],dat$Year[i], sep = "_")
}

## ---- Add info of distance bands ----

# Banda 1: 0-25
# Banda 2: 25-50
# Banda 3: 50-100
# Banda 4: 100-200
# Banda 5: 200-500

dat$distance <- NA # Medium point of each bin

for (i in 1:nrow(dat)){
  if (!is.na(dat$Banda[i])){
    if (dat$Banda[i] == 1) {dat$distance[i] = 12.5}
    else if (dat$Banda[i] == 2) {dat$distance[i] = 37.5}
    else if (dat$Banda[i] == 3) {dat$distance[i] = 75}
    else if (dat$Banda[i] == 4) {dat$distance[i] = 150}
    else if (dat$Banda[i] == 5) {dat$distance[i] = 350}
  }}

## ---- Choose transects from areas where species is present ----

tr <- readOGR("D:/MargSalas/Ganga/Data/FarmdindisDS", "Trans_Ganga_FarmDindis") 
tr@data$Codi[which(tr@data$Pres_Abs == 0)]

# Abs-Pres classification (David Giralt):
#   0: Transectos donde nunca se ha visto ganga (2010-2022)
#   1: Transectos donde se han visto al menos 1 año en el periodo 2010-2022
#   2: Transectos donde se ha visto al menos 1 año entre 2010-2016

## ---- Extract zone of habitat quality where each transect is located ----

pa <- readOGR("D:/MargSalas/Ganga/Data/FarmdindisDS", "zonesGanga")

tr_sf <- st_as_sf(tr)
pa_sf <- st_as_sf(pa)

j <- st_join(tr_sf, pa_sf)

# Check duplicates to keep the one with the majority in the zone X

j <- j[,c(2,3,4,20,23,24)]
jdup <- j[,-c(5)]
check <- jdup[which(duplicated(jdup)), ]
unique(check$Codi)

