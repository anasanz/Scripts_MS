
## -------------------------------------------------
##                 Effort 2017
## ------------------------------------------------- 


rm(list = ls())

library(tidyverse)
library(sf)
library(mapview)


#### FRANCE ####

# Traps positions (Maelis)

setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Oso/PopDyn/SCR/Effort/France")

load("data.RData")

crs_fr <- CRS("+proj=utm +zone=31 +datum=WGS84")

# Plot to see the opportunistic traps

traps_2017_op <- data[["traps"]][["2017"]] %>%
  filter(suivi == "opportunist")

mapview(traps_2017_op)

# Systematic traps

traps_2017 <- data[["traps"]][["2017"]] %>%
  filter(suivi == "systematic")

tdf2017 <- data.frame(trap = as.numeric(traps_2017$trap), # trap id number
                      X = st_coordinates(traps_2017)[,1], # longitude
                      Y = st_coordinates(traps_2017)[,2],# latitude
                      sep = rep("/",length(traps_2017$trap)), # separator
                      site = as.factor(traps_2017$site), # type of trap
                      pays = as.factor(traps_2017$pays), # country
                      suivi = as.factor(traps_2017$suivi)) # type of monitoring

# Effort covariate : factor which account for the country and the number of visits per occasion in France
# 1 = one visit per month in France 
# 2 = two visits per month in France
# 3 = trap is in Spain

nocc <- 7 # number of occasion (month) per session (year) : may to november

for (j in 1 : nocc) {
  eff <- rep(0, length(traps_2017$trap))
  for (i in 1:length(traps_2017$trap)){
    if (traps_2017$pays[i] == "Espagne"){ # if the trap is in Spain
      if(traps_2017$suivi[i] == "systematic"){
        eff[i] <- 3} # site systematic in Spain
      else {
        eff[i] <- 3} # site opportunist in Spain 
    }
    else {
      if (traps_2017$suivi[i] == "systematic"){
        if (traps_2017$effort[i] == "it" ){ # if the trap is on an transect/itinerary
          if (j == 1 | j == 2 | j == 5){ # and if we are in may, june or september
            eff[i] <- 2}  # the site is visited 2 times at this occasion
          else {
            eff[i] <- 1}} # the site is visited once
        else {
          eff[i] <- 1} # the site is visited once
      }
      else{
        eff[i] <- 1} # site opportunist in France 
    }
  } 
  name <- paste("effort", j, sep = ".")
  tdf2017[name] <- as.factor(eff)
}


# Select datos only systematic traps france (because in Spain I will take the year-specific data)
tdf2017_fr <- tdf2017[which(tdf2017$pays == "France" & tdf2017$suivi == "systematic"), ]

coordinates(tdf2017_fr) <- tdf2017_fr[ ,c("X", "Y")]
tdf2017_fr@proj4string <- crs_fr

mapview(tdf2017_fr, zcol = "effort.1")



## ---- Encounter data (ctfc) ----

setwd("D:/MargSalas/Oso/Datos/Tablas_finales")

os <- read.csv("Seguiment_Ossos_Pirineus_1996_2020_taula_final.csv", header = TRUE, row.names = NULL)
det_2017_ctfc <- os[which(os$Year %in% 2017 & 
                       os$Country %in% "France" & 
                       os$Probable_Individual != "Indetermined" &
                       os$Method %in% c("Sampling_station", "Transect") &
                       os$Obs_type %in% c("Hair", "Photo", "Photo/Video", "Video")), ]

coordinates(det_2017_ctfc) <- det_2017_ctfc[,c("x_long","y_lat")]
det_2017_ctfc@proj4string <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
det_2017_ctfc_TRANSF <- spTransform(det_2017_ctfc, crs_fr)

# Spatial join

mapview(tdf2017_fr, col.regions = "red") + mapview(det_2017_ctfc_TRANSF)


det_2017_ctfc
unique(os$Method)
det_2017_ctfc <- det_2017_ctfc %>%
  st_as_sf(coords = c("x_long", "y_lat"), crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

mapview(det_2017_ctfc)








mapview(det_2017_mae)


unique(os$Obs_type)


#  Unir datos espacialmente Francia

#### SPAIN ####

# Datos Spain Maelis: Como puede ser que tenga muestreo sistemático en España?
# Cargar datos España (encounter y traps)

