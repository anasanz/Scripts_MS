## -------------------------------------------------
##      Plot trap locations appendix
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
library(nimbleSCR)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")

load("edf1721.RData")

# As the model is set, there can be only one capture per trap and occasion
# --> The number of trials is 7 (From may to November): So max number of captures per trap = 7
# To fix it in this edf, remove duplicates
edf <- edf[-which(duplicated(edf)), ]

load("tdf2017_effort.RData")
load("tdf2018_effort.RData")
load("tdf2019_effort.RData")
load("tdf2020_effort.RData")
load("tdf2021_effort.RData")

tdf_all <- rbind(tdf2017[,1:3], tdf2018[,1:3], tdf2019[,1:3],
                 tdf2020[,1:3], tdf2021[,1:3]) # Join to define state space
rownames(tdf_all) <- 1:nrow(tdf_all)

# GET TRAPS 
X <- tdf_all[,c(2,3)]
colnames(X) <- c('x', 'y')

# For plot:

#  -> Area
sa <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/clip_pyros2_WGS84_31N_all.shp")
eur <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/esp_fr_2.shp") %>%
  st_transform(dpts, crs = crs(sa))
esp <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/ESP_adm/ESP_adm0.shp") %>%
  st_transform(dpts, crs = crs(sa))
fr <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/FRA_adm/FRA_adm0.shp") %>%
  st_transform(dpts, crs = crs(sa))

#  -> Buffers
setwd("D:/MargSalas/Oso/Datos/GIS/Countries")
Xbuf <- readOGR("Buffer_statespace.shp")
Xbuf2 <- readOGR("Buffer_8500_traps.shp")
Xbuf3 <- readOGR("Buffer_8500_traps_sxyObs.shp")


#  -> Observed individuals outside

# Load z of observed individuals
setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")
load("zObserved_yAgeDeaths.RData")

zdatAGE[is.na(zdatAGE)] <- 0

# Load sxy estimated by the model and unscale

Tt <- 5
M.aug <- 300

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL")
load("myResults_3-3.1_sxy.RData")

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721") # Load original habitat coordinates to unscale
load("habcoord.RData")

sampmat2 <- do.call(rbind, nimOutputSXY) 
s.which <- grep('sxy', colnames(sampmat2)) # ASP: index columns all sxy (sampmat matrix)
sampmat2_sxy <- sampmat2[, s.which]
dim(sampmat2_sxy)
mean_sampmat <- colMeans(sampmat2_sxy) # I will plot the mean location over iterations

sxy <- array(NA, c(300, 2, 5))
for(t in 1:Tt){
  s.which.year <- grep(paste(t,"]", sep = ""), colnames(sampmat2_sxy)) # ASP: index columns all sxy (sampmat matrix)
  sxy[,,t] <- matrix(mean_sampmat[s.which.year] , M.aug, 2) 
}

dimnames(sxy)[[2]] <- c('x','y') 
sxy.uns <- scaleCoordsToHabitatGrid(coordsData = sxy,## this are your sxy
                                    coordsHabitatGridCenter = G,# this is your unscaled habitat (as you used when scaling the habitat/detector to the habitat. G?
                                    scaleToGrid = FALSE)$coordsDataScaled

# Check observed individuals that model places outside
yearnames <- c("2017", "2018", "2019", "2020", "2021")

#setwd("D:/MargSalas/Oso/OPSCR_project/Results/Plots/model3.4_resub")
#pdf("AC_observed_alive.pdf",7,7)
par(mfrow = c(3,2))
for (t in 1:5){
  plot(Xbuf, main = yearnames[t])
  plot(Xbuf2, add = TRUE)
  p <- sxy.uns[which(zdatAGE[,t] == 1),,t]
  sp <- SpatialPoints(p, proj4string=CRS(proj4string(Xbuf2)))
  points(sp, pch= 21, col = "red")
  
  over(sp,Xbuf2)
}
#dev.off()

# Identify which ones are outside
p.outside <- list()
o <- list()
w <- list()
for (t in 1:5){
  which.alive <- which(zdatAGE[,t] == 1)
  p <- sxy.uns[which.alive,,t]
  sp <- SpatialPoints(p, proj4string=CRS(proj4string(Xbuf2)))
  which.out <- which(is.na(over(sp,Xbuf2)))
  w$age <- ageMatAug[which.alive[which.out],t]
  w$sex <- sex[which.alive[which.out]]
  w$ID <- which.alive[which.out]
  o[[t]] <- w # Store info
  p.outside[[t]] <- sp[which.out] # Store points
  
} # All are adult and subadult males

m <- do.call(bind, p.outside)

# Check ID 39
zdatAGE[39,]

par(mfrow = c(3,2))
for (t in 1:5){
  plot(Xbuf, main = yearnames[t])
  plot(Xbuf2, add = TRUE)
  p <- matrix(sxy.uns[39,,t], nrow = 1, ncol = 2)
  sp <- SpatialPoints(p, proj4string=CRS(proj4string(Xbuf2)))
  points(sp, pch= 21, col = "red")
  
}
