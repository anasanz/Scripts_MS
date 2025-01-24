## -------------------------------------------------
##      Define sampling buffer to estimat abundance
##        AIM: Estimate abundance only in trap array
## ------------------------------------------------- 

rm(list = ls())

library(nimble)
library(MCMCvis)
library(nimbleSCR)
library(raster)
library(rgeos)
library(oSCR)
library(terra)
library(sp)
library(dplyr)
library(parallel)
library(sf)
library(lubridate)
library(rgdal)



setwd("D:/MargSalas/Scripts_MS/Stats/Nimble")
#source('dbinomLocal_normalBear.R')
source('dbinomLocal_normalBear_rbinom2.R')

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")

#---- 1. LOAD THE DETECTION DATA ---- 

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


# We remove Nere and Goiat, which are two individuals moving a lot that don't represent
# the rest of the population.

edf <- edf[-which(edf$ind %in% c("Nere", "Goiat")), ]

#---- 2. DEFINE THE TRAP AND THE HABITAT  ---- 
# GET TRAPS 
X <- tdf_all[,c(2,3)]
colnames(X) <- c('x', 'y')
J <- dim(X)[1]

# DEFINE STATE SPACE EXTENT
# State space coordinates
# Buffer: 25000 (used by Maelis, also ~3*sigma in pre-analysis where sig = 6640)
xmin <- min(X[,1]) - 25000
ymin <- min(X[,2]) - 25000
xmax <- max(X[,1]) + 25000
ymax <- max(X[,2]) + 25000
e <- as(raster::extent(xmin, xmax, ymin, ymax), "SpatialPolygons") # Extent of state space

# GET A RASTER FOR THE HABITAT 
# USE A FOREST RASTER TO GET A BASIS FOR THE HABITAT RASTER
# Set up a raster file 
setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
#setwd("~/Data_server/Variables_hrscale")
distcore <- raster("logDistcore_hrbear.tif")
# Crop it to extent of state-space
habitat.r <- crop(distcore, e) 
plot(habitat.r)

# ---- DEFINE THE BUFFER AREA, DECIDE WHICH AREA IS CONSIDERED NOT BUFFER ---- 
# For resubmission, we only do final step
## ---- 6. FINAL definition of sampling buffer: xbuf2 (traps) + observed individuals outside ----

library(nimbleSCR)
library(rgeos)

# Load state space and buffer with traps
setwd("D:/MargSalas/Oso/Datos/GIS/Countries")
Xbuf <- readOGR("Buffer_statespace.shp")
Xbuf2 <- readOGR("Buffer_8500_traps.shp")

## Identify observed individuals that are placed outside the sampling buffer, to join it

# Load z of observed individuals
setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")
load("zObserved_yAgeDeaths_resub.RData")

zdatAGE[is.na(zdatAGE)] <- 0

# Load sxy estimated by the model and unscale

Tt <- 5
M.aug <- 300

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/RESUBMISSION/3-3.4_allparams")
load("myResults_RESUB_3-3.4_sxy.RData")

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

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Plots/model3.4_resub")
pdf("AC_observed_alive.pdf",7,7)
par(mfrow = c(3,2))
for (t in 1:5){
  plot(Xbuf, main = yearnames[t])
  plot(Xbuf2, add = TRUE)
  p <- sxy.uns[which(zdatAGE[,t] == 1),,t]
  sp <- SpatialPoints(p, proj4string=CRS(proj4string(Xbuf2)))
  points(sp, pch= 21, col = "red")
  
  over(sp,Xbuf2)
}
dev.off()

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
  o[[t]] <- w # Store info
  p.outside[[t]] <- sp[which.out] # Store points
  
} # All are adult and subadult males

m <- do.call(bind, p.outside)

# Buffer of 8500 m (as for the traps)
Xbuf3 <- gBuffer(m, width = 8500)
plot(Xbuf3, add = TRUE)

Xbuf3 <- crop(Xbuf3,Xbuf) # Crop with state space

Xbuf3.df <- data.frame( ID=1:length(Xbuf3)) # To spdf to join
Xbuf3 <- SpatialPolygonsDataFrame(Xbuf3, Xbuf3.df) 

XbufFin <- c(Xbuf2, Xbuf3)
XbufFin <- do.call(bind, XbufFin)
XbufFin <- aggregate(XbufFin)

par(mfrow = c(1,1))
plot(Xbuf)
plot(Xbuf2, add = TRUE)
plot(Xbuf3, add = TRUE)
plot(XbufFin, col = "red", add = TRUE)
plot(XbufFin4, col = "blue", add = TRUE)
points(m)


XbufFin.df <- data.frame( ID=1:length(XbufFin)) # To spdf to save
XbufFin <- SpatialPolygonsDataFrame(XbufFin, XbufFin.df) 

### I have edited this one in arcgis as "Buffer_8500_traps_sxyObs_RESUB2" to remove the hole

setwd("D:/MargSalas/Oso/Datos/GIS/Countries")
writeOGR(XbufFin, dsn = "D:/MargSalas/Oso/Datos/GIS/Countries", layer = "Buffer_8500_traps_sxyObs_RESUB", driver = "ESRI Shapefile")  


