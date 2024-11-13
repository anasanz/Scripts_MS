## -------------------------------------------------
##    Plot density evolution all individuals
##          and spatial projections
## ------------------------------------------------- 

rm(list = ls())

library(Rcpp)
library(RcppArmadillo)
library(RcppProgress)
library(rgdal)
library(raster)
library(viridis)

setwd("D:/MargSalas/Scripts_MS/Functions/Nimble")

# Load functions
#sourceCpp("GetSpaceUse_PD.cpp")
sourceCpp("GetDensity_PD.cpp")
source("getDensityInput.R")

# Load political map study area
sa <- readOGR("D:/MargSalas/Oso/Datos/GIS/Countries", layer = "clip_pyros2_WGS84_31N_all")
eur <- readOGR("D:/MargSalas/Oso/Datos/GIS/Countries", "esp_fr_2")
eur <- spTransform(eur, crs(sa))

# Load buffers
setwd("D:/MargSalas/Oso/Datos/GIS/Countries")
Xbuf <- readOGR("Buffer_statespace.shp")
#Xbuf2 <- readOGR("Buffer_8500_traps.shp")
Xbuf2 <- readOGR("Buffer_8500_traps_sxyObs.shp") # This sampling buffer includes AC of observed individuals a bit outside the trapping array

nuc <- readOGR("nuc.shp")
per <- readOGR("per.shp")

# Conver to rasters

setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
r <- raster("logDistcore_hrbear.tif")
raster::values(r) <- NA
r <- crop(r,Xbuf) 

Xbuf_raster <- rasterize(Xbuf, r)
raster::values(Xbuf_raster)[which(is.na(raster::values(Xbuf_raster)))] <- 0 # Raster of 0 and 1


Xbuf2_raster <- rasterize(Xbuf2, r, getCover = TRUE)
raster::values(Xbuf2_raster)[which(raster::values(Xbuf2_raster) > 0.5)] <- 2
raster::values(Xbuf2_raster)[which(raster::values(Xbuf2_raster) < 0.6)] <- 0

ras <- overlay(Xbuf_raster, Xbuf2_raster, fun = max)
f <- rasterize(Xbuf, ras, mask = TRUE)
#values(f)[which(values(f) == 1)] <- 0
#values(f)[which(values(f) == 2)] <- 1
f1 <- f
f1[f1%in%c(1,2)] <-  1
f[] <- as.factor(as.character(f[]))

# Load original habitat coordinates
setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")
load("habcoord.RData")

# Load posterior distribution
library(nimbleSCR) # Load nimbleSCR here, otherwise it gets in conflict with raster package, weird


#### LOAD AND PLOT

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Results_section/Plots")
pdf("2.1.density_evol_allinds3.pdf", 9,6)

par(mfrow = c(3,5),
    mar = c(0,0,0,0),
    oma = c(0.5,1.5,2,2),
    bty = "n")

#palette(c("#d9d9d9", "#bdbdbd", "#969696", "#737373", "#525252", "#252525", "#000000"))
leg <- c(FALSE, FALSE, FALSE, FALSE, TRUE)

## ---- 1. Density evolution all individuals ----

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL")
load("myResults_3-3.1_sxy.RData")

#give dim names to your posteriors sxy
dim(myResultsSXYZ$sims.list$sxy) # YOURPOSTERIORSXY 
dimnames(myResultsSXYZ$sims.list$sxy)[[3]] <- c("x","y")

## first rescale the coordinates to the original scale 
myResultsSXYZ$sims.list$sxy <- scaleCoordsToHabitatGrid(coordsData = myResultsSXYZ$sims.list$sxy,
                                                        coordsHabitatGridCenter = G,
                                                        scaleToGrid = FALSE)$coordsDataScaled

densityInputRegions <- getDensityInput( 
  regions = f,## THIS  A RASTER FILE WITH 0/1 HABITAT VS BUFFER, OR 1/2 fRANCE SPAIN... WHATEVER YOU WANT. 
  habitat = f1,## here put the same than regions argument. 
  s = myResultsSXYZ$sims.list$sxy,
  plot.check = FALSE)

DensityCountriesRegions <- list()
n.years = 5
for(t in 1:n.years){
  DensityCountriesRegions[[t]] <- GetDensity_PD(
    sx = densityInputRegions$sx[,,t],# X COORDINATES
    sy =  densityInputRegions$sy[,,t],# Y COORDINATES
    z = myResultsSXYZ$sims.list$z[,,t],# Z 
    IDmx = densityInputRegions$habitat.id,
    aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
    regionID = densityInputRegions$regions.rgmx,
    returnPosteriorCells = F)
}#t

# Get maximum density of all to scale
maxdens1 <- max(max(DensityCountriesRegions[[1]]$MeanCell), max(DensityCountriesRegions[[2]]$MeanCell), max(DensityCountriesRegions[[3]]$MeanCell),
                  max(DensityCountriesRegions[[4]]$MeanCell), max(DensityCountriesRegions[[5]]$MeanCell))

#plot()
ACDens<-list()
for(t in 1:n.years){
  ACDens[[t]] <- densityInputRegions$regions.r
  ACDens[[t]][] <- NA
  ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  
  ACDens2 <- ACDens[[t]]
  raster::values(ACDens2)[which(raster::values(Xbuf2_raster) == 0)]  <- NA
  
  plot(ACDens2, xaxt = "n", yaxt = "n", bty = "n", zlim=c(0,maxdens1), col = rev(magma(100)), legend = leg[t])
}

## ---- Prediction 1: Constant effect of distcore accross time ----

# Load results from predictions

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL/Predictions/ALLiter")
load("proj_pcr.all.fem.consDist.RData")

# Rescale sxy

dimnames(sxy.proj.all.consDist)[[3]] <- c('x','y')
sxy.proj.all.consDist.uns <- scaleCoordsToHabitatGrid(coordsData = sxy.proj.all.consDist,## this are your sxy
                                                      coordsHabitatGridCenter = G,# this is your unscaled habitat (as you used when scaling the habitat/detector to the habitat. G?
                                                      scaleToGrid = FALSE)$coordsDataScaled

##GET OBJECTS IN SHAPE
densityInputRegions <- getDensityInput( 
  regions = f,## THIS  A RASTER FILE WITH 0/1 HABITAT VS BUFFER, OR 1/2 fRANCE SPAIN... WHATEVER YOU WANT. 
  habitat = f1,## here put the same than regions argument. 
  s = sxy.proj.all.consDist.uns,
  plot.check = FALSE)

DensityCountriesRegions <- list()
n.years = 6
for(t in 1:n.years){
  DensityCountriesRegions[[t]] <- GetDensity_PD(
    sx = densityInputRegions$sx[,,t],# X COORDINATES
    sy =  densityInputRegions$sy[,,t],# Y COORDINATES
    z = z.proj.all.consDist[,,t],# Z 
    IDmx = densityInputRegions$habitat.id,
    aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
    regionID = densityInputRegions$regions.rgmx,
    returnPosteriorCells = F)
}#t

# Maximum density of all years

maxdens2 <- max(max(DensityCountriesRegions[[1]]$MeanCell), max(DensityCountriesRegions[[2]]$MeanCell), max(DensityCountriesRegions[[3]]$MeanCell),
               max(DensityCountriesRegions[[4]]$MeanCell), max(DensityCountriesRegions[[5]]$MeanCell), max(DensityCountriesRegions[[6]]$MeanCell))

leg <- c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)

ACDens<-list()
for(t in 2:n.years){ # Start in 2 because we don't plot year 2021
  ACDens[[t]] <- densityInputRegions$regions.r
  ACDens[[t]][] <- NA
  ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  
  ACDens2 <- ACDens[[t]]
  raster::values(ACDens2)[which(raster::values(Xbuf2_raster) == 0)]  <- NA
  
  plot(ACDens2, xaxt = "n", yaxt = "n", bty = "n", zlim=c(0,maxdens2), col = rev(magma(100)), legend = leg[t])
}

## ---- Prediction 2: Effect of distcore decreases progressively accross time ----

# Load results from predictions

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL/Predictions/ALLiter")
load("proj_pcr.all.fem.decDistProg1.RData")

# Rescale sxy

dimnames(sxy.proj.all.decDist)[[3]] <- c('x','y')
sxy.proj.all.decDist.uns <- scaleCoordsToHabitatGrid(coordsData = sxy.proj.all.decDist,## this are your sxy
                                                     coordsHabitatGridCenter = G,# this is your unscaled habitat (as you used when scaling the habitat/detector to the habitat. G?
                                                     scaleToGrid = FALSE)$coordsDataScaled

##GET OBJECTS IN SHAPE
densityInputRegions <- getDensityInput( 
  regions = f,## THIS  A RASTER FILE WITH 0/1 HABITAT VS BUFFER, OR 1/2 fRANCE SPAIN... WHATEVER YOU WANT. 
  habitat = f1,## here put the same than regions argument. 
  s = sxy.proj.all.decDist.uns,
  plot.check = FALSE)

DensityCountriesRegions <- list()
n.years = 6
for(t in 1:n.years){
  DensityCountriesRegions[[t]] <- GetDensity_PD(
    sx = densityInputRegions$sx[,,t],# X COORDINATES
    sy =  densityInputRegions$sy[,,t],# Y COORDINATES
    z = z.proj.all.decDist[,,t],# Z 
    IDmx = densityInputRegions$habitat.id,
    aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
    regionID = densityInputRegions$regions.rgmx,
    returnPosteriorCells = F)
}#t


# Maximum density of all years

maxdens3 <- max(max(DensityCountriesRegions[[1]]$MeanCell), max(DensityCountriesRegions[[2]]$MeanCell), max(DensityCountriesRegions[[3]]$MeanCell),
               max(DensityCountriesRegions[[4]]$MeanCell), max(DensityCountriesRegions[[5]]$MeanCell), max(DensityCountriesRegions[[6]]$MeanCell))

#plot()
ACDens<-list()
for(t in 2:n.years){
  ACDens[[t]] <- densityInputRegions$regions.r
  ACDens[[t]][] <- NA
  ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell

  ACDens2 <- ACDens[[t]]
  raster::values(ACDens2)[which(raster::values(Xbuf2_raster) == 0)]  <- NA
  
  plot(ACDens2, xaxt = "n", yaxt = "n", bty = "n", zlim=c(0,maxdens3), col = rev(magma(100)), legend = leg[t])
}


dev.off()
