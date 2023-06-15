## -------------------------------------------------
##       Plot predictions DistCore effect
## ------------------------------------------------- 

rm(list = ls())

library(nimbleSCR)
library(nimble)
library(rgdal)
library(raster)
library(viridis)
library(Rcpp)
library(RcppArmadillo)
library(RcppProgress)

## ---- Load necessary stuff ----

setwd("D:/MargSalas/Scripts_MS/Functions/Nimble")
#sourceCpp("GetSpaceUse_PD.cpp")
sourceCpp("GetDensity_PD.cpp")
source("getDensityInput.R")

# Load buffer core area

setwd("D:/MargSalas/Oso/Datos/GIS/Countries")
Xbuf <- readOGR("Buffer_statespace.shp")
Xbuf2 <- readOGR("Buffer_8500_traps.shp")

# Conver buffers state space and core sampling to rasters

setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
r <- raster("logDistcore_hrbear.tif")
raster::values(r) <- NA
r <- crop(r,Xbuf) 

Xbuf_raster <- rasterize(Xbuf, r)
raster::values(Xbuf_raster)[which(is.na(raster::values(Xbuf_raster)))] <- 0 # Raster of 0 and 1

Xbuf2_raster <- rasterize(Xbuf2, r)
raster::values(Xbuf2_raster)[which(raster::values(Xbuf2_raster) == 1)] <- 2
raster::values(Xbuf2_raster)[which(is.na(raster::values(Xbuf2_raster)))] <- 0 # Raster of 0 and 2

ras <- overlay(Xbuf_raster, Xbuf2_raster, fun = max)
f <- rasterize(Xbuf, ras, mask = TRUE)
f1 <- f
f1[f1%in%c(1,2)] <-  1
f[] <- as.factor(as.character(f[]))

# Load habitat coordinates to unscale

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")
load("habcoord.RData")

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
  plot.check = TRUE)

## extract density
#yearnames <- c("2017", "2018", "2019", "2020", "2021", "2022", "2023", "2024", "2025", "2026")
yearnames <- c("2021", "2022", "2023", "2024", "2025", "2026")
leg <- c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE)

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Plots/model3.1/Predictions")
pdf("1_ConstantDistCore_pcr.all_ALLiter.pdf",12,6)


par(mfrow = c(3,2),
    mar = c(0,1,2.5,1),
    oma = c(0.5,1.5,2,2),
    bty = "n")

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

#mean density in each cell 
DensityCountriesRegions[[1]]$MeanCell

# Maximum density of all years

maxdens <- max(max(DensityCountriesRegions[[1]]$MeanCell), max(DensityCountriesRegions[[2]]$MeanCell), max(DensityCountriesRegions[[3]]$MeanCell),
               max(DensityCountriesRegions[[4]]$MeanCell), max(DensityCountriesRegions[[5]]$MeanCell), max(DensityCountriesRegions[[6]]$MeanCell))

#plot()
ACDens<-list()
for(t in 1:n.years){
  ACDens[[t]] <- densityInputRegions$regions.r
  ACDens[[t]][] <- NA
  ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  plot(ACDens[[t]], xaxt = "n", yaxt = "n", bty = "n", zlim=c(0,maxdens), col = rev(magma(100)), legend = leg[t], main = yearnames[t])
}

dev.off()

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
  plot.check = TRUE)

## extract density
#yearnames <- c("2017", "2018", "2019", "2020", "2021", "2022", "2023", "2024", "2025", "2026")
yearnames <- c("2021", "2022", "2023", "2024", "2025", "2026")
leg <- c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE)

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Plots/model3.1/Predictions")
pdf("1_DecProg1DistCore_pcr.all_ALLiter.pdf",12,6)


par(mfrow = c(3,2),
    mar = c(0,1,2.5,1),
    oma = c(0.5,1.5,2,2),
    bty = "n")

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

#mean density in each cell 
DensityCountriesRegions[[1]]$MeanCell

# Maximum density of all years

maxdens <- max(max(DensityCountriesRegions[[1]]$MeanCell), max(DensityCountriesRegions[[2]]$MeanCell), max(DensityCountriesRegions[[3]]$MeanCell),
               max(DensityCountriesRegions[[4]]$MeanCell), max(DensityCountriesRegions[[5]]$MeanCell), max(DensityCountriesRegions[[6]]$MeanCell))

#plot()
ACDens<-list()
for(t in 1:n.years){
  ACDens[[t]] <- densityInputRegions$regions.r
  ACDens[[t]][] <- NA
  ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  plot(ACDens[[t]], xaxt = "n", yaxt = "n", bty = "n", zlim=c(0,maxdens), col = rev(magma(100)), legend = leg[t], main = yearnames[t])
}

dev.off()

## ---- Prediction 3: Effect of distcore decreases to 0 accross time ----

# Load results from predictions

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL/Predictions/ALLiter")
load("proj_pcr.all.fem.decDistProg2.RData")

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
  plot.check = TRUE)

## extract density
#yearnames <- c("2017", "2018", "2019", "2020", "2021", "2022", "2023", "2024", "2025", "2026")
yearnames <- c("2021", "2022", "2023", "2024", "2025", "2026")
leg <- c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE)

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Plots/model3.1/Predictions")
pdf("1_DecProg2DistCore_pcr.all_ALLiter.pdf",12,6)


par(mfrow = c(3,2),
    mar = c(0,1,2.5,1),
    oma = c(0.5,1.5,2,2),
    bty = "n")

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

#mean density in each cell 
DensityCountriesRegions[[1]]$MeanCell

# Maximum density of all years

maxdens <- max(max(DensityCountriesRegions[[1]]$MeanCell), max(DensityCountriesRegions[[2]]$MeanCell), max(DensityCountriesRegions[[3]]$MeanCell),
               max(DensityCountriesRegions[[4]]$MeanCell), max(DensityCountriesRegions[[5]]$MeanCell), max(DensityCountriesRegions[[6]]$MeanCell))

#plot()
ACDens<-list()
for(t in 1:n.years){
  ACDens[[t]] <- densityInputRegions$regions.r
  ACDens[[t]][] <- NA
  ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  plot(ACDens[[t]], xaxt = "n", yaxt = "n", bty = "n", zlim=c(0,maxdens), col = rev(magma(100)), legend = leg[t], main = yearnames[t])
}

dev.off()


