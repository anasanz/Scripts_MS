
## -------------------------------------------------
##               Map density from model 3.1  
## ------------------------------------------------- 

rm(list = ls())

library(Rcpp)
library(RcppArmadillo)
library(RcppProgress)
library(rgdal)
library(raster)

setwd("D:/MargSalas/Scripts_MS/Functions/Nimble")

# Load functions
#sourceCpp("GetSpaceUse_PD.cpp")
sourceCpp("GetDensity_PD.cpp")
source("getDensityInput.R")

# Load buffers
setwd("D:/MargSalas/Oso/Datos/GIS/Countries")
Xbuf <- readOGR("Buffer_statespace.shp")
Xbuf2 <- readOGR("Buffer_8500_traps.shp")

# Conver to rasters

setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
r <- raster("logDistcore_hrbear.tif")
values(r) <- NA
r <- crop(r,Xbuf) 

Xbuf_raster <- rasterize(Xbuf, r)
values(Xbuf_raster)[which(is.na(values(Xbuf_raster)))] <- 0 # Raster of 0 and 1

Xbuf2_raster <- rasterize(Xbuf2, r)
values(Xbuf2_raster)[which(values(Xbuf2_raster) == 1)] <- 2
values(Xbuf2_raster)[which(is.na(values(Xbuf2_raster)))] <- 0 # Raster of 0 and 2

ras <- overlay(Xbuf_raster, Xbuf2_raster, fun = max)
f <- rasterize(Xbuf, ras, mask = TRUE)
#values(f)[which(values(f) == 1)] <- 0
#values(f)[which(values(f) == 2)] <- 1

# Load posterior distribution
library(nimbleSCR) # Load nimbleSCR here, otherwise it gets in conflict with raster package, weird

setwd("D:/MargSalas/Oso/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams")
load("myResults_3-3.1_sxy.RData")

# Load original habitat coordinates
setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")
load("habcoord.RData")

#give dim names to your posteriors sxy
dim(myResultsSXYZ$sims.list$sxy) # YOURPOSTERIORSXY 
dimnames(myResultsSXYZ$sims.list$sxy)[[3]] <- c("x","y")

## first rescale the coordinates to the original scale 
myResultsSXYZ$sims.list$sxy <- scaleCoordsToHabitatGrid(coordsData = myResultsSXYZ$sims.list$sxy,
  coordsHabitatGridCenter = G,
  scaleToGrid = FALSE)$coordsDataScaled

f1 <- f
f1[f1%in%c(1,2)] <-  1
f[] <- as.factor(as.character(f[]))

##GET OBJECTS IN SHAPE
densityInputRegions <- getDensityInput( 
  regions = f,## THIS  A RASTER FILE WITH 0/1 HABITAT VS BUFFER, OR 1/2 fRANCE SPAIN... WHATEVER YOU WANT. 
  habitat = f1,## here put the same than regions argument. 
  s = myResultsSXYZ$sims.list$sxy,
  plot.check = TRUE)

## extract density
yearnames <- c("2017", "2018", "2019", "2020", "2021")

######  ALL INDIVIDUALS  #####

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

#mean density in each cell 
DensityCountriesRegions[[1]]$MeanCell

#plot()
ACDens<-list()
for(t in 1:n.years){
  ACDens[[t]] <- densityInputRegions$regions.r
  ACDens[[t]][] <- NA
  ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  plot(ACDens[[t]], main = paste(yearnames[t], "ALL"))
}

DensityCountriesRegions[[1]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
DensityCountriesRegions[[2]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
DensityCountriesRegions[[3]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
DensityCountriesRegions[[4]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
DensityCountriesRegions[[5]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)

######  CUBS  #####

ZZcubs <- myResultsSXYZ$sims.list$z
ZZcubs[!myResultsSXYZ$sims.list$age.cat %in% c(1,2) ]  <- 3 # Considers all individuals that are not 1 and 2 as dead

DensityCountriesRegions <- list()
n.years = 5
for(t in 1:n.years){
  DensityCountriesRegions[[t]] <- GetDensity_PD(
    sx = densityInputRegions$sx[,,t],# X COORDINATES
    sy =  densityInputRegions$sy[,,t],# Y COORDINATES
    z = ZZcubs[,,t],# Z 
    IDmx = densityInputRegions$habitat.id,
    aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
    regionID = densityInputRegions$regions.rgmx,
    returnPosteriorCells = F)
}#t

#mean density in each cell 
DensityCountriesRegions[[1]]$MeanCell

#plot()
ACDens<-list()
for(t in 1:n.years){
  ACDens[[t]] <- densityInputRegions$regions.r
  ACDens[[t]][] <- NA
  ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  plot(ACDens[[t]], main = paste(yearnames[t], "CUBS"))
}

DensityCountriesRegions[[1]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
DensityCountriesRegions[[2]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
DensityCountriesRegions[[3]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
DensityCountriesRegions[[4]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
DensityCountriesRegions[[5]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)

######  SUBADULTS  #####

ZZsub <- myResultsSXYZ$sims.list$z
ZZsub[!myResultsSXYZ$sims.list$age.cat %in% c(3,4) ]  <- 3 # Considers all individuals that are not 3 and 4 (SUBADULTS) as dead

DensityCountriesRegions <- list()
n.years = 5
for(t in 1:n.years){
  DensityCountriesRegions[[t]] <- GetDensity_PD(
    sx = densityInputRegions$sx[,,t],# X COORDINATES
    sy =  densityInputRegions$sy[,,t],# Y COORDINATES
    z = ZZsub[,,t],# Z 
    IDmx = densityInputRegions$habitat.id,
    aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
    regionID = densityInputRegions$regions.rgmx,
    returnPosteriorCells = F)
}#t

#mean density in each cell 
DensityCountriesRegions[[1]]$MeanCell

#plot()
ACDens<-list()
for(t in 1:n.years){
  ACDens[[t]] <- densityInputRegions$regions.r
  ACDens[[t]][] <- NA
  ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  plot(ACDens[[t]], main = paste(yearnames[t], "SUBADULTS"))
}

DensityCountriesRegions[[1]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
DensityCountriesRegions[[2]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
DensityCountriesRegions[[3]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
DensityCountriesRegions[[4]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
DensityCountriesRegions[[5]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)


######  ADULTS  #####

ZZad <- myResultsSXYZ$sims.list$z
ZZad[!myResultsSXYZ$sims.list$age.cat %in% c(5) ]  <- 3 # Considers all individuals that are not 5 (adults) as dead

DensityCountriesRegions <- list()
n.years = 5
for(t in 1:n.years){
  DensityCountriesRegions[[t]] <- GetDensity_PD(
    sx = densityInputRegions$sx[,,t],# X COORDINATES
    sy =  densityInputRegions$sy[,,t],# Y COORDINATES
    z = ZZad[,,t],# Z 
    IDmx = densityInputRegions$habitat.id,
    aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
    regionID = densityInputRegions$regions.rgmx,
    returnPosteriorCells = F)
}#t

#mean density in each cell 
DensityCountriesRegions[[1]]$MeanCell

#plot()
ACDens<-list()
for(t in 1:n.years){
  ACDens[[t]] <- densityInputRegions$regions.r
  ACDens[[t]][] <- NA
  ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  plot(ACDens[[t]], main = paste(yearnames[t], "ADULTS"))
}

DensityCountriesRegions[[1]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
DensityCountriesRegions[[2]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
DensityCountriesRegions[[3]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
DensityCountriesRegions[[4]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
DensityCountriesRegions[[5]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)

######  FEMALES  #####
# Sex: 0 Females; 1 Males
ZZfem <- myResultsSXYZ$sims.list$z
ZZfem[!myResultsSXYZ$sims.list$sex %in% c(0) ]  <- 3 # Considers all individuals that are not 0 (females) as dead

DensityCountriesRegions <- list()
n.years = 5
for(t in 1:n.years){
  DensityCountriesRegions[[t]] <- GetDensity_PD(
    sx = densityInputRegions$sx[,,t],# X COORDINATES
    sy =  densityInputRegions$sy[,,t],# Y COORDINATES
    z = ZZfem[,,t],# Z 
    IDmx = densityInputRegions$habitat.id,
    aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
    regionID = densityInputRegions$regions.rgmx,
    returnPosteriorCells = F)
}#t

#mean density in each cell 
DensityCountriesRegions[[1]]$MeanCell

#plot()
ACDens<-list()
for(t in 1:n.years){
  ACDens[[t]] <- densityInputRegions$regions.r
  ACDens[[t]][] <- NA
  ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  plot(ACDens[[t]], main = paste(yearnames[t], "FEMALES"))
}

DensityCountriesRegions[[1]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
DensityCountriesRegions[[2]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
DensityCountriesRegions[[3]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
DensityCountriesRegions[[4]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
DensityCountriesRegions[[5]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)

######  MALES  #####
# Sex: 0 Females; 1 Males
ZZmal <- myResultsSXYZ$sims.list$z
ZZmal[!myResultsSXYZ$sims.list$sex %in% c(1) ]  <- 3 # Considers all individuals that are not 0 (females) as dead

DensityCountriesRegions <- list()
n.years = 5
for(t in 1:n.years){
  DensityCountriesRegions[[t]] <- GetDensity_PD(
    sx = densityInputRegions$sx[,,t],# X COORDINATES
    sy =  densityInputRegions$sy[,,t],# Y COORDINATES
    z = ZZmal[,,t],# Z 
    IDmx = densityInputRegions$habitat.id,
    aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
    regionID = densityInputRegions$regions.rgmx,
    returnPosteriorCells = F)
}#t

#mean density in each cell 
DensityCountriesRegions[[1]]$MeanCell

#plot()
ACDens<-list()
for(t in 1:n.years){
  ACDens[[t]] <- densityInputRegions$regions.r
  ACDens[[t]][] <- NA
  ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  plot(ACDens[[t]], main = paste(yearnames[t], "MALES"))
}

DensityCountriesRegions[[1]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
DensityCountriesRegions[[2]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
DensityCountriesRegions[[3]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
DensityCountriesRegions[[4]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
DensityCountriesRegions[[5]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)

######  ONLY OBSERVED INDIVIDUALS  #####

ZZad <- myResultsSXYZ$sims.list$z
ZZad[!myResultsSXYZ$sims.list$age.cat %in% c(5) ]  <- 3 # Considers all individuals that are not 5 (SUBADULTS) as dead

DensityCountriesRegions <- list()
n.years = 5
for(t in 1:n.years){
  DensityCountriesRegions[[t]] <- GetDensity_PD(
    sx = densityInputRegions$sx[,,t],# X COORDINATES
    sy =  densityInputRegions$sy[,,t],# Y COORDINATES
    z = ZZad[,,t],# Z 
    IDmx = densityInputRegions$habitat.id,
    aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
    regionID = densityInputRegions$regions.rgmx,
    returnPosteriorCells = F)
}#t

#mean density in each cell 
DensityCountriesRegions[[1]]$MeanCell

#plot()
ACDens<-list()
for(t in 1:n.years){
  ACDens[[t]] <- densityInputRegions$regions.r
  ACDens[[t]][] <- NA
  ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  plot(ACDens[[t]], main = paste(yearnames[t], "ADULTS"))
}

DensityCountriesRegions[[1]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
DensityCountriesRegions[[2]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
DensityCountriesRegions[[3]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
DensityCountriesRegions[[4]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
DensityCountriesRegions[[5]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
