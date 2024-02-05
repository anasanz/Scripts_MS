## -------------------------------------------------
##    Plot same density values each age category
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

# Load buffers
setwd("D:/MargSalas/Oso/Datos/GIS/Countries")
Xbuf <- readOGR("Buffer_statespace.shp")
Xbuf2 <- readOGR("Buffer_8500_traps_sxyObs.shp") # This sampling buffer includes AC of observed individuals a bit outside the trapping array

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


# Load posterior distribution
library(nimbleSCR) # Load nimbleSCR here, otherwise it gets in conflict with raster package, weird

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL")
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

densityInputRegions <- getDensityInput( 
  regions = f,## THIS  A RASTER FILE WITH 0/1 HABITAT VS BUFFER, OR 1/2 fRANCE SPAIN... WHATEVER YOU WANT. 
  habitat = f1,## here put the same than regions argument. 
  s = myResultsSXYZ$sims.list$sxy,
  plot.check = TRUE)

# Get z of each age class

ZZall <- myResultsSXYZ$sims.list$z

ZZcubs <- myResultsSXYZ$sims.list$z
ZZcubs[!myResultsSXYZ$sims.list$age.cat %in% c(1,2) ]  <- 3 # Considers all individuals that are not 1 and 2 as dead

ZZsub <- myResultsSXYZ$sims.list$z
ZZsub[!myResultsSXYZ$sims.list$age.cat %in% c(3,4) ]  <- 3 # Considers all individuals that are not 3 and 4 (SUBADULTS) as dead

ZZad <- myResultsSXYZ$sims.list$z
ZZad[!myResultsSXYZ$sims.list$age.cat %in% c(5) ]  <- 3 # Considers all individuals that are not 3 and 4 (SUBADULTS) as dead

zs <- list(ZZall, ZZcubs, ZZsub, ZZad)

# Get maximum density of all to scale
maxden <- list()
for(i in 1:length(zs)){
  DensityCountriesRegions <- GetDensity_PD(
    sx = densityInputRegions$sx[,,5],# X COORDINATES
    sy =  densityInputRegions$sy[,,5],# Y COORDINATES
    z = zs[[i]][,,5],# Z 
    IDmx = densityInputRegions$habitat.id,
    aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
    regionID = densityInputRegions$regions.rgmx,
    returnPosteriorCells = F)
  maxden[[i]] <- max(DensityCountriesRegions$MeanCell)
}

maxden.tot <- max(unlist(maxden))


## extract density
yearnames <- c("2017", "2018", "2019", "2020", "2021")
leg <- c(FALSE, FALSE, FALSE, FALSE, TRUE)

palette(c("#d9d9d9", "#bdbdbd", "#969696", "#737373", "#525252", "#252525", "#000000"))

# Load europe
sa <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/clip_pyros2_WGS84_31N_all.shp")
eur <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/esp_fr_2.shp") %>%
  st_transform(dpts, crs = crs(sa))
#and <- sa[which(sa$NAME_0 == "Andorra"),]
#eur2 <- st_union(and,eur) 

esp <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/ESP_adm/ESP_adm0.shp") %>%
  st_transform(dpts, crs = crs(sa))
fr <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/FRA_adm/FRA_adm0.shp") %>%
  st_transform(dpts, crs = crs(sa))


setwd("D:/MargSalas/Oso/OPSCR_project/Results/Results_section/Plots")
pdf("SI_density_scale.pdf",12,6)


######  ALL INDIVIDUALS  #####

par(mfrow = c(4,5),
    mar = c(0,0,0,0),
    oma = c(0.5,4,2,4),
    bty = "n")

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

# Maximum density of all years

maxdens <- max(max(DensityCountriesRegions[[1]]$MeanCell), max(DensityCountriesRegions[[2]]$MeanCell), max(DensityCountriesRegions[[3]]$MeanCell),
               max(DensityCountriesRegions[[4]]$MeanCell), max(DensityCountriesRegions[[5]]$MeanCell))

#plot()
ACDens<-list()
for(t in 1:n.years){
  ACDens[[t]] <- densityInputRegions$regions.r
  ACDens[[t]][] <- NA
  ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  
  ACDens2[[t]] <- ACDens[[t]]
  raster::values(ACDens2[[t]])[which(raster::values(Xbuf2_raster) == 0)]  <- NA
  
  plot(ACDens2[[t]], xaxt = "n", yaxt = "n", bty = "n", zlim=c(0,maxdens), col = palette(), legend = leg[t])
  plot(st_geometry(esp), border = adjustcolor("white", alpha.f = 0.5),  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000), add = TRUE)
  plot(st_geometry(fr), border = adjustcolor("white", alpha.f = 0.5),  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000), add = TRUE)
  
}


######  CUBS  #####

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

# Maximum density of all years

maxdens <- max(max(DensityCountriesRegions[[1]]$MeanCell), max(DensityCountriesRegions[[2]]$MeanCell), max(DensityCountriesRegions[[3]]$MeanCell),
               max(DensityCountriesRegions[[4]]$MeanCell), max(DensityCountriesRegions[[5]]$MeanCell))

ACDens<-list()
for(t in 1:n.years){
  ACDens[[t]] <- densityInputRegions$regions.r
  ACDens[[t]][] <- NA
  ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  
  ACDens2[[t]] <- ACDens[[t]]
  raster::values(ACDens2[[t]])[which(raster::values(Xbuf2_raster) == 0)]  <- NA
  
  plot(ACDens2[[t]], xaxt = "n", yaxt = "n", bty = "n", zlim=c(0,maxdens), col = palette(), legend = leg[t])
  plot(st_geometry(esp), border = adjustcolor("white", alpha.f = 0.5),  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000), add = TRUE)
  plot(st_geometry(fr), border = adjustcolor("white", alpha.f = 0.5),  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000), add = TRUE)
  
}

######  SUBADULTS  #####

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


# Maximum density of all years

maxdens <- max(max(DensityCountriesRegions[[1]]$MeanCell), max(DensityCountriesRegions[[2]]$MeanCell), max(DensityCountriesRegions[[3]]$MeanCell),
               max(DensityCountriesRegions[[4]]$MeanCell), max(DensityCountriesRegions[[5]]$MeanCell))

ACDens<-list()
for(t in 1:n.years){
  ACDens[[t]] <- densityInputRegions$regions.r
  ACDens[[t]][] <- NA
  ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  
  ACDens2[[t]] <- ACDens[[t]]
  raster::values(ACDens2[[t]])[which(raster::values(Xbuf2_raster) == 0)]  <- NA
  
  plot(ACDens2[[t]], xaxt = "n", yaxt = "n", bty = "n", zlim=c(0,maxdens), col = palette(), legend = leg[t])
  plot(st_geometry(esp), border = adjustcolor("white", alpha.f = 0.5),  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000), add = TRUE)
  plot(st_geometry(fr), border = adjustcolor("white", alpha.f = 0.5),  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000), add = TRUE)

  }


######  ADULTS  #####

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


# Maximum density of all years
maxdens <- max(max(DensityCountriesRegions[[1]]$MeanCell), max(DensityCountriesRegions[[2]]$MeanCell), max(DensityCountriesRegions[[3]]$MeanCell),
               max(DensityCountriesRegions[[4]]$MeanCell), max(DensityCountriesRegions[[5]]$MeanCell))

#plot()
ACDens<-list()
for(t in 1:n.years){
  ACDens[[t]] <- densityInputRegions$regions.r
  ACDens[[t]][] <- NA
  ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  
  ACDens2[[t]] <- ACDens[[t]]
  raster::values(ACDens2[[t]])[which(raster::values(Xbuf2_raster) == 0)]  <- NA
  
  plot(ACDens2[[t]], xaxt = "n", yaxt = "n", bty = "n", zlim=c(0,maxdens), col = palette(), legend = leg[t])
  plot(st_geometry(esp), border = adjustcolor("white", alpha.f = 0.5),  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000), add = TRUE)
  plot(st_geometry(fr), border = adjustcolor("white", alpha.f = 0.5),  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000), add = TRUE)

  }

mtext(yearnames, at = c(0.1,0.3,0.5,0.7,0.9), outer = TRUE, line = -1, side = 3)
mtext(c("Adults", "Subadults", "Juveniles", "All"), at = c(0.17,0.45,0.67,0.92), outer = TRUE, line = 0, side = 2, adj = 1)

dev.off()


