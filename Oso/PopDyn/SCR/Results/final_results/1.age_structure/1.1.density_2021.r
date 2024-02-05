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

# Identify cells that are left out buffer when converting to buffer
#Xbuf2_raster_ident <- Xbuf2_raster
#Xbuf2_raster_ident[] <- seq(1:length(Xbuf2_raster_ident))
#
#plot(Xbuf2_raster)
#plot(Xbuf2, border = "red", add = TRUE)
#text(Xbuf2_raster_ident, cex = 0.5)

Xbuf2_raster <- rasterize(Xbuf2, r, getCover = TRUE)
raster::values(Xbuf2_raster)[which(raster::values(Xbuf2_raster) > 0.5)] <- 2
raster::values(Xbuf2_raster)[which(raster::values(Xbuf2_raster) < 0.6)] <- 0

ras <- overlay(Xbuf_raster, Xbuf2_raster, fun = max)
f <- rasterize(Xbuf, ras, mask = TRUE)
#values(f)[which(values(f) == 1)] <- 0
#values(f)[which(values(f) == 2)] <- 1

# Load posterior distribution
library(nimbleSCR) # Load nimbleSCR here, otherwise it gets in conflict with raster package, weird

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL")
load("myResults_3-3.1_sxy.RData")

# Load original habitat coordinates
setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")
load("habcoord.RData")

#palette(c("#d9f0d3", "#e7d4e8", "#c2a5cf", "#9970ab", "#762a83"))
#palette(c("#ffffd9", "#edf8b1", "#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#253494", "#081d58"))
#palette(c("#ffffd9", "#edf8b1", "#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#253494", "#081d58"))
palette(c("#d9d9d9", "#bdbdbd", "#969696", "#737373", "#525252", "#252525", "#000000"))

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

## ---- Plot only in year 2021 ----

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Results_section/Plots")
pdf("1.1.density_scale_2021.pdf", 6,9)

par(mfrow = c(4,1),
    mar = c(0,1,0,3),
    bty = "n")


## CUBS

DensityCountriesRegions <- GetDensity_PD(
  sx = densityInputRegions$sx[,,5],# X COORDINATES
  sy =  densityInputRegions$sy[,,5],# Y COORDINATES
  z = ZZcubs[,,5],# Z 
  IDmx = densityInputRegions$habitat.id,
  aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
  regionID = densityInputRegions$regions.rgmx,
  returnPosteriorCells = F)


ACDens <- densityInputRegions$regions.r
ACDens[] <- NA
ACDens[!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions$MeanCell

#plot(g, bty = "n", xaxt = "n", yaxt = "n", border = adjustcolor("black", alpha.f = 0.1))
#plot(ACDens, xaxt = "n", yaxt = "n", bty = "n", zlim=c(0,maxden.tot), col = palette(), add = TRUE)
#plot(g, bty = "n", xaxt = "n", yaxt = "n", border = adjustcolor("black", alpha.f = 0.1), add = TRUE)
ACDens2 <- ACDens
raster::values(ACDens2)[which(raster::values(Xbuf2_raster) == 0)]  <- NA

plot(ACDens2, xaxt = "n", yaxt = "n", bty = "n", col = palette(), legend = FALSE)
r.range <- c(minValue(ACDens2), maxValue(ACDens2))
plot(ACDens2, legend.only=TRUE, col=palette(),
     legend.width = 0.5,
     axis.args=list(at = seq(r.range[1], r.range[2], by = (r.range[2] - r.range[1])/5),
                    labels = round(seq(r.range[1], r.range[2], by = (r.range[2] - r.range[1])/5),2), 
                    cex.axis = 1))

# Save for other issues (predation)

#setwd('D:/JdC/Predación/Data')
#writeRaster(ACDens2, 'cubDens2021')


## SUBADULTS

DensityCountriesRegions <- GetDensity_PD(
  sx = densityInputRegions$sx[,,5],# X COORDINATES
  sy =  densityInputRegions$sy[,,5],# Y COORDINATES
  z = ZZsub[,,5],# Z 
  IDmx = densityInputRegions$habitat.id,
  aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
  regionID = densityInputRegions$regions.rgmx,
  returnPosteriorCells = F)

ACDens <- densityInputRegions$regions.r
ACDens[] <- NA
ACDens[!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions$MeanCell

#plot(g, bty = "n", xaxt = "n", yaxt = "n", border = adjustcolor("black", alpha.f = 0.1))
#plot(ACDens, xaxt = "n", yaxt = "n", bty = "n", zlim=c(0,maxden.tot), col = palette(), add = TRUE)
#plot(g, bty = "n", xaxt = "n", yaxt = "n", border = adjustcolor("black", alpha.f = 0.1), add = TRUE)
ACDens2 <- ACDens
raster::values(ACDens2)[which(raster::values(Xbuf2_raster) == 0)]  <- NA

r.range <- c(minValue(ACDens2), 0.4)
plot(ACDens2, xaxt = "n", yaxt = "n", bty = "n", zlim=c(r.range[1],r.range[2]), col = palette(), legend = FALSE)
plot(ACDens2, legend.only=TRUE, zlim=c(r.range[1],r.range[2]), col=palette(),
     legend.width = 0.5,
     axis.args = list(at = seq(r.range[1], r.range[2], by = (r.range[2] - r.range[1])/5),
                      labels = round(seq(r.range[1], r.range[2], by = (r.range[2] - r.range[1])/5),2), 
                      cex.axis = 1))
# Save for other issues (predation)

#setwd('D:/JdC/Predación/Data')
#writeRaster(ACDens2, 'subAdDens2021')

## ADULTS

DensityCountriesRegions <- GetDensity_PD(
  sx = densityInputRegions$sx[,,5],# X COORDINATES
  sy =  densityInputRegions$sy[,,5],# Y COORDINATES
  z = ZZad[,,5],# Z 
  IDmx = densityInputRegions$habitat.id,
  aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
  regionID = densityInputRegions$regions.rgmx,
  returnPosteriorCells = F)

ACDens <- densityInputRegions$regions.r
ACDens[] <- NA
ACDens[!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions$MeanCell

#plot(g, bty = "n", xaxt = "n", yaxt = "n", border = adjustcolor("black", alpha.f = 0.1))
#plot(ACDens, xaxt = "n", yaxt = "n", bty = "n", zlim=c(0,maxden.tot), col = palette(), add = TRUE)
#plot(g, bty = "n", xaxt = "n", yaxt = "n", border = adjustcolor("black", alpha.f = 0.1), add = TRUE)
ACDens2 <- ACDens
raster::values(ACDens2)[which(raster::values(Xbuf2_raster) == 0)]  <- NA
plot(ACDens2, xaxt = "n", yaxt = "n", bty = "n", col = palette(), legend = FALSE)
r.range <- c(minValue(ACDens2), maxValue(ACDens2))
plot(ACDens2, legend.only=TRUE, col=palette(),
     legend.width = 0.5,
     axis.args=list(at = seq(r.range[1], r.range[2], by = (r.range[2] - r.range[1])/5),
                    labels = round(seq(r.range[1], r.range[2], by = (r.range[2] - r.range[1])/5),2), 
                    cex.axis = 1))

#mtext(c("Adult", "Subadult", "Cub", "All"), at = c(0.17,0.45,0.67,0.92), outer = TRUE, line = 0, side = 2, adj = 1)

# Save for other issues (predation)

#setwd('D:/JdC/Predación/Data')
#writeRaster(ACDens2, 'AdDens2021')

## ALL INDIVIDUALS

DensityCountriesRegions <- GetDensity_PD(
  sx = densityInputRegions$sx[,,5],# X COORDINATES
  sy =  densityInputRegions$sy[,,5],# Y COORDINATES
  z = ZZall[,,5],# Z 
  IDmx = densityInputRegions$habitat.id,
  aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
  regionID = densityInputRegions$regions.rgmx,
  returnPosteriorCells = F)

ACDens <- densityInputRegions$regions.r
ACDens[] <- NA
ACDens[!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions$MeanCell

#g <- crop(eur, ACDens)
#plot(g, bty = "n", xaxt = "n", yaxt = "n", border = adjustcolor("black", alpha.f = 0.1))
#plot(ACDens, xaxt = "n", yaxt = "n", bty = "n", zlim=c(0,maxden.tot), col = palette())
#plot(g, bty = "n", xaxt = "n", yaxt = "n", border = adjustcolor("black", alpha.f = 0.1), add = TRUE)
ACDens2 <- ACDens
raster::values(ACDens2)[which(raster::values(Xbuf2_raster) == 0)]  <- NA
plot(ACDens2, xaxt = "n", yaxt = "n", bty = "n", col = palette(), legend = FALSE)
r.range <- c(minValue(ACDens2), maxValue(ACDens2))
plot(ACDens2, legend.only=TRUE, col=palette(),
     legend.width = 0.5,
     axis.args=list(at = seq(r.range[1], r.range[2], by = (r.range[2] - r.range[1])/5),
                    labels = round(seq(r.range[1], r.range[2], by = (r.range[2] - r.range[1])/5),2), 
                    cex.axis = 1))

# Save for other issues (predation)

#setwd('D:/JdC/Predación/Data')
#writeRaster(ACDens2, 'AllDens2021')

dev.off()
