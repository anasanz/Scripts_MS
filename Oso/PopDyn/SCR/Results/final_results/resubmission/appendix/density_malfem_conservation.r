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
library(dplyr)
library(sf)

setwd("D:/MargSalas/Scripts_MS/Functions/Nimble")

# Load functions
#sourceCpp("GetSpaceUse_PD.cpp")
sourceCpp("GetDensity_PD.cpp")
source("getDensityInput.R")

# Load buffers
setwd("D:/MargSalas/Oso/Datos/GIS/Countries")
Xbuf <- readOGR("Buffer_statespace.shp")
Xbuf2 <- readOGR("Buffer_8500_traps_sxyObs_RESUB2.shp") # This sampling buffer includes AC of observed individuals a bit outside the trapping array

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

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/RESUBMISSION/3-3.4_allparams")
load("myResults_RESUB_3-3.4_sxy.RData")

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

# ZZall <- myResultsSXYZ$sims.list$z
# 
# ZZcubs <- myResultsSXYZ$sims.list$z
# ZZcubs[!myResultsSXYZ$sims.list$age.cat %in% c(1,2) ]  <- 3 # Considers all individuals that are not 1 and 2 as dead

ZZsub <- myResultsSXYZ$sims.list$z
ZZsub[!myResultsSXYZ$sims.list$age.cat %in% c(3,4) ]  <- 3 # Considers all individuals that are not 3 and 4 (SUBADULTS) as dead

ZZsubFEM <- myResultsSXYZ$sims.list$z
ZZsubFEM[!myResultsSXYZ$sims.list$age.cat %in% c(3,4) ]  <- 3 # Considers all individuals that are not 3 and 4 (SUBADULTS) as dead
ZZsubFEM[!myResultsSXYZ$sims.list$sex %in% c(0) ]  <- 3 # Considers all individuals that dont have a sex 0 (females) as dead

ZZsubMAL <- myResultsSXYZ$sims.list$z
ZZsubMAL[!myResultsSXYZ$sims.list$age.cat %in% c(3,4) ]  <- 3 # Considers all individuals that are not 3 and 4 (SUBADULTS) as dead
ZZsubMAL[!myResultsSXYZ$sims.list$sex %in% c(1) ]  <- 3 # Considers all individuals that dont have a sex 0 (females) as dead

ZZad <- myResultsSXYZ$sims.list$z
ZZad[!myResultsSXYZ$sims.list$age.cat %in% c(5) ]  <- 3 # Considers all individuals that are not 5 (ADULTS) as dead

ZZadFEM <- myResultsSXYZ$sims.list$z
ZZadFEM[!myResultsSXYZ$sims.list$age.cat %in% c(5) ]  <- 3 # Considers all individuals that are not 5 (ADULTS) as dead
ZZadFEM[!myResultsSXYZ$sims.list$sex %in% c(0) ]  <- 3 # Considers all individuals that dont have a sex 0 (females) as dead

ZZadMAL <- myResultsSXYZ$sims.list$z
ZZadMAL[!myResultsSXYZ$sims.list$age.cat %in% c(5) ]  <- 3 # Considers all individuals that are not 5 (ADULTS) as dead
ZZadMAL[!myResultsSXYZ$sims.list$sex %in% c(1) ]  <- 3 # Considers all individuals that dont have a sex 1 (males) as dead

zs <- list(ZZsub, ZZsubFEM , ZZsubMAL, ZZad, ZZadFEM, ZZadMAL)

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
leg <- c(FALSE, FALSE, TRUE)

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


setwd("D:/MargSalas/Oso/OPSCR_project/Results/Results_section/Resubmission/Plots")
pdf("SI_density_adsub_fm.pdf",12,6)

# Plot the 6 categories only the last year


par(mfrow = c(2,3),
    mar = c(0,0,0,0),
    oma = c(0.5,4,2,4),
    bty = "n")

######  SUBADULT INDIVIDUALS  #####

DensityCountriesRegions <- list()

n.cat = 3 # Number of categories in DensityCountriesRegions
t = 5 # Last year

# All subadults
  DensityCountriesRegions[[1]] <- GetDensity_PD(
    sx = densityInputRegions$sx[,,t],# X COORDINATES
    sy =  densityInputRegions$sy[,,t],# Y COORDINATES
    z = ZZsub[,,5],# Z last year
    IDmx = densityInputRegions$habitat.id,
    aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
    regionID = densityInputRegions$regions.rgmx,
    returnPosteriorCells = F)

# SubadultFEM
  DensityCountriesRegions[[2]] <- GetDensity_PD(
    sx = densityInputRegions$sx[,,t],# X COORDINATES
    sy =  densityInputRegions$sy[,,t],# Y COORDINATES
    z = ZZsubFEM[,,5],# Z last year
    IDmx = densityInputRegions$habitat.id,
    aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
    regionID = densityInputRegions$regions.rgmx,
    returnPosteriorCells = F)

  # SubadultMAL
  DensityCountriesRegions[[3]] <- GetDensity_PD(
    sx = densityInputRegions$sx[,,t],# X COORDINATES
    sy =  densityInputRegions$sy[,,t],# Y COORDINATES
    z = ZZsubMAL[,,5],# Z last year
    IDmx = densityInputRegions$habitat.id,
    aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
    regionID = densityInputRegions$regions.rgmx,
    returnPosteriorCells = F)
  
# Maximum density of the three categories

maxdens_all <- max(DensityCountriesRegions[[1]]$MeanCell)
maxdens_femmal <- max(max(DensityCountriesRegions[[2]]$MeanCell), max(DensityCountriesRegions[[3]]$MeanCell))

  
#plot() 
ACDens <- list()
ACDens2 <- list()

c = 1 # For all subadults
  ACDens[[c]] <- densityInputRegions$regions.r
  ACDens[[c]][] <- NA
  ACDens[[c]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[c]]$MeanCell
  
  ACDens2[[c]] <- ACDens[[c]]
  raster::values(ACDens2[[c]])[which(raster::values(Xbuf2_raster) == 0)]  <- NA
  
  plot(ACDens2[[c]], xaxt = "n", yaxt = "n", bty = "n", zlim=c(0,maxdens_all), col = palette(), legend = TRUE)
  plot(st_geometry(esp), border = adjustcolor("white", alpha.f = 0.5),  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000), add = TRUE)
  plot(st_geometry(fr), border = adjustcolor("white", alpha.f = 0.5),  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000), add = TRUE)
  
  mtext(c("Subadults"), line = -1, side = 2)
  

for(c in 2:n.cat){
  ACDens[[c]] <- densityInputRegions$regions.r
  ACDens[[c]][] <- NA
  ACDens[[c]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[c]]$MeanCell
  
  ACDens2[[c]] <- ACDens[[c]]
  raster::values(ACDens2[[c]])[which(raster::values(Xbuf2_raster) == 0)]  <- NA
  
  plot(ACDens2[[c]], xaxt = "n", yaxt = "n", bty = "n", zlim=c(0,maxdens_femmal), col = palette(), legend = leg[[c]])
  plot(st_geometry(esp), border = adjustcolor("white", alpha.f = 0.5),  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000), add = TRUE)
  plot(st_geometry(fr), border = adjustcolor("white", alpha.f = 0.5),  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000), add = TRUE)
}


######  ADULT INDIVIDUALS  #####

DensityCountriesRegions <- list()

n.cat = 3 # Number of categories in DensityCountriesRegions

# All adults
DensityCountriesRegions[[1]] <- GetDensity_PD(
  sx = densityInputRegions$sx[,,t],# X COORDINATES
  sy =  densityInputRegions$sy[,,t],# Y COORDINATES
  z = ZZad[,,5],# Z last year
  IDmx = densityInputRegions$habitat.id,
  aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
  regionID = densityInputRegions$regions.rgmx,
  returnPosteriorCells = F)

# adultFEM
DensityCountriesRegions[[2]] <- GetDensity_PD(
  sx = densityInputRegions$sx[,,t],# X COORDINATES
  sy =  densityInputRegions$sy[,,t],# Y COORDINATES
  z = ZZadFEM[,,5],# Z last year
  IDmx = densityInputRegions$habitat.id,
  aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
  regionID = densityInputRegions$regions.rgmx,
  returnPosteriorCells = F)

# adultMAL
DensityCountriesRegions[[3]] <- GetDensity_PD(
  sx = densityInputRegions$sx[,,t],# X COORDINATES
  sy =  densityInputRegions$sy[,,t],# Y COORDINATES
  z = ZZadMAL[,,5],# Z last year
  IDmx = densityInputRegions$habitat.id,
  aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
  regionID = densityInputRegions$regions.rgmx,
  returnPosteriorCells = F)

# Maximum density of the three categories

maxdens_all <- max(DensityCountriesRegions[[1]]$MeanCell)
maxdens_femmal <- max(max(DensityCountriesRegions[[2]]$MeanCell), max(DensityCountriesRegions[[3]]$MeanCell))


#plot() 
ACDens <- list()
ACDens2 <- list()

c = 1 # For all subadults
ACDens[[c]] <- densityInputRegions$regions.r
ACDens[[c]][] <- NA
ACDens[[c]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[c]]$MeanCell

ACDens2[[c]] <- ACDens[[c]]
raster::values(ACDens2[[c]])[which(raster::values(Xbuf2_raster) == 0)]  <- NA

plot(ACDens2[[c]], xaxt = "n", yaxt = "n", bty = "n", zlim=c(0,maxdens_all), col = palette(), legend = TRUE)
plot(st_geometry(esp), border = adjustcolor("white", alpha.f = 0.5),  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000), add = TRUE)
plot(st_geometry(fr), border = adjustcolor("white", alpha.f = 0.5),  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000), add = TRUE)

mtext(c("Adults"), line = -1, side = 2)


for(c in 2:n.cat){
  ACDens[[c]] <- densityInputRegions$regions.r
  ACDens[[c]][] <- NA
  ACDens[[c]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[c]]$MeanCell
  
  ACDens2[[c]] <- ACDens[[c]]
  raster::values(ACDens2[[c]])[which(raster::values(Xbuf2_raster) == 0)]  <- NA
  
  plot(ACDens2[[c]], xaxt = "n", yaxt = "n", bty = "n", zlim=c(0,maxdens_femmal), col = palette(), legend = leg[[c]])
  plot(st_geometry(esp), border = adjustcolor("white", alpha.f = 0.5),  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000), add = TRUE)
  plot(st_geometry(fr), border = adjustcolor("white", alpha.f = 0.5),  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000), add = TRUE)
}



mtext(c("All", "Females", "Males"), at = c(0.2,0.5, 0.8), outer = TRUE, line = -2, side = 3)

dev.off()


