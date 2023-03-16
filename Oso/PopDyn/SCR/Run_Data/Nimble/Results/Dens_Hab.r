

## -------------------------------------------------
##       Density - Habitats (expansion front)
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

# Load habitats corine
setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Corine_Landcover_2018_100m/CORINE_FINAL_CORRECTED_ElenaPi")
hab <- raster("HABITATS_CORINE_RECLASS_100m_TOT_WGS8431N.tif")

# Load buffers
setwd("D:/MargSalas/Oso/Datos/GIS/Countries")
Xbuf <- readOGR("Buffer_statespace.shp")
Xbuf2 <- readOGR("Buffer_8500_traps.shp")

## ---- 1. Relate density to habitats ----

# Convert density values to raster files

# Raster to store
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
#values(f)[which(values(f) == 1)] <- 0
#values(f)[which(values(f) == 2)] <- 1

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

##GET OBJECTS IN SHAPE
densityInputRegions <- getDensityInput( 
  regions = f,## THIS  A RASTER FILE WITH 0/1 HABITAT VS BUFFER, OR 1/2 fRANCE SPAIN... WHATEVER YOU WANT. 
  habitat = f1,## here put the same than regions argument. 
  s = myResultsSXYZ$sims.list$sxy,
  plot.check = TRUE)

## extract density
yearnames <- c("2017", "2018", "2019", "2020", "2021")

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

#mean density in each cell per year
DensityCountriesRegions[[1]]$MeanCell # YEAR 1

plot(densityInputRegions$regions.r)

# Maximum density of all years

maxdens <- max(max(DensityCountriesRegions[[1]]$MeanCell), max(DensityCountriesRegions[[2]]$MeanCell), max(DensityCountriesRegions[[3]]$MeanCell),
               max(DensityCountriesRegions[[4]]$MeanCell), max(DensityCountriesRegions[[5]]$MeanCell))

#plot()
ACDens<-list()
for(t in 1:n.years){
  ACDens[[t]] <- densityInputRegions$regions.r
  ACDens[[t]][] <- NA
  ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  plot(ACDens[[t]], xaxt = "n", yaxt = "n", bty = "n", zlim=c(0,maxdens), col = rev(magma(100)))
}

#setwd("D:/MargSalas/Oso/OPSCR_project/Results/Density layer")
#for(t in 1:n.years){
#  writeRaster(ACDens[[t]], filename = paste("Densisty_all_", t, ".tif", sep = ""))
#  }

# Layerize habitat (one layer per type)

proj4string(hab) == proj4string(ACDens[[1]])

lcor <- layerize(hab, classes=c(1:12), bylayer=TRUE, suffix='numbers')
hablist <- list()
for (i in 1:12){
  hablist[[i]] <- raster(lcor, layer = i)
}

# Resample to obtain same resolution (habitat to resolution of dens)

hablist_res <- list()
for(i in 1:12){
  hablist_res[[i]] <- resample(hablist[[i]], ACDens[[1]], method = 'bilinear')
}

for(i in 1:12){
  hablist_res[[i]][is.na(ACDens[[1]][])] <- NA
}

# Raster stack of all layers and years
resampld <- c(ACDens, hablist_res)
st <- do.call('stack',resampld)

# Correlation matrix
# Pearson
jnk <- layerStats(st, 'pearson', na.rm=T)
corr_matrix <- jnk$'pearson correlation coefficient'

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Density")
corr_matrix2 <- as.data.frame(corr_matrix[6:17,1:5])
colnames(corr_matrix2) <- yearnames
corr_matrix2$habitat <- c("Artificial surfaces", "Agricultural areas", "Broad-leaved forest",
    "Coniferous forest", "Mixed forest", "Natural grasslands",  "Shrubs",  
    "Other open spaces (beaches, dunes, burnt areas, glaciers)", 
    "Bare rocks", "Wetlands",  "Water bodies",  "Sparse vegetation")
corr_matrix2 <- corr_matrix2[,c(6,1:5)]
#openxlsx::write.xlsx(corr_matrix2, 'corr_matrix.xlsx')

#Spearman (pairwise)
cor1 <- cor(sampleRandom(st, size= ncell(ACDens[[5]]) ), method = "spearman")

st2021_conif <- do.call('stack',c(ACDens[[5]], hablist_res[[4]]))
cor2 <- cor(sampleRandom(st2021_conif, size= ncell(ACDens[[5]]) ), method = "spearman")
df <- corrplot(cor2, method = "number")

# Plots 2020-2021 and forests habitats

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Plots/model3.1")

pdf("Density_Broad-leaved forest.pdf", 9,2)
par(mfrow = c(1,3),
    mar = c(0,2,0,3),
    oma = c(0.5,3,2,1))

plot(ACDens[[4]], xaxt = "n", yaxt = "n", bty = "n", zlim=c(0,maxdens), col = rev(magma(100)), legend = TRUE)
mtext(c("2020"), line = 0.5, side = 3)
plot(hablist_res[[3]], xaxt = "n", yaxt = "n" )
mtext(c("Broad-leaved forest"), line = 0.5, side = 3)
plot(ACDens[[5]], xaxt = "n", yaxt = "n", bty = "n", zlim=c(0,maxdens), col = rev(magma(100)), legend = TRUE)
mtext(c("2021"), line = 0.5, side = 3)

dev.off()


pdf("Density_Coniferous forest.pdf", 9,2)
par(mfrow = c(1,3),
    mar = c(0,2,0,3),
    oma = c(0.5,3,2,1))

plot(ACDens[[4]], xaxt = "n", yaxt = "n", bty = "n", zlim=c(0,maxdens), col = rev(magma(100)), legend = TRUE)
mtext(c("2020"), line = 0.5, side = 3)
plot(hablist_res[[4]], xaxt = "n", yaxt = "n" )
mtext(c("Coniferous forest"), line = 0.5, side = 3)
plot(ACDens[[5]], xaxt = "n", yaxt = "n", bty = "n", zlim=c(0,maxdens), col = rev(magma(100)), legend = TRUE)
mtext(c("2021"), line = 0.5, side = 3)

dev.off()

pdf("Density_Mixed forest.pdf", 9,2)
par(mfrow = c(1,3),
    mar = c(0,2,0,3),
    oma = c(0.5,3,2,1))

plot(ACDens[[4]], xaxt = "n", yaxt = "n", bty = "n", zlim=c(0,maxdens), col = rev(magma(100)), legend = TRUE)
mtext(c("2020"), line = 0.5, side = 3)
plot(hablist_res[[5]], xaxt = "n", yaxt = "n" )
mtext(c("Mixed forest"), line = 0.5, side = 3)
plot(ACDens[[5]], xaxt = "n", yaxt = "n", bty = "n", zlim=c(0,maxdens), col = rev(magma(100)), legend = TRUE)
mtext(c("2021"), line = 0.5, side = 3)

dev.off()

## ---- 2. Relate distcore to habitats and relevant bear variables ----


# Load
setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
distcore <- raster("logDistcore_hrbear.tif")
forest <- raster("forest_hrbear.tif")
dem <- raster("dem_hrbear.tif")
rough <- raster("rough_hrbear.tif")
slope <- raster("slope_hrbear.tif")
roads1 <- raster("roads1_hrbear.tif")
roads4 <- raster("roads4_hrbear.tif")
roads5 <- raster("roads5_hrbear.tif")
roads6 <- raster("roads6_hrbear.tif")

covs <- c(distcore, forest, dem, rough, slope,
          roads1, roads4, roads5, roads6)

# Same extent as state space
covs_crop <- lapply(covs, crop, y = Xbuf)

# Same values as state space
covs_ss <- list()
for (i in 1:length(covs)){
  covs_ss[[i]] <- rasterize(Xbuf, covs_crop[[i]], mask = TRUE)
}

# Resample to obtain same resolution (habitat to resolution of dens)
covs_ss_res <- list()
for(i in 1:length(covs)){
  covs_ss_res[[i]] <- resample(covs_ss[[i]], ACDens[[1]], method = 'bilinear')
}

# Stack
cov_st <- do.call('stack',covs_ss_res)
cor1 <- cor(sampleRandom(cov_st, size= ncell(ACDens[[5]])), method = "spearman") 
# It looks like distcore is slightly correlated with roughness and slope
# Pairwise correlation test:
t1 <- cor.test(x = sampleRandom(cov_st, size= ncell(ACDens[[5]]))[,c(1)], 
         y = sampleRandom(cov_st, size= ncell(ACDens[[5]]))[,c(4)], method = "spearman")
t1$p.value # Significant

# Distcore - Corine forest categories
for_st <- do.call('stack',c(covs_ss_res[[1]], hablist_res[[3]], hablist_res[[4]], hablist_res[[5]], ACDens[[5]]))
cor2 <- cor(sampleRandom(for_st, size= ncell(ACDens[[5]])), method = "spearman") 


plot(covs_ss[[3]])
  crop(r,Xbuf)
distcore <- rasterize(Xbuf, r, mask = TRUE)
compareRaster(distcore, ACDens[[4]])

