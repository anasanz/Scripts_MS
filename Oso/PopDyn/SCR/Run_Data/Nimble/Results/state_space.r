## -------------------------------------------------
##             Save state space without buffer
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
#----   2.1 GET TRAPS---- 
X <- tdf_all[,c(2,3)]
colnames(X) <- c('x', 'y')
J <- dim(X)[1]

#----   2.2 DEFINE STATE SPACE EXTENT ---- 
# State space coordinates
# Buffer: 25000 (used by Maelis, also ~3*sigma in pre-analysis where sig = 6640)
xmin <- min(X[,1]) - 25000
ymin <- min(X[,2]) - 25000
xmax <- max(X[,1]) + 25000
ymax <- max(X[,2]) + 25000
e <- as(raster::extent(xmin, xmax, ymin, ymax), "SpatialPolygons") # Extent of state space

#----   2.3 GET A RASTER FOR THE HABITAT ---- 
# USE A FOREST RASTER TO GET A BASIS FOR THE HABITAT RASTER
# Set up a raster file 
setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
#setwd("~/Data_server/Variables_hrscale")
distcore <- raster("logDistcore_hrbear.tif")
# Crop it to extent of state-space
habitat.r <- crop(distcore, e) 
plot(habitat.r)

#----   2.4 DEFINE THE BUFFER AREA, DECIDE WHICH AREA IS CONSIDERED NOT BUFFER ---- 
# Buffer around traps (5*sigma = 33200)
Xpoints <- X
coordinates(Xpoints) <- Xpoints[,c(1,2)]
Xbuf <- gBuffer(Xpoints, width = 25000)

# Convex hull
hpts <- chull(coordinates(Xpoints))
hpts <- c(hpts, hpts[1])

# Smaller buffer
Xbuf2 <- gBuffer(Xpoints, width = 8500)

# Plots to decide
plot(e)
plot(Xbuf, col = "lightblue", add = TRUE)
points(Xpoints)
lines(X[hpts, ], col = "red")
plot(Xbuf2, col = adjustcolor("pink", alpha = 0.5), add = TRUE)


#RETAIN HABITAT COORDINATES THAT ARE WITHIN THE HABITAT
distcoreMask <- rasterize(Xbuf, habitat.r, mask = TRUE)

G <- coordinates(distcoreMask)[!is.na(distcoreMask[]),]
colnames(G) <- c("x","y")
#PLOT CHECK 
plot(distcoreMask)
points(Xpoints)
points(G[,2]~G[,1],col="red",cex=0.5)

#----   2.5 RESCALE HABITAT AND TRAP COORDINATES ---- 
###scale X and G so that bottom left corner of state space is origin (0,0)
sc.coord <- scaleCoordsToHabitatGrid(coordsData = X,
                                     coordsHabitatGridCenter = G)

##this returns S also in row order (all x for a given y) but starting top left 
## corner 
G.sc <- sc.coord$coordsHabitatGridCenterScaled
X.sc <- sc.coord$coordsDataScaled

#----   2.6 UNSCALE SXY COORDINATES ---- 

# Load results
setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/3.openSCRdenscov_Age/2021/Cyril/3-3.1_sxy_6000")
load("myResults_3-3.1_sxy_6000.RData")


# first add "x", "y" to the dimension of the sxy

dimnames(myResults$sims.list$sxy)[[3]] <- c("x", "y")

myResults$sims.list$sxy <- scaleCoordsToHabitatGrid(coordsData = myResults$sims.list$sxy,## this are your sxy
                                                          
                                                          coordsHabitatGridCenter = G,# this is your unscaled habitat (as you used when scaling the habitat/detector to the habitat. G?
                                                          
                                                          scaleToGrid = FALSE)$coordsDataScaled



