## -------------------------------------------------
##               Check trap location
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


## -------------------------------------------------
##                      2017
## ------------------------------------------------- 

#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Effort/France/Syst_Opport")
setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Oso/PopDyn/SCR/Effort/France/Syst_Opport")
#setwd("C:/Users/cymi/Downloads/Ana")


#---- 1. LOAD THE DETECTION DATA ---- 
load("edf2017_2019_fr.RData")
load("tdf2017_fr.RData")

edf <- edf[which(edf$session == 1), ] # Keep only detections of 2017 (year 1)

# We remove Nere and Goiat, which are two individuals moving a lot that don't represent
# the rest of the population.

edf <- edf[-which(edf$ind %in% c("Néré", "Goiat")), ]

#---- 2. DEFINE THE TRAP AND THE HABITAT  ---- 
#----   2.1 GET TRAPS---- 
X <- tdf2017[,c(2,3)]
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

setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
forest <- raster("forest_hrbear.tif")
habitat.r <- crop(forest, e) 

#----   DEFINE THE BUFFER AREA AND CUT WHAT IS NOT HABITAT---- 
# Buffer around traps (5*sigma = 33200)
Xpoints <- X
coordinates(Xpoints) <- Xpoints[,c(1,2)]
Xbuf <- gBuffer(Xpoints, width = 25000)

forestMask <- rasterize(Xbuf, habitat.r, mask = TRUE)

#RETAIN HABITAT COORDINATES THAT ARE WITHIN THE HABITAT
G <- coordinates(forestMask)[!is.na(forestMask[]),]
colnames(G) <- c("x","y")
#PLOT CHECK 
plot(forestMask)
points(Xpoints, pch = 19)

#----  RESCALE HABITAT AND TRAP COORDINATES ---- 
###scale X and G so that bottom left corner of state space is origin (0,0)
sc.coord <- scaleCoordsToHabitatGrid(coordsData = X,
                                     coordsHabitatGridCenter = G)

##this returns S also in row order (all x for a given y) but starting top left 
## corner 
G.sc <- sc.coord$coordsHabitatGridCenterScaled
X.sc <- sc.coord$coordsDataScaled

###get cell coordinates for G.sc
windowCoords <- getWindowCoords(G.sc)

#----  DETECTION DATA   ---- 
#----  MAKE Y   ---- 

K <- 7 # 7 occasions
n <- length(unique(edf$ind))

Y <- matrix(0, nrow = n, ncol = dim(X.sc)[1])
rownames(Y) <- unique(edf$ind)
xx <- edf[,c(2,4)]

for (obs in 1:nrow(xx)) {
  Y[xx[obs, 1], xx[obs, 2]] <- Y[xx[obs, 1], xx[obs, 2]] + 1
}

# Plot all locations

loc <- list()
for (i in 1:n){
  caps <- which(Y[i, ]>0)
  loc[[i]] <- as.matrix(X.sc[caps,])
  }
all.loc <- do.call("rbind", loc)

# Plot activity centers starting values (average location)
S.in <- matrix(NA, n, 2)
for ( i in 1:n){
  caps <- which(Y[i, ]>0)
  if (length(caps)==1){
    S.in[i,] <- as.matrix(X.sc[caps,])
  }else{
    S.in[i,] <- apply(X.sc[caps,], 2, mean)}
}

# Differenciate opportunistic traps
X.sc.opp <- X.sc[which(tdf2017$suivi == "opportunist"),]
X.sc.sist <- X.sc[which(tdf2017$suivi != "opportunist"),]

points(X.sc.opp, pch = 19, col = "grey")
points(X.sc.sist, pch = 19)
points(all.loc, pch = 19, col = "violet")
points(S.in, pch = 19, col = "yellow")

## -------------------------------------------------
##                      2018
## ------------------------------------------------- 


rm(list = ls())


#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Effort/France/Syst_Opport")
setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Oso/PopDyn/SCR/Effort/France/Syst_Opport")
#setwd("C:/Users/cymi/Downloads/Ana")


#---- 1. LOAD THE DETECTION DATA ---- 
load("edf2017_2019_fr.RData")
load("tdf2018_fr.RData")

edf <- edf[which(edf$session == 2), ] # Keep only detections of 2017 (year 1)

# We remove Nere and Goiat, which are two individuals moving a lot that don't represent
# the rest of the population.

edf <- edf[-which(edf$ind %in% c("Néré", "Goiat")), ]

#---- 2. DEFINE THE TRAP AND THE HABITAT  ---- 
#----   2.1 GET TRAPS---- 
X <- tdf2018[,c(2,3)]
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

setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
forest <- raster("forest_hrbear.tif")
habitat.r <- crop(forest, e) 

#----   DEFINE THE BUFFER AREA AND CUT WHAT IS NOT HABITAT---- 
# Buffer around traps (5*sigma = 33200)
Xpoints <- X
coordinates(Xpoints) <- Xpoints[,c(1,2)]
Xbuf <- gBuffer(Xpoints, width = 25000)

forestMask <- rasterize(Xbuf, habitat.r, mask = TRUE)

#RETAIN HABITAT COORDINATES THAT ARE WITHIN THE HABITAT
G <- coordinates(forestMask)[!is.na(forestMask[]),]
colnames(G) <- c("x","y")
#PLOT CHECK 
#plot(forestMask)
#points(Xpoints, pch = 19)

#----  RESCALE HABITAT AND TRAP COORDINATES ---- 
###scale X and G so that bottom left corner of state space is origin (0,0)
sc.coord <- scaleCoordsToHabitatGrid(coordsData = X,
                                     coordsHabitatGridCenter = G)

##this returns S also in row order (all x for a given y) but starting top left 
## corner 
G.sc <- sc.coord$coordsHabitatGridCenterScaled
X.sc <- sc.coord$coordsDataScaled

###get cell coordinates for G.sc
windowCoords <- getWindowCoords(G.sc)

#----  DETECTION DATA   ---- 
#----  MAKE Y   ---- 

K <- 7 # 7 occasions
n <- length(unique(edf$ind))

Y <- matrix(0, nrow = n, ncol = dim(X.sc)[1])
rownames(Y) <- unique(edf$ind)
xx <- edf[,c(2,4)]

for (obs in 1:nrow(xx)) {
  Y[xx[obs, 1], xx[obs, 2]] <- Y[xx[obs, 1], xx[obs, 2]] + 1
}

# Plot all locations

loc <- list()
for (i in 1:n){
  caps <- which(Y[i, ]>0)
  loc[[i]] <- as.matrix(X.sc[caps,])
}
all.loc <- do.call("rbind", loc)

# Plot activity centers starting values (average location)
S.in <- matrix(NA, n, 2)
for ( i in 1:n){
  caps <- which(Y[i, ]>0)
  if (length(caps)==1){
    S.in[i,] <- as.matrix(X.sc[caps,])
  }else{
    S.in[i,] <- apply(X.sc[caps,], 2, mean)}
}

# Differenciate opportunistic traps
X.sc.opp <- X.sc[which(tdf2018$suivi == "opportunist"),]
X.sc.sist <- X.sc[which(tdf2018$suivi != "opportunist"),]

points(X.sc.opp, pch = 19, col = "grey")
points(X.sc.sist, pch = 19)
points(all.loc, pch = 19, col = "violet")
points(S.in, pch = 19, col = "yellow")


## -------------------------------------------------
##                      2019
## ------------------------------------------------- 


rm(list = ls())


#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Effort/France/Syst_Opport")
setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Oso/PopDyn/SCR/Effort/France/Syst_Opport")
#setwd("C:/Users/cymi/Downloads/Ana")


#---- 1. LOAD THE DETECTION DATA ---- 
load("edf2017_2019_fr.RData")
load("tdf2019_fr.RData")

edf <- edf[which(edf$session == 3), ] # Keep only detections of 2017 (year 1)

# We remove Nere and Goiat, which are two individuals moving a lot that don't represent
# the rest of the population.

edf <- edf[-which(edf$ind %in% c("Néré", "Goiat")), ]

#---- 2. DEFINE THE TRAP AND THE HABITAT  ---- 
#----   2.1 GET TRAPS---- 
X <- tdf2019[,c(2,3)]
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

setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
forest <- raster("forest_hrbear.tif")
habitat.r <- crop(forest, e) 

#----   DEFINE THE BUFFER AREA AND CUT WHAT IS NOT HABITAT---- 
# Buffer around traps (5*sigma = 33200)
Xpoints <- X
coordinates(Xpoints) <- Xpoints[,c(1,2)]
Xbuf <- gBuffer(Xpoints, width = 25000)

forestMask <- rasterize(Xbuf, habitat.r, mask = TRUE)

#RETAIN HABITAT COORDINATES THAT ARE WITHIN THE HABITAT
G <- coordinates(forestMask)[!is.na(forestMask[]),]
colnames(G) <- c("x","y")
#PLOT CHECK 
#plot(forestMask)
#points(Xpoints, pch = 19)

#----  RESCALE HABITAT AND TRAP COORDINATES ---- 
###scale X and G so that bottom left corner of state space is origin (0,0)
sc.coord <- scaleCoordsToHabitatGrid(coordsData = X,
                                     coordsHabitatGridCenter = G)

##this returns S also in row order (all x for a given y) but starting top left 
## corner 
G.sc <- sc.coord$coordsHabitatGridCenterScaled
X.sc <- sc.coord$coordsDataScaled

###get cell coordinates for G.sc
windowCoords <- getWindowCoords(G.sc)

#----  DETECTION DATA   ---- 
#----  MAKE Y   ---- 

K <- 7 # 7 occasions
n <- length(unique(edf$ind))

Y <- matrix(0, nrow = n, ncol = dim(X.sc)[1])
rownames(Y) <- unique(edf$ind)
xx <- edf[,c(2,4)]

for (obs in 1:nrow(xx)) {
  Y[xx[obs, 1], xx[obs, 2]] <- Y[xx[obs, 1], xx[obs, 2]] + 1
}

# Plot all locations

loc <- list()
for (i in 1:n){
  caps <- which(Y[i, ]>0)
  loc[[i]] <- as.matrix(X.sc[caps,])
}
all.loc <- do.call("rbind", loc)

# Plot activity centers starting values (average location)
S.in <- matrix(NA, n, 2)
for ( i in 1:n){
  caps <- which(Y[i, ]>0)
  if (length(caps)==1){
    S.in[i,] <- as.matrix(X.sc[caps,])
  }else{
    S.in[i,] <- apply(X.sc[caps,], 2, mean)}
}

# Differenciate opportunistic traps
X.sc.opp <- X.sc[which(tdf2019$suivi == "opportunist"),]
X.sc.sist <- X.sc[which(tdf2019$suivi != "opportunist"),]

points(X.sc.opp, pch = 19, col = "grey")
points(X.sc.sist, pch = 19)
points(all.loc, pch = 19, col = "violet")
points(S.in, pch = 19, col = "yellow")



