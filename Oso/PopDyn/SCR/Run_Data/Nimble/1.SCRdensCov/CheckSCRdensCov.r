
## -------------------------------------------------
##                 Check model SCRdenscov
## ------------------------------------------------- 

## Check which step is wrong in the SCRhab model with density covariate (it gives wrong estimates of forest effec and N)

#######################################################↕
##############LOAD DATA

rm(list = ls())

library(nimble)
library(MCMCvis)
library(nimbleSCR)
library(raster)
library(rgeos)
library(oSCR)
library(terra)
library(sp)

#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Effort/France/Syst_Opport")
setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Oso/PopDyn/SCR/Effort/France/Syst_Opport")
load("edf2017_2019_fr.RData")
load("tdf2017_fr.RData")

edf <- edf[which(edf$session == 1), ] # Keep only detections of 2017 (year 1)

# We remove Nere and Goiat, which are two individuals moving a lot that don't represent
# the rest of the population.

edf <- edf[-which(edf$ind %in% c("Néré", "Goiat")), ]

## ---- Trap array, state space, habitat covariate ----

##traps
X <- tdf2017[,c(2,3)]
colnames(X) <- c('x', 'y')
J <- dim(X)[1]


# State space coordinates
# Buffer: 25000 (used by Maelis, also ~3*sigma in pre-analysis where sig = 6640)
xmin <- min(X[,1])-25000
ymin <- min(X[,2])-25000
xmax <- max(X[,1])+25000
ymax <- max(X[,2])+25000

##for discrete state space, create 2100 grid cells (center coordinates) of 5x5 km res
#gx <- rep(seq(xmin+2500, xmax-2500,5000),30)
#gy <- rep(seq(ymin+2500, ymax-2500,5000), each=70)

# Set up a raster file 

setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
#setwd("~/Data_server")
forest <- raster("forest_hrbear.tif")

e <- as(raster::extent(xmin, xmax, ymin, ymax), "SpatialPolygons") # Extent of state space
fcr <- crop(forest, e) # Crop it to extent of state-space

# Buffer around traps (5*sigma = 33200)
Xpoints <- X
coordinates(Xpoints) <- Xpoints[,c(1,2)]
Xbuf <- gBuffer(Xpoints,width = 33200)

forestMask <- rasterize(Xbuf, fcr, mask = TRUE)
plot(forestMask)
values(forestMask)

G <- coordinates(forestMask)[!is.na(forestMask[]),]

plot(G)
points(X,col = "red")


###scale X and G so that bottom left corner of state space is origin (0,0)
sc.coord <- scaleCoordsToHabitatGrid(X, G)

##this returns S also in row order (all x for a given y) but starting top left 
## corner 
G.sc <- sc.coord$coordsHabitatGridCenterScaled
X.sc <- sc.coord$coordsDataScaled

###get cell coordinates for G.sc
windowCoords <- getWindowCoords(G.sc)
habitatGrid <- windowCoords$habitatGrid


# Habitat covariate
X.d <- values(forestMask)
X.d <- X.d[complete.cases(X.d)]

#forestMask_matrix <- as.matrix(forestMask)
#X.d <- as.vector(t(forestMask_matrix))

# Habitat mask
habitatMask_rast <- rasterize(Xbuf, fcr)

habitatMask <- as.matrix(habitatMask_rast)
habitatMask[is.na(habitatMask)] <- 0

##determine which traps are within some threshold distance of each habitat grid 
## cell - speeds up computations by only evaluating traps at which an individual
## could have been caught at (rather that traps very far away)

##binary mask that determines which cells in S are suitable (dmax is 5*sigma)
X.sc <- as.matrix(X.sc)

# 3 times log sigma? It would be 26, but I set it to 31 so that it includes all
(5*6640)/5000

localTraps <- getLocalTraps(habitatMask, X.sc, resizeFactor = 1, dmax = 7)

localTrapsIndex <- localTraps$localTrapsIndices
localTrapsNum <- localTraps$numLocalTraps
habitatGridDet <- localTraps$habitatGrid
MaxLocalTraps <- localTraps$numLocalTrapsMax

##name and structure for Nimble model
numHabWindows <- dim(localTrapsIndex)[1] #number of cells in S
numGridRows <- nrow(localTraps$habitatGrid)
numGridCols <- ncol(localTraps$habitatGrid)


#the following two are still passed to function but no longer used by it
lowerHabCoords <- windowCoords$lowerHabCoords
upperHabCoords <- windowCoords$upperHabCoords

## ---- Detection data ----

K <- 7 # 7 occasions
n <- length(unique(edf$ind))


Y <- matrix(0, nrow = n, ncol = dim(X.sc)[1])
rownames(Y) <- unique(edf$ind)
xx <- edf[,c(2,4)]

for (obs in 1:nrow(xx)) {
  Y[xx[obs, 1], xx[obs, 2]] <- Y[xx[obs, 1], xx[obs, 2]] + 1
}

##augment observed data to size M
M <- 400 
y.in <- rbind(Y, 
              matrix(0, nrow = M-n, ncol = J))

#change to 'sparse' format - speeds up computation by reducing file size
y.sparse <- getSparseY(y.in)

##extract pieces to be passed to Nimble

detNums <- y.sparse$detNums # Nº of traps at which each individual was detected
maxDetNums <- y.sparse$maxDetNums # Maximun Nº of traps at which any individual was detected
detIndices <- y.sparse$detIndices # ID of the traps where they were detected
y.sp <- y.sparse$y # Number of detections at each trap


###############################################################################
########## MODEL

## Arguments
beta.dens <- 3.18 # From oSCR
habDens <- X.d
lowerCoords <- lowerHabCoords
upperCoords <- upperHabCoords
habitatGrid

## ---- 1. Generate mu1 (density surface) from the beta estimate obtained in oSCR ----

mu1 <- NULL
mu1[1:numHabWindows] <- exp(beta.dens * habDens[1:numHabWindows]) # Expected dens in each cell

mu1rast <- forestMask
values(mu1rast)[!is.na(values(mu1rast))] <- mu1 # Fill the habitat mask with the values of mu1

# mu1 with scaled coordinates to plot

ex <- matrix(c(0, 70, 0, 30), nrow = 2, ncol = 2, byrow = T)
ex <- extent(ex)
r <- raster(nrows = 30, ncols = 70)
r <- setExtent(r, ex, keepres=F)
crs(r) <- proj4string(forestMask)
values(r)[!is.na(values(mu1rast))] <- mu1

# To calculate logIntensity (transform mu1 into relative probabilities)
logHabIntensity <- NULL 
sumHabIntensity <- sum(mu1[1:numHabWindows]) # Sum of expected densities
logHabIntensity[1:numHabWindows] <- log(mu1[1:numHabWindows]) # Backtransformed expected density
logSumHabIntensity <- log(sumHabIntensity) # Back-transformed sum of expected densities?


## ---- 2. Locate AC with the point process function dbernppAC ----

# 1. Workflow

# Place an AC in a window accordying to the logHabIntensity

# Convert mu1 into relative probabilities per cell
p.cell <- mu1/sum(mu1)

logIntensities <- logHabIntensity/logSumHabIntensity # This the same as pcell in data sim (p.cell <- log(mu1)/log(sum(mu1)))
logIntensities <- log(mu1)/log(sum(mu1))# Used in the model, but it doesn't sum up to one (IT DOESN'T MATTER; rcat contemplates it internally)

# Simulate the allocation of one AC in a window
#####????? In the model it is used exp(logIntensities) within rbern, but it doesn't sum up to one?
windowInd <- rcat(1, p.cell)
numDims <- 2 

# Get exact coordinates of AC. Uniform distribution is used within the target window so that is not right in center
outCoordinates <- lowerCoords[windowInd,] + 
  runif(numDims, 0.0, 1.0) * (upperCoords[windowInd,] - lowerCoords[windowInd,])

dfcoord <- data.frame(x = outCoordinates[1], y = outCoordinates[2])

par(mfrow = c(1,2)) # Plot it
plot(r)
points(dfcoord) 
plot(r, xlim = c(dfcoord[1,1]-3,dfcoord[1,1]+3), ylim = c(dfcoord[1,2]-3,dfcoord[1,2]+3))
points(dfcoord) 

mu1AC <- extract(r,dfcoord, exact = TRUE)

# 2. Check if the model works properly

######## CELLS WITH LOW PROB
# Take a cells with very low density probability (between 1st quantile and median f.e.)
cellsLowp <- which(p.cell >= summary(p.cell)[2] & p.cell <= summary(p.cell)[3])
set.seed(3)
cLow <- sample(cellsLowp,1) 
p.cell[cLow] # Relative probability in this window
mu1[cLow] # Mean expected density in this window

# Get exact coordinates of AC. Uniform distribution is used within the target window so that is not right in center
outCoordinates <- lowerCoords[cLow,] + 
  runif(numDims, 0.0, 1.0) * (upperCoords[cLow,] - lowerCoords[cLow,])
dfcoord <- data.frame(x = outCoordinates[1], y = outCoordinates[2])

par(mfrow = c(1,2)) # Plot it
plot(r)
points(dfcoord) 
plot(r, xlim = c(dfcoord[1,1]-3,dfcoord[1,1]+3), ylim = c(dfcoord[1,2]-3,dfcoord[1,2]+3))
points(dfcoord) 

mu1AC <- extract(r,dfcoord, exact = TRUE)


######## CELLS WITH HIGH PROB
# Take a cells with very low density probability (between 1st quantile and median f.e.)
cellsHighp <- which(p.cell >= summary(p.cell)[5] & p.cell <= summary(p.cell)[6])
set.seed(3)
cHigh <- sample(cellsHighp,1) 
p.cell[cHigh] # Relative probability in this window
mu1[cHigh] # Mean expected density in this window

# Get exact coordinates of AC. Uniform distribution is used within the target window so that is not right in center
outCoordinates <- lowerCoords[cHigh,] + 
  runif(numDims, 0.0, 1.0) * (upperCoords[cHigh,] - lowerCoords[cHigh,])
dfcoord <- data.frame(x = outCoordinates[1], y = outCoordinates[2])

par(mfrow = c(1,2)) # Plot it
plot(r)
points(dfcoord) 
plot(r, xlim = c(dfcoord[1,1]-3,dfcoord[1,1]+3), ylim = c(dfcoord[1,2]-3,dfcoord[1,2]+3))
points(dfcoord) 

mu1AC <- extract(r,dfcoord, exact = TRUE)


c(which(habitatGridDet==cHigh, arr.ind=TRUE)[2], 
  which(habitatGridDet==cHigh, arr.ind=TRUE)[1]) # Problem matching because the lower and upper cords are from the total hab.grid?

###############################################################################
########## LOCATION ACTIVITY CENTERS

##because of local evaluation of possible detectors, activity center initial 
##values have to be specified, eg average capture location
S.in<- matrix(NA, n, 2)
for ( i in 1:n){
  caps<-which(Y[i, ]>0)
  if (length(caps)==1){
    S.in[i,] <- as.matrix(X.sc[caps,])
  }else{
    S.in[i,]<-apply(X.sc[caps,],2,mean)}
}

# First activity center:
S.in[1,1] # Column hab grid
S.in[1,2] # Row hab grid
habitatGrid[trunc(S.in[1,2]),trunc(S.in[1,1])]

# Plot
par(mfrow = c(1,1))
plot(r)
points(S.in)




habitatGrid
windowCoords$lowerHabCoords
windowCoords$upperHabCoords
windowCoords$habitatGrid



