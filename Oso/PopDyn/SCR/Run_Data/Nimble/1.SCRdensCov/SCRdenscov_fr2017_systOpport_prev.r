## -------------------------------------------------
##                      SCR + denscov
##                      French 2017
##                      Syst + Opport
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

################################################################################
#### fit Nimble model ##########################################################

##compile constants
nimConstants <- list(
  M = M, J = J, numHabWindows = numHabWindows, 
  numGridRows = numGridRows, numGridCols = numGridCols, 
  maxDetNums = maxDetNums, MaxLocalTraps = MaxLocalTraps
)

##compile data
nimData <- list(habDens = X.d, y = y.sp, detNums = detNums, 
                lowerHabCoords = lowerHabCoords, upperHabCoords = upperHabCoords,
                habitatGrid = habitatGrid, K = rep(K,J),  X.sc = X.sc,
                habitatGridDet = habitatGridDet, detIndices = detIndices,
                detNums = detNums, localTrapsIndex = localTrapsIndex, 
                localTrapsNum = localTrapsNum 
)
str(nimConstants)
##set up initial values
z.in <- c(rep(1, n), rep(0, M-n))

##because of local evaluation of possible detectors, activity center initial 
##values have to be specified, eg average capture location
S.in<- matrix(NA, M, 2)
for ( i in 1:n){
  caps<-which(Y[i, ]>0)
  if (length(caps)==1){
    S.in[i,] <- as.matrix(X.sc[caps,])
  }else{
    S.in[i,]<-apply(X.sc[caps,],2,mean)}
}


##random ACs for individuals never observed
##simulate them within state space!!
cellACv <- list() # To monitor the cells where there was a problem with initial values
for(i in (n+1) : M){
  cellAC <- sample(c(1:max(habitatGridDet)),1)
  cellACv[[i]] <- cellAC
  S.in[i,] <- c(which(habitatGridDet==cellAC, arr.ind=TRUE)[2] - runif(1, 0.0, 1.0), # The column of the habitatGridDet is the x coordinate
                which(habitatGridDet==cellAC, arr.ind=TRUE)[1] - runif(1, 0.0, 1.0)) # The row of the habitatGridDet is the y coordinate
}

inits<-function(){list(psi=runif(1,0.4, 0.6), 
                       sigma=runif(1,0.5, 1.5),
                       p0=runif(1,0,0.1),
                       z=z.in,
                       sxy=S.in)}

##source model code
##I prefer working on code in a separate script but you can also have everything in
##one script and just execute the code

setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Oso/PopDyn/SCR/Sim/SCRdenscov")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Sim/SCRdenscov")
source('SCR in Nimble.R')

##determine which parameters to monitor
params<-c('N', 'psi', 'sigma', 'p0', 'beta.dens')

#(1) set up model

model <- nimbleModel(SCRhab, constants = nimConstants, 
                     data=nimData, check = FALSE)
##ignore error message, only due to missing initial values at this stage

#(2) Compile model in c++
#     In complex models, this step can take a while (as well as step 5)
#     Much longer than in JAGS, but the model typically runs much faster
cmodel <- compileNimble(model)       

# (3) Configure MCMC - on an uncompiled model - this step allows setting which quantities to monitor
#     Also, nimble allows two sets of monitors, these can be thinned at different rates
#     all of which is more important in complex models but not to start with
conf.mcmc<-configureMCMC(model, monitors = params, thin=1)

# (4) Build the MCMC sampler based on configurations
mcmc <- buildMCMC(conf.mcmc)

# (5) Compile sampler in c++ together with compiled model
cmcmc <- compileNimble(mcmc, project = cmodel, resetFunctions = TRUE)

# (6) Run (monitor time just for fun) [takes 20 seconds on my computer]
system.time(
  (samp <- runMCMC(cmcmc, niter = 15000, nburnin = 5000, nchains=3, inits = inits) )
)

##NOTE: summary command is from MCMCvis package; that also has good plotting options
## summary table for everything in "params" vector
summ <- MCMCsummary(samp)
MCMCtrace(samp)

setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Oso/PopDyn/SCR/Data/Nimble/Results")
save(samp, file = "sampSCRdenscov_habgrid.RData")

