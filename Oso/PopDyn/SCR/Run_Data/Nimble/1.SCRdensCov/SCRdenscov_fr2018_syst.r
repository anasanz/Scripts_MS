## -------------------------------------------------
##                      SCR + denscov
##                  French Data 2018 (fr+Sp)
##                        ONLY Syst
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
##                      FOREST
## ------------------------------------------------- 


setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Effort/France/Data/Systematic")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Effort/France/Data/Systematic")


#---- 1. LOAD THE DETECTION DATA ---- 
load("edf2017_2019_fr.RData")
load("tdf2018_fr.RData")

edf <- edf[which(edf$session == 2), ] # Keep only detections of 2018 (year 2)

# As the model is set, there can be only one capture per trap and occasion
# --> The number of trials is 7 (From may to November): So max number of captures per trap = 7
# To fix it in this edf, remove duplicates
edf <- edf[-which(duplicated(edf)), ]
t <- edf %>% group_by(ind,occ, trap) %>% # All need to be one
  summarise(n()) 

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
# Set up a raster file 
setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
#setwd("~/Data_server/Variables_hrscale")
forest <- raster("forest_hrbear.tif")
# Crop it to extent of state-space
habitat.r <- crop(forest, e) 

#----   2.4 DEFINE THE BUFFER AREA AND CUT WHAT IS NOT HABITAT---- 
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

###get cell coordinates for G.sc
windowCoords <- getWindowCoords(G.sc)
habitatGrid <- windowCoords$habitatGrid

#----   2.6 HABITAT COVARIATES ---- 
#[CM] THIS MIGHT BE WHERE THE ISSUE WAS
# AS YOU DID FOR THE COORDINATES THE HABITAT, YOU NEED TO SELECT HAB COORDS THAT ARE CONSIDERED AS HABITAT
X.d <- values(forestMask)[!is.na(forestMask[])]

# Scale
X.d_mean <- mean(X.d)
X.d_sd <- sd(X.d)
X.d_sc <- (X.d - X.d_mean) / X.d_sd

#----   2.7 GET THE LOCAL DETECTORS OBJECTS  ---- 
# USE THE HABITAT GRID PROVIDED BY GETWINDOWCOORDS
habitatMask <- habitatGrid
habitatMask[habitatMask>0] <- 1 #TURN CELL ID TO 1 TO DEFINE THE HABITAT

##determine which traps are within some threshold distance of each habitat grid 
## cell - speeds up computations by only evaluating traps at which an individual
## could have been caught at (rather that traps very far away)

# 3 times log sigma? It would be 26, but I set it to 31 so that it includes all
(5*6640)/5000
# Im gonna try 3*sigma like in the local evaluation (I have to put 5 because it doesn't allow 4)
(3*6640)/5000
# THERE MIGHT BE SOME TRIAL AND ERROR WHEN DEFINING THE DMAX.
localTraps <- getLocalTraps(habitatMask, X.sc, resizeFactor = 1, dmax = 5)

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

#---- 3. DETECTION DATA   ---- 
#----   3.1 MAKE Y   ---- 

K <- 7 # 7 occasions
n <- length(unique(edf$ind))

Y <- matrix(0, nrow = n, ncol = dim(X.sc)[1])
rownames(Y) <- unique(edf$ind)
xx <- edf[,c(2,4)]

for (obs in 1:nrow(xx)) {
  Y[xx[obs, 1], xx[obs, 2]] <- Y[xx[obs, 1], xx[obs, 2]] + 1
}

max(Y) # Check the number max number of detections (it can't be higher than K)


#----   3.2 AUGMENT Y   ---- 
##augment observed data to size M
M <- 400 
y.in <- rbind(Y, 
              matrix(0, nrow = M-n, ncol = J))

#----   3.3 USE SPARSE FORMAT FOR Y   ---- 
#change to 'sparse' format - speeds up computation by reducing file size
y.sparse <- getSparseY(y.in)

##extract pieces to be passed to Nimble
detNums <- y.sparse$detNums # Nº of traps at which each individual was detected
maxDetNums <- y.sparse$maxDetNums # Maximun Nº of traps at which any individual was detected
detIndices <- y.sparse$detIndices # ID of the traps where they were detected
y.sp <- y.sparse$y # Number of detections at each trap


#---- 4. FIT NIMBLE MODEL    ---- 
#----   4.1 CONSTANT AND DATA    ---- 

##compile constants
nimConstants <- list(
  M = M,
  J = J,
  numHabWindows = numHabWindows, 
  numGridRows = numGridRows,
  numGridCols = numGridCols, 
  maxDetNums = maxDetNums,
  MaxLocalTraps = MaxLocalTraps
)

##compile data
nimData <- list(habDens = X.d_sc,
                y = y.sp,
                detNums = detNums, 
                lowerHabCoords = lowerHabCoords,
                upperHabCoords = upperHabCoords,
                habitatGrid = habitatGrid,
                K = rep(K,J),
                X.sc = X.sc,
                habitatGridDet = habitatGridDet,
                detIndices = detIndices,
                detNums = detNums,
                localTrapsIndex = localTrapsIndex, 
                localTrapsNum = localTrapsNum 
)

#----   4.2 INITIAL VALUES  ---- 
#----     4.2.1 Z   ---- 

##set up initial values
z.in <- c(rep(1, n), rep(0, M-n))

#----     4.2.2 SXY  ---- 
##because of local evaluation of possible detectors, activity center initial 
##values have to be specified, eg average capture location
S.in <- matrix(NA, M, 2)
for ( i in 1:n){
  caps <- which(Y[i, ]>0)
  if (length(caps)==1){
    S.in[i,] <- as.matrix(X.sc[caps,])
  }else{
    S.in[i,] <- apply(X.sc[caps,], 2, mean)}
}


##random ACs for individuals never observed
##simulate them within state space!!
#cellACv <- list() # To monitor the cells where there was a problem with initial values
for(i in (n+1) : M){
  #[CM] THE ALTERNATIVE IS TO USE THE SIMULATION FUNCTIONNALITY OF NIMBLE TO SIMULATE INITIAL AC
  # IT IS THE SAME FUNCTION THAT YOU USE TO FIT THE MODEL, BUT IT STARTS WITH "r"
  S.in[i,] <- rbernppAC(n=1,
                        lowerCoords = nimData$lowerHabCoords,
                        upperCoords = nimData$upperHabCoords,
                        habitatGrid = habitatGrid,
                        numGridRows = nimConstants$numGridRows,
                        numGridCols = nimConstants$numGridCols,
                        logIntensities = rep(1, nimConstants$numHabWindows),#assume equal intensity across habitat cells (all 1)
                        logSumIntensity = sum(rep(1, nimConstants$numHabWindows))
  )
  # 
  # cellAC <- sample(c(1:max(habitatGridDet)),1)
  # cellACv[[i]] <- cellAC
  # S.in[i,] <- c(which(habitatGridDet==cellAC, arr.ind=TRUE)[2] - runif(1, 0.0, 1.0), # The column of the habitatGridDet is the x coordinate
  #               which(habitatGridDet==cellAC, arr.ind=TRUE)[1] - runif(1, 0.0, 1.0)) # The row of the habitatGridDet is the y coordinate
}


#----     4.2.2 COMPILE INITIAL VALUES   ---- 
inits <- function(){list(psi=runif(1,0.4, 0.6), 
                         sigma=runif(1,0.5, 1.5),
                         p0=runif(1,0,0.1),
                         z=z.in,
                         sxy=S.in,
                         beta.dens = runif(1,-0.1, 0.1))}

##source model code
##I prefer working on code in a separate script but you can also have everything in
##one script and just execute the code

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Sim/2.SCRdenscov")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Run_Sim/2.SCRdenscov")

source('SCR in Nimble.R')

##determine which parameters to monitor
params <- c('N', 'psi', 'sigma', 'p0', 'beta.dens')

#---- 5. FIT THE MODEL ---- 
#[CM] INITIAL VALUES ARE IMPORTANT IN NIMBLE AND WITH NIMBLE, PROVIDE THEM
#(1) set up model
model <- nimbleModel(SCRhab,
                     constants = nimConstants, 
                     data=nimData,
                     inits = inits(),
                     check = FALSE)
##ignore error message, only due to missing initial values at this stage
# [CM] I WOULD NOT IGNORE THE ERROR MESSAGE. FROM OUR EXPERIENCE, A SCR MODEL IN NIMBLE THAT IS 
# NOT PROPERLY INITIALIZED MAY NOT RUN CORRECTLY. TO CHECK WHICH PARAMETERS ARE NOT PROPERLY INITIALIZE 
# YOU CAN DO:
model$initializeInfo()
#THEN YOU SHOULD ALWAYS CHECK THAT THE MODEL IS ABLE TO CALCULATE A LIKELIHOOD (RETURN A VALUE) GIVEN THE INITIAL VALUES PROVIDED
model$calculate()#
# IF A -INF OF NA IS RETURNED YOU CAN CHECK WHERE THE PROBLEM COMES FROM WITH (AND THEN TRY TO FIX IT UNTIL THE -INF DISAPEARS):
model$logProb_p0# WILL GIVE YOU LIKELIHOOD OF P0
model$logProb_sxy# THE PROBLEM IS OFTEN WITH SXY OR Y
model$logProb_y

#(2) Compile model in c++
#     In complex models, this step can take a while (as well as step 5)
#     Much longer than in JAGS, but the model typically runs much faster
cmodel <- compileNimble(model)       
#[CM] YOU CAN CALCULATE THE LIKELIHOOD WITH THE COMPILED VERSION, 
#     IT SHOULD GIVE THE SAME VALUE THAN ABOVE (c IS A LOT MORE PICKY SO SOMETIMES IT WORKS IN R BUT NOT IN C...)
cmodel$calculate()
# (3) Configure MCMC - on an uncompiled model - this step allows setting which quantities to monitor
#     Also, nimble allows two sets of monitors, these can be thinned at different rates
#     all of which is more important in complex models but not to start with
conf.mcmc <- configureMCMC(model, monitors = params, thin=1)

# (4) Build the MCMC sampler based on configurations
mcmc <- buildMCMC(conf.mcmc)

# (5) Compile sampler in c++ together with compiled model
cmcmc <- compileNimble(mcmc, project = cmodel, resetFunctions = TRUE)

# (6) Run (monitor time just for fun) [takes 20 seconds on my computer]
system.time(
  (samp <- runMCMC(cmcmc, niter = 25000, nburnin = 1000, nchains=3, inits = inits) )
)

# I DONT HAVE MCMCvis SO HERE IT IS WITH basicMCMCplots
library(basicMCMCplots)
chainsPlot(samp)
# SEEMS TO BE OKAY NOW BETA IS MIXING OKAY. 

##NOTE: summary command is from MCMCvis package; that also has good plotting options
## summary table for everything in "params" vector
summ <- MCMCsummary(samp)
summ
MCMCtrace(samp)

sig <- summ$mean[5]*5000 #meters


setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/SCRdenscov_year/Forest")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/SCRdenscov_year")

save(samp, file = "sampSCRdenscov2018.RData")

load("sampSCRdenscov2018.RData")
summ <- MCMCsummary(samp)

# Density
summ$mean[1]/max(habitatGrid)

## -------------------------------------------------
##                      ELEVATION
## ------------------------------------------------- 

rm(list = ls())

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Effort/France/Data/Systematic")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Effort/France/Data/Systematic")

#---- 1. LOAD THE DETECTION DATA ---- 
load("edf2017_2019_fr.RData")
load("tdf2018_fr.RData")

edf <- edf[which(edf$session == 2), ] 

# As the model is set, there can be only one capture per trap and occasion
# --> The number of trials is 7 (From may to November): So max number of captures per trap = 7
# To fix it in this edf, remove duplicates
edf <- edf[-which(duplicated(edf)), ]

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

# Set up a raster file 
setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
#setwd("~/Data_server/Variables_hrscale")
elevation <- raster("dem_hrbear.tif")

# Crop it to extent of state-space
habitat.r <- crop(elevation, e) 

#----   2.4 DEFINE THE BUFFER AREA AND CUT WHAT IS NOT HABITAT---- 
# Buffer around traps (5*sigma = 33200)
Xpoints <- X
coordinates(Xpoints) <- Xpoints[,c(1,2)]
Xbuf <- gBuffer(Xpoints, width = 25000)

elevMask <- rasterize(Xbuf, habitat.r, mask = TRUE)

#RETAIN HABITAT COORDINATES THAT ARE WITHIN THE HABITAT
G <- coordinates(elevMask)[!is.na(elevMask[]),]
colnames(G) <- c("x","y")
#PLOT CHECK 
plot(elevMask)
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

###get cell coordinates for G.sc
windowCoords <- getWindowCoords(G.sc)
habitatGrid <- windowCoords$habitatGrid

#----   2.6 HABITAT COVARIATES ---- 
#[CM] THIS MIGHT BE WHERE THE ISSUE WAS
# AS YOU DID FOR THE COORDINATES THE HABITAT, YOU NEED TO SELECT HAB COORDS THAT ARE CONSIDERED AS HABITAT
X.d <- values(elevMask)[!is.na(elevMask[])]

# Scale
X.d_mean <- mean(X.d)
X.d_sd <- sd(X.d)
X.d_sc <- (X.d - X.d_mean) / X.d_sd

#----   2.7 GET THE LOCAL DETECTORS OBJECTS  ---- 
# USE THE HABITAT GRID PROVIDED BY GETWINDOWCOORDS
habitatMask <- habitatGrid
habitatMask[habitatMask>0] <- 1 #TURN CELL ID TO 1 TO DEFINE THE HABITAT

# THERE MIGHT BE SOME TRIAL AND ERROR WHEN DEFINING THE DMAX.
localTraps <- getLocalTraps(habitatMask, X.sc, resizeFactor = 1, dmax = 5)

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

#---- 3. DETECTION DATA   ---- 
#----   3.1 MAKE Y   ---- 

K <- 7 # 7 occasions
n <- length(unique(edf$ind))

Y <- matrix(0, nrow = n, ncol = dim(X.sc)[1])
rownames(Y) <- unique(edf$ind)
xx <- edf[,c(2,4)]

for (obs in 1:nrow(xx)) {
  Y[xx[obs, 1], xx[obs, 2]] <- Y[xx[obs, 1], xx[obs, 2]] + 1
}

#----   3.2 AUGMENT Y   ---- 
##augment observed data to size M
M <- 400 
y.in <- rbind(Y, 
              matrix(0, nrow = M-n, ncol = J))

#----   3.3 USE SPARSE FORMAT FOR Y   ---- 
#change to 'sparse' format - speeds up computation by reducing file size
y.sparse <- getSparseY(y.in)

##extract pieces to be passed to Nimble
detNums <- y.sparse$detNums # Nº of traps at which each individual was detected
maxDetNums <- y.sparse$maxDetNums # Maximun Nº of traps at which any individual was detected
detIndices <- y.sparse$detIndices # ID of the traps where they were detected
y.sp <- y.sparse$y # Number of detections at each trap


#---- 4. FIT NIMBLE MODEL    ---- 
#----   4.1 CONSTANT AND DATA    ---- 

##compile constants
nimConstants <- list(
  M = M,
  J = J,
  numHabWindows = numHabWindows, 
  numGridRows = numGridRows,
  numGridCols = numGridCols, 
  maxDetNums = maxDetNums,
  MaxLocalTraps = MaxLocalTraps
)

##compile data
nimData <- list(habDens = X.d_sc,
                y = y.sp,
                detNums = detNums, 
                lowerHabCoords = lowerHabCoords,
                upperHabCoords = upperHabCoords,
                habitatGrid = habitatGrid,
                K = rep(K,J),
                X.sc = X.sc,
                habitatGridDet = habitatGridDet,
                detIndices = detIndices,
                detNums = detNums,
                localTrapsIndex = localTrapsIndex, 
                localTrapsNum = localTrapsNum 
)

#----   4.2 INITIAL VALUES  ---- 
#----     4.2.1 Z   ---- 

##set up initial values
z.in <- c(rep(1, n), rep(0, M-n))

#----     4.2.2 SXY  ---- 
##because of local evaluation of possible detectors, activity center initial 
##values have to be specified, eg average capture location
S.in <- matrix(NA, M, 2)
for ( i in 1:n){
  caps <- which(Y[i, ]>0)
  if (length(caps)==1){
    S.in[i,] <- as.matrix(X.sc[caps,])
  }else{
    S.in[i,] <- apply(X.sc[caps,], 2, mean)}
}


##random ACs for individuals never observed
##simulate them within state space!!
for(i in (n+1) : M){
  #[CM] THE ALTERNATIVE IS TO USE THE SIMULATION FUNCTIONNALITY OF NIMBLE TO SIMULATE INITIAL AC
  # IT IS THE SAME FUNCTION THAT YOU USE TO FIT THE MODEL, BUT IT STARTS WITH "r"
  S.in[i,] <- rbernppAC(n=1,
                        lowerCoords = nimData$lowerHabCoords,
                        upperCoords = nimData$upperHabCoords,
                        habitatGrid = habitatGrid,
                        numGridRows = nimConstants$numGridRows,
                        numGridCols = nimConstants$numGridCols,
                        logIntensities = rep(1, nimConstants$numHabWindows),#assume equal intensity across habitat cells (all 1)
                        logSumIntensity = sum(rep(1, nimConstants$numHabWindows))
  )
}
#----     4.2.2 COMPILE INITIAL VALUES   ---- 
inits <- function(){list(psi=runif(1,0.4, 0.6), 
                         sigma=runif(1,0.5, 1.5),
                         p0=runif(1,0,0.1),
                         z=z.in,
                         sxy=S.in,
                         beta.dens = runif(1,-0.1, 0.1))}

##source model code
setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Sim/2.SCRdenscov")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Run_Sim/2.SCRdenscov")
source('SCR in Nimble.R')

##determine which parameters to monitor
params <- c('N', 'psi', 'sigma', 'p0', 'beta.dens')

#---- 5. FIT THE MODEL ---- 
#[CM] INITIAL VALUES ARE IMPORTANT IN NIMBLE AND WITH NIMBLE, PROVIDE THEM
#(1) set up model
model <- nimbleModel(SCRhab,
                     constants = nimConstants, 
                     data=nimData,
                     inits = inits(),
                     check = FALSE)
##ignore error message, only due to missing initial values at this stage
# [CM] I WOULD NOT IGNORE THE ERROR MESSAGE. FROM OUR EXPERIENCE, A SCR MODEL IN NIMBLE THAT IS 
# NOT PROPERLY INITIALIZED MAY NOT RUN CORRECTLY. TO CHECK WHICH PARAMETERS ARE NOT PROPERLY INITIALIZE 
# YOU CAN DO:
model$initializeInfo()
#THEN YOU SHOULD ALWAYS CHECK THAT THE MODEL IS ABLE TO CALCULATE A LIKELIHOOD (RETURN A VALUE) GIVEN THE INITIAL VALUES PROVIDED
model$calculate()#
# IF A -INF OF NA IS RETURNED YOU CAN CHECK WHERE THE PROBLEM COMES FROM WITH (AND THEN TRY TO FIX IT UNTIL THE -INF DISAPEARS):
model$logProb_p0# WILL GIVE YOU LIKELIHOOD OF P0
model$logProb_sxy# THE PROBLEM IS OFTEN WITH SXY OR Y
model$logProb_y

#(2) Compile model in c++
cmodel <- compileNimble(model)       

# (3) Configure MCMC - on an uncompiled model - this step allows setting which quantities to monitor
conf.mcmc <- configureMCMC(model, monitors = params, thin=1)

# (4) Build the MCMC sampler based on configurations
mcmc <- buildMCMC(conf.mcmc)

# (5) Compile sampler in c++ together with compiled model
cmcmc <- compileNimble(mcmc, project = cmodel, resetFunctions = TRUE)

# (6) Run (monitor time just for fun) [takes 20 seconds on my computer]
system.time(
  (samp <- runMCMC(cmcmc, niter = 25000, nburnin = 1000, nchains=3, inits = inits) )
)


setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/Elevation")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/Elevation")
save(samp, file = "sampSCRdenscov2018.RData")


## -------------------------------------------------
##                      ROUGHNESS
## ------------------------------------------------- 

rm(list = ls())

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Effort/France/Data/Systematic")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Effort/France/Data/Systematic")

#---- 1. LOAD THE DETECTION DATA ---- 
load("edf2017_2019_fr.RData")
load("tdf2018_fr.RData")

edf <- edf[which(edf$session == 2), ] 

# As the model is set, there can be only one capture per trap and occasion
# --> The number of trials is 7 (From may to November): So max number of captures per trap = 7
# To fix it in this edf, remove duplicates
edf <- edf[-which(duplicated(edf)), ]

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

# Set up a raster file 
setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
#setwd("~/Data_server/Variables_hrscale")
roughness <- raster("rough_hrbear.tif")

# Crop it to extent of state-space
habitat.r <- crop(roughness, e) 

#----   2.4 DEFINE THE BUFFER AREA AND CUT WHAT IS NOT HABITAT---- 
# Buffer around traps (5*sigma = 33200)
Xpoints <- X
coordinates(Xpoints) <- Xpoints[,c(1,2)]
Xbuf <- gBuffer(Xpoints, width = 25000)

roughMask <- rasterize(Xbuf, habitat.r, mask = TRUE)

#RETAIN HABITAT COORDINATES THAT ARE WITHIN THE HABITAT
G <- coordinates(roughMask)[!is.na(roughMask[]),]
colnames(G) <- c("x","y")
#PLOT CHECK 
plot(roughMask)
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

###get cell coordinates for G.sc
windowCoords <- getWindowCoords(G.sc)
habitatGrid <- windowCoords$habitatGrid

#----   2.6 HABITAT COVARIATES ---- 
#[CM] THIS MIGHT BE WHERE THE ISSUE WAS
# AS YOU DID FOR THE COORDINATES THE HABITAT, YOU NEED TO SELECT HAB COORDS THAT ARE CONSIDERED AS HABITAT
X.d <- values(roughMask)[!is.na(roughMask[])]

# Scale
X.d_mean <- mean(X.d)
X.d_sd <- sd(X.d)
X.d_sc <- (X.d - X.d_mean) / X.d_sd

#----   2.7 GET THE LOCAL DETECTORS OBJECTS  ---- 
# USE THE HABITAT GRID PROVIDED BY GETWINDOWCOORDS
habitatMask <- habitatGrid
habitatMask[habitatMask>0] <- 1 #TURN CELL ID TO 1 TO DEFINE THE HABITAT

# THERE MIGHT BE SOME TRIAL AND ERROR WHEN DEFINING THE DMAX.
localTraps <- getLocalTraps(habitatMask, X.sc, resizeFactor = 1, dmax = 5)

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

#---- 3. DETECTION DATA   ---- 
#----   3.1 MAKE Y   ---- 

K <- 7 # 7 occasions
n <- length(unique(edf$ind))

Y <- matrix(0, nrow = n, ncol = dim(X.sc)[1])
rownames(Y) <- unique(edf$ind)
xx <- edf[,c(2,4)]

for (obs in 1:nrow(xx)) {
  Y[xx[obs, 1], xx[obs, 2]] <- Y[xx[obs, 1], xx[obs, 2]] + 1
}

#----   3.2 AUGMENT Y   ---- 
##augment observed data to size M
M <- 400 
y.in <- rbind(Y, 
              matrix(0, nrow = M-n, ncol = J))

#----   3.3 USE SPARSE FORMAT FOR Y   ---- 
#change to 'sparse' format - speeds up computation by reducing file size
y.sparse <- getSparseY(y.in)

##extract pieces to be passed to Nimble
detNums <- y.sparse$detNums # Nº of traps at which each individual was detected
maxDetNums <- y.sparse$maxDetNums # Maximun Nº of traps at which any individual was detected
detIndices <- y.sparse$detIndices # ID of the traps where they were detected
y.sp <- y.sparse$y # Number of detections at each trap


#---- 4. FIT NIMBLE MODEL    ---- 
#----   4.1 CONSTANT AND DATA    ---- 

##compile constants
nimConstants <- list(
  M = M,
  J = J,
  numHabWindows = numHabWindows, 
  numGridRows = numGridRows,
  numGridCols = numGridCols, 
  maxDetNums = maxDetNums,
  MaxLocalTraps = MaxLocalTraps
)

##compile data
nimData <- list(habDens = X.d_sc,
                y = y.sp,
                detNums = detNums, 
                lowerHabCoords = lowerHabCoords,
                upperHabCoords = upperHabCoords,
                habitatGrid = habitatGrid,
                K = rep(K,J),
                X.sc = X.sc,
                habitatGridDet = habitatGridDet,
                detIndices = detIndices,
                detNums = detNums,
                localTrapsIndex = localTrapsIndex, 
                localTrapsNum = localTrapsNum 
)

#----   4.2 INITIAL VALUES  ---- 
#----     4.2.1 Z   ---- 

##set up initial values
z.in <- c(rep(1, n), rep(0, M-n))

#----     4.2.2 SXY  ---- 
##because of local evaluation of possible detectors, activity center initial 
##values have to be specified, eg average capture location
S.in <- matrix(NA, M, 2)
for ( i in 1:n){
  caps <- which(Y[i, ]>0)
  if (length(caps)==1){
    S.in[i,] <- as.matrix(X.sc[caps,])
  }else{
    S.in[i,] <- apply(X.sc[caps,], 2, mean)}
}


##random ACs for individuals never observed
##simulate them within state space!!
for(i in (n+1) : M){
  #[CM] THE ALTERNATIVE IS TO USE THE SIMULATION FUNCTIONNALITY OF NIMBLE TO SIMULATE INITIAL AC
  # IT IS THE SAME FUNCTION THAT YOU USE TO FIT THE MODEL, BUT IT STARTS WITH "r"
  S.in[i,] <- rbernppAC(n=1,
                        lowerCoords = nimData$lowerHabCoords,
                        upperCoords = nimData$upperHabCoords,
                        habitatGrid = habitatGrid,
                        numGridRows = nimConstants$numGridRows,
                        numGridCols = nimConstants$numGridCols,
                        logIntensities = rep(1, nimConstants$numHabWindows),#assume equal intensity across habitat cells (all 1)
                        logSumIntensity = sum(rep(1, nimConstants$numHabWindows))
  )
}
#----     4.2.2 COMPILE INITIAL VALUES   ---- 
inits <- function(){list(psi=runif(1,0.4, 0.6), 
                         sigma=runif(1,0.5, 1.5),
                         p0=runif(1,0,0.1),
                         z=z.in,
                         sxy=S.in,
                         beta.dens = runif(1,-0.1, 0.1))}

##source model code
setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Sim/2.SCRdenscov")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Run_Sim/2.SCRdenscov")
source('SCR in Nimble.R')

##determine which parameters to monitor
params <- c('N', 'psi', 'sigma', 'p0', 'beta.dens')

#---- 5. FIT THE MODEL ---- 
#[CM] INITIAL VALUES ARE IMPORTANT IN NIMBLE AND WITH NIMBLE, PROVIDE THEM
#(1) set up model
model <- nimbleModel(SCRhab,
                     constants = nimConstants, 
                     data=nimData,
                     inits = inits(),
                     check = FALSE)
##ignore error message, only due to missing initial values at this stage
# [CM] I WOULD NOT IGNORE THE ERROR MESSAGE. FROM OUR EXPERIENCE, A SCR MODEL IN NIMBLE THAT IS 
# NOT PROPERLY INITIALIZED MAY NOT RUN CORRECTLY. TO CHECK WHICH PARAMETERS ARE NOT PROPERLY INITIALIZE 
# YOU CAN DO:
model$initializeInfo()
#THEN YOU SHOULD ALWAYS CHECK THAT THE MODEL IS ABLE TO CALCULATE A LIKELIHOOD (RETURN A VALUE) GIVEN THE INITIAL VALUES PROVIDED
model$calculate()#
# IF A -INF OF NA IS RETURNED YOU CAN CHECK WHERE THE PROBLEM COMES FROM WITH (AND THEN TRY TO FIX IT UNTIL THE -INF DISAPEARS):
model$logProb_p0# WILL GIVE YOU LIKELIHOOD OF P0
model$logProb_sxy# THE PROBLEM IS OFTEN WITH SXY OR Y
model$logProb_y

#(2) Compile model in c++
cmodel <- compileNimble(model)       

# (3) Configure MCMC - on an uncompiled model - this step allows setting which quantities to monitor
conf.mcmc <- configureMCMC(model, monitors = params, thin=1)

# (4) Build the MCMC sampler based on configurations
mcmc <- buildMCMC(conf.mcmc)

# (5) Compile sampler in c++ together with compiled model
cmcmc <- compileNimble(mcmc, project = cmodel, resetFunctions = TRUE)

# (6) Run (monitor time just for fun) [takes 20 seconds on my computer]
system.time(
  (samp <- runMCMC(cmcmc, niter = 25000, nburnin = 1000, nchains=3, inits = inits) )
)


setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/Roughness")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/Roughness")
save(samp, file = "sampSCRdenscov2018.RData")


## -------------------------------------------------
##                      SLOPE
## ------------------------------------------------- 

rm(list = ls())

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Effort/France/Data/Systematic")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Effort/France/Data/Systematic")

#---- 1. LOAD THE DETECTION DATA ---- 
load("edf2017_2019_fr.RData")
load("tdf2018_fr.RData")

edf <- edf[which(edf$session == 2), ] 

# As the model is set, there can be only one capture per trap and occasion
# --> The number of trials is 7 (From may to November): So max number of captures per trap = 7
# To fix it in this edf, remove duplicates
edf <- edf[-which(duplicated(edf)), ]

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

# Set up a raster file 
setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
#setwd("~/Data_server/Variables_hrscale")
slope <- raster("slope_hrbear.tif")

# Crop it to extent of state-space
habitat.r <- crop(slope, e) 

#----   2.4 DEFINE THE BUFFER AREA AND CUT WHAT IS NOT HABITAT---- 
# Buffer around traps (5*sigma = 33200)
Xpoints <- X
coordinates(Xpoints) <- Xpoints[,c(1,2)]
Xbuf <- gBuffer(Xpoints, width = 25000)

slopeMask <- rasterize(Xbuf, habitat.r, mask = TRUE)

#RETAIN HABITAT COORDINATES THAT ARE WITHIN THE HABITAT
G <- coordinates(slopeMask)[!is.na(slopeMask[]),]
colnames(G) <- c("x","y")
#PLOT CHECK 
plot(slopeMask)
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

###get cell coordinates for G.sc
windowCoords <- getWindowCoords(G.sc)
habitatGrid <- windowCoords$habitatGrid

#----   2.6 HABITAT COVARIATES ---- 
#[CM] THIS MIGHT BE WHERE THE ISSUE WAS
# AS YOU DID FOR THE COORDINATES THE HABITAT, YOU NEED TO SELECT HAB COORDS THAT ARE CONSIDERED AS HABITAT
X.d <- values(slopeMask)[!is.na(slopeMask[])]

# Scale
X.d_mean <- mean(X.d)
X.d_sd <- sd(X.d)
X.d_sc <- (X.d - X.d_mean) / X.d_sd

#----   2.7 GET THE LOCAL DETECTORS OBJECTS  ---- 
# USE THE HABITAT GRID PROVIDED BY GETWINDOWCOORDS
habitatMask <- habitatGrid
habitatMask[habitatMask>0] <- 1 #TURN CELL ID TO 1 TO DEFINE THE HABITAT

# THERE MIGHT BE SOME TRIAL AND ERROR WHEN DEFINING THE DMAX.
localTraps <- getLocalTraps(habitatMask, X.sc, resizeFactor = 1, dmax = 5)

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

#---- 3. DETECTION DATA   ---- 
#----   3.1 MAKE Y   ---- 

K <- 7 # 7 occasions
n <- length(unique(edf$ind))

Y <- matrix(0, nrow = n, ncol = dim(X.sc)[1])
rownames(Y) <- unique(edf$ind)
xx <- edf[,c(2,4)]

for (obs in 1:nrow(xx)) {
  Y[xx[obs, 1], xx[obs, 2]] <- Y[xx[obs, 1], xx[obs, 2]] + 1
}

#----   3.2 AUGMENT Y   ---- 
##augment observed data to size M
M <- 400 
y.in <- rbind(Y, 
              matrix(0, nrow = M-n, ncol = J))

#----   3.3 USE SPARSE FORMAT FOR Y   ---- 
#change to 'sparse' format - speeds up computation by reducing file size
y.sparse <- getSparseY(y.in)

##extract pieces to be passed to Nimble
detNums <- y.sparse$detNums # Nº of traps at which each individual was detected
maxDetNums <- y.sparse$maxDetNums # Maximun Nº of traps at which any individual was detected
detIndices <- y.sparse$detIndices # ID of the traps where they were detected
y.sp <- y.sparse$y # Number of detections at each trap


#---- 4. FIT NIMBLE MODEL    ---- 
#----   4.1 CONSTANT AND DATA    ---- 

##compile constants
nimConstants <- list(
  M = M,
  J = J,
  numHabWindows = numHabWindows, 
  numGridRows = numGridRows,
  numGridCols = numGridCols, 
  maxDetNums = maxDetNums,
  MaxLocalTraps = MaxLocalTraps
)

##compile data
nimData <- list(habDens = X.d_sc,
                y = y.sp,
                detNums = detNums, 
                lowerHabCoords = lowerHabCoords,
                upperHabCoords = upperHabCoords,
                habitatGrid = habitatGrid,
                K = rep(K,J),
                X.sc = X.sc,
                habitatGridDet = habitatGridDet,
                detIndices = detIndices,
                detNums = detNums,
                localTrapsIndex = localTrapsIndex, 
                localTrapsNum = localTrapsNum 
)

#----   4.2 INITIAL VALUES  ---- 
#----     4.2.1 Z   ---- 

##set up initial values
z.in <- c(rep(1, n), rep(0, M-n))

#----     4.2.2 SXY  ---- 
##because of local evaluation of possible detectors, activity center initial 
##values have to be specified, eg average capture location
S.in <- matrix(NA, M, 2)
for ( i in 1:n){
  caps <- which(Y[i, ]>0)
  if (length(caps)==1){
    S.in[i,] <- as.matrix(X.sc[caps,])
  }else{
    S.in[i,] <- apply(X.sc[caps,], 2, mean)}
}


##random ACs for individuals never observed
##simulate them within state space!!
for(i in (n+1) : M){
  #[CM] THE ALTERNATIVE IS TO USE THE SIMULATION FUNCTIONNALITY OF NIMBLE TO SIMULATE INITIAL AC
  # IT IS THE SAME FUNCTION THAT YOU USE TO FIT THE MODEL, BUT IT STARTS WITH "r"
  S.in[i,] <- rbernppAC(n=1,
                        lowerCoords = nimData$lowerHabCoords,
                        upperCoords = nimData$upperHabCoords,
                        habitatGrid = habitatGrid,
                        numGridRows = nimConstants$numGridRows,
                        numGridCols = nimConstants$numGridCols,
                        logIntensities = rep(1, nimConstants$numHabWindows),#assume equal intensity across habitat cells (all 1)
                        logSumIntensity = sum(rep(1, nimConstants$numHabWindows))
  )
}
#----     4.2.2 COMPILE INITIAL VALUES   ---- 
inits <- function(){list(psi=runif(1,0.4, 0.6), 
                         sigma=runif(1,0.5, 1.5),
                         p0=runif(1,0,0.1),
                         z=z.in,
                         sxy=S.in,
                         beta.dens = runif(1,-0.1, 0.1))}

##source model code
setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Sim/2.SCRdenscov")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Run_Sim/2.SCRdenscov")
source('SCR in Nimble.R')

##determine which parameters to monitor
params <- c('N', 'psi', 'sigma', 'p0', 'beta.dens')

#---- 5. FIT THE MODEL ---- 
#[CM] INITIAL VALUES ARE IMPORTANT IN NIMBLE AND WITH NIMBLE, PROVIDE THEM
#(1) set up model
model <- nimbleModel(SCRhab,
                     constants = nimConstants, 
                     data=nimData,
                     inits = inits(),
                     check = FALSE)
##ignore error message, only due to missing initial values at this stage
# [CM] I WOULD NOT IGNORE THE ERROR MESSAGE. FROM OUR EXPERIENCE, A SCR MODEL IN NIMBLE THAT IS 
# NOT PROPERLY INITIALIZED MAY NOT RUN CORRECTLY. TO CHECK WHICH PARAMETERS ARE NOT PROPERLY INITIALIZE 
# YOU CAN DO:
model$initializeInfo()
#THEN YOU SHOULD ALWAYS CHECK THAT THE MODEL IS ABLE TO CALCULATE A LIKELIHOOD (RETURN A VALUE) GIVEN THE INITIAL VALUES PROVIDED
model$calculate()#
# IF A -INF OF NA IS RETURNED YOU CAN CHECK WHERE THE PROBLEM COMES FROM WITH (AND THEN TRY TO FIX IT UNTIL THE -INF DISAPEARS):
model$logProb_p0# WILL GIVE YOU LIKELIHOOD OF P0
model$logProb_sxy# THE PROBLEM IS OFTEN WITH SXY OR Y
model$logProb_y

#(2) Compile model in c++
cmodel <- compileNimble(model)       

# (3) Configure MCMC - on an uncompiled model - this step allows setting which quantities to monitor
conf.mcmc <- configureMCMC(model, monitors = params, thin=1)

# (4) Build the MCMC sampler based on configurations
mcmc <- buildMCMC(conf.mcmc)

# (5) Compile sampler in c++ together with compiled model
cmcmc <- compileNimble(mcmc, project = cmodel, resetFunctions = TRUE)

# (6) Run (monitor time just for fun) [takes 20 seconds on my computer]
system.time(
  (samp <- runMCMC(cmcmc, niter = 25000, nburnin = 1000, nchains=3, inits = inits) )
)


setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/Slope")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/Slope")
save(samp, file = "sampSCRdenscov2018.RData")

## -------------------------------------------------
##              DISTANCE TO CORE (LOG)
## ------------------------------------------------- 

rm(list = ls())

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Effort/France/Data/Systematic")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Effort/France/Data/Systematic")

#---- 1. LOAD THE DETECTION DATA ---- 
load("edf2017_2019_fr.RData")
load("tdf2018_fr.RData")

edf <- edf[which(edf$session == 2), ] # Keep only detections of 2017 (year 1)

# As the model is set, there can be only one capture per trap and occasion
# --> The number of trials is 7 (From may to November): So max number of captures per trap = 7
# To fix it in this edf, remove duplicates
edf <- edf[-which(duplicated(edf)), ]

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

# Set up a raster file 
setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
#setwd("~/Data_server/Variables_hrscale")
distcore <- raster("logDistcore_hrbear.tif")

# Crop it to extent of state-space
habitat.r <- crop(distcore, e) 

#----   2.4 DEFINE THE BUFFER AREA AND CUT WHAT IS NOT HABITAT---- 
# Buffer around traps (5*sigma = 33200)
Xpoints <- X
coordinates(Xpoints) <- Xpoints[,c(1,2)]
Xbuf <- rgeos::gBuffer(Xpoints, width = 25000)

distcoreMask <- rasterize(Xbuf, habitat.r, mask = TRUE)

#RETAIN HABITAT COORDINATES THAT ARE WITHIN THE HABITAT
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

###get cell coordinates for G.sc
windowCoords <- getWindowCoords(G.sc)
habitatGrid <- windowCoords$habitatGrid

#----   2.6 HABITAT COVARIATES ---- 
#[CM] THIS MIGHT BE WHERE THE ISSUE WAS
# AS YOU DID FOR THE COORDINATES THE HABITAT, YOU NEED TO SELECT HAB COORDS THAT ARE CONSIDERED AS HABITAT
X.d <- raster::values(distcoreMask)[!is.na(distcoreMask[])] ## !!! If I load nimbleSCR, the function values gives error

# Scale
X.d_mean <- mean(X.d)
X.d_sd <- sd(X.d)
X.d_sc <- (X.d - X.d_mean) / X.d_sd

#----   2.7 GET THE LOCAL DETECTORS OBJECTS  ---- 
# USE THE HABITAT GRID PROVIDED BY GETWINDOWCOORDS
habitatMask <- habitatGrid
habitatMask[habitatMask>0] <- 1 #TURN CELL ID TO 1 TO DEFINE THE HABITAT

# THERE MIGHT BE SOME TRIAL AND ERROR WHEN DEFINING THE DMAX.
localTraps <- getLocalTraps(habitatMask, X.sc, resizeFactor = 1, dmax = 5)

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

#---- 3. DETECTION DATA   ---- 
#----   3.1 MAKE Y   ---- 

K <- 7 # 7 occasions
n <- length(unique(edf$ind))

Y <- matrix(0, nrow = n, ncol = dim(X.sc)[1])
rownames(Y) <- unique(edf$ind)
xx <- edf[,c(2,4)]

for (obs in 1:nrow(xx)) {
  Y[xx[obs, 1], xx[obs, 2]] <- Y[xx[obs, 1], xx[obs, 2]] + 1
}

#----   3.2 AUGMENT Y   ---- 
##augment observed data to size M
M <- 400 
y.in <- rbind(Y, 
              matrix(0, nrow = M-n, ncol = J))

#----   3.3 USE SPARSE FORMAT FOR Y   ---- 
#change to 'sparse' format - speeds up computation by reducing file size
y.sparse <- getSparseY(y.in)

##extract pieces to be passed to Nimble
detNums <- y.sparse$detNums # Nº of traps at which each individual was detected
maxDetNums <- y.sparse$maxDetNums # Maximun Nº of traps at which any individual was detected
detIndices <- y.sparse$detIndices # ID of the traps where they were detected
y.sp <- y.sparse$y # Number of detections at each trap


#---- 4. FIT NIMBLE MODEL    ---- 
#----   4.1 CONSTANT AND DATA    ---- 

##compile constants
nimConstants <- list(
  M = M,
  J = J,
  numHabWindows = numHabWindows, 
  numGridRows = numGridRows,
  numGridCols = numGridCols, 
  maxDetNums = maxDetNums,
  MaxLocalTraps = MaxLocalTraps
)

##compile data
nimData <- list(habDens = X.d_sc,
                y = y.sp,
                detNums = detNums, 
                lowerHabCoords = lowerHabCoords,
                upperHabCoords = upperHabCoords,
                habitatGrid = habitatGrid,
                K = rep(K,J),
                X.sc = X.sc,
                habitatGridDet = habitatGridDet,
                detIndices = detIndices,
                detNums = detNums,
                localTrapsIndex = localTrapsIndex, 
                localTrapsNum = localTrapsNum 
)

#----   4.2 INITIAL VALUES  ---- 
#----     4.2.1 Z   ---- 

##set up initial values
z.in <- c(rep(1, n), rep(0, M-n))

#----     4.2.2 SXY  ---- 
##because of local evaluation of possible detectors, activity center initial 
##values have to be specified, eg average capture location
S.in <- matrix(NA, M, 2)
for ( i in 1:n){
  caps <- which(Y[i, ]>0)
  if (length(caps)==1){
    S.in[i,] <- as.matrix(X.sc[caps,])
  }else{
    S.in[i,] <- apply(X.sc[caps,], 2, mean)}
}


##random ACs for individuals never observed
##simulate them within state space!!
for(i in (n+1) : M){
  #[CM] THE ALTERNATIVE IS TO USE THE SIMULATION FUNCTIONNALITY OF NIMBLE TO SIMULATE INITIAL AC
  # IT IS THE SAME FUNCTION THAT YOU USE TO FIT THE MODEL, BUT IT STARTS WITH "r"
  S.in[i,] <- rbernppAC(n=1,
                        lowerCoords = nimData$lowerHabCoords,
                        upperCoords = nimData$upperHabCoords,
                        habitatGrid = habitatGrid,
                        numGridRows = nimConstants$numGridRows,
                        numGridCols = nimConstants$numGridCols,
                        logIntensities = rep(1, nimConstants$numHabWindows),#assume equal intensity across habitat cells (all 1)
                        logSumIntensity = sum(rep(1, nimConstants$numHabWindows))
  )
}
#----     4.2.2 COMPILE INITIAL VALUES   ---- 
inits <- function(){list(psi=runif(1,0.4, 0.6), 
                         sigma=runif(1,0.5, 1.5),
                         p0=runif(1,0,0.1),
                         z=z.in,
                         sxy=S.in,
                         beta.dens = runif(1,-0.1, 0.1))}

##source model code
setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Sim/2.SCRdenscov")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Run_Sim/2.SCRdenscov")
source('SCR in Nimble.R')

##determine which parameters to monitor
params <- c('N', 'psi', 'sigma', 'p0', 'beta.dens')

#---- 5. FIT THE MODEL ---- 
#[CM] INITIAL VALUES ARE IMPORTANT IN NIMBLE AND WITH NIMBLE, PROVIDE THEM
#(1) set up model
model <- nimbleModel(SCRhab,
                     constants = nimConstants, 
                     data=nimData,
                     inits = inits(),
                     check = FALSE)
##ignore error message, only due to missing initial values at this stage
# [CM] I WOULD NOT IGNORE THE ERROR MESSAGE. FROM OUR EXPERIENCE, A SCR MODEL IN NIMBLE THAT IS 
# NOT PROPERLY INITIALIZED MAY NOT RUN CORRECTLY. TO CHECK WHICH PARAMETERS ARE NOT PROPERLY INITIALIZE 
# YOU CAN DO:
model$initializeInfo()
#THEN YOU SHOULD ALWAYS CHECK THAT THE MODEL IS ABLE TO CALCULATE A LIKELIHOOD (RETURN A VALUE) GIVEN THE INITIAL VALUES PROVIDED
model$calculate()#
# IF A -INF OF NA IS RETURNED YOU CAN CHECK WHERE THE PROBLEM COMES FROM WITH (AND THEN TRY TO FIX IT UNTIL THE -INF DISAPEARS):
model$logProb_p0# WILL GIVE YOU LIKELIHOOD OF P0
model$logProb_sxy# THE PROBLEM IS OFTEN WITH SXY OR Y
model$logProb_y

#(2) Compile model in c++
cmodel <- compileNimble(model)       

# (3) Configure MCMC - on an uncompiled model - this step allows setting which quantities to monitor
conf.mcmc <- configureMCMC(model, monitors = params, thin=1)

# (4) Build the MCMC sampler based on configurations
mcmc <- buildMCMC(conf.mcmc)

# (5) Compile sampler in c++ together with compiled model
cmcmc <- compileNimble(mcmc, project = cmodel, resetFunctions = TRUE)

# (6) Run (monitor time just for fun) [takes 20 seconds on my computer]
system.time(
  (samp <- runMCMC(cmcmc, niter = 25000, nburnin = 1000, nchains=3, inits = inits) )
)


setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/Distcore")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/Slope")
save(samp, file = "sampSCRdenscov2018.RData")

## -------------------------------------------------
##        DENSITY OF OBSERVATIONS (from 200m res)
## ------------------------------------------------- 

rm(list = ls())

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Effort/France/Data/Systematic")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Effort/France/Data/Systematic")

#---- 1. LOAD THE DETECTION DATA ---- 
load("edf2017_2019_fr.RData")
load("tdf2018_fr.RData")

edf <- edf[which(edf$session == 2), ] # Keep only detections of 2017 (year 1)

# As the model is set, there can be only one capture per trap and occasion
# --> The number of trials is 7 (From may to November): So max number of captures per trap = 7
# To fix it in this edf, remove duplicates
edf <- edf[-which(duplicated(edf)), ]

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

# Set up a raster file 
setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
#setwd("~/Data_server/Variables_hrscale")
obsDens <- raster("obsDens200m_hrbear.tif")

# Crop it to extent of state-space
habitat.r <- crop(obsDens, e) 

#----   2.4 DEFINE THE BUFFER AREA AND CUT WHAT IS NOT HABITAT---- 
# Buffer around traps (5*sigma = 33200)
Xpoints <- X
coordinates(Xpoints) <- Xpoints[,c(1,2)]
Xbuf <- rgeos::gBuffer(Xpoints, width = 25000)

obsDensMask <- rasterize(Xbuf, habitat.r, mask = TRUE)

#RETAIN HABITAT COORDINATES THAT ARE WITHIN THE HABITAT
G <- coordinates(obsDensMask)[!is.na(obsDensMask[]),]
colnames(G) <- c("x","y")
#PLOT CHECK 
plot(obsDensMask)
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

###get cell coordinates for G.sc
windowCoords <- getWindowCoords(G.sc)
habitatGrid <- windowCoords$habitatGrid

#----   2.6 HABITAT COVARIATES ---- 
#[CM] THIS MIGHT BE WHERE THE ISSUE WAS
# AS YOU DID FOR THE COORDINATES THE HABITAT, YOU NEED TO SELECT HAB COORDS THAT ARE CONSIDERED AS HABITAT
X.d <- raster::values(obsDensMask)[!is.na(obsDensMask[])] ## !!! If I load nimbleSCR, the function values gives error

# Scale
X.d_mean <- mean(X.d)
X.d_sd <- sd(X.d)
X.d_sc <- (X.d - X.d_mean) / X.d_sd

#----   2.7 GET THE LOCAL DETECTORS OBJECTS  ---- 
# USE THE HABITAT GRID PROVIDED BY GETWINDOWCOORDS
habitatMask <- habitatGrid
habitatMask[habitatMask>0] <- 1 #TURN CELL ID TO 1 TO DEFINE THE HABITAT

# THERE MIGHT BE SOME TRIAL AND ERROR WHEN DEFINING THE DMAX.
localTraps <- getLocalTraps(habitatMask, X.sc, resizeFactor = 1, dmax = 5)

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

#---- 3. DETECTION DATA   ---- 
#----   3.1 MAKE Y   ---- 

K <- 7 # 7 occasions
n <- length(unique(edf$ind))

Y <- matrix(0, nrow = n, ncol = dim(X.sc)[1])
rownames(Y) <- unique(edf$ind)
xx <- edf[,c(2,4)]

for (obs in 1:nrow(xx)) {
  Y[xx[obs, 1], xx[obs, 2]] <- Y[xx[obs, 1], xx[obs, 2]] + 1
}

#----   3.2 AUGMENT Y   ---- 
##augment observed data to size M
M <- 400 
y.in <- rbind(Y, 
              matrix(0, nrow = M-n, ncol = J))

#----   3.3 USE SPARSE FORMAT FOR Y   ---- 
#change to 'sparse' format - speeds up computation by reducing file size
y.sparse <- getSparseY(y.in)

##extract pieces to be passed to Nimble
detNums <- y.sparse$detNums # Nº of traps at which each individual was detected
maxDetNums <- y.sparse$maxDetNums # Maximun Nº of traps at which any individual was detected
detIndices <- y.sparse$detIndices # ID of the traps where they were detected
y.sp <- y.sparse$y # Number of detections at each trap


#---- 4. FIT NIMBLE MODEL    ---- 
#----   4.1 CONSTANT AND DATA    ---- 

##compile constants
nimConstants <- list(
  M = M,
  J = J,
  numHabWindows = numHabWindows, 
  numGridRows = numGridRows,
  numGridCols = numGridCols, 
  maxDetNums = maxDetNums,
  MaxLocalTraps = MaxLocalTraps
)

##compile data
nimData <- list(habDens = X.d_sc,
                y = y.sp,
                detNums = detNums, 
                lowerHabCoords = lowerHabCoords,
                upperHabCoords = upperHabCoords,
                habitatGrid = habitatGrid,
                K = rep(K,J),
                X.sc = X.sc,
                habitatGridDet = habitatGridDet,
                detIndices = detIndices,
                detNums = detNums,
                localTrapsIndex = localTrapsIndex, 
                localTrapsNum = localTrapsNum 
)

#----   4.2 INITIAL VALUES  ---- 
#----     4.2.1 Z   ---- 

##set up initial values
z.in <- c(rep(1, n), rep(0, M-n))

#----     4.2.2 SXY  ---- 
##because of local evaluation of possible detectors, activity center initial 
##values have to be specified, eg average capture location
S.in <- matrix(NA, M, 2)
for ( i in 1:n){
  caps <- which(Y[i, ]>0)
  if (length(caps)==1){
    S.in[i,] <- as.matrix(X.sc[caps,])
  }else{
    S.in[i,] <- apply(X.sc[caps,], 2, mean)}
}


##random ACs for individuals never observed
##simulate them within state space!!
for(i in (n+1) : M){
  #[CM] THE ALTERNATIVE IS TO USE THE SIMULATION FUNCTIONNALITY OF NIMBLE TO SIMULATE INITIAL AC
  # IT IS THE SAME FUNCTION THAT YOU USE TO FIT THE MODEL, BUT IT STARTS WITH "r"
  S.in[i,] <- rbernppAC(n=1,
                        lowerCoords = nimData$lowerHabCoords,
                        upperCoords = nimData$upperHabCoords,
                        habitatGrid = habitatGrid,
                        numGridRows = nimConstants$numGridRows,
                        numGridCols = nimConstants$numGridCols,
                        logIntensities = rep(1, nimConstants$numHabWindows),#assume equal intensity across habitat cells (all 1)
                        logSumIntensity = sum(rep(1, nimConstants$numHabWindows))
  )
}
#----     4.2.2 COMPILE INITIAL VALUES   ---- 
inits <- function(){list(psi=runif(1,0.4, 0.6), 
                         sigma=runif(1,0.5, 1.5),
                         p0=runif(1,0,0.1),
                         z=z.in,
                         sxy=S.in,
                         beta.dens = runif(1,-0.1, 0.1))}

##source model code
setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Sim/2.SCRdenscov")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Run_Sim/2.SCRdenscov")
source('SCR in Nimble.R')

##determine which parameters to monitor
params <- c('N', 'psi', 'sigma', 'p0', 'beta.dens')

#---- 5. FIT THE MODEL ---- 
#[CM] INITIAL VALUES ARE IMPORTANT IN NIMBLE AND WITH NIMBLE, PROVIDE THEM
#(1) set up model
model <- nimbleModel(SCRhab,
                     constants = nimConstants, 
                     data=nimData,
                     inits = inits(),
                     check = FALSE)
##ignore error message, only due to missing initial values at this stage
# [CM] I WOULD NOT IGNORE THE ERROR MESSAGE. FROM OUR EXPERIENCE, A SCR MODEL IN NIMBLE THAT IS 
# NOT PROPERLY INITIALIZED MAY NOT RUN CORRECTLY. TO CHECK WHICH PARAMETERS ARE NOT PROPERLY INITIALIZE 
# YOU CAN DO:
model$initializeInfo()
#THEN YOU SHOULD ALWAYS CHECK THAT THE MODEL IS ABLE TO CALCULATE A LIKELIHOOD (RETURN A VALUE) GIVEN THE INITIAL VALUES PROVIDED
model$calculate()#
# IF A -INF OF NA IS RETURNED YOU CAN CHECK WHERE THE PROBLEM COMES FROM WITH (AND THEN TRY TO FIX IT UNTIL THE -INF DISAPEARS):
model$logProb_p0# WILL GIVE YOU LIKELIHOOD OF P0
model$logProb_sxy# THE PROBLEM IS OFTEN WITH SXY OR Y
model$logProb_y

#(2) Compile model in c++
cmodel <- compileNimble(model)       

# (3) Configure MCMC - on an uncompiled model - this step allows setting which quantities to monitor
conf.mcmc <- configureMCMC(model, monitors = params, thin=1)

# (4) Build the MCMC sampler based on configurations
mcmc <- buildMCMC(conf.mcmc)

# (5) Compile sampler in c++ together with compiled model
cmcmc <- compileNimble(mcmc, project = cmodel, resetFunctions = TRUE)

# (6) Run (monitor time just for fun) [takes 20 seconds on my computer]
system.time(
  (samp <- runMCMC(cmcmc, niter = 25000, nburnin = 1000, nchains=3, inits = inits) )
)


setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/obsDens")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/obsDens")
save(samp, file = "sampSCRdenscov2018.RData")

## -------------------------------------------------
## DENSITY OF OBSERVATIONS BEFORE STUDY PERIOD 
##          (2010-2016)(from 200m res)
## ------------------------------------------------- 

rm(list = ls())

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Effort/France/Data/Systematic")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Effort/France/Data/Systematic")

#---- 1. LOAD THE DETECTION DATA ---- 
load("edf2017_2019_fr.RData")
load("tdf2018_fr.RData")

edf <- edf[which(edf$session == 2), ] # Keep only detections of 2017 (year 1)

# As the model is set, there can be only one capture per trap and occasion
# --> The number of trials is 7 (From may to November): So max number of captures per trap = 7
# To fix it in this edf, remove duplicates
edf <- edf[-which(duplicated(edf)), ]

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

# Set up a raster file 
setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
#setwd("~/Data_server/Variables_hrscale")
obsDens <- raster("obsDens200m_preST_hrbear.tif")

# Crop it to extent of state-space
habitat.r <- crop(obsDens, e) 

#----   2.4 DEFINE THE BUFFER AREA AND CUT WHAT IS NOT HABITAT---- 
# Buffer around traps (5*sigma = 33200)
Xpoints <- X
coordinates(Xpoints) <- Xpoints[,c(1,2)]
Xbuf <- rgeos::gBuffer(Xpoints, width = 25000)

obsDensMask <- rasterize(Xbuf, habitat.r, mask = TRUE)

#RETAIN HABITAT COORDINATES THAT ARE WITHIN THE HABITAT
G <- coordinates(obsDensMask)[!is.na(obsDensMask[]),]
colnames(G) <- c("x","y")
#PLOT CHECK 
plot(obsDensMask)
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

###get cell coordinates for G.sc
windowCoords <- getWindowCoords(G.sc)
habitatGrid <- windowCoords$habitatGrid

#----   2.6 HABITAT COVARIATES ---- 
#[CM] THIS MIGHT BE WHERE THE ISSUE WAS
# AS YOU DID FOR THE COORDINATES THE HABITAT, YOU NEED TO SELECT HAB COORDS THAT ARE CONSIDERED AS HABITAT
X.d <- raster::values(obsDensMask)[!is.na(obsDensMask[])] ## !!! If I load nimbleSCR, the function values gives error

# Scale
X.d_mean <- mean(X.d)
X.d_sd <- sd(X.d)
X.d_sc <- (X.d - X.d_mean) / X.d_sd

#----   2.7 GET THE LOCAL DETECTORS OBJECTS  ---- 
# USE THE HABITAT GRID PROVIDED BY GETWINDOWCOORDS
habitatMask <- habitatGrid
habitatMask[habitatMask>0] <- 1 #TURN CELL ID TO 1 TO DEFINE THE HABITAT

# THERE MIGHT BE SOME TRIAL AND ERROR WHEN DEFINING THE DMAX.
localTraps <- getLocalTraps(habitatMask, X.sc, resizeFactor = 1, dmax = 5)

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

#---- 3. DETECTION DATA   ---- 
#----   3.1 MAKE Y   ---- 

K <- 7 # 7 occasions
n <- length(unique(edf$ind))

Y <- matrix(0, nrow = n, ncol = dim(X.sc)[1])
rownames(Y) <- unique(edf$ind)
xx <- edf[,c(2,4)]

for (obs in 1:nrow(xx)) {
  Y[xx[obs, 1], xx[obs, 2]] <- Y[xx[obs, 1], xx[obs, 2]] + 1
}

#----   3.2 AUGMENT Y   ---- 
##augment observed data to size M
M <- 400 
y.in <- rbind(Y, 
              matrix(0, nrow = M-n, ncol = J))

#----   3.3 USE SPARSE FORMAT FOR Y   ---- 
#change to 'sparse' format - speeds up computation by reducing file size
y.sparse <- getSparseY(y.in)

##extract pieces to be passed to Nimble
detNums <- y.sparse$detNums # Nº of traps at which each individual was detected
maxDetNums <- y.sparse$maxDetNums # Maximun Nº of traps at which any individual was detected
detIndices <- y.sparse$detIndices # ID of the traps where they were detected
y.sp <- y.sparse$y # Number of detections at each trap


#---- 4. FIT NIMBLE MODEL    ---- 
#----   4.1 CONSTANT AND DATA    ---- 

##compile constants
nimConstants <- list(
  M = M,
  J = J,
  numHabWindows = numHabWindows, 
  numGridRows = numGridRows,
  numGridCols = numGridCols, 
  maxDetNums = maxDetNums,
  MaxLocalTraps = MaxLocalTraps
)

##compile data
nimData <- list(habDens = X.d_sc,
                y = y.sp,
                detNums = detNums, 
                lowerHabCoords = lowerHabCoords,
                upperHabCoords = upperHabCoords,
                habitatGrid = habitatGrid,
                K = rep(K,J),
                X.sc = X.sc,
                habitatGridDet = habitatGridDet,
                detIndices = detIndices,
                detNums = detNums,
                localTrapsIndex = localTrapsIndex, 
                localTrapsNum = localTrapsNum 
)

#----   4.2 INITIAL VALUES  ---- 
#----     4.2.1 Z   ---- 

##set up initial values
z.in <- c(rep(1, n), rep(0, M-n))

#----     4.2.2 SXY  ---- 
##because of local evaluation of possible detectors, activity center initial 
##values have to be specified, eg average capture location
S.in <- matrix(NA, M, 2)
for ( i in 1:n){
  caps <- which(Y[i, ]>0)
  if (length(caps)==1){
    S.in[i,] <- as.matrix(X.sc[caps,])
  }else{
    S.in[i,] <- apply(X.sc[caps,], 2, mean)}
}


##random ACs for individuals never observed
##simulate them within state space!!
for(i in (n+1) : M){
  #[CM] THE ALTERNATIVE IS TO USE THE SIMULATION FUNCTIONNALITY OF NIMBLE TO SIMULATE INITIAL AC
  # IT IS THE SAME FUNCTION THAT YOU USE TO FIT THE MODEL, BUT IT STARTS WITH "r"
  S.in[i,] <- rbernppAC(n=1,
                        lowerCoords = nimData$lowerHabCoords,
                        upperCoords = nimData$upperHabCoords,
                        habitatGrid = habitatGrid,
                        numGridRows = nimConstants$numGridRows,
                        numGridCols = nimConstants$numGridCols,
                        logIntensities = rep(1, nimConstants$numHabWindows),#assume equal intensity across habitat cells (all 1)
                        logSumIntensity = sum(rep(1, nimConstants$numHabWindows))
  )
}
#----     4.2.2 COMPILE INITIAL VALUES   ---- 
inits <- function(){list(psi=runif(1,0.4, 0.6), 
                         sigma=runif(1,0.5, 1.5),
                         p0=runif(1,0,0.1),
                         z=z.in,
                         sxy=S.in,
                         beta.dens = runif(1,-0.1, 0.1))}

##source model code
setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Sim/2.SCRdenscov")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Run_Sim/2.SCRdenscov")
source('SCR in Nimble.R')

##determine which parameters to monitor
params <- c('N', 'psi', 'sigma', 'p0', 'beta.dens')

#---- 5. FIT THE MODEL ---- 
#[CM] INITIAL VALUES ARE IMPORTANT IN NIMBLE AND WITH NIMBLE, PROVIDE THEM
#(1) set up model
model <- nimbleModel(SCRhab,
                     constants = nimConstants, 
                     data=nimData,
                     inits = inits(),
                     check = FALSE)
##ignore error message, only due to missing initial values at this stage
# [CM] I WOULD NOT IGNORE THE ERROR MESSAGE. FROM OUR EXPERIENCE, A SCR MODEL IN NIMBLE THAT IS 
# NOT PROPERLY INITIALIZED MAY NOT RUN CORRECTLY. TO CHECK WHICH PARAMETERS ARE NOT PROPERLY INITIALIZE 
# YOU CAN DO:
model$initializeInfo()
#THEN YOU SHOULD ALWAYS CHECK THAT THE MODEL IS ABLE TO CALCULATE A LIKELIHOOD (RETURN A VALUE) GIVEN THE INITIAL VALUES PROVIDED
model$calculate()#
# IF A -INF OF NA IS RETURNED YOU CAN CHECK WHERE THE PROBLEM COMES FROM WITH (AND THEN TRY TO FIX IT UNTIL THE -INF DISAPEARS):
model$logProb_p0# WILL GIVE YOU LIKELIHOOD OF P0
model$logProb_sxy# THE PROBLEM IS OFTEN WITH SXY OR Y
model$logProb_y

#(2) Compile model in c++
cmodel <- compileNimble(model)       

# (3) Configure MCMC - on an uncompiled model - this step allows setting which quantities to monitor
conf.mcmc <- configureMCMC(model, monitors = params, thin=1)

# (4) Build the MCMC sampler based on configurations
mcmc <- buildMCMC(conf.mcmc)

# (5) Compile sampler in c++ together with compiled model
cmcmc <- compileNimble(mcmc, project = cmodel, resetFunctions = TRUE)

# (6) Run (monitor time just for fun) [takes 20 seconds on my computer]
system.time(
  (samp <- runMCMC(cmcmc, niter = 25000, nburnin = 1000, nchains=3, inits = inits) )
)


setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/obsDensPre")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/obsDensPre")
save(samp, file = "sampSCRdenscov2018.RData")

## -------------------------------------------------
##              ROAD DENSITY CATEGORY 1 
## ------------------------------------------------- 

rm(list = ls())

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Effort/France/Data/Systematic")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Effort/France/Data/Systematic")

#---- 1. LOAD THE DETECTION DATA ---- 
load("edf2017_2019_fr.RData")
load("tdf2018_fr.RData")

edf <- edf[which(edf$session == 2), ] # Keep only detections of 2017 (year 1)

# As the model is set, there can be only one capture per trap and occasion
# --> The number of trials is 7 (From may to November): So max number of captures per trap = 7
# To fix it in this edf, remove duplicates
edf <- edf[-which(duplicated(edf)), ]

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

# Set up a raster file 
setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
#setwd("~/Data_server/Variables_hrscale")
roads1 <- raster("roads1_hrbear.tif")

# Crop it to extent of state-space
habitat.r <- crop(roads1, e) 

#----   2.4 DEFINE THE BUFFER AREA AND CUT WHAT IS NOT HABITAT---- 
# Buffer around traps (5*sigma = 33200)
Xpoints <- X
coordinates(Xpoints) <- Xpoints[,c(1,2)]
Xbuf <- rgeos::gBuffer(Xpoints, width = 25000)

roads1Mask <- rasterize(Xbuf, habitat.r, mask = TRUE)

#RETAIN HABITAT COORDINATES THAT ARE WITHIN THE HABITAT
G <- coordinates(roads1Mask)[!is.na(roads1Mask[]),]
colnames(G) <- c("x","y")
#PLOT CHECK 
plot(roads1Mask)
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

###get cell coordinates for G.sc
windowCoords <- getWindowCoords(G.sc)
habitatGrid <- windowCoords$habitatGrid

#----   2.6 HABITAT COVARIATES ---- 
#[CM] THIS MIGHT BE WHERE THE ISSUE WAS
# AS YOU DID FOR THE COORDINATES THE HABITAT, YOU NEED TO SELECT HAB COORDS THAT ARE CONSIDERED AS HABITAT
X.d <- raster::values(roads1Mask)[!is.na(roads1Mask[])] ## !!! If I load nimbleSCR, the function values gives error

# Scale
X.d_mean <- mean(X.d)
X.d_sd <- sd(X.d)
X.d_sc <- (X.d - X.d_mean) / X.d_sd

#----   2.7 GET THE LOCAL DETECTORS OBJECTS  ---- 
# USE THE HABITAT GRID PROVIDED BY GETWINDOWCOORDS
habitatMask <- habitatGrid
habitatMask[habitatMask>0] <- 1 #TURN CELL ID TO 1 TO DEFINE THE HABITAT

# THERE MIGHT BE SOME TRIAL AND ERROR WHEN DEFINING THE DMAX.
localTraps <- getLocalTraps(habitatMask, X.sc, resizeFactor = 1, dmax = 5)

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

#---- 3. DETECTION DATA   ---- 
#----   3.1 MAKE Y   ---- 

K <- 7 # 7 occasions
n <- length(unique(edf$ind))

Y <- matrix(0, nrow = n, ncol = dim(X.sc)[1])
rownames(Y) <- unique(edf$ind)
xx <- edf[,c(2,4)]

for (obs in 1:nrow(xx)) {
  Y[xx[obs, 1], xx[obs, 2]] <- Y[xx[obs, 1], xx[obs, 2]] + 1
}

#----   3.2 AUGMENT Y   ---- 
##augment observed data to size M
M <- 400 
y.in <- rbind(Y, 
              matrix(0, nrow = M-n, ncol = J))

#----   3.3 USE SPARSE FORMAT FOR Y   ---- 
#change to 'sparse' format - speeds up computation by reducing file size
y.sparse <- getSparseY(y.in)

##extract pieces to be passed to Nimble
detNums <- y.sparse$detNums # Nº of traps at which each individual was detected
maxDetNums <- y.sparse$maxDetNums # Maximun Nº of traps at which any individual was detected
detIndices <- y.sparse$detIndices # ID of the traps where they were detected
y.sp <- y.sparse$y # Number of detections at each trap


#---- 4. FIT NIMBLE MODEL    ---- 
#----   4.1 CONSTANT AND DATA    ---- 

##compile constants
nimConstants <- list(
  M = M,
  J = J,
  numHabWindows = numHabWindows, 
  numGridRows = numGridRows,
  numGridCols = numGridCols, 
  maxDetNums = maxDetNums,
  MaxLocalTraps = MaxLocalTraps
)

##compile data
nimData <- list(habDens = X.d_sc,
                y = y.sp,
                detNums = detNums, 
                lowerHabCoords = lowerHabCoords,
                upperHabCoords = upperHabCoords,
                habitatGrid = habitatGrid,
                K = rep(K,J),
                X.sc = X.sc,
                habitatGridDet = habitatGridDet,
                detIndices = detIndices,
                detNums = detNums,
                localTrapsIndex = localTrapsIndex, 
                localTrapsNum = localTrapsNum 
)

#----   4.2 INITIAL VALUES  ---- 
#----     4.2.1 Z   ---- 

##set up initial values
z.in <- c(rep(1, n), rep(0, M-n))

#----     4.2.2 SXY  ---- 
##because of local evaluation of possible detectors, activity center initial 
##values have to be specified, eg average capture location
S.in <- matrix(NA, M, 2)
for ( i in 1:n){
  caps <- which(Y[i, ]>0)
  if (length(caps)==1){
    S.in[i,] <- as.matrix(X.sc[caps,])
  }else{
    S.in[i,] <- apply(X.sc[caps,], 2, mean)}
}


##random ACs for individuals never observed
##simulate them within state space!!
for(i in (n+1) : M){
  #[CM] THE ALTERNATIVE IS TO USE THE SIMULATION FUNCTIONNALITY OF NIMBLE TO SIMULATE INITIAL AC
  # IT IS THE SAME FUNCTION THAT YOU USE TO FIT THE MODEL, BUT IT STARTS WITH "r"
  S.in[i,] <- rbernppAC(n=1,
                        lowerCoords = nimData$lowerHabCoords,
                        upperCoords = nimData$upperHabCoords,
                        habitatGrid = habitatGrid,
                        numGridRows = nimConstants$numGridRows,
                        numGridCols = nimConstants$numGridCols,
                        logIntensities = rep(1, nimConstants$numHabWindows),#assume equal intensity across habitat cells (all 1)
                        logSumIntensity = sum(rep(1, nimConstants$numHabWindows))
  )
}
#----     4.2.2 COMPILE INITIAL VALUES   ---- 
inits <- function(){list(psi=runif(1,0.4, 0.6), 
                         sigma=runif(1,0.5, 1.5),
                         p0=runif(1,0,0.1),
                         z=z.in,
                         sxy=S.in,
                         beta.dens = runif(1,-0.1, 0.1))}

##source model code
setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Model")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Run_Sim/2.SCRdenscov")
source('2.SCRdenscov in Nimble.R')

##determine which parameters to monitor
params <- c('N', 'psi', 'sigma', 'p0', 'beta.dens')

#---- 5. FIT THE MODEL ---- 
#[CM] INITIAL VALUES ARE IMPORTANT IN NIMBLE AND WITH NIMBLE, PROVIDE THEM
#(1) set up model
model <- nimbleModel(SCRhab,
                     constants = nimConstants, 
                     data=nimData,
                     inits = inits(),
                     check = FALSE)
##ignore error message, only due to missing initial values at this stage
# [CM] I WOULD NOT IGNORE THE ERROR MESSAGE. FROM OUR EXPERIENCE, A SCR MODEL IN NIMBLE THAT IS 
# NOT PROPERLY INITIALIZED MAY NOT RUN CORRECTLY. TO CHECK WHICH PARAMETERS ARE NOT PROPERLY INITIALIZE 
# YOU CAN DO:
model$initializeInfo()
#THEN YOU SHOULD ALWAYS CHECK THAT THE MODEL IS ABLE TO CALCULATE A LIKELIHOOD (RETURN A VALUE) GIVEN THE INITIAL VALUES PROVIDED
model$calculate()#
# IF A -INF OF NA IS RETURNED YOU CAN CHECK WHERE THE PROBLEM COMES FROM WITH (AND THEN TRY TO FIX IT UNTIL THE -INF DISAPEARS):
model$logProb_p0# WILL GIVE YOU LIKELIHOOD OF P0
model$logProb_sxy# THE PROBLEM IS OFTEN WITH SXY OR Y
model$logProb_y

#(2) Compile model in c++
cmodel <- compileNimble(model)       

# (3) Configure MCMC - on an uncompiled model - this step allows setting which quantities to monitor
conf.mcmc <- configureMCMC(model, monitors = params, thin=1)

# (4) Build the MCMC sampler based on configurations
mcmc <- buildMCMC(conf.mcmc)

# (5) Compile sampler in c++ together with compiled model
cmcmc <- compileNimble(mcmc, project = cmodel, resetFunctions = TRUE)

# (6) Run (monitor time just for fun) [takes 20 seconds on my computer]
system.time(
  (samp <- runMCMC(cmcmc, niter = 25000, nburnin = 1000, nchains=3, inits = inits) )
)


setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/roads1")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/obsDensPre")
save(samp, file = "sampSCRdenscov2018.RData")

## -------------------------------------------------
##              ROAD DENSITY CATEGORY 4 
## ------------------------------------------------- 

rm(list = ls())

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Effort/France/Data/Systematic")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Effort/France/Data/Systematic")

#---- 1. LOAD THE DETECTION DATA ---- 
load("edf2017_2019_fr.RData")
load("tdf2018_fr.RData")

edf <- edf[which(edf$session == 2), ] # Keep only detections of 2017 (year 1)

# As the model is set, there can be only one capture per trap and occasion
# --> The number of trials is 7 (From may to November): So max number of captures per trap = 7
# To fix it in this edf, remove duplicates
edf <- edf[-which(duplicated(edf)), ]

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

# Set up a raster file 
setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
#setwd("~/Data_server/Variables_hrscale")
roads4 <- raster("roads4_hrbear.tif")

# Crop it to extent of state-space
habitat.r <- crop(roads4, e) 

#----   2.4 DEFINE THE BUFFER AREA AND CUT WHAT IS NOT HABITAT---- 
# Buffer around traps (5*sigma = 33200)
Xpoints <- X
coordinates(Xpoints) <- Xpoints[,c(1,2)]
Xbuf <- rgeos::gBuffer(Xpoints, width = 25000)

roads4Mask <- rasterize(Xbuf, habitat.r, mask = TRUE)

#RETAIN HABITAT COORDINATES THAT ARE WITHIN THE HABITAT
G <- coordinates(roads4Mask)[!is.na(roads4Mask[]),]
colnames(G) <- c("x","y")
#PLOT CHECK 
plot(roads4Mask)
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

###get cell coordinates for G.sc
windowCoords <- getWindowCoords(G.sc)
habitatGrid <- windowCoords$habitatGrid

#----   2.6 HABITAT COVARIATES ---- 
#[CM] THIS MIGHT BE WHERE THE ISSUE WAS
# AS YOU DID FOR THE COORDINATES THE HABITAT, YOU NEED TO SELECT HAB COORDS THAT ARE CONSIDERED AS HABITAT
X.d <- raster::values(roads4Mask)[!is.na(roads4Mask[])] ## !!! If I load nimbleSCR, the function values gives error

# Scale
X.d_mean <- mean(X.d)
X.d_sd <- sd(X.d)
X.d_sc <- (X.d - X.d_mean) / X.d_sd

#----   2.7 GET THE LOCAL DETECTORS OBJECTS  ---- 
# USE THE HABITAT GRID PROVIDED BY GETWINDOWCOORDS
habitatMask <- habitatGrid
habitatMask[habitatMask>0] <- 1 #TURN CELL ID TO 1 TO DEFINE THE HABITAT

# THERE MIGHT BE SOME TRIAL AND ERROR WHEN DEFINING THE DMAX.
localTraps <- getLocalTraps(habitatMask, X.sc, resizeFactor = 1, dmax = 5)

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

#---- 3. DETECTION DATA   ---- 
#----   3.1 MAKE Y   ---- 

K <- 7 # 7 occasions
n <- length(unique(edf$ind))

Y <- matrix(0, nrow = n, ncol = dim(X.sc)[1])
rownames(Y) <- unique(edf$ind)
xx <- edf[,c(2,4)]

for (obs in 1:nrow(xx)) {
  Y[xx[obs, 1], xx[obs, 2]] <- Y[xx[obs, 1], xx[obs, 2]] + 1
}

#----   3.2 AUGMENT Y   ---- 
##augment observed data to size M
M <- 400 
y.in <- rbind(Y, 
              matrix(0, nrow = M-n, ncol = J))

#----   3.3 USE SPARSE FORMAT FOR Y   ---- 
#change to 'sparse' format - speeds up computation by reducing file size
y.sparse <- getSparseY(y.in)

##extract pieces to be passed to Nimble
detNums <- y.sparse$detNums # Nº of traps at which each individual was detected
maxDetNums <- y.sparse$maxDetNums # Maximun Nº of traps at which any individual was detected
detIndices <- y.sparse$detIndices # ID of the traps where they were detected
y.sp <- y.sparse$y # Number of detections at each trap


#---- 4. FIT NIMBLE MODEL    ---- 
#----   4.1 CONSTANT AND DATA    ---- 

##compile constants
nimConstants <- list(
  M = M,
  J = J,
  numHabWindows = numHabWindows, 
  numGridRows = numGridRows,
  numGridCols = numGridCols, 
  maxDetNums = maxDetNums,
  MaxLocalTraps = MaxLocalTraps
)

##compile data
nimData <- list(habDens = X.d_sc,
                y = y.sp,
                detNums = detNums, 
                lowerHabCoords = lowerHabCoords,
                upperHabCoords = upperHabCoords,
                habitatGrid = habitatGrid,
                K = rep(K,J),
                X.sc = X.sc,
                habitatGridDet = habitatGridDet,
                detIndices = detIndices,
                detNums = detNums,
                localTrapsIndex = localTrapsIndex, 
                localTrapsNum = localTrapsNum 
)

#----   4.2 INITIAL VALUES  ---- 
#----     4.2.1 Z   ---- 

##set up initial values
z.in <- c(rep(1, n), rep(0, M-n))

#----     4.2.2 SXY  ---- 
##because of local evaluation of possible detectors, activity center initial 
##values have to be specified, eg average capture location
S.in <- matrix(NA, M, 2)
for ( i in 1:n){
  caps <- which(Y[i, ]>0)
  if (length(caps)==1){
    S.in[i,] <- as.matrix(X.sc[caps,])
  }else{
    S.in[i,] <- apply(X.sc[caps,], 2, mean)}
}


##random ACs for individuals never observed
##simulate them within state space!!
for(i in (n+1) : M){
  #[CM] THE ALTERNATIVE IS TO USE THE SIMULATION FUNCTIONNALITY OF NIMBLE TO SIMULATE INITIAL AC
  # IT IS THE SAME FUNCTION THAT YOU USE TO FIT THE MODEL, BUT IT STARTS WITH "r"
  S.in[i,] <- rbernppAC(n=1,
                        lowerCoords = nimData$lowerHabCoords,
                        upperCoords = nimData$upperHabCoords,
                        habitatGrid = habitatGrid,
                        numGridRows = nimConstants$numGridRows,
                        numGridCols = nimConstants$numGridCols,
                        logIntensities = rep(1, nimConstants$numHabWindows),#assume equal intensity across habitat cells (all 1)
                        logSumIntensity = sum(rep(1, nimConstants$numHabWindows))
  )
}
#----     4.2.2 COMPILE INITIAL VALUES   ---- 
inits <- function(){list(psi=runif(1,0.4, 0.6), 
                         sigma=runif(1,0.5, 1.5),
                         p0=runif(1,0,0.1),
                         z=z.in,
                         sxy=S.in,
                         beta.dens = runif(1,-0.1, 0.1))}

##source model code
setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Model")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Run_Sim/2.SCRdenscov")
source('2.SCRdenscov in Nimble.R')

##determine which parameters to monitor
params <- c('N', 'psi', 'sigma', 'p0', 'beta.dens')

#---- 5. FIT THE MODEL ---- 
#[CM] INITIAL VALUES ARE IMPORTANT IN NIMBLE AND WITH NIMBLE, PROVIDE THEM
#(1) set up model
model <- nimbleModel(SCRhab,
                     constants = nimConstants, 
                     data=nimData,
                     inits = inits(),
                     check = FALSE)
##ignore error message, only due to missing initial values at this stage
# [CM] I WOULD NOT IGNORE THE ERROR MESSAGE. FROM OUR EXPERIENCE, A SCR MODEL IN NIMBLE THAT IS 
# NOT PROPERLY INITIALIZED MAY NOT RUN CORRECTLY. TO CHECK WHICH PARAMETERS ARE NOT PROPERLY INITIALIZE 
# YOU CAN DO:
model$initializeInfo()
#THEN YOU SHOULD ALWAYS CHECK THAT THE MODEL IS ABLE TO CALCULATE A LIKELIHOOD (RETURN A VALUE) GIVEN THE INITIAL VALUES PROVIDED
model$calculate()#
# IF A -INF OF NA IS RETURNED YOU CAN CHECK WHERE THE PROBLEM COMES FROM WITH (AND THEN TRY TO FIX IT UNTIL THE -INF DISAPEARS):
model$logProb_p0# WILL GIVE YOU LIKELIHOOD OF P0
model$logProb_sxy# THE PROBLEM IS OFTEN WITH SXY OR Y
model$logProb_y

#(2) Compile model in c++
cmodel <- compileNimble(model)       

# (3) Configure MCMC - on an uncompiled model - this step allows setting which quantities to monitor
conf.mcmc <- configureMCMC(model, monitors = params, thin=1)

# (4) Build the MCMC sampler based on configurations
mcmc <- buildMCMC(conf.mcmc)

# (5) Compile sampler in c++ together with compiled model
cmcmc <- compileNimble(mcmc, project = cmodel, resetFunctions = TRUE)

# (6) Run (monitor time just for fun) [takes 20 seconds on my computer]
system.time(
  (samp <- runMCMC(cmcmc, niter = 25000, nburnin = 1000, nchains=3, inits = inits) )
)


setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/roads4")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/obsDensPre")
save(samp, file = "sampSCRdenscov2018.RData")

## -------------------------------------------------
##              ROAD DENSITY CATEGORY 5 
## ------------------------------------------------- 

rm(list = ls())

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Effort/France/Data/Systematic")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Effort/France/Data/Systematic")

#---- 1. LOAD THE DETECTION DATA ---- 
load("edf2017_2019_fr.RData")
load("tdf2018_fr.RData")

edf <- edf[which(edf$session == 2), ] # Keep only detections of 2017 (year 1)

# As the model is set, there can be only one capture per trap and occasion
# --> The number of trials is 7 (From may to November): So max number of captures per trap = 7
# To fix it in this edf, remove duplicates
edf <- edf[-which(duplicated(edf)), ]

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

# Set up a raster file 
setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
#setwd("~/Data_server/Variables_hrscale")
roads5 <- raster("roads5_hrbear.tif")

# Crop it to extent of state-space
habitat.r <- crop(roads5, e) 

#----   2.4 DEFINE THE BUFFER AREA AND CUT WHAT IS NOT HABITAT---- 
# Buffer around traps (5*sigma = 33200)
Xpoints <- X
coordinates(Xpoints) <- Xpoints[,c(1,2)]
Xbuf <- rgeos::gBuffer(Xpoints, width = 25000)

roads5Mask <- rasterize(Xbuf, habitat.r, mask = TRUE)

#RETAIN HABITAT COORDINATES THAT ARE WITHIN THE HABITAT
G <- coordinates(roads5Mask)[!is.na(roads5Mask[]),]
colnames(G) <- c("x","y")
#PLOT CHECK 
plot(roads5Mask)
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

###get cell coordinates for G.sc
windowCoords <- getWindowCoords(G.sc)
habitatGrid <- windowCoords$habitatGrid

#----   2.6 HABITAT COVARIATES ---- 
#[CM] THIS MIGHT BE WHERE THE ISSUE WAS
# AS YOU DID FOR THE COORDINATES THE HABITAT, YOU NEED TO SELECT HAB COORDS THAT ARE CONSIDERED AS HABITAT
X.d <- raster::values(roads5Mask)[!is.na(roads5Mask[])] ## !!! If I load nimbleSCR, the function values gives error

# Scale
X.d_mean <- mean(X.d)
X.d_sd <- sd(X.d)
X.d_sc <- (X.d - X.d_mean) / X.d_sd

#----   2.7 GET THE LOCAL DETECTORS OBJECTS  ---- 
# USE THE HABITAT GRID PROVIDED BY GETWINDOWCOORDS
habitatMask <- habitatGrid
habitatMask[habitatMask>0] <- 1 #TURN CELL ID TO 1 TO DEFINE THE HABITAT

# THERE MIGHT BE SOME TRIAL AND ERROR WHEN DEFINING THE DMAX.
localTraps <- getLocalTraps(habitatMask, X.sc, resizeFactor = 1, dmax = 5)

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

#---- 3. DETECTION DATA   ---- 
#----   3.1 MAKE Y   ---- 

K <- 7 # 7 occasions
n <- length(unique(edf$ind))

Y <- matrix(0, nrow = n, ncol = dim(X.sc)[1])
rownames(Y) <- unique(edf$ind)
xx <- edf[,c(2,4)]

for (obs in 1:nrow(xx)) {
  Y[xx[obs, 1], xx[obs, 2]] <- Y[xx[obs, 1], xx[obs, 2]] + 1
}

#----   3.2 AUGMENT Y   ---- 
##augment observed data to size M
M <- 400 
y.in <- rbind(Y, 
              matrix(0, nrow = M-n, ncol = J))

#----   3.3 USE SPARSE FORMAT FOR Y   ---- 
#change to 'sparse' format - speeds up computation by reducing file size
y.sparse <- getSparseY(y.in)

##extract pieces to be passed to Nimble
detNums <- y.sparse$detNums # Nº of traps at which each individual was detected
maxDetNums <- y.sparse$maxDetNums # Maximun Nº of traps at which any individual was detected
detIndices <- y.sparse$detIndices # ID of the traps where they were detected
y.sp <- y.sparse$y # Number of detections at each trap


#---- 4. FIT NIMBLE MODEL    ---- 
#----   4.1 CONSTANT AND DATA    ---- 

##compile constants
nimConstants <- list(
  M = M,
  J = J,
  numHabWindows = numHabWindows, 
  numGridRows = numGridRows,
  numGridCols = numGridCols, 
  maxDetNums = maxDetNums,
  MaxLocalTraps = MaxLocalTraps
)

##compile data
nimData <- list(habDens = X.d_sc,
                y = y.sp,
                detNums = detNums, 
                lowerHabCoords = lowerHabCoords,
                upperHabCoords = upperHabCoords,
                habitatGrid = habitatGrid,
                K = rep(K,J),
                X.sc = X.sc,
                habitatGridDet = habitatGridDet,
                detIndices = detIndices,
                detNums = detNums,
                localTrapsIndex = localTrapsIndex, 
                localTrapsNum = localTrapsNum 
)

#----   4.2 INITIAL VALUES  ---- 
#----     4.2.1 Z   ---- 

##set up initial values
z.in <- c(rep(1, n), rep(0, M-n))

#----     4.2.2 SXY  ---- 
##because of local evaluation of possible detectors, activity center initial 
##values have to be specified, eg average capture location
S.in <- matrix(NA, M, 2)
for ( i in 1:n){
  caps <- which(Y[i, ]>0)
  if (length(caps)==1){
    S.in[i,] <- as.matrix(X.sc[caps,])
  }else{
    S.in[i,] <- apply(X.sc[caps,], 2, mean)}
}


##random ACs for individuals never observed
##simulate them within state space!!
for(i in (n+1) : M){
  #[CM] THE ALTERNATIVE IS TO USE THE SIMULATION FUNCTIONNALITY OF NIMBLE TO SIMULATE INITIAL AC
  # IT IS THE SAME FUNCTION THAT YOU USE TO FIT THE MODEL, BUT IT STARTS WITH "r"
  S.in[i,] <- rbernppAC(n=1,
                        lowerCoords = nimData$lowerHabCoords,
                        upperCoords = nimData$upperHabCoords,
                        habitatGrid = habitatGrid,
                        numGridRows = nimConstants$numGridRows,
                        numGridCols = nimConstants$numGridCols,
                        logIntensities = rep(1, nimConstants$numHabWindows),#assume equal intensity across habitat cells (all 1)
                        logSumIntensity = sum(rep(1, nimConstants$numHabWindows))
  )
}
#----     4.2.2 COMPILE INITIAL VALUES   ---- 
inits <- function(){list(psi=runif(1,0.4, 0.6), 
                         sigma=runif(1,0.5, 1.5),
                         p0=runif(1,0,0.1),
                         z=z.in,
                         sxy=S.in,
                         beta.dens = runif(1,-0.1, 0.1))}

##source model code
setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Model")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Run_Sim/2.SCRdenscov")
source('2.SCRdenscov in Nimble.R')

##determine which parameters to monitor
params <- c('N', 'psi', 'sigma', 'p0', 'beta.dens')

#---- 5. FIT THE MODEL ---- 
#[CM] INITIAL VALUES ARE IMPORTANT IN NIMBLE AND WITH NIMBLE, PROVIDE THEM
#(1) set up model
model <- nimbleModel(SCRhab,
                     constants = nimConstants, 
                     data=nimData,
                     inits = inits(),
                     check = FALSE)
##ignore error message, only due to missing initial values at this stage
# [CM] I WOULD NOT IGNORE THE ERROR MESSAGE. FROM OUR EXPERIENCE, A SCR MODEL IN NIMBLE THAT IS 
# NOT PROPERLY INITIALIZED MAY NOT RUN CORRECTLY. TO CHECK WHICH PARAMETERS ARE NOT PROPERLY INITIALIZE 
# YOU CAN DO:
model$initializeInfo()
#THEN YOU SHOULD ALWAYS CHECK THAT THE MODEL IS ABLE TO CALCULATE A LIKELIHOOD (RETURN A VALUE) GIVEN THE INITIAL VALUES PROVIDED
model$calculate()#
# IF A -INF OF NA IS RETURNED YOU CAN CHECK WHERE THE PROBLEM COMES FROM WITH (AND THEN TRY TO FIX IT UNTIL THE -INF DISAPEARS):
model$logProb_p0# WILL GIVE YOU LIKELIHOOD OF P0
model$logProb_sxy# THE PROBLEM IS OFTEN WITH SXY OR Y
model$logProb_y

#(2) Compile model in c++
cmodel <- compileNimble(model)       

# (3) Configure MCMC - on an uncompiled model - this step allows setting which quantities to monitor
conf.mcmc <- configureMCMC(model, monitors = params, thin=1)

# (4) Build the MCMC sampler based on configurations
mcmc <- buildMCMC(conf.mcmc)

# (5) Compile sampler in c++ together with compiled model
cmcmc <- compileNimble(mcmc, project = cmodel, resetFunctions = TRUE)

# (6) Run (monitor time just for fun) [takes 20 seconds on my computer]
system.time(
  (samp <- runMCMC(cmcmc, niter = 25000, nburnin = 1000, nchains=3, inits = inits) )
)


setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/roads5")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/obsDensPre")
save(samp, file = "sampSCRdenscov2018.RData")

