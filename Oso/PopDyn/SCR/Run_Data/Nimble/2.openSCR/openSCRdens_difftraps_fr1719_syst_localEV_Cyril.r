## -------------------------------------------------
##                      openSCR + denscov
##                      French 2017-2019
##                        ONLY Syst
##                  Different trap arrays per year
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

##NOTE: load two functions modified by Cyril to allow yearly efficient local evaluation

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Sim/3.openSCR")
source('getLocalObjects.R')
source('dbinomLocal_normal.R')

#---- 1. LOAD THE DETECTION DATA ---- 

#setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Effort/France/Data/Systematic")
setwd("~/Scripts_MS/Oso/PopDyn/SCR/Effort/France/Data/Systematic")

load("edf2017_2019_fr.RData")

# As the model is set, there can be only one capture per trap and occasion
# --> The number of trials is 7 (From may to November): So max number of captures per trap = 7
# To fix it in this edf, remove duplicates
edf <- edf[-which(duplicated(edf)), ]
t <- edf %>% group_by(ind,occ, session,trap) %>% # All need to be one
  summarise(n()) 

load("tdf2017_fr.RData")
load("tdf2018_fr.RData")
load("tdf2019_fr.RData")

tdf_all <- rbind(tdf2017[,1:3], tdf2018[,1:3], tdf2019[,1:3]) # Join to define state space
rownames(tdf_all) <- 1:nrow(tdf_all)


# We remove Nere and Goiat, which are two individuals moving a lot that don't represent
# the rest of the population.

edf <- edf[-which(edf$ind %in% c("Néré", "Goiat")), ]

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
#setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
setwd("~/Data_server/Variables_hrscale")

distcore <- raster("logDistcore_hrbear.tif")

# Crop it to extent of state-space
habitat.r <- crop(distcore, e) 

#----   2.4 DEFINE THE BUFFER AREA AND CUT WHAT IS NOT HABITAT---- 
# Buffer around traps (5*sigma = 33200)
Xpoints <- X
coordinates(Xpoints) <- Xpoints[,c(1,2)]
Xbuf <- gBuffer(Xpoints, width = 25000)

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
X.d <- values(distcoreMask)[!is.na(distcoreMask[])]

# Scale
X.d_mean <- mean(X.d)
X.d_sd <- sd(X.d)
X.d_sc <- (X.d - X.d_mean) / X.d_sd

#----   2.7 GET THE LOCAL DETECTORS OBJECTS  ---- 
# USE THE HABITAT GRID PROVIDED BY GETWINDOWCOORDS
habitatMask <- habitatGrid
habitatMask[habitatMask>0] <- 1 #TURN CELL ID TO 1 TO DEFINE THE HABITAT

# Format trap array for nimble: ARRAY WITH TRAP MATRIX PER YEAR

Jyear <- c(nrow(tdf2017), nrow(tdf2018), nrow(tdf2019)) # Traps indexes
idx <- 1:nrow(tdf_all)
J.year <- list(idx[1]:Jyear[1], 
               (idx[Jyear[1]]+1):(idx[Jyear[1]] + idx[Jyear[2]]), 
               (idx[Jyear[1]] + idx[Jyear[2]]+1): length(idx))
Tt <- 3
Yrs <- seq(1:Tt)

Xt <- Xt.sc <- list() 
for (t in 1:Tt){
  Xt.sc[[t]] <- X.sc[J.year[[t]],] # For the getLocalTraps function
  Xt[[t]] <- X[J.year[[t]],] # To simulate yearly detection data later on
}

Xt.sc.array <- array(NA, c(max(Jyear), 2, Tt))
for (t in 1:Tt){
  Xt.sc.array[1:Jyear[t],,t] <- as.matrix(Xt.sc[[t]])
}

Xt.sc.array[1:Jyear[1], 1:2, 1] # Example: to get traps from year 1

# Visualize traps per year before local evaluation
for (t in 1:Tt){
  plot(distcoreMask, main = t)
  points(Xpoints[J.year[[t]], ])
  points(G[,2]~G[,1], col="red", cex=0.5)
} # We will need to use a different dmax per year (in year 1 there are no traps in a region)


# Get one localtraps per year

##determine which traps are within some threshold distance of each habitat grid 
## cell - speeds up computations by only evaluating traps at which an individual
## could have been caught at (rather that traps very far away)
# THERE MIGHT BE SOME TRIAL AND ERROR WHEN DEFINING THE DMAX.

# 5 times log sigma
(5*6640)/5000

dmax <- c(10,10,10) ## The three dmax should be the same with the new function

localTraps <- localTrapsNum.l <- MaxLocalTraps.l <- list()
for (t in 1:Tt){
  localTraps[[t]] <- getLocalObjects(habitatMask, Xt.sc[[t]], resizeFactor = 1, dmax = dmax[[t]])
  localTrapsNum.l[[t]]  <- localTraps[[t]]$numLocalIndices
  MaxLocalTraps.l[[t]]  <- localTraps[[t]]$numLocalIndicesMax
}

##check that there are always local traps available
lapply(localTrapsNum.l, function(x){min(x)})
##not in year 1 but there is no problem now, the function has been modified

# Store objects of localtraps as numeric objects because nimble doesn't support lists

localTrapsNum <- do.call(cbind, localTrapsNum.l) # Matrix with cols = years
habitatGridDet <- localTraps[[1]]$habitatGrid # The habitat grid is the same for all
MaxLocalTraps <- unlist(MaxLocalTraps.l) # This will also act as an index to map the dimensions of localTrapsIndex in the model

localTrapsIndex <- array(NA, c(max(habitatGridDet), max(MaxLocalTraps), Tt)) # Array with dimensions of maximun maxlocaltraps, will need to subset in the model
for (t in 1:Tt){
  localTrapsIndex[,1:MaxLocalTraps[t],t] <- localTraps[[t]]$localIndices
}

##some characteristics of S, not affected by changing trap array
numHabWindows <- sum(habitatGrid !=0) #number of cells in S
numGridRows <- nrow(localTraps[[1]]$habitatGrid) # I take the first year of local traps but it doesn't matter, all the same
numGridCols <- ncol(localTraps[[1]]$habitatGrid)


#the following two are still passed to function but no longer used by it
lowerHabCoords <- windowCoords$lowerHabCoords
upperHabCoords <- windowCoords$upperHabCoords

#---- 3. DETECTION DATA   ---- 
#----   3.1 MAKE Y   ---- 

K <- 7 # 7 occasions
n <- length(unique(edf$ind))


Y <- array(0, c(n, max(Jyear), Tt))
rownames(Y) <- unique(edf$ind)
xx <- edf[,c(1,2,4)]

for (t in 1:Tt){
  xxt <- xx[which(xx$session == t), ]
  for (obs in 1:nrow(xxt)) {
    Y[xxt[obs, 2], xxt[obs, 3], t] <-  Y[xxt[obs, 2], xxt[obs, 3], t] + 1
  }
}

max(Y) # Check the number max number of detections (it can't be higher than K)

#----   3.2 AUGMENT Y   ---- 
##augment observed data to size M
M <- 400 

y.in <- array(0, c(M, max(Jyear), Tt))
y.in[1:n,,] <- Y
c <- Y[,,2]

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
  J = Jyear,
  numHabWindows = numHabWindows, 
  numGridRows = numGridRows,
  numGridCols = numGridCols, 
  maxDetNums = maxDetNums,
  MaxLocalTraps = MaxLocalTraps,
  nobs=n, 
  Nyr=Tt
)

##compile data
nimData <- list(habDens = X.d_sc,
                y = y.sp,
                detNums = detNums, 
                lowerHabCoords = lowerHabCoords,
                upperHabCoords = upperHabCoords,
                habitatGrid = habitatGrid,
                K = rep(K,max(Jyear)),
                X.sc = Xt.sc.array,
                habitatGridDet = habitatGridDet,
                detIndices = detIndices,
                detNums = detNums,
                localTrapsIndex = localTrapsIndex, 
                localTrapsNum = localTrapsNum 
)


#----   4.2 INITIAL VALUES  ---- 
#----     4.2.1 Z   ---- 

##set up initial values

# Alive state
# Get first year detected
z <- matrix(NA,nrow = n, ncol = Tt)
for (t in 1:Tt){
  Ymat <- Y[,,t]
  z[,t] <- apply(Ymat,1,sum)
}
z[z>1] <- 1
first <- apply(z,1,function(x)min(which(x==1)))

## ASP: From the first year detected fill as alive
z.in <- matrix(0, M, Tt)
for(i in 1:n){
  z.in[i,first[i]:Tt]<-1
}  


#----     4.2.2 SXY  ---- 

##because of local evaluation of possible detectors, activity center initial 
##values have to be specified, eg average capture location in 'model' space

S.in <- array(NA, c(M, 2, Tt))

for ( i in 1:n){
  for (t in 1:Tt){
    caps <- which(Y[i,,t] > 0) ## ASP: Get in which traps the ind was captured at year t
    
    if (length(caps)==0) next #fill in missing ACs with reasonable values later
    if (length(caps)==1){ ## ASP: If its only in 1 trap, a put the trap location as AC
      S.in[i,,t]<- as.numeric(Xt[[t]][caps,])
    }else{ 
      S.in[i,,t]<- as.numeric(apply(Xt[[t]][caps,],2,mean))}  ## ASP: If its > 1 trap average location
  }
}

##fill in missing ACs as average of 'observed' ACs in nearby time step

for (i in 1:n){
  #which ACs unobserved
  nac <- which(is.na(S.in[i,1,]))
  wac <- (1:Tt)[-nac] ## ASP: Only the years that are observed
  #for those, use mean observed
  S.in[i,1,nac]<-mean(S.in[i,1,wac]) ## ASP: Use the mean of the observed to fill unobserved
  S.in[i,2,nac]<-mean(S.in[i,2,wac]) 
}


##random ACs for individuals never observed
for(i in (n+1) : M){
  for (t in 1:Tt){
    ssg<-sample(1:length(X.d), 1)
    S.in[i,,t]<-G[ssg,]
  }
}

colnames(S.in) <- c('x', 'y')
S.in.sc <- scaleCoordsToHabitatGrid(S.in, G)

###random ACs for individuals never observed
###simulate them within state space!!
##cellACv <- list() # To monitor the cells where there was a problem with initial values
#for(i in (n+1) : M){
#  #[CM] THE ALTERNATIVE IS TO USE THE SIMULATION FUNCTIONNALITY OF NIMBLE TO SIMULATE INITIAL AC
#  # IT IS THE SAME FUNCTION THAT YOU USE TO FIT THE MODEL, BUT IT STARTS WITH "r"
#  S.in[i,] <- rbernppAC(n=1,
#                        lowerCoords = nimData$lowerHabCoords,
#                        upperCoords = nimData$upperHabCoords,
#                        habitatGrid = habitatGrid,
#                        numGridRows = nimConstants$numGridRows,
#                        numGridCols = nimConstants$numGridCols,
#                        logIntensities = rep(1, nimConstants$numHabWindows),#assume equal intensity across habitat cells (all 1)
#                        logSumIntensity = sum(rep(1, nimConstants$numHabWindows))
#  )
#  # 
#  # cellAC <- sample(c(1:max(habitatGridDet)),1)
#  # cellACv[[i]] <- cellAC
#  # S.in[i,] <- c(which(habitatGridDet==cellAC, arr.ind=TRUE)[2] - runif(1, 0.0, 1.0), # The column of the habitatGridDet is the x coordinate
#  #               which(habitatGridDet==cellAC, arr.ind=TRUE)[1] - runif(1, 0.0, 1.0)) # The row of the habitatGridDet is the y coordinate
#}


#----     4.2.2 COMPILE INITIAL VALUES   ---- 

inits<-function(){list(gamma=c(0.5, rep(0.1, (Tt-1))), 
                       sigma=runif(1,0.5, 1.5),
                       p0=runif(1,0,1),
                       phi=runif(1,0.5,1),
                       beta.dens = runif(1,-0.1, 0.1), # Added
                       z=z.in,
                       sxy=S.in.sc$coordsDataScaled,
                       sigD=runif(1, 1.5, 2.5))}

##source model code
##I prefer working on code in a separate script but you can also have everything in
##one script and just execute the code

#setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Model")
setwd("~/Scripts_MS/Oso/PopDyn/SCR/Model")

source('3.2.SCRopen_diftraps in Nimble.R')

##determine which parameters to monitor
params<-c('N', 'gamma', 'sigma', 'p0', 'phi', 'beta.dens', 'sigD','R', 'pc.gam', 'Nsuper')

#(1) set up model

model <- nimbleModel(SCRhab.Open.diftraps, constants = nimConstants, 
                     data=nimData, inits=inits(), check = FALSE)
##ignore error message, only due to missing initial values at this stage

model$initializeInfo()
#THEN YOU SHOULD ALWAYS CHECK THAT THE MODEL IS ABLE TO CALCULATE A LIKELIHOOD (RETURN A VALUE) GIVEN THE INITIAL VALUES PROVIDED
model$calculate()#
# IF A -INF OF NA IS RETURNED YOU CAN CHECK WHERE THE PROBLEM COMES FROM WITH (AND THEN TRY TO FIX IT UNTIL THE -INF DISAPEARS):
model$logProb_p0# WILL GIVE YOU LIKELIHOOD OF P0
model$logProb_sxy# THE PROBLEM IS OFTEN WITH SXY OR Y
model$logProb_y[,,3]



#(2) Compile model in c++
#     In complex models, this step can take a while (as well as step 5)
#     Much longer than in JAGS, but the model typically runs much faster
cmodel <- compileNimble(model)       

# (3) Configure MCMC - on an uncompiled model - this step allows setting which quantities to monitor
#     Also, nimble allows two sets of monitors, these can be thinned at different rates
#     all of which is more important in complex models but not to start with
conf.mcmc<-configureMCMC(model, monitors = params, thin=10)

# (4) Build the MCMC sampler based on configurations
mcmc <- buildMCMC(conf.mcmc)

# (5) Compile sampler in c++ together with compiled model
cmcmc <- compileNimble(mcmc, project = cmodel, resetFunctions = TRUE)

# (6) Run (monitor time just for fun) [takes 20 seconds on my computer]
system.time(
  (samp <- runMCMC(cmcmc, niter = 150000, nburnin = 100000, nchains=3, inits = inits) )
)

#setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/2.openSCRdenscov")
setwd("~/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/2.openSCRdenscov")

save(samp, file = "sampOpenSCR_diftraps.RData")

##remove NAs
inn<-colnames(samp[[1]])
remm<-pmatch(c("R[1]", "pc.gam[1]"), inn)
samp2<-lapply(samp, function(x)x[,-remm])

##NOTE: summary command is from MCMCvis package; that also has good plotting options
## summary table for everything in "params" vector
summ<-MCMCsummary(samp2)
MCMCtrace(samp)


# Density
summ$mean[1]/max(habitatGrid) # 65/605 = 0.1074799