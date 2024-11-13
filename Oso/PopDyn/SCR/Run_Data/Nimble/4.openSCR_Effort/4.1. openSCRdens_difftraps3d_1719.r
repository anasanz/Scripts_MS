## -------------------------------------------------
##                      openSCR + denscov
##                      French 2017-2019
##                        ONLY Syst
##                  Different trap arrays per year
##          Added 3d structure in effort but NOT covariate
## -------------------------------------------------


rm(list = ls())

library(nimble)
library(nimbleSCR)
library(raster)
library(rgeos)
library(terra)
library(sp)
library(dplyr)
library(parallel)

#setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")
setwd("~/data/data/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")


#---- 1. LOAD THE DETECTION DATA ---- 

# DETECTIONS 

load("edf1721.RData")

edf <- edf[-which(duplicated(edf)), ] # As the model is set, there can be only one capture per trap and occasion
                                      # --> The number of trials is 7 (From may to November): So max number of captures per trap = 7
                                      # To fix it in this edf, remove duplicates

# For the moment only 2017-2019, because I am not sure of the effort the rest of years

edf <- edf[which(edf$session %in% c(1,2,3)),] 

# TRAPS

load("tdf2017.RData")
load("tdf2018.RData")
load("tdf2019.RData")
#load("tdf2020.RData")
#load("tdf2021.RData")

tdf_all <- rbind(tdf2017[,1:4], tdf2018[,1:4], tdf2019[,1:4]
                 #,tdf2020[,1:4], tdf2021[,1:4] # Join to define state space
                    )
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
#setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
setwd("~/data/data/Data_server/Variables_hrscale")
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

dmax <- c(20,10,10
          #,10,8
          ) 

localTraps <- localTrapsNum.l <- MaxLocalTraps.l <- list()
for (t in 1:Tt){
  localTraps[[t]] <- getLocalObjects(habitatMask, Xt.sc[[t]], resizeFactor = 1, dmax = dmax[[t]])
  localTrapsNum.l[[t]]  <- localTraps[[t]]$numLocalIndices
  MaxLocalTraps.l[[t]]  <- localTraps[[t]]$numLocalIndicesMax
}

##check that there are always local traps available
lapply(localTrapsNum.l, function(x){min(x)})
##yes, every cell has at least some local traps

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

# Now Y needs to be filled by occasion

Y <- array(0, c(n, max(Jyear), K, Tt))
rownames(Y) <- unique(edf$ind)
xx <- edf[,c(1,2,3,4)]

for (t in 1:Tt){
  xxt <- xx[which(xx$session == t), ]
  
  for (k in 1:K){
    xxtk <- xxt[which(xxt$occ == k), ]
    
    if(nrow(xxtk) == 0) next # If there were no detections in that occasion
    
    for (obs in 1:nrow(xxtk)) {
      Y[xxtk[obs, 2], xxtk[obs, 4], xxtk[obs, 3], t] <- 1 # ASP: 1 because it can only be detected once per trap and occasion
    }
  }
}

max(Y) # Check the number max number of detections (it can't be higher than 1)


#----   3.2 AUGMENT Y   ---- 
##augment observed data to size M
M <- 400 

y.in <- array(0, c(M, max(Jyear), K, Tt))
y.in[1:n,,,] <- Y

#----   3.3 USE SPARSE FORMAT FOR Y   ---- 

##change to 'sparse' format - speeds up computation by reducing file size
## getSparseY cannot handle 4d arrays, so loop over years to get a 3d array per year
y.sparse <- list()
for (t in 1:Tt){
  y.sparse[[t]] <- getSparseY(y.in[,,,t]) 
}

##extract pieces to be passed to Nimble - this is now a little more complex
##because each year's data has different dimensions
#get max(maxDetNums) over years
max.max <- max(sapply(y.sparse, function(x)x$maxDetNums)) # ASP: Maximun of the maxDetNums

y.sp <- detIndices <- array(NA, c(M, max.max, K, Tt))

for (t in 1:Tt){
  y.sp[,1:y.sparse[[t]]$maxDetNums,,t] <- y.sparse[[t]]$y
  detIndices[,1:y.sparse[[t]]$maxDetNums,,t] <- y.sparse[[t]]$detIndices
}

detNums<-array(NA, c(M, K, Tt))
for (t in 1:Tt){
  detNums[,,t] <- y.sparse[[t]]$detNums 
}

maxDetNums <- sapply(y.sparse, function(x)x$maxDetNums)

##number of trials per trap - now always 1 but needs to be passed to Nimble
# always 1, because we model each occasion separately
ones <- rep(1, max(Jyear))

#---- 4. FIT NIMBLE MODEL    ----

############################################################################################
### running a model in parallel ############################################################

##source code to run model in parallel 
#setwd("D:/MargSalas/Scripts_MS/Stats/Nimble")
setwd("~/data/data/Scripts_MS/Stats/Nimble")
source("Parallel Nimble function FOR aNA2_model3.2.r") 

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
  Nyr=Tt,
  K=K
)

##compile data
nimData <- list(habDens = X.d_sc,
                y = y.sp,
                detNums = detNums, 
                lowerHabCoords = lowerHabCoords,
                upperHabCoords = upperHabCoords,
                habitatGrid = habitatGrid,
                ones=ones,           # ASP: ones substitutes K=rep(K, n.max.traps)
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
  Ymat <- Y[,,,t]
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
    caps <- which(apply(Y[i,,,t],1,sum) > 0) ## ASP: Get in which traps the ind was captured at year t
                                              # With 3d data Need to sum the rows (all occasions) to get idx right
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

#----     4.2.2 COMPILE INITIAL VALUES   ---- 

S.in.sc_coords <- S.in.sc$coordsDataScaled

inits<-function(){list(gamma = c(0.5, rep(0.1, (Tt-1))), 
                       sigma = runif(1,0.5, 1.5),
                       p0 = runif(1,0,1),
                       phi = runif(1,0.5,1),
                       beta.dens = runif(1,-0.1, 0.1), # Added
                       z = z.in,
                       sxy = S.in.sc_coords,
                       sigD  =runif(1, 1.5, 2.5))}

##source model code
#setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Model")
setwd("~/data/data/Scripts_MS/Oso/PopDyn/SCR/Model")
source('SCR in Nimble_diftraps_difeff.R')

##determine which parameters to monitor
params<-c('N', 'gamma', 'sigma', 'p0', 'phi', 'beta.dens', 'sigD','R', 'pc.gam', 'Nsuper')


#### OPTION 1: PARALLEL ####
detectCores()

##start cluster w/ 3 cores (for 3 chains)
this_cluster <- makeCluster(3)


##run wrapper function in parallel - cl and X need to be this way
## cl defines which cluster to use, X provides random number seeds to each core

old <- Sys.time()

chain_output <- parLapply(cl = this_cluster, X = 1:3, 
                          fun = run_MCMC_allcode,      ##function in "Parallel Nimble function.R"
                          data = nimData,              ##your data list
                          code = SCRhab.Open.diftraps.3d,   ##your model code
                          inits = inits,                 ##your inits function
                          constants = nimConstants,      ##your list of constants
                          params = params,               ##your vector with params to monitor
                          niter = 150000,                  ##iterations per chain
                          nburnin = 100000,                ##burn-in
                          nthin = 10,                  ##thinning, main parameters
                          Tt = Tt,                     ##additional objects needed within inits
                          z.in = z.in,
                          S.in.sc_coords = S.in.sc_coords 
)
new <- Sys.time() - old

## ALWAYS close cluster when model is done
stopCluster(this_cluster)

### output is a list 

#setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/2.openSCRdenscov")
setwd("~/data/data/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/2.openSCRdenscov")
save(chain_output, file = "sampOpenSCR_diftraps_1719_FINALDATA_3d.RData")



#### OPTION 2: NO PARALLEL (TO TRY INITIAL VALUES AND SEE IF MODEL WORKS) ####

#(1) set up model

model <- nimbleModel(SCRhab.Open.diftraps.3d, constants = nimConstants, 
                     data=nimData, inits=inits(), check = FALSE)
##ignore error message, only due to missing initial values at this stage

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

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/2.openSCRdenscov")
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