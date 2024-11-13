## -------------------------------------------------
##                      openSCR + denscov
##                      ALL DATA 2017-2021
##                  Different trap arrays per year
##          Cov in p0: Effort + Type of trap + Behavioral response
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

setwd("D:/MargSalas/Scripts_MS/Stats/Nimble")
#source('dbinomLocal_normalBear.R')
source('dbinomLocal_normalBear_rbinom2.R')


setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")
#setwd("~/data/data/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")


#---- 1. LOAD THE DETECTION DATA ---- 

# DETECTIONS 

load("edf1721.RData")

edf <- edf[-which(duplicated(edf)), ] # As the model is set, there can be only one capture per trap and occasion
# --> The number of trials is 7 (From may to November): So max number of captures per trap = 7
# To fix it in this edf, remove duplicates

# TRAPS

load("tdf2017_effort.RData")
load("tdf2018_effort.RData")
load("tdf2019_effort.RData")
load("tdf2020_effort.RData")
load("tdf2021_effort.RData")

tdf_all <- rbind(tdf2017[,1:3], tdf2018[,1:3], tdf2019[,1:3]
                 ,tdf2020[,1:3], tdf2021[,1:3] # Join to define state space
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

# USE A RASTER TO GET A BASIS FOR THE HABITAT RASTER
# Set up a raster file 
setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
#setwd("~/data/data/Data_server/Variables_hrscale")
coreRest <- raster("core_hrbear_rest.tif")
distcore <- raster("logDistcore_hrbear.tif")

# Crop both to extent of state-space
habitat.coreRest <- crop(coreRest, e) 
habitat.distcore <- crop(distcore, e) 


habitat.coreRest.resampled <- resample(habitat.coreRest, habitat.distcore, method = 'bilinear') # Resample so that both cov have the same extent
values(habitat.coreRest.resampled)[values(habitat.coreRest.resampled)>0.49] <- 1
values(habitat.coreRest.resampled)[values(habitat.coreRest.resampled)<0.5] <- 0
habitat.coreRest.resampled[which(is.na(habitat.coreRest.resampled[]))] <- 0

#----   2.4 DEFINE THE BUFFER AREA AND CUT WHAT IS NOT HABITAT---- 
# Buffer around traps (5*sigma = 33200)
Xpoints <- X
coordinates(Xpoints) <- Xpoints[,c(1,2)]
Xbuf <- gBuffer(Xpoints, width = 25000)


distcoreMask <- rasterize(Xbuf, habitat.distcore, mask = TRUE)
coreRestMask <- rasterize(Xbuf, habitat.coreRest.resampled, mask = TRUE)

which(is.na(distcoreMask[])) == which(is.na(coreRestMask[]))


#RETAIN HABITAT COORDINATES THAT ARE WITHIN THE HABITAT
G <- coordinates(distcoreMask)[!is.na(distcoreMask[]),] #I get the habitat coordinates from the dist.core variable (but it doesnt matter, same extents)
colnames(G) <- c("x","y")

#PLOT CHECK 
plot(distcoreMask)
plot(coreRestMask)
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
windowCoords <- getWindowCoords(G.sc, plot.check = FALSE)
habitatGrid <- windowCoords$habitatGrid

#----   2.6 HABITAT COVARIATES ---- 
# AS YOU DID FOR THE COORDINATES THE HABITAT, YOU NEED TO SELECT HAB COORDS THAT ARE CONSIDERED AS HABITAT

# ---- 2.6.1. Density covariate ----

X.d <- values(distcoreMask)[!is.na(distcoreMask[])]

# Scale
X.d_mean <- mean(X.d)
X.d_sd <- sd(X.d)
X.d_sc <- (X.d - X.d_mean) / X.d_sd

# ---- 2.6.2. Survival covariate ----

X.surv <- values(coreRestMask)[!is.na(coreRestMask[])]
# I don't scale it, as is categorical


#----   2.7 GET THE LOCAL DETECTORS OBJECTS  ---- 
# USE THE HABITAT GRID PROVIDED BY GETWINDOWCOORDS
habitatMask <- habitatGrid
habitatMask[habitatMask>0] <- 1 #TURN CELL ID TO 1 TO DEFINE THE HABITAT

# Format trap array for nimble: ARRAY WITH TRAP MATRIX PER YEAR

Jyear <- c(nrow(tdf2017), nrow(tdf2018), nrow(tdf2019), nrow(tdf2020), nrow(tdf2021)) # Traps indexes
idx <- 1:nrow(tdf_all)
J.year <- list(idx[1]:Jyear[1], 
               (idx[Jyear[1]]+1):(idx[Jyear[1]] + idx[Jyear[2]]),
               (idx[Jyear[1]] + idx[Jyear[2]]+1): (idx[Jyear[1]] + idx[Jyear[2]] + idx[Jyear[3]]),
               (idx[Jyear[1]] + idx[Jyear[2]] + idx[Jyear[3]] + 1) : (idx[Jyear[1]] + idx[Jyear[2]] + idx[Jyear[3]] + idx[Jyear[4]]),
               (idx[Jyear[1]] + idx[Jyear[2]] + idx[Jyear[3]] + idx[Jyear[4]] + 1) : (idx[Jyear[1]] + idx[Jyear[2]] + idx[Jyear[3]] + idx[Jyear[4]] + idx[Jyear[5]]))


Tt <- 5
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


# Get one localtraps per year

##determine which traps are within some threshold distance of each habitat grid 
## cell - speeds up computations by only evaluating traps at which an individual
## could have been caught at (rather that traps very far away)
# THERE MIGHT BE SOME TRIAL AND ERROR WHEN DEFINING THE DMAX.

# 5 times log sigma
(5*6640)/5000

dmax <- c(20,10,10,10,8) ## PROBLEM: First year the local evaluation is almost useless


localTraps <- localTrapsNum.l <- MaxLocalTraps.l <- list()
for (t in 1:Tt){
  localTraps[[t]] <- getLocalObjects(habitatMask, Xt.sc[[t]], resizeFactor = 1, dmax = dmax[[t]], plot.check = FALSE)
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

#----   2.8. EFFORT COVARIATE ---- 

K <- 7 # 7 occasions

effort<-array(NA, c(max(Jyear), K, Tt))
tdf_list <- list(tdf2017, tdf2018, tdf2019, tdf2020, tdf2021)

for (t in 1:Tt){
  effort[1:Jyear[t],,t] <- as.numeric(as.matrix(tdf_list[[t]][1:Jyear[t], 6:12]))
}

# Create dummy variable

effort.dummy <- array(1, c(max(Jyear), K, Tt, 2)) # 4th dimension includes 2 arrays: one per level (intercept doesnt count)

# effort.dummy[,,,1] =0  effort.dummy[,,,2] =0 ==> 1 visit in France (cat 1): Intercept, no need to add
# effort.dummy[,,,1] =1  effort.dummy[,,,2] =0 ==> 2 visit in France (cat 2): Multiply b.effort1*array #1 
# effort.dummy[,,,1] =0  effort.dummy[,,,2] =1 ==> Spain (cat 3): Multiply b.effort2*array #2

for (t in 1:Tt){
  tmp <-tmp2 <- tmp3 <- effort[,,t]
  
  # Dummy variable 2 visits in France (only de 2 appear as 1)
  tmp2[tmp2[] %in% c(1,3)] <- 0
  tmp2[tmp2[] %in% c(2)] <- 1
  effort.dummy[,,t,1] <- tmp2
  
  # Dummy variable trap in Spain (only de 3 appear as 1)
  tmp3[tmp3[] %in% c(1,2)] <- 0
  tmp3[tmp3[] %in% c(3)] <- 1
  effort.dummy[,,t,2] <- tmp3
}

#----   2.9. TRAP COVARIATE ---- 

# 3D array to store it together with effort (fast run)
trap <- array(NA, c(max(Jyear), K, Tt)) 

for (t in 1:Tt){
  for(k in 1:K){
    trap[1:Jyear[t],k,t] <- as.numeric(as.factor(as.matrix(tdf_list[[t]][1:Jyear[t], 4])))
  }
}
trap[trap[] %in% c(2)] <- 0

#---- TRAP + EFFORT AS ARRAY FOR FAST FUNCTION ----#

effortarray <- array(1, c(max(Jyear), K, Tt, 3))
effortarray[,,,1] <- effort.dummy[,,,1]
effortarray[,,,2] <- effort.dummy[,,,2]
effortarray[,,,3] <- trap


#----   2.10. SEX COVARIATE ---- 

setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
info <- readxl::read_xlsx("D:/MargSalas/Oso/Datos/Tablas_finales/2022/info_individuals_2021.xlsx", sheet = 1)

#setwd("~/data/data/Data_server")
#info <- readxl::read_xlsx("~/data/data/Data_server/info_individuals_2021.xlsx", sheet = 1)

info <- info[,c(4,5)]
colnames(info)[1] <- "ind"

# Arrange in the same order than detections data frame below
indsex <- data.frame(ind = unique(edf$ind))
indsex <- left_join(indsex,info, by = "ind")

indsex$Sex[indsex$Sex == "F"] <- 0
indsex$Sex[indsex$Sex == "M"] <- 1

sex <- as.numeric(indsex$Sex)

#---- 3. DETECTION DATA   ---- 
#----   3.1 MAKE Y   ---- 

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

####
dim(Y)

#----   3.1.1 BEHAVIOURAL RESPONSE COVARIATE FROM Y   ---- 

prevcap <- list()
for (s in 1:dim(Y)[4]) {
  Ys <- Y[,,,s]
  prevcap[[s]] <- array(0, dim = c(dim(Ys)[1], dim(Ys)[2], 
                                   dim(Ys)[3]))
  first <- matrix(0, dim(Ys)[1], dim(Ys)[2])
  for (i in 1:dim(Ys)[1]) {
    for (j in 1:dim(Ys)[2]) {
      if (sum(Ys[i, j, ]) > 0) {
        first[i, j] <- min((1:(dim(Ys)[3]))[Ys[i, j, 
        ] > 0])
        prevcap[[s]][i, j, 1:first[i, j]] <- 0
        if (first[i, j] < dim(Ys)[3]) 
          prevcap[[s]][i, j, (first[i, j] + 1):(dim(Ys)[3])] <- 1
      }
    }
  }
  # zeros <- array(0, c(1, dim(prevcap[[s]])[2], dim(prevcap[[s]])[3]))
  # prevcap[[s]] <- abind(prevcap[[s]], zeros, along = 1)
}

##
prevcapArray <- array(0,dim(Y))
for(t in 1:dim(Y)[4]){
  prevcapArray[,,,t] <- prevcap[[t]]
}

#----   3.2 AUGMENT Y   ---- 
##augment observed data to size M
M <- 400 

y.in <- array(0, c(M, max(Jyear), K, Tt))
y.in[1:n,,,] <- Y

prevcapArray.in <- array(0, c(M, max(Jyear), K, Tt))
prevcapArray.in[1:n,,,] <- prevcapArray

##augment observed sex variable to size M (it becomes latent)
sex <- c(sex,rep(NA,length((n+1):M)))


#----   3.3 USE SPARSE FORMAT FOR Y   ---- 
##change to 'sparse' format - speeds up computation by reducing file size
## getSparseY cannot handle 4d arrays, so loop over years to get a 3d array per year
y.sparse <- list()
for (t in 1:Tt){
  y.sparse[[t]] <- getSparseY(y.in[,,,t]) 
}

##extract pieces to be passed to Nimble 
## Cyril changed this part of the code, so that we can use the function dbinomLocal_normal that speeds up the
# model. Now the detIndices and detNums are also contained within y.sp (lengthYCombined).
# This formulation is also necesary for making predictions


max.max <- max(sapply(y.sparse, function(x)x$lengthYCombined)) # ASP: Maximun of the maxDetNums

y.sp  <- array(NA, c(M, max.max, K, Tt))

for (t in 1:Tt){
  y.sp[,1:y.sparse[[t]]$lengthYCombined,,t] <- y.sparse[[t]]$yCombined
  #detIndices[,1:y.sparse[[t]]$maxDetNums,,t] <- y.sparse[[t]]$detIndices
}

# detNums<-array(NA, c(M, K, Tt))
# for (t in 1:Tt){
#   detNums[,,t] <- y.sparse[[t]]$detNums 
# }
# 
lengthYCombined <- sapply(y.sparse, function(x)x$lengthYCombined)


##number of trials per trap - now always 1 but needs to be passed to Nimble
# always 1, because we model each occasion separately
ones <- rep(1, max(Jyear))

#---- 4. FIT NIMBLE MODEL    ---- 

############################################################################################
### running a model in parallel ############################################################

##source code to run model in parallel 
setwd("D:/MargSalas/Scripts_MS/Stats/Nimble")
#setwd("~/data/data/Scripts_MS/Stats/Nimble")
#source("Parallel Nimble function FOR aNA2_model5-2.2.r")
source("Parallel Nimble function FOR aNA2_model6-1.r")



#----   4.1 CONSTANT AND DATA    ---- 

##compile constants
nimConstants <- list(
  M = M,
  J = Jyear,
  numHabWindows = numHabWindows, 
  numGridRows = numGridRows,
  numGridCols = numGridCols, 
  #maxDetNums = maxDetNums,
  lengthYCombined = lengthYCombined,
  MaxLocalTraps = MaxLocalTraps,
  nobs = n, 
  Nyr = Tt,
  K = K,
  effort = effortarray,
  nTrapCovs = dim(effortarray)[4],
  prevcap = prevcapArray.in
)

##compile data
nimData <- list(habDens = X.d_sc,
                habSurv = X.surv,
                y = y.sp,
                #detNums = detNums, 
                lowerHabCoords = lowerHabCoords,
                upperHabCoords = upperHabCoords,
                habitatGrid = habitatGrid,
                ones=ones,           # ASP: ones substitutes K=rep(K, n.max.traps)
                X.sc = Xt.sc.array,
                habitatGridDet = habitatGridDet,
                #detIndices = detIndices,
                localTrapsIndex = localTrapsIndex, 
                localTrapsNum = localTrapsNum,
                sex = sex
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

#----     4.2.3 SEX  ---- 

sex.in <- c(rep(NA,n), rep(0,length((n+1):M))) # NA for known values and 
sex.in[sample((n+1):M, 30, replace = FALSE)] <- 1 # random sex for augmented individuals

#----     4.2.2 COMPILE INITIAL VALUES   ---- 

S.in.sc_coords <- S.in.sc$coordsDataScaled

inits<-function(){list(gamma =c(0.5, rep(0.1, (Tt-1))), 
                       sigma = runif(2,0.5, 1.5),
                       p0=runif(1,0,0.1), # Value for p0 on the probability scale (0-1)
                       trapBetas = runif(3, 0.5,1),
                       # b.effort2 = runif(1, 0.5,1),
                       #b.trap = runif(1, 0.5,1),
                       b.bh = runif(1, 0.5,1),
                       omega = runif(1, 0.5,1),
                       sex = sex.in,
                       mu.phi=runif(1,-0.2,0.2),
                       beta.cov=runif(1,-0.2,0.2),
                       beta.dens = runif(1,-0.1, 0.1), # Added
                       z = z.in,
                       sxy = S.in.sc_coords,
                       sigD = runif(1, 1.5, 2.5))}

##source model code
setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Model")
#setwd("~/data/data/Scripts_MS/Oso/PopDyn/SCR/Model")
source('7.0.5. JS_OPSCR in Nimble_CatSurvCov_effortTrapBhCov_fastNewY.r')

##determine which parameters to monitor
params<-c('N', 'gamma', 'sigma', 'p0', 'trapBetas', 'b.bh', 'omega', 'mu.phi', 'beta.dens', 'sigD','R', 'pc.gam', 'Nsuper', 'beta.cov')


###### SAVE FOR RUNNING #####

modelcode = JS.SCRhab.Open.diftraps.sex.effortTrapBhCov.fast.surv

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/6.OPSCR_survCov/Data_server")
save(nimData, nimConstants, 
     inits, Tt, sex.in, z.in, S.in.sc_coords, 
     params, run_MCMC_allcode, modelcode, file = "Data_Model6-3.1.RData")


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
                          code = SCRhab.Open.diftraps.3d.effortTrapBhCov.sigSex.fast,   ##your model code
                          inits = inits,                 ##your inits function
                          constants = nimConstants,      ##your list of constants
                          params = params,               ##your vector with params to monitor
                          niter = 5,                  ##iterations per chain
                          nburnin = 1,                ##burn-in
                          nthin = 1,                  ##thinning, main parameters
                          Tt = Tt,                     ##additional objects needed within inits
                          z.in = z.in,
                          S.in.sc_coords = S.in.sc_coords,
                          sex.in = sex.in
)
new <- Sys.time() - old

## ALWAYS close cluster when model is done
stopCluster(this_cluster)

### output is a list 

#setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/4.openSCRdenscov_Effort")
setwd("~/data/data/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/4.openSCRdenscov_Effort")
save(chain_output, file = "sampOpenSCR_diftraps_effortTrapBhCov_1721_FINALDATA_3d.RData")





#### OPTION 2: NO PARALLEL (TO TRY INITIAL VALUES AND SEE IF MODEL WORKS) ####

#(1) set up model

model <- nimbleModel(JS.SCRhab.Open.diftraps.sex.effortTrapBhCov.fast.surv, constants = nimConstants, 
                     data=nimData, inits=inits(), check = FALSE,calculate = F)
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
  (samp <- runMCMC(cmcmc, niter = 5, nburnin = 0, nchains=3, inits = inits) )
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