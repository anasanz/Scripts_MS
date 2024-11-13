## -------------------------------------------------
##                      SCR + denscov
##                      French 2017
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


setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Effort/France/Data/forOpenPop")

#---- 1. LOAD THE DETECTION DATA ---- 
load("tdf_all.RData")
load("edf_all.RData")

# We remove Nere and Goiat, which are two individuals moving a lot that don't represent
# the rest of the population.

edf <- new_edf[-which(new_edf$ind %in% c("Néré", "Goiat")), ]

edf[which(edf$ind == "Cannellito" & edf$trap == 210 & edf$session == 3),]

#---- 2. DEFINE THE TRAP AND THE HABITAT  ---- 
#----   2.1 GET TRAPS---- 
X <- tdf[,c(2,3)]
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
#setwd("~/Data_server")
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

# 5 times log sigma
(5*6640)/5000

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
Tt <- 3 #number of years


Y <- array(0, c(n, dim(X.sc)[1], Tt))
rownames(Y) <- unique(edf$ind)
xx <- edf[,c(1,2,4)]

for (t in 1:Tt){
  xxt <- xx[which(xx$session == t), ]
  for (obs in 1:nrow(xxt)) {
    Y[xxt[obs, 2], xxt[obs, 3], t] <-  Y[xxt[obs, 2], xxt[obs, 3], t] + 1
  }
}


#----   3.2 AUGMENT Y   ---- 
##augment observed data to size M
M <- 400 

y.in <- array(0, c(M, J, Tt))
y.in[1:n,,] <- Y

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
      S.in[i,,t]<- as.numeric(X[caps,])
    }else{ 
      S.in[i,,t]<- as.numeric(apply(X[caps,],2,mean))}  ## ASP: If its > 1 trap average location
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
                       z=z.in,
                       sxy=S.in.sc$coordsDataScaled,
                       sigD=runif(1, 1.5, 2.5))}

##source model code
##I prefer working on code in a separate script but you can also have everything in
##one script and just execute the code

setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Oso/PopDyn/SCR/Sim/openSCR")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Sim/SCRdenscov")
source('SCR in Nimble.R')

##determine which parameters to monitor
params<-c('N', 'gamma', 'sigma', 'p0', 'phi', 'beta.dens', 'sigD','R', 'pc.gam', 'Nsuper')

#(1) set up model

model <- nimbleModel(SCRhab.Open, constants = nimConstants, 
                     data=nimData, inits=inits(), check = FALSE)
##ignore error message, only due to missing initial values at this stage

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


################
#warning: logProb of data node y[7, 1:12, 3]: logProb is -Inf.

##remove NAs
inn<-colnames(samp[[1]])
remm<-pmatch(c("R[1]", "pc.gam[1]"), inn)
samp2<-lapply(samp, function(x)x[,-remm])

##NOTE: summary command is from MCMCvis package; that also has good plotting options
## summary table for everything in "params" vector
summ<-MCMCsummary(samp2)
MCMCtrace(samp)




setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Oso/PopDyn/SCR/Data/Nimble/Results/SCRdenscov_year")
save(samp, file = "sampSCRdenscov2017.RData")

# Density
summ$mean[1]/max(habitatGrid) # 65/605 = 0.1074799

