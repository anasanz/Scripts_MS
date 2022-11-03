## -------------------------------------------------
##                      SCR + denscov
##                  Final Systematic Data 2020 (fr+Sp)
## ------------------------------------------------- 
rm(list = ls())

library(nimble)
library(MCMCvis)
library(nimbleSCR)
library(raster)
library(rgeos)
library(terra)
library(sp)
library(dplyr)


setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_prefinal_1721")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_final_1719")

#---- 1. LOAD THE DETECTION DATA ---- 
load("edf1721.RData")
load("tdf2021.RData")

edf <- edf[which(edf$session == 5), ] # Keep only detections of 2017 (year 1)

# As the model is set, there can be only one capture per trap and occasion
# --> The number of trials is 7 (From may to November): So max number of captures per trap = 7
# To fix it in this edf, remove duplicates
edf <- edf[-which(duplicated(edf)), ]
t <- edf %>% group_by(ind,occ, trap) %>% # All need to be one
  summarise(n()) 


# We remove Nere and Goiat, which are two individuals moving a lot that don't represent
# the rest of the population.

edf <- edf[-which(edf$ind %in% c("Nere", "Goiat")), ]

#---- 2. DEFINE THE TRAP AND THE HABITAT  ---- 
#----   2.1 GET TRAPS---- 
X <- tdf2021[,c(2,3)]
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

## ---- *** START LOOP TO RUN SINGLE-COVARIATE MODELS *** ----

covs <- c("forest", "dem", "rough", "slope", "logDistcore", "obsDens200m", "obsDens200m_preST",
          "roads1", "roads4", "roads5", "roads6")

for (xxx in 1:length(covs)) {
  
  #----   2.3 GET A RASTER FOR THE HABITAT ---- 
  # USE A FOREST RASTER TO GET A BASIS FOR THE HABITAT RASTER
  
  # Set up a raster file 
  setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
  #setwd("~/Data_server/Variables_hrscale")
  cov <- raster(paste(covs[xxx],"_hrbear.tif", sep = ""))
  
  # Crop it to extent of state-space
  habitat.r <- crop(cov, e) 
  
  #----   2.4 DEFINE THE BUFFER AREA AND CUT WHAT IS NOT HABITAT---- 
  # Buffer around traps (5*sigma = 33200)
  Xpoints <- X
  coordinates(Xpoints) <- Xpoints[,c(1,2)]
  Xbuf <- rgeos::gBuffer(Xpoints, width = 25000)
  
  covMask <- rasterize(Xbuf, habitat.r, mask = TRUE)
  
  #RETAIN HABITAT COORDINATES THAT ARE WITHIN THE HABITAT
  G <- coordinates(covMask)[!is.na(covMask[]),]
  colnames(G) <- c("x","y")
  #PLOT CHECK 
  plot(covMask)
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
  X.d <- raster::values(covMask)[!is.na(covMask[])] ## !!! If I load nimbleSCR, the function values gives error
  
  # Scale
  X.d_mean <- mean(X.d)
  X.d_sd <- sd(X.d)
  X.d_sc <- (X.d - X.d_mean) / X.d_sd
  
  #----   2.7 GET THE LOCAL DETECTORS OBJECTS  ---- 
  # USE THE HABITAT GRID PROVIDED BY GETWINDOWCOORDS
  habitatMask <- habitatGrid
  habitatMask[habitatMask>0] <- 1 #TURN CELL ID TO 1 TO DEFINE THE HABITAT
  
  # THERE MIGHT BE SOME TRIAL AND ERROR WHEN DEFINING THE DMAX.
  localTraps <- getLocalTraps(habitatMask, X.sc, resizeFactor = 1, dmax = 10)
  
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
  #$logProb_p0# WILL GIVE YOU LIKELIHOOD OF P0
  #model$logProb_sxy# THE PROBLEM IS OFTEN WITH SXY OR Y
  #model$logProb_y
  
  # Error in y, individual 4 (Rodri)
  # warning: logProb of data node y[4, 1:10, 1]: logProb is -Inf.
  # Y is exactly the same as in the openSCR model that works perfectly, with detections in traps 16,385 and 453
  # TE LOCAL EVALUATION WAS TOO SMALL! Increasing to 10 it worked
  
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
  
  setwd(paste("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_17_19/", covs[xxx], sep = ""))
  #setwd("~/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/obsDensPre")
  save(samp, file = "sampSCRdenscov2021.RData")
  
}