## -------------------------------------------------
##                      SCR + denscov
##                  Final Systematic Data 2020 (fr+Sp)
##                      Covs on p0
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
library(parallel)



setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")

#---- 1. LOAD THE DETECTION DATA ---- 
load("edf1721.RData")

edf <- edf[which(edf$session == 4), ] # Keep only detections of 2017 (year 1)

# As the model is set, there can be only one capture per trap and occasion
# --> The number of trials is 7 (From may to November): So max number of captures per trap = 7
# To fix it in this edf, remove duplicates
edf <- edf[-which(duplicated(edf)), ]
t <- edf %>% group_by(ind,occ, trap) %>% # All need to be one
  summarise(n()) 


# We remove Nere and Goiat, which are two individuals moving a lot that don't represent
# the rest of the population.

edf <- edf[-which(edf$ind %in% c("Nere", "Goiat")), ]

## Traps
load("tdf2020_effort.RData")
tdf <- tdf2020


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

## ---- *** START LOOP TO RUN SINGLE-COVARIATE MODELS *** ----

covs <- c("forest", "slope", "logDistcore", "roads1")
#,"dem", "rough", "obsDens200m", "obsDens200m_preST",
#, "roads4", "roads5", "roads6")

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
  #plot(covMask)
  #points(Xpoints)
  #points(G[,2]~G[,1],col="red",cex=0.5)
  
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
  localTraps <- getLocalObjects(habitatMask, X.sc, resizeFactor = 1, dmax = 10)
  
  localTrapsIndex <- localTraps$localIndices
  localTrapsNum <- localTraps$numLocalIndices
  habitatGridDet <- localTraps$habitatGrid
  MaxLocalTraps <- localTraps$numLocalIndicesMax
  
  ##name and structure for Nimble model
  numHabWindows <- dim(localTrapsIndex)[1] #number of cells in S
  numGridRows <- nrow(localTraps$habitatGrid)
  numGridCols <- ncol(localTraps$habitatGrid)
  
  
  #the following two are still passed to function but no longer used by it
  lowerHabCoords <- windowCoords$lowerHabCoords
  upperHabCoords <- windowCoords$upperHabCoords
  
  #----   2.8. EFFORT COVARIATE ---- 
  
  K <- 7 # 7 occasions
  
  effort<-matrix(as.numeric(as.matrix(tdf[ ,6:12])), nrow = nrow(tdf), ncol = K)
  
  # Create dummy variable
  
  effort.dummy <- array(1, c(nrow(tdf), K, 2)) # 4th dimension includes 2 arrays: one per level (intercept doesnt count)
  
  # effort.dummy[,,,1] =0  effort.dummy[,,,2] =0 ==> 1 visit in France (cat 1): Intercept, no need to add
  # effort.dummy[,,,1] =1  effort.dummy[,,,2] =0 ==> 2 visit in France (cat 2): Multiply b.effort1*array #1 
  # effort.dummy[,,,1] =0  effort.dummy[,,,2] =1 ==> Spain (cat 3): Multiply b.effort2*array #2
  
  tmp <-tmp2 <- tmp3 <- effort
  
  # Dummy variable 2 visits in France (only de 2 appear as 1)
  tmp2[tmp2[] %in% c(1,3)] <- 0
  tmp2[tmp2[] %in% c(2)] <- 1
  effort.dummy[,,1] <- tmp2
  
  # Dummy variable trap in Spain (only de 3 appear as 1)
  tmp3[tmp3[] %in% c(1,2)] <- 0
  tmp3[tmp3[] %in% c(3)] <- 1
  effort.dummy[,,2] <- tmp3
  
  
  #----   2.9. TRAP COVARIATE ---- 
  
  trap <- as.numeric(as.factor(as.matrix(tdf[,4])))
  trap[trap[] %in% c(2)] <- 0
  
  
  #---- 3. DETECTION DATA   ---- 
  #----   3.1 MAKE Y   ---- 
  
  n <- length(unique(edf$ind))
  
  
  # Now Y needs to be filled by occasion
  
  Y <- array(0, c(n, nrow(tdf), K))
  rownames(Y) <- unique(edf$ind)
  xx <- edf[,c(1,2,3,4)]
  
  for (k in 1:K){
    xxtk <- xx[which(xx$occ == k), ]
    
    if(nrow(xxtk) == 0) next # If there were no detections in that occasion
    
    for (obs in 1:nrow(xxtk)) {
      Y[xxtk[obs, 2], xxtk[obs, 4], xxtk[obs, 3]] <- 1 # ASP: 1 because it can only be detected once per trap and occasion
    }
  }
  
  
  
  #----   3.1.1 BEHAVIOURAL RESPONSE COVARIATE FROM Y   ---- 
  
  dim(Y)
  
  Ys <- Y
  prevcap <- array(0, dim = c(dim(Ys)[1], dim(Ys)[2], 
                              dim(Ys)[3]))
  first <- matrix(0, dim(Ys)[1], dim(Ys)[2])
  for (i in 1:dim(Ys)[1]) {
    for (j in 1:dim(Ys)[2]) {
      if (sum(Ys[i, j, ]) > 0) {
        first[i, j] <- min((1:(dim(Ys)[3]))[Ys[i, j, 
        ] > 0])
        prevcap[i, j, 1:first[i, j]] <- 0
        if (first[i, j] < dim(Ys)[3]) 
          prevcap[i, j, (first[i, j] + 1):(dim(Ys)[3])] <- 1
      }
    }
  }
  
  dim(prevcap)
  
  
  #----   3.2 AUGMENT Y   ---- 
  ##augment observed data to size M
  M <- 400 
  
  y.in <- array(0, c(M, nrow(tdf), K))
  y.in[1:n,,] <- Y
  
  prevcap.in <- array(0, c(M, nrow(tdf), K))
  prevcap.in[1:n,,] <- prevcap
  
  
  #----   3.3 USE SPARSE FORMAT FOR Y   ---- 
  #change to 'sparse' format - speeds up computation by reducing file size
  y.sparse <- getSparseY(y.in)
  
  ##extract pieces to be passed to Nimble
  detNums <- y.sparse$detNums # Nº of traps at which each individual was detected
  maxDetNums <- y.sparse$maxDetNums # Maximun Nº of traps at which any individual was detected
  detIndices <- y.sparse$detIndices # ID of the traps where they were detected
  y.sp <- y.sparse$y # Number of detections at each trap
  
  ##number of trials per trap - now always 1 but needs to be passed to Nimble
  # always 1, because we model each occasion separately
  ones <- rep(1, nrow(tdf))
  
  
  #---- 4. FIT NIMBLE MODEL    ---- 
  
  ############################################################################################
  ### running a model in parallel ############################################################
  
  ##source code to run model in parallel 
  setwd("D:/MargSalas/Scripts_MS/Stats/Nimble")
  #setwd("~/data/data/Scripts_MS/Stats/Nimble")
  source("Parallel Nimble function FOR aNA2_model2.1.r")  
  
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
    K = K,
    effort = effort.dummy,
    trap = trap,
    prevcap = prevcap.in
  )
  
  ##compile data
  nimData <- list(habDens = X.d_sc,
                  y = y.sp,
                  detNums = detNums, 
                  lowerHabCoords = lowerHabCoords,
                  upperHabCoords = upperHabCoords,
                  habitatGrid = habitatGrid,
                  ones=ones,           # ASP: ones substitutes K=rep(K, n.max.traps)
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
    caps <- which(apply(Y[i,,],1,sum) > 0)
    if (length(caps)==0) next #fill in missing ACs with reasonable values later
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
                           b.effort1 = runif(1, 0.5,1),
                           b.effort2 = runif(1, 0.5,1),
                           b.trap = runif(1, 0.5,1),
                           b.bh = runif(1, 0.5,1),
                           z=z.in,
                           sxy=S.in,
                           beta.dens = runif(1,-0.1, 0.1))}
  
  ##source model code
  setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Model")
  #setwd("~/Scripts_MS/Oso/PopDyn/SCR/Run_Sim/2.SCRdenscov")
  source('2.2.SCRdenscov_p0Covs in Nimble.r')
  
  ##determine which parameters to monitor
  params <- c('N', 'psi', 'sigma', 'p0', 'beta.dens', 'b.effort1', 'b.effort2', 'b.trap', 'b.bh')
  
  ###### SAVE FOR RUNNING #####
  
  model = SCRhab_covsp0
  
  setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/1.SCRdensCov/Data_server_Cyril_covsp0")
  save(nimData, nimConstants, 
       inits, z.in, S.in, 
       params, model, file = paste("Data_Model1-1_2020_", covs[xxx], ".RData", sep = ""))
  
  
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
                            code = SCRhab_covsp0_sigSex,   ##your model code
                            inits = inits,                 ##your inits function
                            constants = nimConstants,      ##your list of constants
                            params = params,               ##your vector with params to monitor
                            niter = 150000,                  ##iterations per chain
                            nburnin = 100000,                ##burn-in
                            nthin = 10,                  ##thinning, main parameters
                            z.in = z.in,
                            S.in = S.in,
                            sex.in = sex.in
  )
  new <- Sys.time() - old
  
  ## ALWAYS close cluster when model is done
  stopCluster(this_cluster)
  
  ### output is a list 
  
  setwd(paste("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_covsp0/", covs[xxx], sep = ""))
  #setwd(paste("~/data/data/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_covsp0/", covs[xxx], sep = ""))
  save(chain_output, file = paste("Results_Model1-1_", covs[xxx], ".RData", sep = ""))
  
  
  
  
  #### OPTION 2: NO PARALLEL (TO TRY INITIAL VALUES AND SEE IF MODEL WORKS) ####
  model <- nimbleModel(SCRhab_covsp0,
                       constants = nimConstants, 
                       data=nimData,
                       inits = inits(),
                       check = FALSE)
  model$initializeInfo()
  model$calculate()
  
  
}