## -------------------------------------------------
##                openSCR + denscov + Age structure
##                    Our data 2017-2021
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
library(parallel)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")

#---- 1. LOAD THE DETECTION DATA ---- 

load("edf1721.RData")

# As the model is set, there can be only one capture per trap and occasion
# --> The number of trials is 7 (From may to November): So max number of captures per trap = 7
# To fix it in this edf, remove duplicates
edf <- edf[-which(duplicated(edf)), ]
t <- edf %>% group_by(ind,occ, session,trap) %>% # All need to be one
  summarise(n()) 
p <- t[which(t$session == 4 & t$ind == "Cannellito"),]

load("tdf2017.RData")
load("tdf2018.RData")
load("tdf2019.RData")
load("tdf2020.RData")
load("tdf2021.RData")

tdf_all <- rbind(tdf2017[,1:3], tdf2018[,1:3], tdf2019[,1:3],
                 tdf2020[,1:3], tdf2021[,1:3]) # Join to define state space
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
setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
#setwd("~/Data_server/Variables_hrscale")
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

dmax <- c(20,10,10,10,8) ## PROBLEM: First year the local evaluation is almost useless

localTraps <- localTrapsNum.l <- MaxLocalTraps.l <- list()
for (t in 1:Tt){
  localTraps[[t]] <- getLocalTraps(habitatMask, Xt.sc[[t]], resizeFactor = 1, dmax = dmax[[t]])
  localTrapsNum.l[[t]]  <- localTraps[[t]]$numLocalTraps
  MaxLocalTraps.l[[t]]  <- localTraps[[t]]$numLocalTrapsMax
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
  localTrapsIndex[,1:MaxLocalTraps[t],t] <- localTraps[[t]]$localTrapsIndices
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
M <- 250
nz <- M-n

y.in <- array(0, c(M, max(Jyear), Tt))
y.in[1:n,,] <- Y

#----   3.3 USE SPARSE FORMAT FOR Y   ---- 
#change to 'sparse' format - speeds up computation by reducing file size
y.sparse <- getSparseY(y.in)

##extract pieces to be passed to Nimble
detNums <- y.sparse$detNums # Nº of traps at which each individual was detected
maxDetNums <- y.sparse$maxDetNums # Maximun Nº of traps at which any individual was detected
detIndices <- y.sparse$detIndices # ID of the traps where they were detected
y.sp <- y.sparse$y # Number of detections at each trap


## ---- 3.4 INCLUDE AGE DATA ----

setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
#setwd("~/Data_server")

info <- readxl::read_excel("info_individuals_2021.xlsx")
#info <- readxl::read_xlsx("D:/MargSalas/Oso/Datos/Tablas_finales/2022/Info_individuals_2021.xlsx", sheet = 1)

info <- info[,c(4,8,10)]

# Create matrix with age data

# Matrix with exact ages
x <- matrix(NA, nrow = nrow(Y), ncol = Tt) # Matrix to fill
colnames(x) <- c(2017,2018,2019, 2020, 2021)
rownames(x) <- rownames(Y)

# While filling, 0 counts as an age (first year)

for (i in 1:nrow(x)){
  birth <- as.numeric(info$Year_birth[which(info$ID %in% rownames(x)[i])])
  if (birth >= 2017){
    x[i,which(colnames(x) %in% birth):Tt] <- 0:(Tt-which(colnames(x) %in% birth))
  } else { x[i,] <- (2017-birth):((2017-birth)+(Tt-1)) }
}

# Sum 1 so that the first year of life is age = 1
age <- x+1
#set known unborn individuals to 0
age[is.na(age)]<-0

###convert to 'raw' age categories: 
### age = 0: not yet recruited (category 1)
### age = 1: cub1 ( category 2), 
### age = 2: cub2 (category 3), 
### age = 3: subadult 1 (category 4)
### age = 4: subadult 2 (category 5)
### age >=5: adult (category 6)
### category 1: dead

age.cat<-age
age.cat[age==0]<-1
age.cat[age==1]<-2
age.cat[age==2]<-3
age.cat[age==3]<-4
age.cat[age==4]<-5
age.cat[age>=5]<-6

ageMatAug <- matrix(NA,M,Tt)
ageMatAug[1:n,] <- age.cat #Matrix with age categories augmented

# here, this is highest age category excluding dead and unrecruited, so 5
max.age <- 5

## ---- Known z ----

# To include the death information into the z, I create an age matrix including the deaths (value 0)
# The function of this matrix is ONLY to then create the last_alive vector
age.cat.z <- age.cat
for (i in 1:nrow(age.cat)){
  death <- info$Confirmed_death[which(info$ID %in% rownames(age.cat.z)[i])]
  if (death == "Alive") next
  if (death == 2021 | death == 2022) next #This has to be modified depending on the years of study
  death <- as.numeric(death) + 1 # We don't know when it died exactly, so we set it as death the year after
  # I could dig into this and get the specific date for some individuals
  age.cat.z[i,which(colnames(age.cat.z) %in% death):Tt] <- 0 
}

##some stuff for known z state below
#year of recruitment for observed individuals
r <- apply(age,1,function(x) min(which(x!=0)))
Y2d <- apply(Y, c(1,3), sum, na.rm=T)

# Combined vector to inform of the death or last capture (if death not available) of individuals
last_cap <- apply(Y2d,1,function(x) max(which(x>0)))  # last year with detection
last_alive <- apply(age.cat.z,1,function(x) min(which(x == 0))-1) # last occasion alive

last <- last_cap
last[which(last_alive != Inf)] <- last_alive[which(last_alive != Inf)] 

# Known z: state process is informed by both y, age data AND death recoveries

zdatAGE <- matrix(NA, M, Tt)
for(i in 1:n){
  #zdatNoAGE[i, first[i]:last[i]] <- 1         # when ignoring age data, we known an individuals is alive between first and last observation
  #if age is used to inform z
  zdatAGE[i, r[i]:last[i]] <- 1               # alive between known recruitment year and last capture or last occasion alive when known
  if(r[i]>1)  zdatAGE[i, 1:(r[i]-1)] <- 0     # Not entered prior to age==1
  if(last_alive[i] != Inf) zdatAGE[i, (last_alive[i]+1):Tt] <- 0 # Death after the last occasion alive (when available, so when last_alive !=Inf)
}


#---- 4. FIT NIMBLE MODEL    ---- 

############################################################################################
### running a model in parallel ############################################################

##source code to run model in parallel 
setwd("D:/MargSalas/Scripts_MS/Stats/Nimble")
source("Parallel Nimble function FOR aNA2.r") ##sorry, caps lock...

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
  nobs = n, 
  Nyr = Tt,
  max.age = max.age
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
                localTrapsNum = localTrapsNum,
                agePlusOne = ageMatAug[,1],       #known age in yr 1
                w = c(rep(1, n), rep(NA, nz)),    #membership in superpopulation
                u = zdatAGE,                      #known alive states
                b = rep(1,K), a = rep(1, max.age) #prior params for recruitment, age distribution
)


#----   4.2 INITIAL VALUES  ---- 
#----     4.2.1 Z AND AGE COMPONENT   ---- 

piAGE.in <- c(0.23, 0.09, 0.04, 0.14, 0.5)  # estimates from bear model (we will see, our data does not contain that many cubs)

zstAGE <- zdatAGE # ASP: u (alive state = zstAGE)
for(i in 1:n){
  if(last[i]<Tt) zstAGE[i,(last[i]+1):Tt] <- 0 ## ASP: Start as dead after the last observation
}
zstAGE[(n+1):M,] <- 0  # start augmented individuals as not entered to prevent error messages.
zstAGE[!is.na(zdatAGE)] <- NA #set known z values to NA

#because setting all augmented guys to 0 did not work:
#randomly generate augmented individuals as part of superpop
w.in <- c(rep(NA, n), rbinom(nz, 1, 0.5))

##generate entrance occasion for all individuals
ent.occ.aug <- sample(1:Tt, nz, replace=TRUE)

## ASP: GEnerate age category of entry individuals
age.cat.in <- c(rep(NA, n), rep(1, nz)) #1=not yet recruited

#randomly assign starting age category in year 1 to all individuals 
#alive in year 1
age.cat.in[(n+1):M][ent.occ.aug==1] <- sample(2:(max.age+1),sum(ent.occ.aug==1) , 
                                              prob=piAGE.in,
                                              replace=TRUE)
## ASP: Only assign age category for year 1, because it is what is used in the model to construct the 
## probabilities

##adjust starting values for u accordingly
for(i in (n+1):M){
  zstAGE[i,ent.occ.aug[i-n]:Tt]<-1
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

S.in.sc_coords <- S.in.sc$coordsDataScaled

inits<-function(){list(beta=c(0.15,rep(0.85/(Tt-1), Tt-1)), 
                       sigma=runif(1,0.5, 1.5),
                       p.ad=runif(1,0,0.5),
                       p.sub=runif(1,0,0.5),
                       p.cub=runif(1,0,0.5),
                       phi.ad=runif(1,0.5,1),
                       phi.sub=runif(1,0.5,1),
                       phi.cub=runif(1,0.5,1),
                       beta.dens = runif(1,-0.1, 0.1), # Added
                       piAGE=piAGE.in,
                       u = zstAGE,
                       w=w.in,
                       agePlusOne=age.cat.in,
                       psi=0.7,
                       sxy=S.in.sc_coords,
                       sigD=runif(1, 1.5, 2.5))}


##source model code
##I prefer working on code in a separate script but you can also have everything in
##one script and just execute the code

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Model")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Model")
source('3.3.SCRopen_diftraps_Age in Nimble.r')

##determine which parameters to monitor
params<-c("p.ad", "p.sub","p.cub","phi.ad","phi.sub","phi.cub", 
          "beta", "psi", "piAGE", "Nsuper", "N", "B", "N.cub", "N.sub", "N.ad", 'sigma', 'beta.dens', 'sigD')


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
                          code = SCRhab.Open.diftraps.age,   ##your model code
                          inits = inits,                 ##your inits function
                          constants = nimConstants,      ##your list of constants
                          params = params,               ##your vector with params to monitor
                          niter = 150000,                  ##iterations per chain
                          nburnin = 100000,                ##burn-in
                          nthin = 10,                  ##thinning, main parameters
                          Tt = Tt,                     ##additional objects needed within inits
                          piAGE.in = piAGE.in,
                          zstAGE = zstAGE,
                          w.in = w.in,
                          age.cat.in = age.cat.in,
                          S.in.sc_coords = S.in.sc_coords 
)
new <- Sys.time() - old

## ALWAYS close cluster when model is done
stopCluster(this_cluster)

### output is a list 

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/3.openSCRdenscov_Age")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/3.openSCRdenscov_Age")
save(chain_output, file = "sampOpenSCR_diftraps_age_2021.RData")




#### OPTION 2: NO PARALLEL (TO TRY INITIAL VALUES AND SEE IF MODEL WORKS) ####

#(1) set up model

model <- nimbleModel(SCRhab.Open.diftraps.age, constants = nimConstants, 
                     data=nimData, inits=inits(), check = FALSE)
##ignore error message, only due to missing initial values at this stage

model$initializeInfo()
#THEN YOU SHOULD ALWAYS CHECK THAT THE MODEL IS ABLE TO CALCULATE A LIKELIHOOD (RETURN A VALUE) GIVEN THE INITIAL VALUES PROVIDED
model$calculate()#
# IF A -INF OF NA IS RETURNED YOU CAN CHECK WHERE THE PROBLEM COMES FROM WITH (AND THEN TRY TO FIX IT UNTIL THE -INF DISAPEARS):
#model$logProb_p0# WILL GIVE YOU LIKELIHOOD OF P0
model$logProb_p0# WILL GIVE YOU LIKELIHOOD OF P0
model$logProb_sxy# THE PROBLEM IS OFTEN WITH SXY OR Y
model$logProb_y

model$logProb_y[,,4] # Año 2020, Canellito was not working because local evaluation was too low!


#(2) Compile model in c++
#     In complex models, this step can take a while (as well as step 5)
#     Much longer than in JAGS, but the model typically runs much faster
cmodel <- compileNimble(model)       

# (3) Configure MCMC - on an uncompiled model - this step allows setting which quantities to monitor
#     Also, nimble allows two sets of monitors, these can be thinned at different rates
#     all of which is more important in complex models but not to start with
conf.mcmc <- configureMCMC(model, monitors = params, thin=10)

# (4) Build the MCMC sampler based on configurations
mcmc <- buildMCMC(conf.mcmc)

# (5) Compile sampler in c++ together with compiled model
cmcmc <- compileNimble(mcmc, project = cmodel, resetFunctions = TRUE)

# (6) Run (monitor time just for fun) [takes 20 seconds on my computer]
system.time(
  (samp <- runMCMC(cmcmc, niter = 10, nburnin = 5, nchains=3, inits = inits) )
)

#setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/2.openSCRdenscov")
setwd("~/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/3.openSCRdenscov_Age")

save(samp, file = "sampOpenSCR_diftraps_age2.RData")

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