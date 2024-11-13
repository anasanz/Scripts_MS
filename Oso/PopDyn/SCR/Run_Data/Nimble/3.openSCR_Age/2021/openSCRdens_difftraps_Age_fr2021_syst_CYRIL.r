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

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_prefinal_1721/Cyril")


#---- 1. LOAD THE DETECTION DATA ---- 

load("edf1721.RData")

# As the model is set, there can be only one capture per trap and occasion
# --> The number of trials is 7 (From may to November): So max number of captures per trap = 7
# To fix it in this edf, remove duplicates
edf <- edf[-which(duplicated(edf)), ]
t <- edf %>% group_by(ind,occ, session,trap) %>% # All need to be one
  summarise(n()) 
p <- t[which(t$session == 4 & t$ind == "Canellito"),]

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

edf <- edf[-which(edf$ind %in% c("Nere", "Goiat", "Canellito")), ]

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

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_prefinal_1721/Cyril")

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


# Get one localtraps per year

dmax <- c(20,10,10,9,8) ## PROBLEM: First year the local evaluation is almost useless

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

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_prefinal_1721/Cyril")

info <- readxl::read_excel("info_individuals_2021.xlsx")

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

#S.in[7,,4] <- as.numeric(Xt[[4]][288,])

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

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_prefinal_1721/Cyril")

save(piAGE.in, file = "piAGE.in.RData")
save(zstAGE, file = "zstAGE.RData")
save(w.in, file = "w.in.RData")
save(age.cat.in, file = "age.cat.in.RData")
save(S.in.sc_coords, file = "S.in.sc_coords.RData")


inits<-function(){list(beta=c(0.15,rep(0.85/(5-1), 5-1)), 
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


## MODEL 
# This model is a combination of the age-structured JS model (code_JS_AGEcatV3.1)
# and the open SCR model with different trap arrays/year (SCR in Nimble_diftraps)

SCRhab.Open.diftraps.age<-nimbleCode({
  
  ### PRIORS ###
  psi~ dbeta(1,1)           # data augmentation
  
  # recruitment prob at k 
  beta[1:Nyr] ~ ddirch(b[1:Nyr])
  
  eta[1] <- beta[1]
  for(k in 2:Nyr){
    eta[k] <- beta[k]/(1-sum(beta[1:(k-1)]))
  }
  
  # starting age distribution is not yet recruited (age = 0), or age if recruited
  piAGE[1:max.age] ~ ddirch(a[1:max.age])
  piAGEuncond[1:(max.age+1)] <- c( (1-eta[1] ), eta[1]*piAGE[1:max.age] )  
  
  ##movement parameters: detection model and between-year AC movement model
  sigma~dunif(0,5) #adjust to units of trap array and space use of species
  sigD~dunif(0,5)  #dispersal Kernel SD, adjust to units of trap array
  
  ##detection parameter - p0 (baseline detection probability), per age category
  p.ad ~ dbeta(1,1)            # detection per age class
  p.sub ~ dbeta(1,1) 
  p.cub ~ dbeta(1,1) 
  
  p0[1]<-0  #not recruited yet, placeholder to make indexing work
  p0[2]<-p.cub
  p0[3]<-p.cub
  p0[4]<-p.sub
  p0[5]<-p.sub
  p0[6]<-p.ad
  
  ##survival probability, per age category
  phi.ad ~ dbeta(1,1)          # survival adults
  phi.sub ~ dbeta(1,1)          # survival subadults
  phi.cub ~ dbeta(1,1)          # survival cubs
  
  phi[1]<-0  #not recruited yet, placeholder to make indexing work
  phi[2]<-phi.cub
  phi[3]<-phi.cub
  phi[4]<-phi.sub
  phi[5]<-phi.sub
  phi[6]<-phi.ad
  
  ##effect of habitat cov on density
  beta.dens ~ dnorm(0, 0.01)
  
  ##model for density surface - vectorized, plus some other stuff for the function
  ##for initial AC 
  mu1[1:numHabWindows] <- exp(beta.dens * habDens[1:numHabWindows])
  sumHabIntensity <- sum(mu1[1:numHabWindows])
  logHabIntensity[1:numHabWindows] <- log(mu1[1:numHabWindows])
  logSumHabIntensity <- log(sumHabIntensity)
  
  ##  FIRST YEAR
  ####activity centers, alive state, age yr 1
  for (i in 1:M){
    
    w[i] ~ dbern(psi) # part of superpopulation?
    
    # Age process
    agePlusOne[i] ~ dcat(piAGEuncond[1:(max.age+1)]) 
    age[i,1] <- agePlusOne[i]-1 #age 0 = not yet entered, 5 = adult
    age.cat[i,1]<-age[i,1] #age category, everything above 5 == 5
    
    # State process
    u[i,1] <- step(agePlusOne[i]-1.1)  # alive if age[i,1] >0
    z[i,1] <- u[i,1]*w[i]  
    
    # derived stuff
    avail[i,1] <- 1- u[i,1]            # still available for recruitment. 
    recruit[i,1] <- z[i,1]             # recruited at k
    
    
    #activity centers according to density surface
    sxy[i, 1:2,1] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],#not used; from getWindowCoords()
      upperCoords = upperHabCoords[1:numHabWindows, 1:2],#not used; getWindowCoords()
      logIntensities = logHabIntensity[1:numHabWindows], #from model
      logSumIntensity = logSumHabIntensity, #from model
      habitatGrid = habitatGrid[1:numGridRows,1:numGridCols],#from getWindowCoords()
      numGridRows =  numGridRows, #calculated from habitatGrid (data)
      numGridCols = numGridCols
    )
  } # end M loop for yr 1 ACs
  
  ## FOLLOWING YEARS
  ###activity centers, demographic model, t>1
  for (t in 2:Nyr){
    ##for observed inds, model movement of ACs between years
    for (i in 1:nobs){
      sxy[i, 1:2, t] ~ dbernppACmovement_normal(
        lowerCoords            = lowerHabCoords[1:numHabWindows, 1:2],#data getWindowCoords()
        upperCoords            = upperHabCoords[1:numHabWindows, 1:2],#data getWindowCoords()
        s                      = sxy[i, 1:2, t-1], #model parameter
        sd                     = sigD, #model parameter
        baseIntensities        = mu1[1:numHabWindows], 
        habitatGrid            = habitatGrid[1:numGridRows,1:numGridCols],
        numGridRows            = numGridRows,
        numGridCols            = numGridCols,
        numWindows             = numHabWindows
      )
    }
    
    #for never observed, always generate random ACs each year
    for (i in (nobs+1):M){
      sxy[i, 1:2, t] ~ dbernppAC(
        lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],#not used; from getWindowCoords()
        upperCoords = upperHabCoords[1:numHabWindows, 1:2],#not used; getWindowCoords()
        logIntensities = logHabIntensity[1:numHabWindows], #from model
        logSumIntensity = logSumHabIntensity, #from model
        habitatGrid = habitatGrid[1:numGridRows,1:numGridCols],#from getWindowCoords()
        numGridRows =  numGridRows, #calculated from habitatGrid (data)
        numGridCols = numGridCols
      )
    }
    
    ###demographic model
    for (i in 1:M){
      # State process
      u[i,t] ~ dbern( u[i,t-1]*phi[ (age.cat[i,t-1]+1) ] + avail[i,t-1]*eta[t] )   #
      z[i,t] <- u[i,t]*w[i]  
      
      # Age process
      age[i,t] <- age[i,t-1] + max(u[i,1:t]) # ages by one year after recruitment
      ##make sure age>5 get converted to age class 5 (adult)
      age.cat[i,t]<-min(age[i,t], max.age)
      
      # derived stuff
      avail[i,t] <- 1- max(u[i,1:t])       
      recruit[i,t] <- equals(z[i,t]-z[i,t-1],1) # recruited at k
      
    }
    
  }#end yr loop for ACs, demographic model
  
  #derived population level stuff
  for (t in 1:Nyr){
    N[t] <- sum(z[1:M,t])               # Annual abundance
    B[t] <- sum(recruit[1:M,t])         # Number of entries
    
    for(c in 1:max.age){ # Abundance per age class.
      N.age[c,t] <- sum(age.cat[1:M,t]==c & z[1:M,t]== 1 )
    } # c
    
    N.cub[t] <- sum(N.age[1:2,t])
    N.sub[t] <- sum(N.age[3:4,t])
    N.ad[t] <- N.age[5,t]
    
  } #t
  
  # In case it doesn't work the sum
  #for(c in 1:max.age){
  #  for (i in 1:M){
  #  isAgeAlive[i,t,c] <- (age.cat[i,t]==c) & (z[i,t]==1)  
  #  }
  #  N.age[c,t] <- sum(isAgeAlive[1:M,t,c])
  #}
  
  Nsuper <- sum(w[1:M])            # Superpopulation size
  
  ##################################################################################################################################
  ##detection model
  for (t in 1:Nyr){
    for(i in 1:M){
      #detection model; data are from getSparseY()$y
      #maxDetNums = getSparseY()$maxDetNums
      #Note that in this example, the detection and the habitat grid have the same dimensions
      #otherwise numHabWindows would need to be provided separately for the detection grid
      y[i,1:maxDetNums,t]~dbinomLocal_normal(detNums = detNums[i,t],#getSparseY()$detNums
                                             detIndices = detIndices[i,1:maxDetNums,t],#getSparseY()$detIndices; ASP: Links with trapID
                                             size = K[1:J[t]], ##number of trials per trap (data); ASP: different each year
                                             p0 = p0[age.cat[i,t]+1], #model parameter, age-specific
                                             sigma = sigma, #model parameter
                                             s = sxy[i,1:2,t], #model parameter
                                             trapCoords = X.sc[1:J[t],1:2,t], #trap coordinates (data); ASP: Year specific trap array
                                             localTrapsIndices = localTrapsIndex[1:numHabWindows,1:MaxLocalTraps[t],t], #from getLocalTraps()
                                             localTrapsNum = localTrapsNum[1:numHabWindows,t], #from getLocalTraps()
                                             resizeFactor = 1, #no resizing
                                             habitatGrid = habitatGridDet[1:numGridRows,1:numGridCols],#from getLocalTraps()
                                             indicator = z[i,t]) #model parameter
    }#end ind loop
    
  }
  
  
})



##determine which parameters to monitor
params<-c("p.ad", "p.sub","p.cub","phi.ad","phi.sub","phi.cub", 
          "beta", "psi", "piAGE", "Nsuper", "N", "B", "N.cub", "N.sub", "N.ad", 'sigma', 'beta.dens', 'sigD')


#### OPTION 2: NO PARALLEL (TO TRY INITIAL VALUES AND SEE IF MODEL WORKS) ####

# Save for Cyril:

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_prefinal_1721/Cyril")
save(SCRhab.Open.diftraps.age, file = "SCRhab.Open.diftraps.age.RData")
save(nimConstants, file = "nimConstants.RData")
save(nimData, file = "nimData.RData")
save(params, file = "params.RData")
save(inits, file = "inits.RData")

#(1) set up model

model <- nimbleModel(SCRhab.Open.diftraps.age, constants = nimConstants, 
                     data=nimData, inits=inits(), check = FALSE)
##ignore error message, only due to missing initial values at this stage

model$initializeInfo()
model$calculate()#


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
  (samp <- runMCMC(cmcmc, niter = 150000, nburnin = 100000, nchains=3, inits = inits) )
)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_prefinal_1721/Cyril")

save(samp, file = "sampOpenSCR_diftraps_age2.RData")


