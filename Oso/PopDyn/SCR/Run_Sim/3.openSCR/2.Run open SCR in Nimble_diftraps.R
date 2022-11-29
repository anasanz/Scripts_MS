rm(list = ls())

library(nimble)
library(MCMCvis)
library(nimbleSCR)

################################################################################
#### helper functions ##########################################################

##distances function between two sets of locations
e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}


### function to generate correlated covariate 
spcov2 <- function(D, alpha=2, standardize=TRUE) {
  V <- exp(-D/alpha)
  cov1 <- t(chol(V)) %*% rnorm(nrow(D))
  if(standardize)
    z <- as.numeric((cov1 - mean(cov1)) / sd(cov1))
  else
    z <- as.numeric(cov1)
  return(z)
}

################################################################################
###############################################################################

### generate some basic SCR data with density covariate

## ASP: I generate a larger state space with larger trap array, as I am going to 
## subset a different number of traps each year

##detection parameters
#sigma=movement
sigma <- 1
#p0=baseline detection
p0 <- 0.2

##trap array (ASP: All years)
X<-as.matrix(expand.grid(seq(-6,6,1), seq(-6,7,1))) 
colnames(X)<-c('x', 'y')
J<-dim(X)[1]

##state space coordinates, max trap coords (below) + 3*sigma
xmin<-min(X[,1])-3*sigma
ymin<-min(X[,2])-3*sigma
xmax<-max(X[,1])+3*sigma
ymax<-max(X[,2])+3*sigma

##for discrete state space, create 342 grid cells (center coordinates)
gx <- rep(seq(xmin+0.5, xmax-0.5,1), 19)
gy <- rep(seq(ymax-0.5, ymin+0.5, -1), each=18)
##note to self: this is order by row (all x for a given y) and starts with 
##lowest y (ie, bottom left corner of grid)
G <- cbind(gx, gy)
colnames(G) <- c('x','y')

##remove some as non-habitat
### ASP: Make sure that there are no traps in non-habitat!!
rem <- which(G[,2]>2 & G[,1]<(-6.5))
G <- G[-rem,]

###scale X and G so that bottom left corner of state space is origin (0,0)
sc.coord <- scaleCoordsToHabitatGrid(X, G)

##this returns S also in row order (all x for a given y) but starting top left 
## corner 
G.sc <- sc.coord$coordsHabitatGridCenterScaled
X.sc <- sc.coord$coordsDataScaled

###get cell coordinates for G.sc
windowCoords <- getWindowCoords(G.sc)
habitatGrid <- windowCoords$habitatGrid


##determine which traps are within some threshold distance of each habitat cell 

##binary mask that determines which cells in S are suitable 
habitatMask <-matrix(1, nrow(habitatGrid), ncol(habitatGrid))
habitatMask[habitatGrid==0] <- 0

# Trap array 
# Subset the traps sampled every year from scaled vector X.sc
Tt <- 5 #number of years
Yrs <- seq(1:Tt)

# Set it up so that there are a few number of traps not sampled every year 
# (otherwise local evaluation does not work because it doesn't find traps, need a very high dmax param)
set.seed(2022)
J.year <- list() # Which traps from X.sc are sampled every year
  for (t in 1:Tt){
    trap.year <- rbinom(nrow(X.sc),1, 0.75)
    J.year[[t]] <- which(trap.year>0)
  }
Jyear <- unlist(lapply(J.year,length)) # Number of traps per year


# Format trap array for nimble: ARRAY WITH TRAP MATRIX PER YEAR

Xt <- Xt.sc <- list() 
for (t in 1:Tt){
  Xt.sc[[t]] <- X.sc[J.year[[t]],] # For the getLocalTraps function
  Xt[[t]] <- X[J.year[[t]],] # To simulate yearly detection data later on
}

Xt.sc.array <- array(NA, c(max(Jyear), 2, Tt))
for (t in 1:Tt){
  Xt.sc.array[1:Jyear[t],,t] <- Xt.sc[[t]]
}
Xt.sc.array[1:Jyear[1], 1:2, 1] # Example: to get traps from year 1

# Get one localtraps per year

localTraps <- localTrapsNum.l <- MaxLocalTraps.l <- list()
for (t in 1:Tt){
  localTraps[[t]] <- getLocalTraps(habitatMask, Xt.sc[[t]], resizeFactor = 1, dmax = 5*sigma)
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
#are for entire S so not affected by changing trap array
lowerHabCoords<-windowCoords$lowerHabCoords
upperHabCoords<-windowCoords$upperHabCoords


################################################################################
## generate data in 'real' space

##calculate distances between all cells
dg <- e2dist(G, G)

##generate spatially correlated covariate
X.d <- spcov2(dg)

##effect of covariate on density
beta.d <- 1

## calculate cell probabilities
p.cell <- exp(beta.d*X.d)/sum(exp(beta.d*X.d))


################################################################################
### data generation ############################################################

##true abundance within state space
N <- 50 #target population size each year
M <- 150 #total potentially ever alive

phi <- rep(0.8, Tt) #survival
sigD <- 1 # dispersal sigma

##recruitment
##ASP: Probability of entering the population (be born or immigrate)
gamma <- NULL
gamma[1] <- N/M 

## ASP: Number of individuals that were not recruited (or not recorded as survivors) 
##  / Number of individuals available for recruitment from the augmented ones
## ASP: The first one is Pop size/Augmented

###demographic model
## ASP: Important assumption: an individual cannot be recruited twice!
## It can't leave, and come back, the only possibility of leaving is dying

z <- r <- al <- matrix(0, nrow = M, ncol = Tt)
r[,1] <- rbinom(M, 1, gamma[1]) ## ASP: Recruited/non-recruited, depends on prob.recruitment gamma 
z[,1] <- r[,1]  ## ASP: Alive state

for (t in 2:Tt){
  # survival
  surv <- rbinom(M, 1, z[,t-1] * phi[t]) 
  ## ASP: For each individual, whether it survives (surv = 1) in t=2 depends on 
  # whether it was alive in t=1 and the survival probability in t=2
  
  # recruitment
  al[,t] <- apply(matrix(z[,1:(t-1)], nrow = M, byrow = FALSE), 1, sum) > 0
  ## ASP: Sum alive states previous years to know if it has ever been alive?
  ## 0 = FALSE (never alive or in the population) -> Available to be recruited
  idx <- 1 - as.numeric(al[,t]) 
  ## ASP: Contrary to al[,t] so that 1 means available to be recruited
  gamma[t] <- (N - sum(surv)) / sum(idx)
  ## ASP: Probability of recruitment at t=2 
  ## Number of individuals that were not recruited (or not recorded as survivors) 
  ##  / Number of individuals available for recruitment from the augmented ones
  ## In the case t=2 -> 5/91: 5 out of 91 are recruited, so very low probability
  if (gamma[t] < 0) gamma[t] <- 0
  r[,t] <- rbinom(M, idx, gamma[t]) 
  ## ASP: Recruited/non recruited in t=2
  ## Number of trials (idx) is 0 if is not available
  z[,t] <- surv + r[,t]
  ## ASP: Alive state in t=2 depends on whether it survived from the last year
  ## OR it was recruited (born or immigrated). 
  ## It can only be recruited if it is available (i.e. not in survived category),
  ## so it can never be >1
}
n.all <- sum(apply(z,1,sum)>0) ## ASP: Number of individuals ever detected (all years)
z.all <- z[apply(z,1,sum)>0,] ## ASP: Capture histories of individuals ever detected

Sx<-Sy<-matrix(NA,n.all,Tt)


##first activity center is random

##get first time alive
first <- apply(z.all,1,function(x)min(which(x==1)))

s.multi <- rmultinom(n.all, 1, p.cell) 
## ASP: Place each individual AC in each of the n grids-categories according to p.cell

s.g <- apply(s.multi,2, function(x){which(x==1)})
## ASP: In which habitat window is the AC?

for (i in 1:n.all){
  Sx[i,first[i]]<-G[s.g[i],1] ## ASP: Coordinates of the habitat window where the AC is placed
  Sy[i,first[i]]<-G[s.g[i],2]       # Put into the first year when the individual was detected 
}                                 # (because it's random, the rest of AC will depend on where it was at t-1)

#####generate activity centers for t>1

for (t in 2:Tt){
  for (i in 1:n.all){
    if (z.all[i,t]==0) next ## ASP: Only for alive individuals
    if (first[i]== t) next  ## ASP: Only if it's not the first occasion (already placed randomly)
    Dd <- e2dist(cbind(Sx[i,t-1], Sy[i,t-1]), G) 
    ## ASP: Distance from where it was in the previous time step to each grid center
    pd <- exp(-Dd^2/(2*sigD^2))
    ## ASP: Probability of detection in each cell according to how far is from the AC
    ## ??? is this really a probability of detection??
    ## Detection function: The probability of detection in each cells decreases 
    ## with the distance to the AC
    ## sigD determines how much it moves in the space across years
    pcomb <- pd*p.cell
    ## ASP: Probability of an individual being detected in a cell is a combination of (????being detected or actually present in year t???)
    ## 1) the probability of detection in that cell pd (depends on where it was the individual at t-1)
    ## 2) The density probability in that cell p.cell (related to habitat)
    pcomb <- pcomb/sum(pcomb)
    ## ASP: Convert to relative probability
    ssg <- rmultinom(1,1,pcomb)
    ## ASP: Place the AC in one of the cells according to combinated probabilities
    ## The probability that an individual is in a cell at year t depends on the pd
    ## (which depends on where it was at t-1) and p.cell (expected density and habitat)
    Sx[i,t]<-G[which(ssg==1),1]
    Sy[i,t]<- G[which(ssg==1),2] 
    ## ASP: Fill in the activity center coordinates at year t
  }
}


##generate detection data

K <-5 #number of sampling occasions

## Because each year there is a different number of traps, the column dimension of the array is the maximum number of traps
## the rest is filled as NA 

n.max.traps <- max(Jyear)
obs <- array(NA, c(n.all, n.max.traps, Tt))

for (t in 1:Tt){
  D <- e2dist(cbind(Sx[,t], Sy[,t]), Xt[[t]])
  ## ASP: Distance from each AC to each traps
  
  for (i in 1:n.all){
    if(z.all[i,t]==0){obs[i,,t] <- 0} else{
      
      p.eff <- p0*exp(-D[i,]^2/(2*sigma^2)) ## ASP: The p is a function of the distance to the traps
      
      obs[i,1:nrow(Xt[[t]]),t] <- rbinom(nrow(Xt[[t]]), K, p.eff)
      }
  }
}

##how many individuals detected?
n <- sum(apply(obs, 1, sum, na.rm = TRUE) >0)

##subset obs to include only individuals detected at least once

seen <- which(apply(obs, 1, sum, na.rm = TRUE)>0)
Y <- obs[seen,,] ## ASP: Only detected individuals (capture histories)

##this now has NAs - does that matter?? Should prob all be 0,
##they won't enter the model anyway


##augment observed data to size Maug
Maug <- 150 
y.in <- array(0, c(Maug, n.max.traps, Tt)) # One matrix of nxJ every year (5 years)
y.in[1:n,,] <- Y


##change to 'sparse' format - speeds up computation by reducing file size
y.sparse <- getSparseY(y.in) 

## ASP: I have checked the function and there is no problem of having NA in the traps not sampled a given year

##extract pieces to be passed to Nimble
y.sp <- y.sparse$y
detIndices <- y.sparse$detIndices
detNums <- y.sparse$detNums
maxDetNums <- y.sparse$maxDetNums


################################################################################
#### fit Nimble model ##########################################################

##compile constants
nimConstants <- list(
  M=Maug, J=Jyear, numHabWindows=numHabWindows, 
  numGridRows=numGridRows, numGridCols=numGridCols, 
  maxDetNums=maxDetNums, MaxLocalTraps=MaxLocalTraps,
  nobs=n, Nyr=Tt)

##compile data
nimData <- list(habDens=X.d, y=y.sp,detNums=detNums, 
                lowerHabCoords=lowerHabCoords, upperHabCoords=upperHabCoords,
                habitatGrid=habitatGrid, K=rep(K, n.max.traps),  X.sc = Xt.sc.array, # ASP: Vector of k for bern trials with length = max number of traps
                habitatGridDet=habitatGridDet,detIndices=detIndices,
                detNums=detNums, localTrapsIndex=localTrapsIndex, 
                localTrapsNum=localTrapsNum
)

##set up initial values
##in real data, get yr of first detection
f.in <- first[seen] 
z.in <- matrix(0, Maug, Tt)
for(i in 1:n){
  z.in[i,f.in[i]:Tt]<-1
}  ## ASP: From the first year detected fill as alive


##because of local evaluation of possible detectors, activity center initial 
##values have to be specified, eg average capture location in 'model' space
##RS CORRECTION: Use random capture location to avoid ACs in non-habitat cells

S.in<- array(NA, c(Maug, 2, Tt))

for ( i in 1:n){
  for (t in 1:Tt){
    caps <- which(Y[i, ,t] > 0) ## ASP: Get in which traps the ind was captured at year t
    
    if (length(caps)==0) next #fill in missing ACs with reasonable values later
    if (length(caps)==1){ ## ASP: If its only in 1 trap, a put the trap location as AC
      S.in[i,,t]<-Xt[[t]][caps,]
    }else{ 
      #ran.cap<-sample(caps,1)
      S.in[i,,t]<-apply(Xt[[t]][caps,],2,mean)}  ## ASP: If its > 1 trap average location
  }                                               # The average location should always be habitat (as all traps are located in habitat?)
  }


##fill in missing ACs as average of 'observed' ACs in nearby time step
##RS CORRECTION: use random nearby capture occasion

for (i in 1:n){
  #which ACs unobserved
  nac <- which(is.na(S.in[i,1,]))
  wac <- (1:Tt)[-nac] ## ASP: Only the years that are observed
  #for those, use mean observed
  
  if(length(wac)==1) {
    S.in[i,1,nac]<-S.in[i,1,wac] 
    S.in[i,2,nac]<-S.in[i,2,wac]   
  } else {
    ran.ac<-sample(wac,length(nac), replace=TRUE)
    S.in[i,1,nac]<-S.in[i,1,ran.ac] ## ASP: Sample as many random locations as missing values
    S.in[i,2,nac]<-S.in[i,2,ran.ac] ## (NOT average location, as this locates them out of the habitat grid!)
  }

}


##random ACs for individuals never observed
for(i in (n+1) : Maug){
  for (t in 1:Tt){
    ssg<-sample(1:length(X.d), 1)
    S.in[i,,t]<-G[ssg,]
  }
}

colnames(S.in) <- c('x', 'y')
S.in.sc <- scaleCoordsToHabitatGrid(S.in, G)

inits<-function(){list(gamma=c(0.5, rep(0.1, (Tt-1))), 
                       sigma=runif(1,0.5, 1.5),
                       p0=runif(1,0,1),
                       phi=runif(1,0.5,1),
                       beta.dens = runif(1,-0.1, 0.1), # Added
                       z=z.in,
                       sxy=S.in.sc$coordsDataScaled,
                       sigD=runif(1, 1.5, 2.5))}

##source model code
setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Sim/3.openSCR")
source('SCR in Nimble_diftraps.R')

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
  (samp <- runMCMC(cmcmc, niter = 20, nburnin = 10, nchains=1, inits = inits) )
)

##remove NAs
inn<-colnames(samp[[1]])
remm<-pmatch(c("R[1]", "pc.gam[1]"), inn)
samp2<-lapply(samp, function(x)x[,-remm])

##NOTE: summary command is from MCMCvis package; that also has good plotting options
## summary table for everything in "params" vector
summ<-MCMCsummary(samp2)
MCMCtrace(samp)



