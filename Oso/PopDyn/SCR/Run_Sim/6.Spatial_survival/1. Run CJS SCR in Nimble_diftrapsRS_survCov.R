rm(list = ls())

library(nimble)
library(MCMCvis)
library(nimbleSCR)
##NOTE: load two helper functions at bottom of script for script to run

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
###YOU HAD TRAPS IN NON-HABITAT
rem <- which(G[,2]>2 & G[,1]<(-6.5))
G <- G[-rem,]

###scale X and G so that bottom left corner of state space is origin (0,0)
sc.coord <- scaleCoordsToHabitatGrid(X, G)

##this returns S also in row order (all x for a given y) but starting top left 
## corner 
G.sc <- sc.coord$coordsHabitatGridCenterScaled
X.sc <- sc.coord$coordsDataScaled

###get cell coordinates for G.sc
windowCoords <- getWindowCoords(G.sc, plot.check = FALSE)
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


# Trap array for nimble

## OPTION 3: ARRAY WITH TRAP MATRIX PER YEAR, WORKS
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
  localTraps[[t]] <- getLocalObjects(habitatMask, Xt.sc[[t]], resizeFactor = 1, dmax = 5*sigma, plot.check = FALSE)
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

#c <- localTrapsIndex[,,1]


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

##generate spatially correlated covariates, one for density, one for survival
X.d <- spcov2(dg)
X.surv<-spcov2(dg)

##effect of covariate on density and survival
beta.d <- 1
beta.surv <- 0.5

## calculate cell probabilities
p.cell <- exp(beta.d*X.d)/sum(exp(beta.d*X.d))


################################################################################
### data generation ############################################################

##true abundance within state space
N <- 50 #target population size each year
M <- 200 #total potentially ever alive

mu.phi <- rep(1, Tt) #baseline survival, logit scale
sigD <- 1 # dispersal sigma

##recruitment
##ASP: Probability of entering the population (be born or immigrate)
gamma <- NULL
gamma[1] <- N/M 

## ASP: Number of individuals that were not recruited (or not recorded as survivors) 
##  / Number of individuals available for recruitment from the augmented ones
## ASP: The first one is Pop size/Augmented

###demographic model

z <- r <- al <-age <- matrix(0, nrow = M, ncol = Tt)
r[,1] <- rbinom(M, 1, gamma[1]) ## ASP: Recruited/non-recruited, depends on prob.recruitment gamma 
z[,1] <- r[,1]  ## ASP: Alive state

## ASP: Individuals enter the model when they are first detected, so all individuals enter in year 1
Sx<-Sy<-s.g<-matrix(NA,M,Tt) #activity centers
#firs AC is random
s.multi <- rmultinom(M, 1, p.cell) 
s.g[,1] <- apply(s.multi,2, function(x){which(x==1)})
for (i in 1:M){
  Sx[i,1]<-G[s.g[i,1],1] ## ASP: Coordinates of the habitat window where the AC is placed
  Sy[i,1]<-G[s.g[i,1],2]       # Put into the first year when the individual was detected 
}    

## ASP: To integrate spatial covariate: 
# integrate the activity center generation into the simulation of the alive states, because we need
# to know where was the AC the previous year
for (t in 2:Tt){
  # survival
  phi.eff<-plogis(mu.phi[t]+beta.surv*X.surv[s.g[,t-1]])
  surv <- rbinom(M, 1, z[,t-1] * phi.eff) 
  
  ## ASP: For each individual, whether it survives (surv = 1) in t=2 depends on 
  # whether it was alive in t=1 and the survival probability in t=2. 
  # This survival probability depends on the habitat covariates of the cells of the
  # activity centers the previous year
  
  # recruitment
  al[,t] <- apply(matrix(z[,1:(t-1)], nrow = M, byrow = FALSE), 1, sum) > 0
  idx <- 1 - as.numeric(al[,t]) 
  gamma[t] <- (N - sum(surv)) / sum(idx)
  if (gamma[t] < 0) gamma[t] <- 0
  r[,t] <- rbinom(M, idx, gamma[t]) 
  z[,t] <- surv + r[,t]

  # activity center movement
  for (i in 1:M){
  if (z[i,t-1]==1 & z[i,t]==1) { # if not alive last year or not alive this year, random
    Dd <- e2dist(cbind(Sx[i,t-1], Sy[i,t-1]), G) 
    ## ASP: Distance from where it was in the previous time step to each grid center
    pd <- exp(-Dd^2/(2*sigD^2))
    ## sigD determines how much it moves in the space across years
    pcomb <- pd*p.cell
    pcomb <- pcomb/sum(pcomb)
    ssg <- rmultinom(1,1,pcomb)
    s.g[i,t]<-apply(ssg,2, function(x){which(x==1)})
    Sx[i,t]<-G[s.g[i,t],1]
    Sy[i,t]<- G[s.g[i,t],2] 
  }else{
    
    s.multi <- rmultinom(1, 1, p.cell) 
    s.g[i,t] <- apply(s.multi,2, function(x){which(x==1)})
    
    Sx[i,t]<-G[s.g[i,t],1] ## ASP: Coordinates of the habitat window where the AC is placed
    Sy[i,t]<-G[s.g[i,t],2]       # Put into the first year when the individual was detected 
}
}

}
n.all <- sum(apply(z,1,sum)>0) ##
z.all <- z[apply(z,1,sum)>0,] ## 
is.in<-which(apply(z,1,sum)>0)

###if this throws an error, M may have to be increased

##generate detection data

K <-5 #number of sampling occasions

## Because each year there is a different number of traps, the column dimension of the array is the maximum number of traps
## the rest is filled as NA 

n.max.traps <- max(Jyear)
obs <- array(NA, c(n.all, n.max.traps, Tt))

for (t in 1:Tt){
  D <- e2dist(cbind(Sx[is.in,t], Sy[is.in,t]), Xt[[t]])
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
##for CJS also exclude individuals seen first in occasion 5
##get first year an individual was detected
y2d<-apply(obs,c(1,3),sum, na.rm=T )
first<-apply(y2d,1, function(x){min(which(x>0))} )
##Inf means never detected

seen <- which(apply(obs, 1, sum, na.rm = TRUE)>0)
early<-which(first %in% (1:4)) # to exclude individuals seen first in occasion 5
keep<-which( (1:M) %in% seen & (1:M) %in% early)
Y <- obs[keep,,] ## ASP: Only detected individuals (capture histories)
n.obs<-dim(Y)[1]
##this now has NAs - does that matter?? Should prob all be 0,
##they won't enter the model anyway, tehy refer to inactive traps



##no augmentation for CJS
##change to 'sparse' format - speeds up computation by reducing file size
y.sparse <- getSparseY(Y) 

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
  J=Jyear, numHabWindows=numHabWindows, 
  numGridRows=numGridRows, numGridCols=numGridCols, 
  maxDetNums=maxDetNums, MaxLocalTraps=MaxLocalTraps,
  nobs=n.obs, Nyr=Tt,
  first=first[keep],
  first.po=first[keep]+1 
  #whichtrap.year = whichtrap.year, Jyear.index = Jyear.index # To map the year specific trap arrays
)

##compile data
nimData <- list(habDens=X.d, habSurv=X.surv, 
                y=y.sp,detNums=detNums, 
                lowerHabCoords=lowerHabCoords, upperHabCoords=upperHabCoords,
                habitatGrid=habitatGrid, K=rep(K, n.max.traps),  X.sc = Xt.sc.array, # ASP: Vector of k for bern trials with length = max number of traps
                habitatGridDet=habitatGridDet,detIndices=detIndices,
                detNums=detNums, localTrapsIndex=localTrapsIndex, 
                localTrapsNum=localTrapsNum
)

##set up initial values
##in real data, get yr of first detection
f.in <- first[keep]
z.in <- matrix(NA, n.obs, Tt)
for(i in 1:n.obs){
  z.in[i,(f.in[i]+1):Tt]<-1
}  ## ASP: From the first year detected fill as alive


##because of local evaluation of possible detectors, activity center initial 
##values have to be specified, eg average capture location in 'model' space
##CORRECTION: Use random capture location to avoid ACs in non-habitat cells

S.in<- array(NA, c(n.obs, 2, Tt))

for ( i in 1:n.obs){
  for (t in 1:Tt){
    caps <- which(Y[i, ,t] > 0) ## ASP: Get in which traps the ind was captured at year t
    
    if (length(caps)==0) next #fill in missing ACs with reasonable values later
    if (length(caps)==1){ ## ASP: If its only in 1 trap, a put the trap location as AC
      S.in[i,,t]<-Xt[[t]][caps,]
    }else{ 
      #ran.cap<-sample(caps,1)
      S.in[i,,t]<-apply(Xt[[t]][caps,],2,mean)}  ## ASP: If its > 1 trap average location
  }
  }


##fill in missing ACs as average of 'observed' ACs in nearby time step
##CORRECTION: use random nearby capture occasion

for (i in 1:n.obs){
  #which ACs unobserved
  nac <- which(is.na(S.in[i,1,]))
  wac <- (1:Tt)[-nac] ## ASP: Only the years that are observed
  #for those, use mean observed
  
  if(length(wac)==1) {
    S.in[i,1,nac]<-S.in[i,1,wac] ## ASP: Use the mean of the observed to fill unobserved
    S.in[i,2,nac]<-S.in[i,2,wac]   
  } else {
    ran.ac<-sample(wac,length(nac), replace=TRUE)
    S.in[i,1,nac]<-S.in[i,1,ran.ac] ## ASP: Use the mean of the observed to fill unobserved
    S.in[i,2,nac]<-S.in[i,2,ran.ac] 
  }

}

colnames(S.in) <- c('x', 'y')
S.in.sc <- scaleCoordsToHabitatGrid(S.in, G)

inits<-function(){list(sigma=runif(1,0.5, 1.5),
                       p.ad=runif(1,0,1),
                       mu.phi=runif(1,-0.2,0.2),
                       beta.cov=runif(1,-0.2,0.2),
                       z=z.in,
                       sxy=S.in.sc$coordsDataScaled,
                       sigD=runif(1, 1.5, 2.5),
                       beta.dens=runif(1,-0.5, 0.5))}

##source model code
##I prefer working on code in a separate script but you can also have everything in
##one script and just execute the code
setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Model")

source('7.CJS_SCR in Nimble_diftraps_SpatialCov.R')

##determine which parameters to monitor
params<-c('sigma', 'p.ad', 'mu.phi', 'beta.dens', 'sigD','beta.cov')

#(1) set up model

model <- nimbleModel(CJS.SCRhab.Open.diftraps, constants = nimConstants, 
                     data=nimData, inits=inits(), check = FALSE)
model$initializeInfo
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
  (samp <- runMCMC(cmcmc, niter = 5000, nburnin = 2000, nchains=3, inits = inits) )
)

# ##remove NAs
# inn<-colnames(samp[[1]])
# remm<-pmatch(c("R[1]", "pc.gam[1]"), inn)
# samp2<-lapply(samp, function(x)x[,-remm])

##NOTE: summary command is from MCMCvis package; that also has good plotting options
## summary table for everything in "params" vector
summ<-MCMCsummary(samp)
MCMCtrace(samp)



