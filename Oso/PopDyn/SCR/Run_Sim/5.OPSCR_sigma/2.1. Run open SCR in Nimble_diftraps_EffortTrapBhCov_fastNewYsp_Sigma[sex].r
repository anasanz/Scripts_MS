rm(list = ls())

library(nimble)
library(MCMCvis)
library(nimbleSCR)
##NOTE: load two helper functions at bottom of script for script to run

setwd("D:/MargSalas/Scripts_MS/Stats/Nimble")
source('dbinomLocal_normalBear.R')

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
sigma <- c(0.5,1)

#p0=baseline detection, intercept on logit scale
p0 <- -2.534601

#   # Covariate effects on p (b.effort1, b.effort2, b.trap)
trapBetas <- c(0.8,-0.4, 0.3)

# 1 beta for behavioral response (effect of being already captured in that trap)
b.bh <- 0.5

##trap array (ASP: All years)
X<-as.matrix(expand.grid(seq(-6,6,1), seq(-6,7,1))) 
colnames(X)<-c('x', 'y')
J<-dim(X)[1]

##state space coordinates, max trap coords (below) + 3*sigma
xmin<-min(X[,1])-3*sigma[2]
ymin<-min(X[,2])-3*sigma[2]
xmax<-max(X[,1])+3*sigma[2]
ymax<-max(X[,2])+3*sigma[2]

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
#Xt.sc.array[1:Jyear[1], 1:2, 1] # Example: to get traps from year 1

# Get one localtraps per year

localTraps <- localTrapsNum.l <- MaxLocalTraps.l <- list()
for (t in 1:Tt){
  localTraps[[t]] <- getLocalObjects(habitatMask, Xt.sc[[t]], resizeFactor = 1, dmax = 5*sigma[2], plot.check = FALSE)
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

K <-7 #number of sampling occasions, here, constant for all traps in all years

## COVARIATES IN DETECTION
# 1. effort covariate
## create an effort array (trap by occasion per year), leaving most at 1, but setting some to 2
## here, randomly set ca. 30% of active occasions to 2
## Add third level, set 20% of active occasions to 3 (so set level 2 to 40%, as some of these will be 3)
n.max.traps <- max(Jyear)

effort<-array(1, c(n.max.traps, K, Tt))
n2<-round(prod(dim(effort))*0.4, dig=0) #ASP: To get how many numbers we will set as effort 2 (more or less, some will be 3)
n3<-round(prod(dim(effort))*0.2, dig=0) #ASP: To get how many numbers we will set as effort 3
effort[sample(1:prod(dim(effort)),n2, replace=FALSE)] <- 2
effort[sample(1:prod(dim(effort)),n3, replace=FALSE)] <- 3

# Create dummy variable

effort.dummy <- array(1, c(n.max.traps, K, Tt, 2)) # 4th dimension includes 2 arrays: one per level (intercept doesnt count)

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

# 2. trap covariate (Type of trap: 0 -> Hair, 1 -> Both)

trap <- array(0, c(n.max.traps, K, Tt)) 
n4 <- round(prod(dim(trap))*0.3, dig=0) #ASP: 30% will be traps of type "both
trap[sample(1:prod(dim(trap)),n4, replace=FALSE)] <- 1

#---- TRAP + EFFORT AS ARRAY FOR FAST FUNCTION ----#

effortarray <- array(1, c(n.max.traps, K, Tt, 3))
effortarray[,,,1] <- effort.dummy[,,,1]
effortarray[,,,2] <- effort.dummy[,,,2]
effortarray[,,,3] <- trap

#3. Sex covariate (Female = 0; Male = 1)

sex0 <- rep(0,n.all)
sex0[sample(1:length(sex0), 50, replace = FALSE)] <- 1 # Sex ratio around 1:1

# 4. Behavioural response
# 1 if an individual i has been already captured in that year in any occasion in all traps

b <- array(0, c(n.all, n.max.traps, K, Tt)) # Default needs to be 0, otherwise NA enter in model
# when not observed that year, and gives NA logliklhood

## ARRAY FOR DETECTION DATA
## Because each year there is a different number of traps, the column dimension of the array is the maximum number of traps
## the rest is filled as NA 
##add 4th dimension, which is occasion 
obs <- array(NA, c(n.all, n.max.traps, K, Tt))
cap <- array(0, c(n.all, n.max.traps, K, Tt)) # FOr behavioral response (default is 0, not captured)

#x <- model$logProb_y[,1,,1] # Individuals 1:4, 6:8 in all occasions year 1

# All years same NA until id 54 (not augmented individuals). Year 5 all are NA


for (t in 1:Tt){
  D <- e2dist(cbind(Sx[,t], Sy[,t]), Xt[[t]])
  ## ASP: Distance from each AC to each traps
  
  for (i in 1:n.all){
    if(z.all[i,t]==0){obs[i,,,t] <- 0} else {
      
      ##calculate distance component of detection, which is constant over occasions
      p.d <- exp(-D[i,]^2/(2*sigma[sex0[i]+1]^2)) ## 
      
      for (j in 1:nrow(Xt[[t]])){
        
        for (k in 1:K){
          if(k==1){
            b[i,j,k,t] <- 0 # At k = 1 is always 0
            p <- plogis(p0+trapBetas[1]*effortarray[j,k,t,1] + trapBetas[2]*effortarray[j,k,t,2] + trapBetas[3]*effortarray[j,k,t,3]+
                          b.bh*b[i,j,k,t])
            
            p.eff <- p*p.d[j] 
            obs[i,j,k,t] <- rbinom(1, 1, p.eff)
            
            if(obs[i,j,k,t] > 0){ # Check if has been detected in that trap, occasion and year
              cap[i,j,k,t] <- 1                       
            } 
          } else { # if k>1
            b[i,j,k,t] <- ifelse(sum(cap[i,j,1:(k-1),t]) > 0, 1, 0) # If it has been already captured before that k we put a 1
            p <- plogis(p0+trapBetas[1]*effortarray[j,k,t,1] + trapBetas[2]*effortarray[j,k,t,2] + trapBetas[3]*effortarray[j,k,t,3]+
                          b.bh*b[i,j,k,t])
            
            p.eff <- p*p.d[j] 
            obs[i,j,k,t] <- rbinom(1, 1, p.eff)
            
            if(obs[i,j,k,t] > 0){ # Check if has been detected in that trap, occasion and year
              cap[i,j,k,t] <- 1 } 
          }} # Loop k
      } # Loop j
    }} # Loop i
} # Loop t


##how many individuals detected?
n <- sum(apply(obs, 1, sum, na.rm = TRUE) >0)

##subset obs to include only individuals detected at least once

seen <- which(apply(obs, 1, sum, na.rm = TRUE)>0)
Y <- obs[seen,,,] ## ASP: Only detected individuals (capture histories)

# Sex of seen individuals
sex1 <- sex0[seen]

# b is an individual covariate, should have the i dimensions of Y
b2 <- b[seen,,,]


##augment observed data to size Maug
Maug <- 160 
y.in <- array(0, c(Maug, n.max.traps,K, Tt)) # One array of nxJxK every year (5 years)
y.in[1:n,,,] <- Y

##augment sex with NA
sex <- c(sex1,rep(NA,length((n+1):Maug)))


##augment behavioral response (latent variable)
b.aug <- array(0, c(Maug, n.max.traps, K, Tt)) 
b.aug[1:n,,,] <- b2

##change to 'sparse' format - speeds up computation by reducing file size
y.sparse <- list()
for (t in 1:Tt){
  y.sparse[[t]] <- getSparseY(y.in[,,,t]) 
}

##extract pieces to be passed to Nimble 
## Cyril changed this part of the code, so that we can use the function dbinomLocal_normal that speeds up the
# model. Now the detIndices and detNums are also contained within y.sp (lengthYCombined).
# This formulation is also necesary for making predictions

max.max <- max(sapply(y.sparse, function(x)x$lengthYCombined)) # ASP: Maximun of the maxDetNums

y.sp <- array(NA, c(Maug, max.max, K, Tt))

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
ones <- rep(1, n.max.traps)

################################################################################
#### fit Nimble model ##########################################################

##compile constants
nimConstants <- list(
  M = Maug, 
  J = Jyear, 
  numHabWindows = numHabWindows, 
  numGridRows = numGridRows, 
  numGridCols = numGridCols, 
  #maxDetNums=maxDetNums, 
  lengthYCombined = lengthYCombined,
  MaxLocalTraps = MaxLocalTraps,
  nobs = n, 
  Nyr = Tt, 
  K = K , 
  effort = effortarray, 
  nTrapCovs = dim(effortarray)[4],
  prevcap = b.aug
)

##compile data
nimData <- list(habDens = X.d, 
                y = y.sp,
                #detNums = detNums, 
                lowerHabCoords = lowerHabCoords, 
                upperHabCoords = upperHabCoords,
                habitatGrid = habitatGrid, 
                ones = ones,  
                X.sc = Xt.sc.array, 
                habitatGridDet = habitatGridDet,
                #detIndices=detIndices,     
                localTrapsIndex = localTrapsIndex, 
                localTrapsNum = localTrapsNum, 
                sex = sex
)

##set up initial values
##in real data, get yr of first detection
f.in <- first[seen] 

z.in <- matrix(0, Maug, Tt)
for(i in 1:n){
  z.in[i,f.in[i]:Tt]<-1
}  ## ASP: From the first year detected fill as alive

# LIKELIHOOD PROBLEM: There seems to be a problem when they are all started as alive

##because of local evaluation of possible detectors, activity center initial 
##values have to be specified, eg average capture location in 'model' space
##CORRECTION: Use random capture location to avoid ACs in non-habitat cells

S.in<- array(NA, c(Maug, 2, Tt))

#sum captures within year over occasions
Yt<-apply(Y, c(1,2,4), sum)

for ( i in 1:n){
  for (t in 1:Tt){
    caps <- which(Yt[i, ,t] > 0) ## ASP: Get in which traps the ind was captured at year t
    
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

for (i in 1:n){
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

#for (i in 1:n){
#  #which ACs unobserved
#  nac <- which(is.na(S.in[i,1,]))
#  wac <- (1:Tt)[-nac] ## ASP: Only the years that are observed
#  #for those, use mean observed
#  S.in[i,1,nac]<-mean(S.in[i,1,wac]) ## ASP: Use the mean of the observed to fill unobserved
#  S.in[i,2,nac]<-mean(S.in[i,2,wac]) 
#}


##random ACs for individuals never observed
for(i in (n+1) : Maug){
  for (t in 1:Tt){
    ssg<-sample(1:length(X.d), 1)
    S.in[i,,t]<-G[ssg,]
  }
}

## initial values for sex
sex.in <- c(rep(NA,n),rep(0,length((n+1):Maug)))
sex.in[sample((n+1):Maug, 30, replace = FALSE)] <- 1 # Random sex for augmented individuals

##INITIAL VALUES

colnames(S.in) <- c('x', 'y')
S.in.sc <- scaleCoordsToHabitatGrid(S.in, G)

inits<-function(){list(gamma=c(0.5, rep(0.1, (Tt-1))), 
                       sigma=runif(2,0.5, 1.5),
                       p0=runif(1,0,0.1), # Value for p0 on the probability scale (0-1)
                       trapBetas = runif(3, 0.5,1),
                       #b.effort1=runif(1, 0.5,1),
                       #b.effort2=runif(1, 0.5,1),
                       #b.trap=runif(1, 0.5,1),
                       b.bh=runif(1, 0.5,1),
                       omega = runif(1, 0.5,1),
                       sex = sex.in,
                       phi=runif(1,0.5,1),
                       z=z.in,
                       sxy=S.in.sc$coordsDataScaled,
                       sigD=runif(1, 1.5, 2.5),
                       beta.dens=runif(1,-0.5, 0.5))}

##source model code
setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Model")
#setwd("~/data/data/Scripts_MS/Oso/PopDyn/SCR/Model")
source('6.1.SCRopen_diftraps_difeff_effortTrapBhCov_fast_Sigma[sex] in Nimble.r')

##determine which parameters to monitor
params<-c('N', 'gamma', 'sigma', 'p0', 'trapBetas', 'b.bh', 'omega', 'phi', 'beta.dens', 'sigD','R', 'pc.gam', 'Nsuper')

#(1) set up model

model <- nimbleModel(SCRhab.Open.diftraps.3d.effortTrapBhCov.sigSex.fast, constants = nimConstants, 
                     data=nimData, inits=inits(), check = FALSE)
##ignore error message, only due to missing initial values at this stage
model$calculate()

model$initializeInfo()



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


setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Sim/Results/5.OPSCR_sigma")
save(samp, file = "samp_2.1.RData")

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Sim/Results/5.OPSCR_sigma")
load("samp_2.1.RData")

library(secr)
invlogit(-2.534601) # p0 was on the logit scale, but you model it in the non-logit scale

library(MCMCvis)
##remove NAs
inn<-colnames(samp[[1]])
remm<-pmatch(c("R[1]", "pc.gam[1]"), inn)
samp2<-lapply(samp, function(x)x[,-remm])

##NOTE: summary command is from MCMCvis package; that also has good plotting options
## summary table for everything in "params" vector
summ<-MCMCsummary(samp2)
summ
#MCMCtrace(samp)