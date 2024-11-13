


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


##detection parameters
#sigma=movement
sigma <- 1
#p0=baseline detection
p0 <- 0.2

##trap array
X<-as.matrix(expand.grid(seq(-3,3,1), seq(-3,4,1)))
colnames(X)<-c('x', 'y')
J<-dim(X)[1]

##state space coordinates, max trap coords (below) + 3*sigma
xmin<-min(X[,1])-3*sigma
ymin<-min(X[,2])-3*sigma
xmax<-max(X[,1])+3*sigma
ymax<-max(X[,2])+3*sigma

##for discrete state space, create 144 grid cells (center coordinates)
gx<-rep(seq(xmin+0.5, xmax-0.5,1), 13)
gy<-rep(seq(ymax-0.5,ymin+0.5, -1), each=12)
##note to self: this is order by row (all x for a given y) and starts with 
##lowest y (ie, bottom left corner of grid)
G<-cbind(gx, gy)
colnames(G)<-c('x','y')

##remove some as non-habitat
rem<-which(G[,2]>2 & G[,1]<(-4))
G<-G[-rem,]

###scale X and G so that bottom left corner of state space is origin (0,0)
sc.coord<-scaleCoordsToHabitatGrid(X, G)

##this returns S also in row order (all x for a given y) but starting top left 
## corner 
G.sc<-sc.coord$coordsHabitatGridCenterScaled
X.sc<-sc.coord$coordsDataScaled

###get cell coordinates for G.sc
windowCoords<-getWindowCoords(G.sc)
habitatGrid<-windowCoords$habitatGrid


##determine which traps are within some threshold distance of each habitat cell 

##binary mask that determines which cells in S are suitable 
habitatMask <-matrix(1, nrow(habitatGrid), ncol(habitatGrid))
habitatMask[habitatGrid==0]<-0
localTraps<-getLocalTraps(habitatMask, X.sc, resizeFactor = 1, dmax = 5*sigma)

localTrapsIndex<-localTraps$localTrapsIndices
localTrapsNum<-localTraps$numLocalTraps
habitatGridDet<-localTraps$habitatGrid
MaxLocalTraps<-localTraps$numLocalTrapsMax

##name and structure for Nimble model
numHabWindows<-sum(habitatGrid !=0) #number of cells in S
numGridRows<-nrow(localTraps$habitatGrid)
numGridCols<-ncol(localTraps$habitatGrid)
  
#the following two are still passed to function but no longer used by it
lowerHabCoords<-windowCoords$lowerHabCoords
upperHabCoords<-windowCoords$upperHabCoords


################################################################################
## generate data in 'real' space

##calculate distances between all cells
dg<-e2dist(G, G)

##generate spatially correlated covariate
X.d<-spcov2(dg)

##effect of covariate on density
beta.d<-1

## calculate cell probabilities
p.cell<-exp(beta.d*X.d)/sum(exp(beta.d*X.d))


################################################################################
### data generation ############################################################

##true abundance within state space
N <- 50 #target population size each year
M <- 150 #total potentially ever alive

Tt <- 5 #number of years
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
  ## ASP: Place each individual AC in each of the 146 grids-categories according to p.cell

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
    ## Detection function: The probability of detection in each cells decreases 
    ## with the distance to the AC
    ## sigD determines how much it moves in the space across years
  pcomb <- pd*p.cell
    ## ASP: Probability of an individual being detected in a cell is a combination of
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
obs <- array(NA, c(n.all, J, Tt))

for (t in 1:Tt){
  D <- e2dist(cbind(Sx[,t], Sy[,t]), X)
    ## ASP: Distance from each AC to each traps
for (i in 1:n.all){
  if(z.all[i,t]==0){obs[i,,t]<-0} else{
  p.eff <- p0*exp(-D[i,]^2/(2*sigma^2)) ## ASP: The p is a function of the distance to the traps
  obs[i,,t]<-rbinom(J, K, p.eff)}
}
}

##how many individuals detected?
n <- sum(apply(obs, 1, sum)>0)

##subset obs to include only individuals detected at least once
seen <- which(apply(obs, 1, sum)>0)
Y <- obs[seen,,] ## ASP: Only detected individuals (capture histories)

##augment observed data to size Maug
Maug <- 150 
y.in <- array(0, c(Maug, J, Tt)) # One matrix of nxJ every year (5 years)
y.in[1:n,,] <- Y


##change to 'sparse' format - speeds up computation by reducing file size
y.sparse <- getSparseY(y.in)

##extract pieces to be passed to Nimble
y.sp <- y.sparse$y
detIndices <- y.sparse$detIndices
detNums <- y.sparse$detNums
maxDetNums <- y.sparse$maxDetNums


################################################################################
#### fit Nimble model ##########################################################

##compile constants
nimConstants <- list(
  M=Maug,J=J, numHabWindows=numHabWindows, 
  numGridRows=numGridRows, numGridCols=numGridCols, 
  maxDetNums=maxDetNums, MaxLocalTraps=MaxLocalTraps,
  nobs=n, Nyr=Tt
)

##compile data
nimData <- list(habDens=X.d, y=y.sp,detNums=detNums, 
              lowerHabCoords=lowerHabCoords, upperHabCoords=upperHabCoords,
              habitatGrid=habitatGrid,K=rep(K,J),  X.sc=X.sc,
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

S.in<- array(NA, c(Maug, 2, Tt))

for ( i in 1:n){
  for (t in 1:Tt){
  caps <- which(Y[i, ,t] > 0) ## ASP: Get in which traps the ind was captured at year t
  
  if (length(caps)==0) next #fill in missing ACs with reasonable values later
  if (length(caps)==1){ ## ASP: If its only in 1 trap, a put the trap location as AC
    S.in[i,,t]<-X[caps,]
  }else{ 
  S.in[i,,t]<-apply(X[caps,],2,mean)}  ## ASP: If its > 1 trap average location
  }}


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
                       z=z.in,
                       sxy=S.in.sc$coordsDataScaled,
                       sigD=runif(1, 1.5, 2.5))}

##source model code
##I prefer working on code in a separate script but you can also have everything in
##one script and just execute the code
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

##remove NAs
inn<-colnames(samp[[1]])
remm<-pmatch(c("R[1]", "pc.gam[1]"), inn)
samp2<-lapply(samp, function(x)x[,-remm])

##NOTE: summary command is from MCMCvis package; that also has good plotting options
## summary table for everything in "params" vector
summ<-MCMCsummary(samp2)
MCMCtrace(samp)



