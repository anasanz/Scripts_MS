

rm(list = ls())

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
#### packages ##########################################################

library(nimble)
library(MCMCvis)
library(nimbleSCR)
##NOTE: load two helper functions at bottom of script for script to run

### generate some basic SCR data with density covariate


##detection parameters
#sigma=movement
sigma <- 1
#p0=baseline detection
p0 <- 0.2

##trap array
X <- as.matrix(expand.grid(seq(-3,3,1), seq(-3,3,1)))
colnames(X) <- c('x', 'y')
J <- dim(X)[1]

##state space coordinates, max trap coords (below) + 3*sigma
xmin <- min(X[,1])-3*sigma
ymin <- min(X[,2])-3*sigma
xmax <- max(X[,1])+3*sigma
ymax <- max(X[,2])+3*sigma

##for discrete state space, create 144 grid cells (center coordinates)
gx <- rep(seq(xmin+0.5, xmax-0.5,1), 12)
gy <- rep(seq(ymin+0.5, ymax-0.5,1), each=12)
##note to self: this is order by row (all x for a given y) and starts with 
##lowest y (ie, bottom left corner of grid)
G <- cbind(gx, gy)
colnames(G) <- c('x','y')

###scale X and G so that bottom left corner of state space is origin (0,0)
sc.coord <- scaleCoordsToHabitatGrid(X, G)

##this returns S also in row order (all x for a given y) but starting top left 
## corner 
G.sc <- sc.coord$coordsHabitatGridCenterScaled
X.sc <- sc.coord$coordsDataScaled

plot(G.sc)

###get cell coordinates for G.sc
windowCoords <- getWindowCoords(G.sc)
habitatGrid <- windowCoords$habitatGrid
##coordinates start in upper left corner, go by row 
##habitatGrid appears to be in a different order, as the number 1 is in the 
##bottom left - but that's correct as habitatGrid is not a 1:1 representation
##of space - see comment below
##Note: any cells in habitatGrid not habitat should be coded as 0!!


##determine which traps are within some threshold distance of each habitat grid 
## cell - speeds up computations by only evaluating traps at which an individual
## could have been caught at (rather that traps very far away)
## not important for this toy example but for large areas

##binary mask that determines which cells in S are suitable - in this toy 
##example, all are suitable
habitatMask <- matrix(1, 12, 12)
localTraps <- getLocalTraps(habitatMask, X.sc, resizeFactor = 1, dmax = 5*sigma)


##this returns traps for each habitat grid cell in column order
##this is to be used in the detection model function, dbinomLocal_normal
##as for the habitatGrid above, the returned objects appear not to match 
##(grid has value 1 in top left, local traps index starts bottom left)
##they do work together correctly as only one of them is a 1:1 representation
##of space and the other one is upside down so that coordinates can be matched
##with rows and columns in the other

##IMPORTANT: even though in this example the detection and the habitat grid 
##           have the same dimensions, they cannot be used interchangeably
##           because one goes by row and one by column, so use objects in 
##           'localTraps' for detection model and objects in 'windowCoords'
##           in activity center/density model

localTrapsIndex <- localTraps$localTrapsIndices
localTrapsNum <- localTraps$numLocalTraps
habitatGridDet <- localTraps$habitatGrid
MaxLocalTraps <- localTraps$numLocalTrapsMax

##name and structure for Nimble model
numHabWindows <- prod(dim(localTraps$habitatGrid)) #number of cells in S
numGridRows <- nrow(localTraps$habitatGrid)
numGridCols <- ncol(localTraps$habitatGrid)
  
#the following two are still passed to function but no longer used by it
lowerHabCoords <- windowCoords$lowerHabCoords
upperHabCoords <- windowCoords$upperHabCoords

################################################################################
### data generation ############################################################

##effect of covariate on density
beta.d <- 1

##true abundance within state space
N <- 50

##calculate distances between all cells
dg <- e2dist(G.sc, G.sc)

##generate spatially correlated covariate
X.d <- spcov2(dg)

##Order has to match that represented in habitatGrid - but it's a little 
## complicated:
## habitatGrid is NOT a representation of space: as you look at it, it's NOT
## a proxy for a map. It's upside down. Because coordinates are used to index
## rows and columns in habitatGrid. Low y coordinates in true space are at the
## bottom, but in habitatGrid they correspond to low row numbers and are thus
## at the top!!
##So: the correct order of the habitat covariate is starting in the top left and 
##    going down by row. If you land in the first cell (top left), you get a 
##    y coordinate of (here) 11.5, and an x coordinate of 0.5
##    In habitatGrid that puts you in row 12, column 1 (bottom left) which
##    correctly references back to the first habWindow/value of the covariate

## calculate cell probabilities
p.cell <- exp(beta.d*X.d)/sum(exp(beta.d*X.d))

##draw activity center cells from multinomial
s.multi <- rmultinom(N, 1, p.cell)
s.g <- apply(s.multi,2, function(x){which(x==1)}) # Where are the activity centers??

##match cell ID to coordinates
Sx<-G.sc[s.g,1]
Sy<-G.sc[s.g,2]

###############################################################################
##plot with image:
## Has nothing to do with model, just visualizes data

##image has x axis = rows and y axis = columns, with first column at the bottom

##looking at this matrix it represents space, ie, top corresponds to high y
## coordinates
Xd.spatial<-matrix(X.d, 12, 12, byrow=TRUE)

##I think this is the correct plot... 
image(t(Xd.spatial[12:1,]), 
      x=seq(0.5,11.5, 1), y=seq(0.5,11.5, 1))
##random activity centers, based on density model
#calculate probability of s falling in each cell
points(Sx, Sy, pch=19)
im <- t(Xd.spatial[12:1,])
################################################################################
# If I want to go from im (order covariate so that it looks good in space) to the matrix again:
f <- t(im) # First t 
d <- f[12:1,] # Then from 12:1
as.vector(t(d))

##calculate distance from activity centers to detectors
D <- e2dist(cbind(Sx, Sy), X.sc)

##generate detection data
K <- 5 #number of sampling occasions
obs <- matrix(NA, N, J)

for (i in 1:N){
  p.eff <- p0*exp(-D[i,]^2/(2*sigma^2))
  obs[i,] <- rbinom(J, K, p.eff)
}

##how many individuals detected?
n <- sum(apply(obs, 1, sum) > 0)

##subset obs to include only individuals detected at least once
seen <- which(apply(obs, 1, sum)>0)
Y <- obs[seen,]

##augment observed data to size M
M <- 80 
y.in <- rbind(Y, 
            matrix(0, nrow = M-n, ncol = J))

##change to 'sparse' format - speeds up computation by reducing file size
y.sparse <- getSparseY(y.in)

##extract pieces to be passed to Nimble

detNums <- y.sparse$detNums # Nº of traps at which each individual was detected
maxDetNums <- y.sparse$maxDetNums # Maximun Nº of traps at which any individual was detected
detIndices <- y.sparse$detIndices # ID of the traps where they were detected
y.sp <- y.sparse$y # Number of detections at each trap

################################################################################
#### fit Nimble model ##########################################################

##compile constants
nimConstants <- list(
  M = M, J = J, numHabWindows = numHabWindows, 
  numGridRows = numGridRows, numGridCols = numGridCols, 
  maxDetNums = maxDetNums, MaxLocalTraps = MaxLocalTraps
)

##compile data
nimData <- list(habDens = X.d, y = y.sp, detNums = detNums, 
              lowerHabCoords = lowerHabCoords, upperHabCoords = upperHabCoords,
              habitatGrid = habitatGrid, K = rep(K,J),  X.sc = X.sc,
              habitatGridDet = habitatGridDet, detIndices = detIndices,
              detNums = detNums, localTrapsIndex = localTrapsIndex, 
              localTrapsNum = localTrapsNum 
              )

##set up initial values
z.in <- c(rep(1, n), rep(0, M-n))

##because of local evaluation of possible detectors, activity center initial 
##values have to be specified, eg average capture location
S.in<- matrix(NA, 80, 2)
for ( i in 1:length(seen)){
  caps<-which(Y[i, ]>0)
  if (length(caps)==1){
    S.in[i,]<-X.sc[caps,]
  }else{
  S.in[i,]<-apply(X.sc[caps,],2,mean)}
}

##random ACs for individuals never observed
for(i in (length(seen)+1) : M){
  S.in[i,]<-c(runif(1,0, 12),runif(1,0, 12))
}


inits<-function(){list(psi=runif(1,0.4, 0.6), 
                       sigma=runif(1,0.5, 1.5),
                       p0=runif(1,0,1),
                       z=z.in,
                       sxy=S.in)}

##source model code
##I prefer working on code in a separate script but you can also have everything in
##one script and just execute the code

setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Oso/PopDyn/SCR/Sim/SCRdenscov")
source('SCR in Nimble.R')

##determine which parameters to monitor
params<-c('N', 'psi', 'sigma', 'p0', 'beta.dens')

#(1) set up model

model <- nimbleModel(SCRhab, constants = nimConstants, 
                     data=nimData, check = FALSE)
##ignore error message, only due to missing initial values at this stage

#(2) Compile model in c++
#     In complex models, this step can take a while (as well as step 5)
#     Much longer than in JAGS, but the model typically runs much faster
cmodel <- compileNimble(model)       

# (3) Configure MCMC - on an uncompiled model - this step allows setting which quantities to monitor
#     Also, nimble allows two sets of monitors, these can be thinned at different rates
#     all of which is more important in complex models but not to start with
conf.mcmc<-configureMCMC(model, monitors = params, thin=1)

# (4) Build the MCMC sampler based on configurations
mcmc <- buildMCMC(conf.mcmc)

# (5) Compile sampler in c++ together with compiled model
cmcmc <- compileNimble(mcmc, project = cmodel, resetFunctions = TRUE)

# (6) Run (monitor time just for fun) [takes 20 seconds on my computer]
system.time(
  (samp <- runMCMC(cmcmc, niter = 15000, nburnin = 5000, nchains=3, inits = inits) )
)

##NOTE: summary command is from MCMCvis package; that also has good plotting options
## summary table for everything in "params" vector
summ <- MCMCsummary(samp)
MCMCtrace(samp)



