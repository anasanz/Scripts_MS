## -------------------------------------------------
##                 Predict future population
##                      Data script 3.1
##                         Model 3.7
## ------------------------------------------------- 

## Our abundance estimates are based in subsetting our population to our trapping area (core buffer)
## To predict abundance spatially we will estimate abundance choosing as starting point 
# all individuals from state-space, and subset afterwards (as we do in the model).
#   1. Choose the individuals in all state space as starting point
#   2. Load estimated pcr from individuals: In whole state space: R/NºAdults and in core buffer only: R core/NºAdults core
#   3. Subset afterwards to the individuals located in the core in all study years


rm(list = ls())

## Load packages

library(nimble)
library(nimbleSCR)
library(basicMCMCplots)
library(R.utils)
library(abind)
library(coda)
library(rgdal)


source("D:/MargSalas/Scripts_MS/Functions/ProcessCodaOutput.R")
source("D:/MargSalas/Scripts_MS/Functions/plot.violins2.r")
source("D:/MargSalas/Scripts_MS/Functions/DoScale.r")

setwd("D:/MargSalas/Scripts_MS/Functions/Nimble")
source('dbinomLocal_normalBear_rbinom2.R')

# Load buffer core area

setwd("D:/MargSalas/Oso/Datos/GIS/Countries")
Xbuf2 <- readOGR("Buffer_8500_traps.shp")


# Load data and model results

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/3.openSCR_Age/Data_server")
load("Data_Model3-3.1_CYRIL_allparams.RData")

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL")
load("myResults_3-3.1_param.RData")
sampmat1 <- do.call(rbind, nimOutput)

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL")
load("myResults_3-3.1_sxy.RData")
sampmat2 <- do.call(rbind, nimOutputSXY)

sampmat <- cbind(sampmat1, sampmat2)

#save(sampmat, file = "sampmat_3.1.RData")

#setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL")
#load("sampmat_3.1.RData")

M.aug <- nimConstants$M
Tt <- nimConstants$Nyr

## ---- 1. Re-Build model with extra dimensions ----

### I re-build model structure by taking values of one iteration

## ---- 1.1. Create starting z's, age's and sxy's ----

# These will be based on model output (last year), new M and new time frame
t.new <- 5 # sim 5 addl years
M.new <- 800 # new augmentation limit

# Use estimates of z, age etc from last year of estimation model only

# Z

z.start <- matrix(NA, M.new, 1+t.new)

##fill z.start with estimates of z form fitted model; set addl individuals in years
##covered by model to 0 - known NOT to be alive
#first, make matrix of z estimates
z.which1 <- pmatch(paste('z[', 1:M.aug, ', ', Tt,']', sep=''), colnames(sampmat)) # ASP: index columns z for LAST YEAR (sampmat2 matrix)
z.which <- grep('z\\[', colnames(sampmat)) # ASP: index columns all z (sampmat2 matrix)
##ordered: all individuals for yr 1, then all for yr2, etc
z.est <- matrix(sampmat[1, z.which], M.aug, Tt) # ASP: All z from estimation model
z.start[1:M.aug, 1] <- sampmat[1,z.which1] # ASP: Starting z for projection: z of year 5 
z.start[(M.aug+1):M.new, 1] <- 0 # ASP: Augmented as not yet recruited

# Repeat this for: AGE

age.start <- matrix(NA, M.new, 1+t.new)

##NOTE: age.start starts at 0 and thus works for age[i,1] AND age.cat[i,1]
##      age in later years acn exceed 5, but in age.cat we always have a max of 5
##      but we need to add 1 for agePlusOne (the categorical variable)
age.which <- pmatch(paste('age.cat[', 1:M.aug, ', ', Tt,']', sep=''), colnames(sampmat)) # ASP: index columns AGE for LAST YEAR
age.est <- sampmat[1,age.which] # ASP: Ages last year
##set everyone not alive at all to all-0, ie, can be recruited in the NEW recruitment 
##model (they were never part of the superpopulation -> ASP: because of how the model is formulated z=u*w)
z.nosuper <- which(apply(z.est,1,sum)==0)
age.est[z.nosuper] <- 0
#age.est[age.est>5] <- 5 #set to max age category. ASP: NO NEED OF THIS HERE, already age.cat
age.start[1:M.aug, 1] <- age.est # ASP: Starting AGE for projection: AGE of year 5 
age.start[(M.aug+1):M.new, 1] <- 0 # ASP: Augmented as not yet recruited

# Repeat this for: SXY

sxy.start <- array(NA, c(M.new, 2, 1+t.new))

##make array of activity centers from model output
s.which <- grep('sxy', colnames(sampmat)) # ASP: index columns all sxy (sampmat matrix)
##all individuals X Yr1, then all inds Y Yr1; all inds X yr 2, then y yr 2 etc ## ASP: This is the order the sxy vector is placed (s.which)
indx <- (Tt*M.aug*2-(M.aug*2-1)):(Tt*M.aug*2) # ASP: Index the LAST year of sxy. *2 because there are the double of columns: x and y. (M.aug*2-1) corresponds to 1 year of data
sxy.start[1:M.aug,,1] <- matrix(sampmat[ 1, s.which[indx] ], M.aug, 2) # ASP: Starting SXY for projection: SXY of year 5 

# Repeat this for: SEX
sex.start <- rep(NA, M.new)
sex.which <- grep('sex', colnames(sampmat))
sex.start[1:M.aug] <- sampmat[1,sex.which] 

## ---- 1.2. Per capita recruitment ----

# Load results from function in script 0.PCR_all_core.r

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL")
#load("pcr_corebuf.RData")
#load("pcr_all.RData")
#load("pcr_range.RData")
load("pcr_corebuf_fem.RData")
load("pcr_all_fem.RData")

pcr1_corebuf <- apply(pcr_corebuf_fem[[3]],1,mean) # ASP: Mean per capita recruitment per iteration
pcr2_all <- apply(pcr_all_fem[[3]],1,mean) # ASP: Mean per capita recruitment per iteration


# Take pcr value of one iteration to build the model
pcr <- pcr2_all[1]
  
#B <- matrix(sampmat[1:nrow(sampmat), grep('B\\[', colnames(sampmat))], nrow(sampmat), Tt) 
#pcr_all[[1]] == B[,c(2:5)] # ASP: B is the same than R in all state space! (except column 1, which is abundance)

## ---- 1.3. Other parameters ----
# We don't include the parameters related with detection (not in projection model)

phi.ad <- sampmat[1,'phi.ad']
phi.cub <- sampmat[1,'phi.cub']
phi.sub <- sampmat[1,'phi.sub']
beta.dens <- sampmat[1,'beta.dens']
omega <- sampmat[1,'omega']
sigD <- sampmat[1,'sigD']

##calculate psi and unconditional age structure at t=1, which don't really matter 
## as they  don't apply to the projection years
##but need to be in the model for proper structure
##won't be changed in simulations
psi <- mean(sampmat[,'N[5]'])/M.new # ASP: This is the data augmentation parameter (N last year/AUG), not really understand but it doesn't matter?
##unconditional age structure at t=1 (old t=5, from single iteration)
pi.uncond <- table(age.start[,1])/M.new # ASP: Reminder -> it is unconditional because it adds the prob. of not being recruited for augmented individuals

##psi and pi.uncond don't match because one is based on posterior mean, 
##one based on single iteration, but again, does not matter

## ---- 1.4. Set un constants and build model ----

##set up constants
nimConstants <- list(
  M=(M.new), numHabWindows = nimConstants$numHabWindows, 
  numGridRows = nimConstants$numGridRows, numGridCols = nimConstants$numGridCols, 
  Nyr=(1+t.new),
  max.age = nimConstants$max.age
)

##set up data
##new piece: agePlusOne; also note - age and age.cat the same 
##because they only have values for yr 1

names(nimData)
nimData<-list(habDens = nimData$habDens, 
              lowerHabCoords = nimData$lowerHabCoords, upperHabCoords = nimData$upperHabCoords,
              habitatGrid = nimData$habitatGrid,
              z = z.start,
              sxy = sxy.start,
              agePlusOne=age.start[,1]+1,
              age = age.start,
              age.cat = age.start,
              u = z.start,
              pcr = pcr, ##transformed per capita recruitment
              phi.ad = phi.ad,
              phi.cub = phi.cub,
              phi.sub = phi.sub,
              beta.dens = beta.dens,
              sigD = sigD,
              psi = psi,
              pi.uncond = pi.uncond,
              sex = sex.start,
              omega = omega
              
)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Model")
source('3.7.SCRopen_PREDICT in Nimble RS.r')

# Create model using projection code (per capita recruitment, no observations)
model <- nimbleModel( code = SCRhab.Open.diftraps.age.effortTrapBhCov.sigsexage.PR.pcr.fem,
                      constants = nimConstants,
                      data = nimData,
                      #inits = nimInits,
                      check = F,       
                      calculate = T)  

nodesToSim <- model$getDependencies(c("sigD","phi.ad","phi.cub","phi.sub",
                                      "pcr", "beta.dens", 'pi.uncond', 'psi', 'omega'),
                                    self = F,
                                    downstream = T,
                                    returnScalarComponents = TRUE)
model$simulate(nodesToSim, includeData = FALSE)
N <- apply(model$z,2,sum)

model$sxy[,,2] # Here there are no NA in the AC
sxy.start[1:M.aug,,1] == model$sxy[1:M.aug,,1] # And values in year 1 are = as data

# Are there individuals available for recruitment?
apply(model$recruit,1,sum) 
sum(apply(model$recruit,1,sum) ) 

##compile for faster simulation
cmodelSims <- compileNimble(model)



## ---- 2. Run projection for several iterations and subset abundance ----

# WE FIND THE S AND Z NODES THAT SHOULDNT BE PREDICTED FOR. 
# added more age related nodes
samplerConfList <- model$getNodeNames()
zNodes <- samplerConfList[grep("z\\[",samplerConfList)]
sNodes <- samplerConfList[grep("sxy",samplerConfList)]
uNodes <- samplerConfList[grep("u\\[",samplerConfList)]
ageNodes <- samplerConfList[grep("age\\[",samplerConfList)]
age.cat.Nodes <- samplerConfList[grep("age.cat\\[",samplerConfList)]
agePO.Nodes <-samplerConfList[grep("agePlusOne",samplerConfList)]
sexNodes <- samplerConfList[grep("sex\\[",samplerConfList)]


length(values(cmodelSims, sNodes)) #800*2*6 (sNodes are 4800 because dim 2 is 1:2)
sxy.start[1:M.aug,,1] == model$sxy[1:M.aug,,1]
values(cmodelSims, sNodes)[1:600] == c(t(sxy.start[1:M.aug,,1]))

## ---- 2.1. PCR inside core buffer ----
## ------ 2.1.1. Projection ----

pcr.core <- pcr1_corebuf

##get some random iterations from posterior
itera <- sample(1:nrow(sampmat), 20)
nimData1 <- nimData
Nmat.core <- Rmat.core <- matrix(NA, length(itera), 1+t.new)

# To store as simlist
sxy.proj.core <- array(NA, c(length(itera), M.new, 2, t.new+1))
z.proj.core <- age.cat.proj.core <- array(NA, c(length(itera), M.new, t.new+1))
sex.proj.core <- matrix(NA, length(itera), M.new)

for(ite in 1:length(itera)){
  # WE UPDATE THE Z VALUES USING THE POSTERIORS PREDICTED Z FOR THE FIVE FIRST YEARS AND THEN 
  # USE NA FOR THE years to predict#
  
  # WE SET NA FOR Z FOR YEARS TO PREDICT - from posterior
  z.est <- matrix(sampmat[itera[ite],z.which], M.aug, Tt)
  z.start[1:M.aug, 1] <- sampmat[itera[ite],z.which1]
  z.start[(M.aug+1):M.new, 1]<-0
  
  nimData1$z <- z.start
  nimData1$u <- z.start
  
  # we set the values in the model - yr 1 aren't nodes in model (only data)
  values(cmodelSims, zNodes) <- nimData1$z#[,2:10] # Fill the values predicted by the model by new values where the extra years are NA
  values(cmodelSims, uNodes) <- nimData1$z#[,2:10] # Fill the values predicted by the model by new values where the extra years are NA
  
  # WE SET SXY
  sxy.start[1:M.aug,,1] <- matrix(sampmat[ itera[ite], s.which[indx] ], M.aug, 2)
  
  nimData1$sxy <- sxy.start
  #length(c(t(sxy.start[1:M.aug,,1])))
  
  # we set the values in the model
  #values(cmodelSims, sNodes) <- nimData1$sxy
  
  ## FIND s NODES FROM THE FIRST YEAR (SHOULD BE REPLACED)
  sNodes <- sNodes[grep( ", 1]",sNodes)]
  for(i in 1:length(sNodes)){
  values(cmodelSims, sNodes[i]) <-sxy.start[i,,1]
  }
  
  ###set age
  age.est <- sampmat[itera[ite],age.which]
  ##set everyone not alive at all to all-0, ie, can be recruited in the NEW recruitment 
  ##model (they were never part of the superpopulation)
  z.nosuper <- which(apply(z.est,1,sum)==0)
  age.est[z.nosuper] <- 0
  #age.est[age.est>5] <- 5 #set to max age category
  age.start[1:M.aug, 1] <- age.est
  age.start[(M.aug+1):M.new, 1] <- 0
  
  nimData1$age <- age.start 
  values(cmodelSims, ageNodes) <- age.start #[,2:(Tt+t.new)] 
  
  ##also replace age.cat and agePlusOne
  nimData1$agePlusOne <- age.start[,1]+1 
  values(cmodelSims, agePO.Nodes) <- age.start[,1]+1
  
  nimData1$age.cat <- age.start 
  values(cmodelSims, age.cat.Nodes) <- age.start
  
  ## Also replace sex, because for the augmented individuals it changes each iteration
  
  sex.start[1:M.aug] <- sampmat[itera[ite],sex.which]
  nimData1$sex <- sex.start
  values(cmodelSims, sexNodes) <- sex.start
  
  ##set pcr, phi, sigD, beta.dens
  ##not sure if needed to set in data and cmodelSims...
  ##simplified since we calculated pcr for all iterations above
  nimData1$pcr <- pcr.core[itera[ite]]
  values(cmodelSims,"pcr") <- nimData1$pcr
  
  nimData1$phi.ad<-sampmat[itera[ite],'phi.ad']
  values(cmodelSims,"phi.ad") <- nimData1$phi.ad
  nimData1$phi.cub<-sampmat[itera[ite],'phi.cub']
  values(cmodelSims,"phi.cub") <- nimData1$phi.cub  
  nimData1$phi.sub<-sampmat[itera[ite],'phi.sub']
  values(cmodelSims,"phi.sub") <- nimData1$phi.sub
  
  nimData1$sigD<-sampmat[itera[ite],'sigD']
  values(cmodelSims,"sigD") <- nimData1$sigD
  
  nimData1$beta.dens<-sampmat[itera[ite],'beta.dens']
  values(cmodelSims,"beta.dens") <- nimData1$beta.dens
  
  #now we simulate 
  cmodelSims$simulate(nodes = nodesToSim,#c(zNodes,sNodes,yNodes),
                      includeData = F)#---if TRUE: want to simulate new values also for nodes considered as data
  
  # Store results from iteration
  #t=2
  #cmodelSims$pcr * cmodelSims$N.ad[t - 1]
  #cmodelSims$R[t]/sum(cmodelSims$avail[, t - 1])
  #cmodelSims$u[i, t - 1] * cmodelSims$phi[(cmodelSims$age.cat[i, t -1] + 1)] +
  #  cmodelSims$avail[i, t - 1] * cmodelSims$gamma[t]
  
  
  Nmat.core[ite,] <- apply(cmodelSims$z,2,sum)
  Rmat.core[ite,] <- cmodelSims$R
  
  sxy.proj.core[ite,,,] <- c(cmodelSims$sxy)
  z.proj.core[ite,,] <- cmodelSims$z
  age.cat.proj.core[ite,,] <- cmodelSims$age.cat
  
  sex.proj.core[ite,] <- cmodelSims$sex
  
  #cmodelSims$sxy[,,2]
  
}

##plot trajectory
plot(1:(1+t.new), apply(Nmat.core,2,mean), type='l', ylim=range(Nmat.core, na.rm=TRUE))
for (ite in 1:20){
  points(1:(1+t.new), Nmat.core[ite,], type='l', col='lightgrey')
}
points(1:(1+t.new), apply(Nmat.core,2,mean), type='l')

# There were NA in z, but it was solved by increasing the augmented individuals 
# I think it was not enough

## ------ 2.1.2. Subset abundance ----
## Unscale sxy coordinates to subset only individuals within core buffer

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")
load("habcoord.RData")

dimnames(sxy.proj.core)[[3]] <- c('x','y')
sxy.proj.core.uns <- scaleCoordsToHabitatGrid(coordsData = sxy.proj.core,## this are your sxy
                                              coordsHabitatGridCenter = G,# this is your unscaled habitat (as you used when scaling the habitat/detector to the habitat. G?
                                              scaleToGrid = FALSE)$coordsDataScaled

# Matrix to store abundance in the buffer each iteration and year

NIn_coreBuf <- matrix(NA,nrow = dim(z.proj.core)[1], ncol=dim(z.proj.core)[3]) # nrow = iterations, ncol = year
sp.check <- sp.t.check <- list()


for(ite in 1:dim(z.proj.core)[1]){
  for(t in 1:dim(z.proj.core)[3]){
    
    which.alive <- which(z.proj.core[ite,,t]==1) # Select only the individuals alive (z=1)
    
    which.aliveSXY <- sxy.proj.core.uns[ite,which.alive,,t] # Retrieve the activity center for those individuals
    
    sp <- SpatialPoints(which.aliveSXY, proj4string=CRS(proj4string(Xbuf2))) # CONVERT SXY TO SPATIAL POINTS 
    
    which.In <- over(sp, Xbuf2) # Check which ones are in the buffer
    
    NIn_coreBuf[ite,t] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
    
    sp.t.check[[t]] <- sp # To check where points fall
    
  }
  
  sp.check[[ite]] <-  sp.t.check # To check where points fall
}

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL/Predictions/20iter")
save(Nmat.core, Rmat.core, sxy.proj.core, z.proj.core, age.cat.proj.core, sxy.proj.core.uns, NIn_coreBuf, file = "proj_pcr.core.RData")

#average number of individuals without the buffer for each year 
colMeans(NIn_coreBuf)
NIn_coreBuf[,1]#posterior distrib for the first year

par(mfrow = c(2,3))
for (t in 1:dim(z.proj.core)[3]){
  plot(density(NIn_coreBuf[,t]), main = t)
  abline(v = colMeans(NIn_coreBuf)[t], col = "blue")
}

# Sum of individuals alive in total each year (without buffer)
NALL <- apply(z.proj.core,c(1,3),function(x) sum(x==1))
colMeans(NALL)

#library(raster)
## Load distcore variable 
#setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
#distcore <- raster("logDistcore_hrbear.tif")

par(mfrow = c(1,1))
plot(distcore)
plot(Xbuf2, add = TRUE)
points(sp.check[[8]][[6]]) 


## ---- 2.2. PCR in whole state space ----
## ------ 2.2.1. Projection ----

pcr.all <- pcr2_all

##get some random iterations from posterior
itera <- sample(1:nrow(sampmat), 20)
nimData1 <- nimData
Nmat.all <- Rmat.all <- matrix(NA, length(itera), 1+t.new)

# To store as simlist
sxy.proj.all <- array(NA, c(length(itera), M.new, 2, t.new+1))
z.proj.all <- age.cat.proj.all <- array(NA, c(length(itera), M.new, t.new+1))
sex.proj.all <- matrix(NA, length(itera), M.new)


for(ite in 1:length(itera)){
  # WE UPDATE THE Z VALUES USING THE POSTERIORS PREDICTED Z FOR THE FIVE FIRST YEARS AND THEN 
  # USE NA FOR THE years to predict#
  
  # WE SET NA FOR Z FOR YEARS TO PREDICT - from posterior
  z.est <- matrix(sampmat[itera[ite],z.which], M.aug, Tt)
  z.start[1:M.aug, 1] <- sampmat[itera[ite],z.which1]
  z.start[(M.aug+1):M.new, 1]<-0
  
  nimData1$z <- z.start
  nimData1$u <- z.start
  
  # we set the values in the model - yr 1 aren't nodes in model (only data)
  values(cmodelSims, zNodes) <- nimData1$z#[,2:10] # Fill the values predicted by the model by new values where the extra years are NA
  values(cmodelSims, uNodes) <- nimData1$z#[,2:10] # Fill the values predicted by the model by new values where the extra years are NA
  
  # WE SET SXY
  sxy.start[1:M.aug,,1]<-matrix(sampmat[ itera[ite], s.which[indx] ], M.aug, 2)
  nimData1$sxy <- sxy.start
  
  ## FIND s NODES FROM THE FIRST YEAR (SHOULD BE REPLACED)
  sNodes <- sNodes[grep( ", 1]",sNodes)]
  for(i in 1:length(sNodes)){
    values(cmodelSims, sNodes[i]) <-sxy.start[i,,1]
  }
  
  ###set age
  age.est <- sampmat[itera[ite],age.which]
  ##set everyone not alive at all to all-0, ie, can be recruited in the NEW recruitment 
  ##model (they were never part of the superpopulation)
  z.nosuper <- which(apply(z.est,1,sum)==0)
  age.est[z.nosuper] <- 0
  #age.est[age.est>5] <- 5 #set to max age category
  age.start[1:M.aug, 1] <- age.est
  age.start[(M.aug+1):M.new, 1] <- 0
  
  nimData1$age <- age.start 
  values(cmodelSims, ageNodes) <- age.start #[,2:(Tt+t.new)] 
  
  ##also replace age.cat and agePlusOne
  nimData1$agePlusOne <- age.start[,1]+1 
  values(cmodelSims, agePO.Nodes) <- age.start[,1]+1
  
  nimData1$age.cat <- age.start 
  values(cmodelSims, age.cat.Nodes) <- age.start
  
  ## Also replace sex, because for the augmented individuals it changes each iteration
  
  sex.start[1:M.aug] <- sampmat[itera[ite],sex.which]
  nimData1$sex <- sex.start
  values(cmodelSims, sexNodes) <- sex.start
  
  ##set pcr, phi, sigD, beta.dens
  ##not sure if needed to set in data and cmodelSims...
  ##simplified since we calculated pcr for all iterations above
  nimData1$pcr <- pcr.all[itera[ite]]
  values(cmodelSims,"pcr") <- nimData1$pcr
  
  nimData1$phi.ad<-sampmat[itera[ite],'phi.ad']
  values(cmodelSims,"phi.ad") <- nimData1$phi.ad
  nimData1$phi.cub<-sampmat[itera[ite],'phi.cub']
  values(cmodelSims,"phi.cub") <- nimData1$phi.cub  
  nimData1$phi.sub<-sampmat[itera[ite],'phi.sub']
  values(cmodelSims,"phi.sub") <- nimData1$phi.sub
  
  nimData1$sigD<-sampmat[itera[ite],'sigD']
  values(cmodelSims,"sigD") <- nimData1$sigD
  
  nimData1$beta.dens<-sampmat[itera[ite],'beta.dens']
  values(cmodelSims,"beta.dens") <- nimData1$beta.dens
  
  #now we simulate 
  cmodelSims$simulate(nodes = nodesToSim,#c(zNodes,sNodes,yNodes),
                      includeData = F)#---if TRUE: want to simulate new values also for nodes considered as data
  
  # Store results from iteration
  
  Nmat.all[ite,] <- apply(cmodelSims$z,2,sum)
  Rmat.all[ite,] <- cmodelSims$R
  
  sxy.proj.all[ite,,,] <- c(cmodelSims$sxy)
  z.proj.all[ite,,] <- cmodelSims$z
  age.cat.proj.all[ite,,] <- cmodelSims$age.cat
  sex.proj.all[ite,] <- cmodelSims$sex
  

}



##plot trajectory
plot(1:(1+t.new), apply(Nmat.all,2,mean, na.rm = TRUE), type='l', ylim=range(Nmat.all, na.rm=TRUE))
for (ite in 1:20){
  points(1:(1+t.new), Nmat.all[ite,], type='l', col='lightgrey')
}
points(1:(1+t.new), apply(Nmat.all,2,mean), type='l')


## ------ 2.1.2. Subset abundance ----
## Unscale sxy coordinates to subset only individuals within core buffer

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")
load("habcoord.RData")

dimnames(sxy.proj.all)[[3]] <- c('x','y')
sxy.proj.all.uns <- scaleCoordsToHabitatGrid(coordsData = sxy.proj.all,## this are your sxy
                                              coordsHabitatGridCenter = G,# this is your unscaled habitat (as you used when scaling the habitat/detector to the habitat. G?
                                              scaleToGrid = FALSE)$coordsDataScaled

# Matrix to store abundance in the buffer each iteration and year

NIn_coreBuf.all <- matrix(NA,nrow = dim(z.proj.all)[1], ncol=dim(z.proj.all)[3]) # nrow = iterations, ncol = year
sp.check <- sp.t.check <- list()
ite=1
t=1

for(ite in 1:dim(z.proj.all)[1]){
  for(t in 1:dim(z.proj.all)[3]){
    
    which.alive <- which(z.proj.all[ite,,t]==1) # Select only the individuals alive (z=1)
    
    which.aliveSXY <- sxy.proj.all.uns[ite,which.alive,,t] # Retrieve the activity center for those individuals
    
    sp <- SpatialPoints(which.aliveSXY, proj4string=CRS(proj4string(Xbuf2))) # CONVERT SXY TO SPATIAL POINTS 
    
    which.In <- over(sp, Xbuf2) # Check which ones are in the buffer
    
    NIn_coreBuf.all[ite,t] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
    
    sp.t.check[[t]] <- sp # To check where points fall
    
  }
  
  sp.check[[ite]] <-  sp.t.check # To check where points fall
}

# Save
setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL/Predictions/20iter")
save(Nmat.all, Rmat.all, sxy.proj.all, z.proj.all, age.cat.proj.all, sxy.proj.all.uns, NIn_coreBuf.all, file = "proj_pcr.all.RData")

#average number of individuals without the buffer for each year 
colMeans(NIn_coreBuf.all)
NIn_coreBuf.all[,1]#posterior distrib for the first year

par(mfrow = c(2,3))
for (t in 1:dim(z.proj.all)[3]){
  plot(density(NIn_coreBuf.all[,t]), main = t)
  abline(v = colMeans(NIn_coreBuf.all)[t], col = "blue")
}

# Sum of individuals alive in total each year (without buffer)
NALL <- apply(z.proj.all,c(1,3),function(x) sum(x==1))
colMeans(NALL)

## ---- 2.3. PCR whole range (from min to max) ----
## ------ 2.2.1. Projection ----

pcr.all <- pcr_range 

# pcr_range contains joint all minimum and all maximum values of pcr (ordered by iteration, 1st minimums, 2nd maximums)
# Multiply sampmat so that iterations are repeated twice and they fit with pcr_range
# I don't know if this matters but in case. Then I can subet after anyway

sampmatx2 <- rbind(sampmat,sampmat)
##get some random iterations from posterior
itera <- sample(1:nrow(sampmatx2), 40)
nimData1 <- nimData
Nmat.range <- Rmat.range <- matrix(NA, length(itera), 1+t.new)

# To store as simlist
sxy.proj.range <- array(NA, c(length(itera), M.new, 2, t.new+1))
z.proj.range <- age.cat.proj.range <- array(NA, c(length(itera), M.new, t.new+1))

for(ite in 1:length(itera)){
  # WE UPDATE THE Z VALUES USING THE POSTERIORS PREDICTED Z FOR THE FIVE FIRST YEARS AND THEN 
  # USE NA FOR THE years to predict#
  
  # WE SET NA FOR Z FOR YEARS TO PREDICT - from posterior
  z.est <- matrix(sampmatx2[itera[ite],z.which], M.aug, Tt)
  z.start[1:M.aug, 1] <- sampmatx2[itera[ite],z.which1]
  z.start[(M.aug+1):M.new, 1]<-0
  
  nimData1$z <- z.start
  nimData1$u <- z.start
  
  # we set the values in the model - yr 1 aren't nodes in model (only data)
  values(cmodelSims, zNodes) <- nimData1$z#[,2:10] # Fill the values predicted by the model by new values where the extra years are NA
  values(cmodelSims, uNodes) <- nimData1$z#[,2:10] # Fill the values predicted by the model by new values where the extra years are NA
  
  # WE SET SXY
  sxy.start[1:M.aug,,1]<-matrix(sampmatx2[ itera[ite], s.which[indx] ], M.aug, 2)
  nimData1$sxy <- sxy.start
  
  ## FIND s NODES FROM THE FIRST YEAR (SHOULD BE REPLACED)
  sNodes <- sNodes[grep( ", 1]",sNodes)]
  for(i in 1:length(sNodes)){
    values(cmodelSims, sNodes[i]) <-sxy.start[i,,1]
  }
  
  ###set age
  age.est <- sampmatx2[itera[ite],age.which]
  ##set everyone not alive at all to all-0, ie, can be recruited in the NEW recruitment 
  ##model (they were never part of the superpopulation)
  z.nosuper <- which(apply(z.est,1,sum)==0)
  age.est[z.nosuper] <- 0
  #age.est[age.est>5] <- 5 #set to max age category
  age.start[1:M.aug, 1] <- age.est
  age.start[(M.aug+1):M.new, 1] <- 0
  
  nimData1$age <- age.start 
  values(cmodelSims, ageNodes) <- age.start #[,2:(Tt+t.new)] 
  
  ##also replace age.cat and agePlusOne
  nimData1$agePlusOne <- age.start[,1]+1 
  values(cmodelSims, agePO.Nodes) <- age.start[,1]+1
  
  nimData1$age.cat <- age.start 
  values(cmodelSims, age.cat.Nodes) <- age.start
  
  ##set pcr, phi, sigD, beta.dens
  ##not sure if needed to set in data and cmodelSims...
  ##simplified since we calculated pcr for all iterations above
  nimData1$pcr <- pcr.all[itera[ite]]
  values(cmodelSims,"pcr") <- nimData1$pcr
  
  nimData1$phi.ad<-sampmatx2[itera[ite],'phi.ad']
  values(cmodelSims,"phi.ad") <- nimData1$phi.ad
  nimData1$phi.cub<-sampmatx2[itera[ite],'phi.cub']
  values(cmodelSims,"phi.cub") <- nimData1$phi.cub  
  nimData1$phi.sub<-sampmatx2[itera[ite],'phi.sub']
  values(cmodelSims,"phi.sub") <- nimData1$phi.sub
  
  nimData1$sigD<-sampmatx2[itera[ite],'sigD']
  values(cmodelSims,"sigD") <- nimData1$sigD
  
  nimData1$beta.dens<-sampmatx2[itera[ite],'beta.dens']
  values(cmodelSims,"beta.dens") <- nimData1$beta.dens
  
  #now we simulate 
  cmodelSims$simulate(nodes = nodesToSim,#c(zNodes,sNodes,yNodes),
                      includeData = F)#---if TRUE: want to simulate new values also for nodes considered as data
  
  # Store results from iteration
  
  Nmat.range[ite,] <- apply(cmodelSims$z,2,sum)
  Rmat.range[ite,] <- cmodelSims$R
  
  sxy.proj.range[ite,,,] <- c(cmodelSims$sxy)
  z.proj.range[ite,,] <- cmodelSims$z
  age.cat.proj.range[ite,,] <- cmodelSims$age.cat
  
}


##plot trajectory
plot(1:(1+t.new), apply(Nmat.range,2,mean, na.rm = TRUE), type='l', ylim=range(Nmat.range, na.rm=TRUE))
for (ite in 1:20){
  points(1:(1+t.new), Nmat.range[ite,], type='l', col='lightgrey')
}
points(1:(1+t.new), apply(Nmat.range,2,mean), type='l')


# I think this doesn't make sense. We can maybe play with the results afterwards? I can run in anyway but it will reflect
# the same results as if we would join the abundance projection from both pcr (core and all). No need to reproject

