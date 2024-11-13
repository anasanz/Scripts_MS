## -------------------------------------------------
##                 Predict future population
##                      Data script 3.1
##                         Model 3.7
##              DECREASING effect of distance to core
## ------------------------------------------------- 


rm(list = ls())

## Load packages

library(Rcpp)
library(RcppArmadillo)
library(RcppProgress)
library(nimble)
library(nimbleSCR)
library(basicMCMCplots)
library(R.utils)
library(abind)
library(coda)
library(rgdal)
library(viridis)

detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
  
}

#source("D:/MargSalas/Scripts_MS/Functions/ProcessCodaOutput.R")
#source("D:/MargSalas/Scripts_MS/Functions/plot.violins2.r")
#source("D:/MargSalas/Scripts_MS/Functions/DoScale.r")

#setwd("D:/MargSalas/Scripts_MS/Functions/Nimble")
##sourceCpp("GetSpaceUse_PD.cpp")
#sourceCpp("GetDensity_PD.cpp")
#source("getDensityInput.R")

#setwd("D:/MargSalas/Scripts_MS/Functions/Nimble")
setwd("~/Scripts_MS/Functions/Nimble")
source('dbinomLocal_normalBear_rbinom2.R')

# Load buffer core area
setwd("~/Data_server/Oso/Data_for_Prediction")
#setwd("D:/MargSalas/Oso/Datos/GIS/Countries")
Xbuf <- readOGR("Buffer_statespace.shp")
Xbuf2 <- readOGR("Buffer_8500_traps.shp")


# Load data and model results

#setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/3.openSCR_Age/Data_server")
setwd("~/Data_server/Oso/Data_for_Prediction")
load("Data_Model3-3.1_CYRIL_allparams.RData")

#setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL")
setwd("~/Data_server/Oso/Data_for_Prediction")
load("myResults_3-3.1_param.RData")
sampmat1 <- do.call(rbind, nimOutput)

#setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL")
setwd("~/Data_server/Oso/Data_for_Prediction")
load("myResults_3-3.1_sxy.RData")
sampmat2 <- do.call(rbind, nimOutputSXY)

sampmat <- cbind(sampmat1, sampmat2)

#save(sampmat, file = "sampmat_3.1.RData")

#setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL")
#load("sampmat_3.1.RData")

M.aug <- nimConstants$M
Tt <- nimConstants$Nyr

## ---- 1. Re-Build model with extra dimensions ----
## ---- 1.1. Create starting z's, age's and sxy's ----

# These will be based on model output (last year), new M and new time frame
t.new <- 5 # sim 5 addl years
M.new <- 700 # new augmentation limit

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

#setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL")
#load("pcr_corebuf.RData")
#load("pcr_all.RData")
#load("pcr_range.RData")
setwd("~/Data_server/Oso/Data_for_Prediction")
load("pcr_corebuf_fem.RData")
load("pcr_all_fem.RData")

pcr1_corebuf <- apply(pcr_corebuf_fem[[3]],1,mean) # ASP: Mean per capita recruitment per iteration
pcr2_all <- apply(pcr_all_fem[[3]],1,mean) # ASP: Mean per capita recruitment per iteration


# Take pcr value of one iteration to build the model
pcr <- pcr1_corebuf[1]

#B <- matrix(sampmat[1:nrow(sampmat), grep('B\\[', colnames(sampmat))], nrow(sampmat), Tt) 
#pcr_all[[1]] == B[,c(2:5)] # ASP: B is the same than R in all state space! (except column 1, which is abundance)

## ---- 1.3. Other parameters ----
# We don't include the parameters related with detection (not in projection model)

phi.ad <- sampmat[1,'phi.ad']
phi.cub <- sampmat[1,'phi.cub']
phi.sub <- sampmat[1,'phi.sub']
omega <- sampmat[1,'omega']
sigD <- sampmat[1,'sigD']

# Increasing beta dens
beta.dens <- c(sampmat[1,'beta.dens'], -0.30, -0.25, -0.20, -0.10, 0)


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

#setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Model")
setwd("~/Scripts_MS/Oso/PopDyn/SCR/Model")
source('3.7.SCRopen_PREDICT in Nimble RS.r')

# Create model using projection code (per capita recruitment, no observations)
model <- nimbleModel( code = SCROpen.PR.distcore.t.pcrFem,
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
betaDensNodes <- samplerConfList[grep("beta.dens\\[",samplerConfList)]
sexNodes <- samplerConfList[grep("sex\\[",samplerConfList)]

##get some random iterations from posterior and save 
#itera <- sample(1:nrow(sampmat), 5000)
setwd("~/Model_results/Oso/ALLiter")
load("itera.RData")


## ---- 2.1. PCR inside core buffer ----
## ------ 2.1.1. Projection ----

pcr.core <- pcr1_corebuf

## ---- SCENARIO 1: CONSTANT DISTCORE ----

beta.dens.predict1 <- rep(mean(sampmat[,'beta.dens']),6)

nimData1 <- nimData
#Nmat.core <- Rmat.core <- matrix(NA, length(itera), 1+t.new)

# To store as simlist
sxy.proj.core.consDist <- array(NA, c(length(itera), M.new, 2, t.new+1))
z.proj.core.consDist <- age.cat.proj.core.consDist <- array(NA, c(length(itera), M.new, t.new+1))
sex.proj.core.consDist <- matrix(NA, length(itera), M.new)


detachAllPackages()
library(nimbleSCR)

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
  
  #nimData1$beta.dens<-sampmat[itera[ite],'beta.dens']
  # Beta dens: Always the same incresing values: -0.3538166 -0.3000000 -0.2500000 -0.2000000 -0.1000000  0.0000000 
  values(cmodelSims,"beta.dens") <- beta.dens.predict1
  
  #now we simulate 
  cmodelSims$simulate(nodes = nodesToSim,#c(zNodes,sNodes,yNodes),
                      includeData = F)#---if TRUE: want to simulate new values also for nodes considered as data
  
  #Nmat.core[ite,] <- apply(cmodelSims$z,2,sum)
  #Rmat.core[ite,] <- cmodelSims$R
  
  sxy.proj.core.consDist[ite,,,] <- c(cmodelSims$sxy)
  z.proj.core.consDist[ite,,] <- cmodelSims$z
  age.cat.proj.core.consDist[ite,,] <- cmodelSims$age.cat
  
  sex.proj.core.consDist[ite,] <- cmodelSims$sex
  
}

#setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL/Predictions/20iter")
setwd("~/Model_results/Oso/ALLiter")
save(sxy.proj.core.consDist, z.proj.core.consDist, age.cat.proj.core.consDist,sex.proj.core.consDist, file = "proj_pcr.core.fem.consDist.RData")

## ---- SCENARIO 2: DISTCORE DECREASES PROGRESSIVELY ----

#beta.dens.predict2 <- c(mean(sampmat[,'beta.dens']), -0.30, -0.25, -0.20, -0.10, 0)
beta.dens.predict2 <- c(sampmat[1,'beta.dens'], -0.27, -0.24, -0.21, -0.18, -0.15) # More proggresive than before


nimData1 <- nimData
#Nmat.core <- Rmat.core <- matrix(NA, length(itera), 1+t.new)

# To store as simlist
sxy.proj.core.decDist <- array(NA, c(length(itera), M.new, 2, t.new+1))
z.proj.core.decDist <- age.cat.proj.core.decDist <- array(NA, c(length(itera), M.new, t.new+1))
sex.proj.core.decDist <- matrix(NA, length(itera), M.new)


detachAllPackages()
library(nimbleSCR)

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
  
  #nimData1$beta.dens<-sampmat[itera[ite],'beta.dens']
  # Beta dens: Always the same incresing values: -0.3538166 -0.3000000 -0.2500000 -0.2000000 -0.1000000  0.0000000 
  values(cmodelSims,"beta.dens") <- beta.dens.predict2
  
  #now we simulate 
  cmodelSims$simulate(nodes = nodesToSim,#c(zNodes,sNodes,yNodes),
                      includeData = F)#---if TRUE: want to simulate new values also for nodes considered as data
  
  # Nmat.core[ite,] <- apply(cmodelSims$z,2,sum)
  # Rmat.core[ite,] <- cmodelSims$R
  
  sxy.proj.core.decDist[ite,,,] <- c(cmodelSims$sxy)
  z.proj.core.decDist[ite,,] <- cmodelSims$z
  age.cat.proj.core.decDist[ite,,] <- cmodelSims$age.cat
  
  sex.proj.core.decDist[ite,] <- cmodelSims$sex
  
  
}

#setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL/Predictions/20iter")
setwd("~/Model_results/Oso/ALLiter")
save(sxy.proj.core.decDist, z.proj.core.decDist, age.cat.proj.core.decDist, sex.proj.core.decDist, file = "proj_pcr.core.fem.decDist.RData")


## ---- 2.2. PCR in whole state space ----
## ------ 2.2.1. Projection ----

pcr.all <- pcr2_all

## ---- SCENARIO 1: CONSTANT DISTCORE ----

beta.dens.predict1 <- rep(mean(sampmat[,'beta.dens']),6)

nimData1 <- nimData
#Nmat.all <- Rmat.all <- matrix(NA, length(itera), 1+t.new)

# To store as simlist
sxy.proj.all.consDist <- array(NA, c(length(itera), M.new, 2, t.new+1))
z.proj.all.consDist <- age.cat.proj.all.consDist <- array(NA, c(length(itera), M.new, t.new+1))
sex.proj.all.consDist <- matrix(NA, length(itera), M.new)


detachAllPackages()
library(nimbleSCR)

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
  
  nimData1$beta.dens<-beta.dens.predict1
  #values(cmodelSims,"beta.dens") <- nimData1$beta.dens
  values(cmodelSims,"beta.dens") <- beta.dens.predict1
  
  
  #now we simulate 
  cmodelSims$simulate(nodes = nodesToSim,#c(zNodes,sNodes,yNodes),
                      includeData = F)#---if TRUE: want to simulate new values also for nodes considered as data
  
  # Store results from iteration
  
  Nmat.all[ite,] <- apply(cmodelSims$z,2,sum)
  Rmat.all[ite,] <- cmodelSims$R
  
  sxy.proj.all.consDist[ite,,,] <- c(cmodelSims$sxy)
  z.proj.all.consDist[ite,,] <- cmodelSims$z
  age.cat.proj.all.consDist[ite,,] <- cmodelSims$age.cat
  
  sex.proj.all.consDist[ite,] <- cmodelSims$sex
  
}

# Save
#setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL/Predictions/20iter")
setwd("~/Model_results/Oso/ALLiter")
save(z.proj.all.consDist, age.cat.proj.all.consDist, sxy.proj.all.consDist, sex.proj.all.consDist, file = "proj_pcr.all.fem.consDist.RData")

## ---- SCENARIO 2: DISTCORE DECREASES PROGRESSIVELY ----

#beta.dens.predict2 <- c(mean(sampmat[,'beta.dens']), -0.30, -0.25, -0.20, -0.10, 0)
beta.dens.predict2 <- c(sampmat[1,'beta.dens'], -0.27, -0.24, -0.21, -0.18, -0.15) # More proggresive than before


nimData1 <- nimData
#Nmat.all <- Rmat.all <- matrix(NA, length(itera), 1+t.new)

# To store as simlist
sxy.proj.all.decDist <- array(NA, c(length(itera), M.new, 2, t.new+1))
z.proj.all.decDist <- age.cat.proj.all.decDist <- array(NA, c(length(itera), M.new, t.new+1))
sex.proj.all.decDist <- matrix(NA, length(itera), M.new)


detachAllPackages()
library(nimbleSCR)

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
  
  nimData1$beta.dens<-beta.dens.predict1
  #values(cmodelSims,"beta.dens") <- nimData1$beta.dens
  values(cmodelSims,"beta.dens") <- beta.dens.predict2
  
  
  #now we simulate 
  cmodelSims$simulate(nodes = nodesToSim,#c(zNodes,sNodes,yNodes),
                      includeData = F)#---if TRUE: want to simulate new values also for nodes considered as data
  
  # Store results from iteration
  
  Nmat.all[ite,] <- apply(cmodelSims$z,2,sum)
  Rmat.all[ite,] <- cmodelSims$R
  
  sxy.proj.all.decDist[ite,,,] <- c(cmodelSims$sxy)
  z.proj.all.decDist[ite,,] <- cmodelSims$z
  age.cat.proj.all.decDist[ite,,] <- cmodelSims$age.cat
  
  sex.proj.all.decDist[ite,] <- cmodelSims$sex
  
  
}

# Save
#setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL/Predictions/20iter")
setwd("~/Model_results/Oso/ALLiter")
save(z.proj.all.decDist, age.cat.proj.all.decDist, sxy.proj.all.decDist, sex.proj.all.decDist, file = "proj_pcr.all.fem.decDist.RData")

