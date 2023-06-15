## -------------------------------------------------
##                 Predict future population
##                      Data script 3.1
##                         Model 3.7
## ------------------------------------------------- 

## Our abundance estimates are based in subsetting our population to our trapping area (core buffer)
## To predict abundance spatially we will estimate abundance choosing as starting point 
# all individuals from state-space, and subset afterwards (as we do in the model).
#   1. Choose the individuals in all state space as starting point
#   2. Load estimated pcr from individuals: WE CHOSE THE ONE IN WHOLE STATE SPACE, as we use ind in all state space and then subset
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
Xbuf3 <- gBuffer(Xbuf2, width = -8000) # To ensure that females that I add are simulated in and don't go out so that they get inside the pop.estimates


setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")
load("habcoord.RData")


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

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL")
load("pcr_all_fem.RData")

pcr2_all <- apply(pcr_all_fem[[3]],1,mean) # ASP: Mean per capita recruitment per iteration

# Take pcr value of one iteration to build the model
pcr <- pcr2_all[1]

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

##get some random iterations from posterior
itera <- sample(1:nrow(sampmat), 20)
setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL/Predictions/20iter_sc3in")
#save(itera, file = "itera.RData")
load("itera.RData")

## ---- 2.0. Scenario 0: Normal projection ----

pcr.all <- pcr2_all # We have decided to run the projection with the PCR in the whole state space

nimData1 <- nimData
Nmat.all.sc0 <- Rmat.all.sc0 <- matrix(NA, length(itera), 1+t.new)

# To store as simlist
sxy.proj.all.sc0 <- array(NA, c(length(itera), M.new, 2, t.new+1))
z.proj.all.sc0 <- age.cat.proj.all.sc0 <- array(NA, c(length(itera), M.new, t.new+1))
sex.proj.all.sc0 <- matrix(NA, length(itera), M.new)


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
  
  Nmat.all.sc0[ite,] <- apply(cmodelSims$z,2,sum)
  Rmat.all.sc0[ite,] <- cmodelSims$R
  
  sxy.proj.all.sc0[ite,,,] <- c(cmodelSims$sxy)
  z.proj.all.sc0[ite,,] <- cmodelSims$z
  age.cat.proj.all.sc0[ite,,] <- cmodelSims$age.cat
  sex.proj.all.sc0[ite,] <- cmodelSims$sex
  
}

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL/Predictions/20iter_sc3in")
#setwd("~/Model_results/Oso/ALLiter")
save(sxy.proj.all.sc0, z.proj.all.sc0, age.cat.proj.all.sc0,sex.proj.all.sc0, file = "proj_pcr.all.fem.sc0.RData")


## ---- 2.1. Scenario 1: Death of 10% of the females in year 2021 ----

pcr.all <- pcr2_all # We have decided to run the projection with the PCR in the whole state space

##get some random iterations from posterior
#itera <- sample(1:nrow(sampmat), 20)
nimData1 <- nimData
Nmat.all.sc1 <- Rmat.all.sc1 <- matrix(NA, length(itera), 1+t.new)

# To store as simlist
sxy.proj.all.sc1 <- array(NA, c(length(itera), M.new, 2, t.new+1))
z.proj.all.sc1 <- age.cat.proj.all.sc1 <- array(NA, c(length(itera), M.new, t.new+1))
sex.proj.all.sc1 <- matrix(NA, length(itera), M.new)


for(ite in 1:length(itera)){
  # WE UPDATE THE Z VALUES USING THE POSTERIORS PREDICTED Z FOR THE FIVE FIRST YEARS AND THEN 
  # USE NA FOR THE years to predict#
  
  ## We first replace the sex. We need to replace because for the augmented individuals it changes each iteration
  
  sex.start[1:M.aug] <- sampmat[itera[ite],sex.which]
  nimData1$sex <- sex.start
  values(cmodelSims, sexNodes) <- sex.start
  
  # WE SET NA FOR Z FOR YEARS TO PREDICT - from posterior
  z.est <- matrix(sampmat[itera[ite],z.which], M.aug, Tt)
  z.start[1:M.aug, 1] <- sampmat[itera[ite],z.which1]
  z.start[(M.aug+1):M.new, 1]<-0
  
  ## APPLY SCENARIO 1: We turn 10% of the females alive to death (z=0)
  which.females <- which(sex.start[1:M.aug] == 0) 
  fem.alive <- sum(z.start[which.females, 1]) # How many females are alive?
  kill <- round(fem.alive*10/100,0) # How many do we kill? 10%
  which.kill <- sample(which(z.start[which.females, 1] == 1), kill) # Randomly sample the FEMALES we will kill
  z.start[which.females[which.kill],1] <- 0 # Set to 0 the females that we wanted to kill
  
  fem.alive - kill == sum(z.start[which.females, 1]) # Check I did it right
  
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
  
  Nmat.all.sc1[ite,] <- apply(cmodelSims$z,2,sum)
  Rmat.all.sc1[ite,] <- cmodelSims$R
  
  sxy.proj.all.sc1[ite,,,] <- c(cmodelSims$sxy)
  z.proj.all.sc1[ite,,] <- cmodelSims$z
  age.cat.proj.all.sc1[ite,,] <- cmodelSims$age.cat
  sex.proj.all.sc1[ite,] <- cmodelSims$sex
  
}

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL/Predictions/20iter_sc3in")
#setwd("~/Model_results/Oso/ALLiter")
save(sxy.proj.all.sc1, z.proj.all.sc1, age.cat.proj.all.sc1,sex.proj.all.sc1, file = "proj_pcr.all.fem.sc1.RData")

## ---- 2.2. Scenario 2: Death of 20% of the females in year 2021 ----

pcr.all <- pcr2_all # We have decided to run the projection with the PCR in the whole state space

##get some random iterations from posterior
#itera <- sample(1:nrow(sampmat), 20)
nimData1 <- nimData
Nmat.all.sc2 <- Rmat.all.sc2 <- matrix(NA, length(itera), 1+t.new)

# To store as simlist
sxy.proj.all.sc2 <- array(NA, c(length(itera), M.new, 2, t.new+1))
z.proj.all.sc2 <- age.cat.proj.all.sc2 <- array(NA, c(length(itera), M.new, t.new+1))
sex.proj.all.sc2 <- matrix(NA, length(itera), M.new)


for(ite in 1:length(itera)){
  # WE UPDATE THE Z VALUES USING THE POSTERIORS PREDICTED Z FOR THE FIVE FIRST YEARS AND THEN 
  # USE NA FOR THE years to predict#
  
  ## We first replace the sex. We need to replace because for the augmented individuals it changes each iteration
  
  sex.start[1:M.aug] <- sampmat[itera[ite],sex.which]
  nimData1$sex <- sex.start
  values(cmodelSims, sexNodes) <- sex.start
  
  # WE SET NA FOR Z FOR YEARS TO PREDICT - from posterior
  z.est <- matrix(sampmat[itera[ite],z.which], M.aug, Tt)
  z.start[1:M.aug, 1] <- sampmat[itera[ite],z.which1]
  z.start[(M.aug+1):M.new, 1]<-0
  
  ## APPLY SCENARIO 2: We turn 20% of the females alive to death (z=0)
  which.females <- which(sex.start[1:M.aug] == 0) 
  fem.alive <- sum(z.start[which.females, 1]) # How many females are alive?
  kill <- round(fem.alive*20/100,0) # How many do we kill? 20%
  which.kill <- sample(which(z.start[which.females, 1] == 1), kill) # Randomly sample the FEMALES we will kill
  z.start[which.females[which.kill],1] <- 0 # Set to 0 the females that we wanted to kill
  
  fem.alive - kill == sum(z.start[which.females, 1]) # Check I did it right
  
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
  
  Nmat.all.sc2[ite,] <- apply(cmodelSims$z,2,sum)
  Rmat.all.sc2[ite,] <- cmodelSims$R
  
  sxy.proj.all.sc2[ite,,,] <- c(cmodelSims$sxy)
  z.proj.all.sc2[ite,,] <- cmodelSims$z
  age.cat.proj.all.sc2[ite,,] <- cmodelSims$age.cat
  sex.proj.all.sc2[ite,] <- cmodelSims$sex
  
}

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL/Predictions/20iter_sc3in")
#setwd("~/Model_results/Oso/ALLiter")
save(sxy.proj.all.sc2, z.proj.all.sc2, age.cat.proj.all.sc2,sex.proj.all.sc2, file = "proj_pcr.all.fem.sc2.RData")

## ---- 2.3. Scenario 3: Death of 20% of the females in year 2021, reintroduction of 5 females in year 2024----

# I need to include extra augmented individuals and start the model in 2024, having only 2 extra years
# Dimensions of this model change, so I need to set up the model again

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL/Predictions/20iter_sc3in")
load("proj_pcr.all.fem.sc2.RData")

## ------ 2.3.1. Create starting z's, age's and sxy's AGAIN ----
# In this projection the starting point is year 2024, so I need to set model structure again with new augmented individuals

t.new2 <- 2 # Only project for 2 years extra (2025, 2026)
M.aug2 <- 700 # Augmented individuals from previous model
M.new2 <- 1000 # All augmented individuals (700 previous + 300 new)

# Z

z.start2 <- matrix(NA, M.new2, 1+t.new2)
z.start2[1:M.aug2, 1] <- z.proj.all.sc2[1,,4] # Year 2024 as starting point, which is the 4th time point of scenario 2
z.start2[(M.aug2+1):M.new2, 1] <- 0 # ASP: Augmented as not yet recruited

# AGE
## I need to join all z from previous years to know the ind that have never been alive 
# from 2017 - 2024, to modify the age variable

zAUGM <- matrix(0,nrow = 400, ncol = Tt)
z.est2.0 <- matrix(sampmat[1, z.which], M.aug, Tt) # Years 2017 - 2021 from data
z.est2.0 <- rbind(z.est2.0, zAUGM)

z.est2.1 <- matrix(z.proj.all.sc2[1,,c(2:4)], M.aug2,3) # Take years 2022, 2023 and 2024 for predictions
z.est2 <- cbind(z.est2.0, z.est2.1)

age.start2 <- matrix(NA, M.new2, 1+t.new2)

age.est2 <- age.cat.proj.all.sc2[1,,4] # ASP: Ages last year (year 2024)

##set everyone not alive at all to all-0, ie, can be recruited in the NEW recruitment 
##model (they were never part of the superpopulation -> ASP: because of how the model is formulated z=u*w)
z.nosuper2 <- which(apply(z.est2,1,sum)==0)
age.est2[z.nosuper2] <- 0

age.start2[1:M.aug2, 1] <- age.est2 # ASP: Starting AGE for projection: AGE of year 5 
age.start2[(M.aug2+1):M.new2, 1] <- 0 # ASP: Augmented as not yet recruited

# SXY

sxy.start2 <- array(NA, c(M.new2, 2, 1+t.new2))
sxy.start2[1:M.aug2,,1] <- sxy.proj.all.sc2[1,,,4] # ASP: Starting SXY for projection: SXY of year 2024

# SEX
sex.start2 <- rep(NA, M.new2)
sex.start2[1:M.aug2] <- sex.proj.all.sc2[1,]


## ------ 2.3.2. Per capita recruitment ----
# Pcr is estimated for each iteration, and constitute the mean of PCR in all years. So in this
# new prediction model, I take the same pcr of each iteration (pcr.all)

pcr.all <- pcr2_all # We have decided to run the projection with the PCR in the whole state space
pcr # To build model, we only take value of one iteration

## ------ 2.3.3. Other parameters ----
# Remain the same as for previous model

## ------ 2.3.4. Set un constants and build model ----

##set up constants
nimConstants2 <- list(
  M=(M.new2), numHabWindows = nimConstants$numHabWindows, 
  numGridRows = nimConstants$numGridRows, numGridCols = nimConstants$numGridCols, 
  Nyr=(1+t.new2),
  max.age = nimConstants$max.age
)

##set up data
nimData2 <- list(habDens = nimData$habDens, 
                 lowerHabCoords = nimData$lowerHabCoords, upperHabCoords = nimData$upperHabCoords,
                 habitatGrid = nimData$habitatGrid,
                 z = z.start2,
                 sxy = sxy.start2,
                 agePlusOne=age.start2[,1]+1,
                 age = age.start2,
                 age.cat = age.start2,
                 u = z.start2,
                 pcr = pcr, ##transformed per capita recruitment
                 phi.ad = phi.ad,
                 phi.cub = phi.cub,
                 phi.sub = phi.sub,
                 beta.dens = beta.dens,
                 sigD = sigD,
                 psi = psi,
                 pi.uncond = pi.uncond,
                 sex = sex.start2,
                 omega = omega
                 
)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Model")
source('3.7.SCRopen_PREDICT in Nimble RS.r')

# Create model using projection code (per capita recruitment, no observations)
model2 <- nimbleModel( code = SCRhab.Open.diftraps.age.effortTrapBhCov.sigsexage.PR.pcr.fem,
                       constants = nimConstants2,
                       data = nimData2,
                       #inits = nimInits,
                       check = F,       
                       calculate = T)  

nodesToSim2 <- model2$getDependencies(c("sigD","phi.ad","phi.cub","phi.sub",
                                        "pcr", "beta.dens", 'pi.uncond', 'psi', 'omega'),
                                      self = F,
                                      downstream = T,
                                      returnScalarComponents = TRUE)
model2$simulate(nodesToSim2, includeData = FALSE)

##compile for faster simulation
cmodelSims2 <- compileNimble(model2)

## ------ 2.3.5. Run projection for several iterations and subset abundance ----

# WE FIND THE S AND Z NODES THAT SHOULDNT BE PREDICTED FOR. 
samplerConfList2 <- model2$getNodeNames()
zNodes2 <- samplerConfList2[grep("z\\[",samplerConfList2)]
sNodes2 <- samplerConfList2[grep("sxy",samplerConfList2)]
uNodes2 <- samplerConfList2[grep("u\\[",samplerConfList2)]
ageNodes2 <- samplerConfList2[grep("age\\[",samplerConfList2)]
age.cat.Nodes2 <- samplerConfList2[grep("age.cat\\[",samplerConfList2)]
agePO.Nodes2 <-samplerConfList2[grep("agePlusOne",samplerConfList2)]
sexNodes2 <- samplerConfList2[grep("sex\\[",samplerConfList2)]

#length(values(cmodelSims, sNodes)) #800*2*6 (sNodes are 4800 because dim 2 is 1:2)
#sxy.start2[1:M.aug2,,1] == model2$sxy[1:M.aug2,,1]
#values(cmodelSims2, sNodes2)[1:1400] == c(t(sxy.start2[1:M.aug2,,1]))

##get some random iterations from posterior
#itera <- sample(1:nrow(sampmat), 20)
#setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL/Predictions/20iter_sc3in")
#save(itera, file = "itera.RData")
#load("itera.RData")

# Store

nimData1 <- nimData2
Nmat.all.sc3 <- Rmat.all.sc3 <- matrix(NA, length(itera), 1+t.new2)

sxy.proj.all.sc3 <- array(NA, c(length(itera), M.new2, 2, t.new2+1)) # To store as simlist
z.proj.all.sc3 <- age.cat.proj.all.sc3 <- array(NA, c(length(itera), M.new2, t.new2+1))
sex.proj.all.sc3 <- matrix(NA, length(itera), M.new2)


for(ite in 1:length(itera)){
  # WE UPDATE THE Z VALUES USING THE POSTERIORS PREDICTED Z FOR THE FIVE FIRST YEARS AND THEN 
  # USE NA FOR THE years to predict#
  
  ## We first replace the sex. We need to replace because for the augmented individuals it changes each iteration
  
  sex.start2[1:M.aug2] <- sex.proj.all.sc2[ite,] # itera[ite] when I have the full projection from sc2, and the same in whole code
  nimData1$sex <- sex.start2
  values(cmodelSims2, sexNodes2) <- sex.start2
  
  # WE SET SXY (We do it first to know if the females that we turn alive fall within the sampling buffer)
  sxy.start2[1:M.aug2,,1] <- sxy.proj.all.sc2[ite,,,4]
  nimData1$sxy <- sxy.start2
  
  ## FIND s NODES FROM THE FIRST YEAR (SHOULD BE REPLACED)
  sNodes2 <- sNodes2[grep( ", 1]",sNodes2)]
  for(i in 1:length(sNodes2)){
    values(cmodelSims2, sNodes2[i]) <-sxy.start2[i,,1]
  }
  
  # WE SET NA FOR Z FOR YEARS TO PREDICT - from posterior
  
  z.start2[1:M.aug2, 1] <- z.proj.all.sc2[ite,,4] # Year 2024 as starting point, which is the 4th time point of scenario 2
  z.start2[(M.aug2+1):M.new2, 1] <- 0 # ASP: Augmented as not yet recruited
  
  ## APPLY SCENARIO 3: We turn 7 of the females death (z=0) to alive (z = 1)
  ## WARNING: We turn any of the previous/new augmented individuals. I think it doesn't matter, but keep in mind in case it doesn't work
  
  which.females <- which(sex.start2[1:M.aug2] == 0) 
  fem.alive <- sum(z.start2[which.females, 1]) # How many females are alive? Just to check, not useful in this scenario
  turnAlive <- 7 # How many do we turn alive?
  
  # Activate females only that fall in ss: Unscale ss
  dimnames(sxy.start2)[[2]] <- c('x','y') 
  sxy.start2.uns <- scaleCoordsToHabitatGrid(coordsData = sxy.start2,## this are your sxy
                                             coordsHabitatGridCenter = G,# this is your unscaled habitat (as you used when scaling the habitat/detector to the habitat. G?
                                             scaleToGrid = FALSE)$coordsDataScaled
  allwhich.turnAlive <- list()
  for (a in 1:turnAlive){ # Sample points until they fall in the buffer
    repeat{
    which.turnAlive <- sample(which(z.start2[which.females, 1] == 0),1)
    which.turnAliveSXY <- sxy.start2.uns[which.turnAlive,,1]
    which.turnAliveSXY <- data.frame(x = which.turnAliveSXY[1], y = which.turnAliveSXY[2])
    sp <- SpatialPoints(which.turnAliveSXY, proj4string=CRS(proj4string(Xbuf2)))
    ovr <- is.na(over(sp, Xbuf3)) # ovr is TRUE when there is NA (not overlap)
    if(ovr == FALSE) break
    }
    allwhich.turnAlive[[a]] <- which.turnAlive
  }
  
  which.turnAliveIn <- unlist(allwhich.turnAlive)

  z.start2[which.females[which.turnAliveIn],1] <- 1 # Set to 1 the females that we wanted to turn alive
  
  #fem.alive + turnAlive == sum(z.start2[which.females, 1]) # Check I did it right
  
  nimData1$z <- z.start2
  nimData1$u <- z.start2
  
  # we set the values in the model - yr 1 aren't nodes in model (only data)
  values(cmodelSims2, zNodes2) <- nimData1$z#[,2:10] # Fill the values predicted by the model by new values where the extra years are NA
  values(cmodelSims2, uNodes2) <- nimData1$z#[,2:10] # Fill the values predicted by the model by new values where the extra years are NA
  
  # AGE
  ## I need to join all z from previous years to know the ind that have never been alive 
  # from 2017 - 2024, to modify the age variable
  
  zAUGM <- matrix(0,nrow = 400, ncol = Tt)
  z.est2.0 <- matrix(sampmat[ite, z.which], M.aug, Tt)  # Take years 2017 - 2021 from data
  z.est2.0 <- rbind(z.est2.0, zAUGM)
  
  z.est2.1 <- matrix(z.proj.all.sc2[ite,,c(2,3)], M.aug2, 2) # Take years 2022, 2023 from predictions
  z.est2.2 <- matrix(z.start2[1:M.aug2,1], M.aug2, 1) # Take year 2024 from modified reintroduced individuals (otherwise age can turn 0 and dead again!)
  
  z.est2 <- cbind(z.est2.0, z.est2.1, z.est2.2)
  
  age.est2 <- age.cat.proj.all.sc2[ite,,4] # ASP: Ages last year (year 2024)
  
  ##set everyone not alive at all to all-0, ie, can be recruited in the NEW recruitment 
  ##model (they were never part of the superpopulation -> ASP: because of how the model is formulated z=u*w)
  z.nosuper2 <- which(apply(z.est2,1,sum)==0)
  age.est2[z.nosuper2] <- 0
  
  age.start2[1:M.aug2, 1] <- age.est2 # ASP: Starting AGE for projection: AGE of year 5 
  age.start2[(M.aug2+1):M.new2, 1] <- 0 # ASP: Augmented as not yet recruited
  
  nimData1$age <- age.start2 
  values(cmodelSims2, ageNodes2) <- age.start2 #[,2:(Tt+t.new)] 
  
  ##also replace age.cat and agePlusOne
  nimData1$agePlusOne <- age.start2[,1]+1 
  values(cmodelSims2, agePO.Nodes2) <- age.start2[,1]+1
  
  nimData1$age.cat <- age.start2 
  values(cmodelSims2, age.cat.Nodes2) <- age.start2
  
  
  ##set pcr, phi, sigD, beta.dens
  ##not sure if needed to set in data and cmodelSims...
  ##simplified since we calculated pcr for all iterations above
  nimData1$pcr <- pcr.all[itera[ite]]
  values(cmodelSims2,"pcr") <- nimData1$pcr
  
  nimData1$phi.ad<-sampmat[itera[ite],'phi.ad']
  values(cmodelSims2,"phi.ad") <- nimData1$phi.ad
  nimData1$phi.cub<-sampmat[itera[ite],'phi.cub']
  values(cmodelSims2,"phi.cub") <- nimData1$phi.cub  
  nimData1$phi.sub<-sampmat[itera[ite],'phi.sub']
  values(cmodelSims2,"phi.sub") <- nimData1$phi.sub
  
  nimData1$sigD <- sampmat[itera[ite],'sigD']
  values(cmodelSims2,"sigD") <- nimData1$sigD
  
  nimData1$beta.dens<-sampmat[itera[ite],'beta.dens']
  values(cmodelSims2,"beta.dens") <- nimData1$beta.dens
  
  #now we simulate 
  cmodelSims2$simulate(nodes = nodesToSim2,#c(zNodes,sNodes,yNodes),
                       includeData = F)#---if TRUE: want to simulate new values also for nodes considered as data
  
  # Store results from iteration
  
  Nmat.all.sc3[ite,] <- apply(cmodelSims2$z,2,sum)
  Rmat.all.sc3[ite,] <- cmodelSims2$R
  
  sxy.proj.all.sc3[ite,,,] <- c(cmodelSims2$sxy)
  z.proj.all.sc3[ite,,] <- cmodelSims2$z
  age.cat.proj.all.sc3[ite,,] <- cmodelSims2$age.cat
  sex.proj.all.sc3[ite,] <- cmodelSims2$sex
  
}

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL/Predictions/20iter_sc3in")
#setwd("~/Model_results/Oso/ALLiter")
save(sxy.proj.all.sc3, z.proj.all.sc3, age.cat.proj.all.sc3,sex.proj.all.sc3, file = "proj_pcr.all.fem.sc3.RData")
