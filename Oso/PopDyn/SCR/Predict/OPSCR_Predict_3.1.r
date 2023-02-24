
## -------------------------------------------------
##                 Predict future population
##                      Data script 3.1
##                         Model 3.7
## ------------------------------------------------- 

rm(list = ls())

## Load packages

library(nimble)
library(nimbleSCR)
library(basicMCMCplots)
library(R.utils)
library(abind)
library(coda)

source("D:/MargSalas/Scripts_MS/Functions/ProcessCodaOutput.R")
source("D:/MargSalas/Scripts_MS/Functions/plot.violins2.r")
source("D:/MargSalas/Scripts_MS/Functions/DoScale.r")

setwd("D:/MargSalas/Scripts_MS/Functions/Nimble")
source('dbinomLocal_normalBear_rbinom2.R')

# Load data and model results

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/3.openSCR_Age/Data_server")
load("Data_Model3-3.1_CYRIL_allparams.RData")

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams")
load("myResults_3-3.1_param.RData")

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1")
load("myResults_3-3.1_sxy.RData")

#### PREDICT FROM THE MODEL POSTERIORS ####
# WE TAKE THE SAME MODEL, BUT WE ADD MORE YEARS.

## ---- 1. Re-Build model with extra dimensions ----

# HERE WE PREDICT FOR 6 EXTRA YEARS #

NExtraYears <- 6
nimConstants$Nyr <- nimConstants$Nyr+NExtraYears

nimData1 <- nimData

## WE NEED TO ADD EXTRA TIME DIMENSIONS TO ALL OBJECTS
yPredict <- array(0,c(dim(nimData$y)[1:3],NExtraYears ))
nimData1$y <- abind(nimData$y,yPredict,along=4) 

##SINCE WE ALSO ADD MORE YEARS, WE INCREASE AUGMENTATION TO MAKE SURE WE DONT RUN OUT OF IDS AVAILABLE. 
NEXTRA <- 500

yAUGM <- array(0,c(NEXTRA,dim(nimData1$y)[2:4]))
nimData1$y <- abind(nimData1$y,yAUGM, along = 1)

nimConstants$M <- nimConstants$M+NEXTRA

nimData1$sex <- c(nimData1$sex, rep(NA,NEXTRA))
nimData1$agePlusOne <- c(nimData1$agePlusOne, rep(NA,NEXTRA))
nimData1$w <- c(nimData1$w, rep(NA,NEXTRA))

# These values dont matter because detections are anyway 0? Only important to increase dimensions
# Dimensions are sometimes specific to each year. In those cases I augment with the maximun number times 2

nimConstants$lengthYCombined <- c(nimConstants$lengthYCombined, rep(max(nimConstants$MaxLocalTraps)*2,NExtraYears))

nimConstants$MaxLocalTraps <- c(nimConstants$MaxLocalTraps, rep(max(nimConstants$MaxLocalTraps)*2,NExtraYears))

effortAUG <- abind(nimConstants$effort,nimConstants$effort[,,1,], along = 3) # Effort data to augment (6 years, so 5+1, every year is the same)
nimConstants$effort <- abind(nimConstants$effort, effortAUG, along = 3)

prevcapAUG1 <- array(0,c(NEXTRA, dim(nimConstants$prevcap)[2:4])) # Extra individuals in current data-years
nimConstants$prevcap <- abind(nimConstants$prevcap,prevcapAUG1, along = 1)

prevcapAUG2 <- array(0,c(dim(nimData1$y)[1], dim(nimConstants$prevcap)[2:3], NExtraYears)) # Extra years + extra individuals
nimConstants$prevcap <- abind(nimConstants$prevcap,prevcapAUG2, along = 4)

localTrapsIndexAUG <- array(0,c(dim(nimData1$localTrapsIndex)[1:2],NExtraYears))
nimData1$localTrapsIndex <- abind(nimData1$localTrapsIndex,localTrapsIndexAUG, along = 3)

localTrapsNumAUG <- matrix(max(nimData1$localTrapsNum)*2, nrow = dim(nimData1$localTrapsNum)[1], ncol = NExtraYears)
nimData1$localTrapsNum <- abind(nimData1$localTrapsNum,localTrapsNumAUG, along = 2)

X.scAUG <- abind(nimData1$X.sc,nimData1$X.sc[,,1], along = 3)
nimData1$X.sc <- abind(nimData1$X.sc, X.scAUG, along = 3)


##NOW WE USE THE POSTERIOR DISTRIBUTION OBTAINED FROM THE FITTED MODEL AS PARAMETER VALUES#
# THEY WILL BE USED FOR SIMULATING THE TRAJECTORY OF THE POPULATION

# We just want to build a model object, so we take the posterior of only one iteration

ite <- 1

# WE SET NA FOR Z FOR YEARS TO PREDICT (Z is called U in our data)
zPredict <- array(NA,c(dim(nimData$u)[1], NExtraYears))
nimData1$u <- abind(myResultsSXYZ$sims.list$z[ite,,], zPredict, along=2)
# WE AUGMENT 
zAUGM <- array(1,c(NEXTRA,dim(nimData1$u)[2]))
zAUGM[,(nimConstants$Nyr-NExtraYears + 1):nimConstants$Nyr] <- NA

nimData1$u <- abind(nimData1$u, zAUGM, along=1)

# WE SET NA FOR s FOR YEARS TO PREDICT
sPredict <- array(NA, c(dim(myResultsSXYZ$sims.list$sxy)[2:3],NExtraYears ))
nimData1$s <- abind(myResultsSXYZ$sims.list$sxy[ite,,,], sPredict, along=3)
# WE AUGMENT 
sAUGM <- array(0, c(NEXTRA, dim(nimData1$s)[2:3]))
nimData1$s <- abind(nimData1$s, sAUGM, along=1)

##NOW WE BUILD THE MODEL AGAIN WITH THE NEW DATA AND THE YEARS TO PREDICT (EXTRA DIMENSIONS)###

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Model")
source('3.7.SCRopen_diftraps_Age_trapsBhCov_sigsexAge in Nimble.r')
modelCode <- SCRhab.Open.diftraps.age.effortTrapBhCov.sigsexage



modelSims <- nimbleModel( code = modelCode,
                          constants = nimConstants,
                          data = nimData1,
                          #inits = nimInits,
                          check = F,
                          calculate = F)

# Error: Maybe we need to make assumptions because we don't predict initial values???


##WE MAKE SOME ASSUMPTIONS FOR THE OTHER PARAMETERS#
# Now it doesn't matter, this will be substituted when we simulate the predictions

names(myResults$sims.list)
params

nimData1$phi.ad <- myResults$sims.list$phi.ad[ite]
nimData1$phi.sub <- myResults$sims.list$phi.sub[ite]
nimData1$phi.ad <- myResults$sims.list$phi.cub[ite]

nimData1$p.ad <- myResults$sims.list$p.ad[ite]
nimData1$p.sub <- myResults$sims.list$p.sub[ite]
nimData1$p.cub <- myResults$sims.list$p.cub[ite]

nimData1$psi <- myResults$sims.list$psi[ite]
nimData1$beta.dens <- myResults$sims.list$beta.dens[ite]
nimData1$omega <- myResults$sims.list$omega[ite]

nimData1$b.bh <- myResults$sims.list$b.bh[ite]
nimData1$trapBetas <- myResults$sims.list$trapBetas[ite,]

nimData1$sigma <- myResults$sims.list$sigma[ite,]
nimData1$sigD <- myResults$sims.list$sigD[ite]

### This two I need to augment them
nimData1$beta <- myResults$sims.list$beta[ite,] # This needs to sum to one
nimData1$piAGE <- myResults$sims.list$piAGE[ite,]


# HERE WE ASSUME THAT GAMMA IS 0.2 FOR THE YEARS TO PREDICT
nimData1$beta <- c(myResults$sims.list$beta[ite,], rep(0.2,NExtraYears))

sum(myResults$sims.list$beta[ite,])

beta <- myResults$sims.list$beta[ite,]
eta[1] <- beta[1]
for(k in 2:5){
  eta[k] <- beta[k]/(1-sum(beta[1:(k-1)]))
}

myResultsSXYZ$sims.list$z[ite,,4]

# We simulate the trajectory for each version of the reality (iteration)
##****Need to decide: Which is survival for the years to predict?


# we take random posteriors to use for simulations as parameter values for the first 5 years.  
itera <- sample(1:dim(myResultsSXYZ$sims.list$z)[1], 15) # We start by taking 15

