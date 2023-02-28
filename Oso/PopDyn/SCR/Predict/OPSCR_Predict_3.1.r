
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
sampmat1 <- do.call(rbind, nimOutput)

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1")
load("myResults_3-3.1_sxy.RData")
sampmat2 <- do.call(rbind, nimOutputSXY)


## ---- 1. Re-Build model with extra dimensions ----

### I re-build model structure by taking values of one iteration

## ---- 1.1. Create starting z's, age's and sxy's ----

# These will be based on model output (last year), new M and new time frame
t.new <- 5 # sim 5 addl years
M.new <- 400 # new augmentation limit

# Use estimates of z, age etc from last year of estimation model only

# Z

z.start <- matrix(NA, M.new, 1+t.new)

##fill z.start with estimates of z form fitted model; set addl individuals in years
##covered by model to 0 - known NOT to be alive
#first, make matrix of z estimates
z.which1 <- pmatch(paste('z[', 1:nimConstants$M, ', ', nimConstants$Nyr,']', sep=''), colnames(sampmat2)) # ASP: index columns z for LAST YEAR (sampmat2 matrix)
z.which <- grep('z\\[', colnames(sampmat2)) # ASP: index columns all z (sampmat2 matrix)
##ordered: all individuals for yr 1, then all for yr2, etc
z.est <- matrix(sampmat2[1, z.which], nimConstants$M, nimConstants$Nyr) # ASP: All z from estimation model
z.start[1:nimConstants$M, 1] <- sampmat2[1,z.which1] # ASP: Starting z for projection: z of year 5 
z.start[(nimConstants$M+1):M.new, 1] <- 0 # ASP: Augmented as not yet recruited

# Repeat this for: AGE

age.start <- matrix(NA, M.new, 1+t.new)

age.which <- pmatch(paste('age.cat[', 1:nimConstants$M, ', ', nimConstants$Nyr,']', sep=''), colnames(sampmat2)) # ASP: index columns AGE for LAST YEAR
age.est <- sampmat2[1,age.which] # ASP: Ages last year
##set everyone not alive at all to all-0, ie, can be recruited in the NEW recruitment 
##model (they were never part of the superpopulation -> ASP: because of how the model is formulated z=u*w)
z.nosuper <- which(apply(z.est,1,sum)==0)
age.est[z.nosuper] <- 0
#age.est[age.est>5] <- 5 #set to max age category. ASP: NO NEED OF THIS HERE, already age.cat
age.start[1:nimConstants$M, 1] <- age.est # ASP: Starting AGE for projection: AGE of year 5 
age.start[(nimConstants$M+1):M.new, 1] <- 0 # ASP: Augmented as not yet recruited

# Repeat this for: SXY

sxy.start <- array(NA, c(M.new, 2, 1+t.new))

##make array of activity centers from model output
s.which <- grep('sxy', colnames(sampmat2)) # ASP: index columns all sxy (sampmat matrix)
##all individuals X Yr1, then all inds Y Yr1; all inds X yr 2, then y yr 2 etc ## ASP: This is the order the sxy vector is placed (s.which)
indx <- (nimConstants$Nyr*nimConstants$M*2-(nimConstants$M*2-1)):(nimConstants$Nyr*nimConstants$M*2) # ASP: Index the LAST year of sxy. *2 because there are the double of columns: x and y. (nimConstants$M*2-1) corresponds to 1 year of data
sxy.start[1:nimConstants$M,,1]<-matrix(sampmat2[ 1, s.which[indx] ], nimConstants$M, 2) # ASP: Starting SXY for projection: SXY of year 5 

## ---- 1.2. Calculate per capita recruitment ----





























# HERE WE PREDICT FOR 6 EXTRA YEARS #

NExtraYears <- 6
nimConstants$Nyr <- nimConstants$Nyr+NExtraYears

nimData1 <- nimData

## WE NEED TO ADD EXTRA TIME DIMENSIONS TO ALL OBJECTS
#yPredict <- array(0,c(dim(nimData$y)[1:3],NExtraYears ))
#nimData1$y <- abind(nimData$y,yPredict,along=4) 

##SINCE WE ALSO ADD MORE YEARS, WE INCREASE AUGMENTATION TO MAKE SURE WE DONT RUN OUT OF IDS AVAILABLE. 
NEXTRA <- 500

#yAUGM <- array(0,c(NEXTRA,dim(nimData1$y)[2:4]))
#nimData1$y <- abind(nimData1$y,yAUGM, along = 1)

nimConstants$M <- nimConstants$M+NEXTRA

nimData1$sex <- c(nimData1$sex, rep(NA,NEXTRA))
nimData1$agePlusOne <- c(nimData1$agePlusOne, rep(NA,NEXTRA))
nimData1$w <- c(nimData1$w, rep(NA,NEXTRA))
nimData1$b <- c(nimData1$b, rep(1,NExtraYears))

# These values dont matter because detections are anyway 0? Only important to increase dimensions
# Dimensions are sometimes specific to each year. In those cases I augment with the maximun number times 2

#nimConstants$lengthYCombined <- c(nimConstants$lengthYCombined, rep(max(nimConstants$MaxLocalTraps)*2,NExtraYears))
#
#nimConstants$MaxLocalTraps <- c(nimConstants$MaxLocalTraps, rep(max(nimConstants$MaxLocalTraps)*2,NExtraYears))
#
#effortAUG <- abind(nimConstants$effort,nimConstants$effort[,,1,], along = 3) # Effort data to augment (6 years, so 5+1, every year is the same)
#nimConstants$effort <- abind(nimConstants$effort, effortAUG, along = 3)
#
#prevcapAUG1 <- array(0,c(NEXTRA, dim(nimConstants$prevcap)[2:4])) # Extra individuals in current data-years
#nimConstants$prevcap <- abind(nimConstants$prevcap,prevcapAUG1, along = 1)
#
#prevcapAUG2 <- array(0,c(dim(nimData1$y)[1], dim(nimConstants$prevcap)[2:3], NExtraYears)) # Extra years + extra individuals
#nimConstants$prevcap <- abind(nimConstants$prevcap,prevcapAUG2, along = 4)
#
#localTrapsIndexAUG <- array(0,c(dim(nimData1$localTrapsIndex)[1:2],NExtraYears))
#nimData1$localTrapsIndex <- abind(nimData1$localTrapsIndex,localTrapsIndexAUG, along = 3)
#
#localTrapsNumAUG <- matrix(max(nimData1$localTrapsNum)*2, nrow = dim(nimData1$localTrapsNum)[1], ncol = NExtraYears)
#nimData1$localTrapsNum <- abind(nimData1$localTrapsNum,localTrapsNumAUG, along = 2)
#
#X.scAUG <- abind(nimData1$X.sc,nimData1$X.sc[,,1], along = 3)
#nimData1$X.sc <- abind(nimData1$X.sc, X.scAUG, along = 3)


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
nimData1$sxy <- abind(myResultsSXYZ$sims.list$sxy[ite,,,], sPredict, along=3)
# WE AUGMENT 
sAUGM <- array(0, c(NEXTRA, dim(nimData1$sxy)[2:3]))
nimData1$sxy <- abind(nimData1$sxy, sAUGM, along=1)

##NOW WE BUILD THE MODEL AGAIN WITH THE NEW DATA AND THE YEARS TO PREDICT (EXTRA DIMENSIONS)###

nimData1[['y']] <- NULL


setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Model")
source('3.7.SCRopen_PREDICT in Nimble.r')
modelCode <- SCRhab.Open.diftraps.age.effortTrapBhCov.sigsexage

nimData1[['y']] <- NULL
nimData1[['ones']] <- NULL
nimData1[['X.sc']] <- NULL
nimData1[['habitatGridDet']] <- NULL
nimData1[['localTrapsIndex']] <- NULL
nimData1[['localTrapsNum']] <- NULL


modelSims <- nimbleModel( code = modelCode,
                          constants = nimConstants,
                          data = nimData1,
                          #inits = nimInits,
                          check = F,
                          calculate = F)

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

