

## -------------------------------------------------
##          ProcessingBites from output Cyril 
##          1.1. Model parameters (params1)
##          1.2. Sxy, z, age.cat (params2)
## ------------------------------------------------- 

rm(list = ls())

library(coda)
library(MCMCvis)


#setwd("C:/Users/cymi/Downloads/AnaNimbleSCR")
setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/3.openSCRdenscov_Age/2021/Cyril/3-3.2")

# Number of iter ran: Bite size = 100: 1000 iterations per bite, thinned by 10
# 39 bites * 1000 iter / 10 = 3900


### ====  1.GET AND COMPILE BITES ====
# COMPILE CHARACTERISTICS 
bitesize <- 100
burnin <- 2500
NSkipBites <- burnin/bitesize

### ====    1.1 MODEL PARAMETERS ====

outDirectories <- list.files()[grep("NimbleOut", list.files())]
path.list <-outDirectories# file.path(myVars$WD, myVars$modelNameF, outDirectories)

# Retrieve the minimum number of bites per chain
numBites <- unlist(lapply(path.list, function(x){
  files <- list.files(x)
  files <- files[grep(".RData", files)]
  length(files)/2
}))
minBites <- floor(min(numBites))

## GO TROUGH THE BITES AND GET THEM
nimOutput <- RUNTIME <- list()
for(p in 1:length(path.list)){
  print(path.list[p])
  outfiles <- list.files(path.list[p])
  out  <- runtime <- list()#[CM]
  for(x in NSkipBites:minBites){
    print(x)
    load(file.path(path.list[p], paste("bite_", x, ".RData", sep = "")))
    runtime[[x]] <- RunTime[3]
    out[[x]] <- this.sample
  }#x
  RUNTIME[[p]] <- unlist(runtime)#[CM]
  out.mx <- do.call(rbind, out)
  
  nimOutput[[p]] <- as.mcmc(out.mx)
}#p


## COMPILE THE RESULTS
nimOutput <- as.mcmc.list(nimOutput)
summary(nimOutput)

setwd("D:/PhD/MyScripts_PhD/Ch. 2-3/Ch. 3/Results/Functions")
source("ProcessCodaOutput.R")

myResults <- ProcessCodaOutput(nimOutput,params.omit = c("sxy","z"))

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/3.openSCRdenscov_Age/2021/Cyril/3-3.2") 
save(nimOutput,myResults, file = "myResults_3-3.2_param.RData")

MCMCtrace(nimOutput,   
          ind = TRUE)

library(boot)
inv.logit(-6.995e-04) #p.sex
inv.logit(0.0029) #p.cub


### ====    1.2 SXY,Z,AGE.CAT ====

## GO TROUGH THE BITES AND GET THEM (SXY LOCATIONS)
nimOutputSXY <- RUNTIME <- list()
for(p in 1:length(path.list)){
  print(path.list[p])
  outfiles <- list.files(path.list[p])
  out  <- runtime <- list()#[CM]
  for(x in NSkipBites:minBites){
    print(x)
    load(file.path(path.list[p], paste("biteSxyZ_", x, ".RData", sep = "")))
    runtime[[x]] <- RunTime[3]
    out[[x]] <- this.sample2
  }#x
  RUNTIME[[p]] <- unlist(runtime)#[CM]
  out.mx <- do.call(rbind, out)
  
  nimOutputSXY[[p]] <- as.mcmc(out.mx)
}#p

## COMPILE THE RESULTS
nimOutputSXY <- as.mcmc.list(nimOutputSXY)

setwd("D:/PhD/MyScripts_PhD/Ch. 2-3/Ch. 3/Results/Functions")
source("ProcessCodaOutput.R")

myResultsSXYZ <- ProcessCodaOutput(nimOutputSXY)

# Export to process
setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/3.openSCRdenscov_Age/2021/Cyril/3-3.2")
save(myResultsSXYZ, nimOutputSXY, file = "myResults_3-3.2_sxy.RData")

