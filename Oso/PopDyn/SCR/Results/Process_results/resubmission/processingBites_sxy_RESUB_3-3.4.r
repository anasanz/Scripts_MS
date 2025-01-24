## -------------------------------------------------
##          ProcessingBites from output Cyril 
##          1.1. Model parameters (params1)
##          1.2. Sxy, z, age.cat (params2)
## ------------------------------------------------- 

rm(list = ls())

library(coda)
library(MCMCvis)


#setwd("C:/Users/cymi/Downloads/AnaNimbleSCR")
setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/RESUBMISSION/3-3.4_allparams")

# Number of iter ran: Bite size = 100: 1000 iterations per bite, thinned by 10
# 39 bites * 1000 iter / 10 = 3900


### ====  1.GET AND COMPILE BITES ====
# COMPILE CHARACTERISTICS 
bitesize <- 100
burnin <- 1000
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
  for(x in NSkipBites:minBites){ # ASP: Here Im skiping the first 1000 iterations (10 bites * 100 iter = 1000 iter)
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
sum <- summary(nimOutput)
sum$statistics

setwd("D:/PhD/MyScripts_PhD/Ch. 2-3/Ch. 3/Results/Functions")
source("ProcessCodaOutput.R")

myResults <- ProcessCodaOutput(nimOutput,params.omit = c("sxy","z"))

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/RESUBMISSION/3-3.4_allparams")
save(nimOutput,myResults, file = "myResults_RESUB_3-3.4_param.RData")

MCMCtrace(nimOutput,   
          ind = TRUE,
          pdf = FALSE)

library(boot)
inv.logit(0.004) # Adults
inv.logit(0.0009198) # Subadults

### ====    1.2 SXY,Z,AGE.CAT ====

# COMPILE CHARACTERISTICS 
#bitesize <- 50 # ASP: Atención, bites are already thinned :( :()
#burnin <- 1000
#NSkipBites <- burnin/bitesize

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
setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/RESUBMISSION/3-3.4_allparams")
save(myResultsSXYZ, nimOutputSXY, file = "myResults_RESUB_3-3.4_sxy.RData")

