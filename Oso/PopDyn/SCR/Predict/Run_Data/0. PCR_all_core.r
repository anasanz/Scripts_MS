## -------------------------------------------------
##                       PCR
## ------------------------------------------------- 

# For projecting the population we need to estimate recruitment based on PCR
# We obtain PCR from our model results, but it can be calculated in two ways:
# 1. From all individuals: R/NºAdults
# 2. From individuals in core buffer: R in core/NºAdults in core

# Recruitment is the same all over the state space, but the number of adults is lower 
# in core area (buffer), therefore PCR will be higher in core area. 
# BEcause of how we subset abundance (keeping adults only in core), this PCR (core) would make more sense
# At the end we use that of the whole state space

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


M.aug <- nimConstants$M
Tt <- nimConstants$Nyr

# To constrain individuals inside, we need to unscale the sxy first.

# Load coordinates habitat grid scaled (G) -> Needed to unscale sxy
setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")
load("habcoord.RData")

dimnames(myResultsSXYZ$sims.list$sxy)[[3]] <- c('x','y')
myResultsSXYZ$sims.list$sxy <- scaleCoordsToHabitatGrid(coordsData = myResultsSXYZ$sims.list$sxy,## this are your sxy
                                                        coordsHabitatGridCenter = G,# this is your unscaled habitat (as you used when scaling the habitat/detector to the habitat. G?
                                                        scaleToGrid = FALSE)$coordsDataScaled

# Function to estimate pcr in the whole state space (pcr_core = FALSE) 
# or only in core buffer (pcr_core = TRUE)

calc.pcr <- function(pcr_core = TRUE){
  
  R <- matrix(NA, nrow(sampmat), Tt-1) # Number of recruits
  N.ad <- matrix(NA, nrow(sampmat), Tt-1)# Number of adults
  pcrmat <- matrix(NA, nrow(sampmat), Tt-1) # PCR matrix
  
  ZZad <- myResultsSXYZ$sims.list$z # To calculate recruitment per nº of adults: set z adults
  ZZad[!myResultsSXYZ$sims.list$age.cat %in% c(5) ]  <- 3 # Considers all individuals that are not 5 (ADULTS) as dead
  
  for (ite in 1:nrow(sampmat)){
    
    for (t in 2:Tt){
      
      zwt <- paste('z[', 1:M.aug,', ', t, ']', sep='' )
      zwtm <- paste('z[', 1:M.aug,', ', t-1, ']', sep='' )
      
      if(pcr_core == FALSE){ # If pcr is calculated as Nº recruits in all state space/Nº Adults in all state space
        
        # N recruits: From year t-1 to t
        R[ite,t-1] <- sum(ifelse((sampmat[ite,zwt]-sampmat[ite,zwtm])==1, 1, 0)) # Number of (ALL) recruits in t-1
        
        # N adults: Important, on year t-1!
        N.ad[ite,t-1] <- length(which(ZZad[ite,,t-1]==1)) # Number of (ALL) adults alive in t-1
        
        # PCR
        pcrmat[ite,t-1] <- R[ite,t-1]/N.ad[ite,t-1]
        
      } else { # If pcr is calculated as Nº recruits in core buffer/Nº Adults in core buffer
        
        # N recruits: From year t-1 to t
        
        sp.t <- SpatialPoints(myResultsSXYZ$sims.list$sxy[ite,,,t],proj4string=CRS(proj4string(Xbuf2))) # Activity centers of t. Convert to spatial
        sp.tm <- SpatialPoints(myResultsSXYZ$sims.list$sxy[ite,,,t-1],proj4string=CRS(proj4string(Xbuf2))) # Activity centers of t-1. Convert to spatial
        
        which.In.t <- over(sp.t, Xbuf2) # AC inside core buffer
        which.In.tm <- over(sp.tm, Xbuf2)
        
        t.in <- sampmat[ite,zwt] # z of t and t - 1
        tm.in <- sampmat[ite,zwtm]
        
        t.in[is.na(which.In.t)] <- 0  # All outside buffer(NA) considered dead
        tm.in[is.na(which.In.tm)] <- 0 
        
        R[ite,t-1]<-sum(ifelse((t.in-tm.in)==1, 1, 0))
        
        # N adults: Important, on year t-1!
        
        which.alive <- which(ZZad[ite,,t-1]==1) # Which adults are alive
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t-1] # Retrieve the activity center for those individuals
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(Xbuf2))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, Xbuf2) # Which inside core buffer
        
        N.ad[ite,t-1] <- sum(which.In,na.rm = T) # Number of (ALL) adults in t-1
        
        # PCR
        pcrmat[ite,t-1] <- R[ite,t-1]/N.ad[ite,t-1]
        
      } # pcr_core
    } # t
  }# ite
  
  results <- list(R,N.ad,pcrmat)
  names(results) <- c("R", "N.ad", "pcrmat")
  
  return(results)
  
}

# Estimate pcr 

pcr_corebuf <- calc.pcr(pcr_core = TRUE)
pcr_all <- calc.pcr(pcr_core = FALSE)

# Save pcr because it takes forever to calculate
setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL")
save(pcr_corebuf, file = "pcr_corebuf.RData")
save(pcr_all, file = "pcr_all.RData")

# Create new matrix of pcr to show all possibilities: Containts absolute maximum, absolute minimum, mean pcr in core

pcr1_corebuf <- apply(pcr_corebuf[[3]],1,mean) # ASP: Mean per capita recruitment per iteration
pcr2_all <- apply(pcr_all[[3]],1,mean) # ASP: Mean per capita recruitment per iteration

pcr_range_min <- NULL
for (i in 1:length(pcr1_corebuf)) {
  pcr_range_min <- c(pcr_range_min,min(pcr1_corebuf[i], pcr2_all[i]))
}

pcr_range_max <- NULL
for (i in 1:length(pcr1_corebuf)) {
  pcr_range_max <- c(pcr_range_max,max(pcr1_corebuf[i], pcr2_all[i]))
}

pcr_range <- c(pcr_range_min, pcr_range_max)

summary(pcr1_corebuf)
summary(pcr2_all)

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL")
save(pcr_range, file = "pcr_range.RData")

## -------------------------------------------------
##                 PCR females
## ------------------------------------------------- 
# Modify function to estimate PCR only for females

calc.pcr.fem <- function(pcr_core = TRUE){
  
  R <- matrix(NA, nrow(sampmat), Tt-1) # Number of recruits
  N.ad <- matrix(NA, nrow(sampmat), Tt-1)# Number of adults
  pcrmat <- matrix(NA, nrow(sampmat), Tt-1) # PCR matrix
  
  ZZad <- myResultsSXYZ$sims.list$z # To calculate recruitment per nº of adults: set z adults
  ZZad[!myResultsSXYZ$sims.list$age.cat %in% c(5) ]  <- 3 # Considers all individuals that are not 5 (ADULTS) as dead
  ZZad[!myResultsSXYZ$sims.list$sex %in% c(0) ]  <- 3 # Considers all individuals that are not 0 (females) as dead
  
  for (ite in 1:nrow(sampmat)){
    
    for (t in 2:Tt){
      
      zwt <- paste('z[', 1:M.aug,', ', t, ']', sep='' )
      zwtm <- paste('z[', 1:M.aug,', ', t-1, ']', sep='' )
      
      if(pcr_core == FALSE){ # If pcr is calculated as Nº recruits in all state space/Nº Adults in all state space
        
        # N recruits: From year t-1 to t
        R[ite,t-1] <- sum(ifelse((sampmat[ite,zwt]-sampmat[ite,zwtm])==1, 1, 0)) # Number of (ALL) recruits in t-1
        
        # N adults: Important, on year t-1!
        N.ad[ite,t-1] <- length(which(ZZad[ite,,t-1]==1)) # Number of (ALL) adult FEMALES alive in t-1
        
        # PCR
        pcrmat[ite,t-1] <- R[ite,t-1]/N.ad[ite,t-1]
        
      } else { # If pcr is calculated as Nº recruits in core buffer/Nº Adults in core buffer
        
        # N recruits: From year t-1 to t
        
        sp.t <- SpatialPoints(myResultsSXYZ$sims.list$sxy[ite,,,t],proj4string=CRS(proj4string(Xbuf2))) # Activity centers of t. Convert to spatial
        sp.tm <- SpatialPoints(myResultsSXYZ$sims.list$sxy[ite,,,t-1],proj4string=CRS(proj4string(Xbuf2))) # Activity centers of t-1. Convert to spatial
        
        which.In.t <- over(sp.t, Xbuf2) # AC inside core buffer
        which.In.tm <- over(sp.tm, Xbuf2)
        
        t.in <- sampmat[ite,zwt] # z of t and t - 1
        tm.in <- sampmat[ite,zwtm]
        
        t.in[is.na(which.In.t)] <- 0  # All outside buffer(NA) considered dead
        tm.in[is.na(which.In.tm)] <- 0 
        
        R[ite,t-1]<-sum(ifelse((t.in-tm.in)==1, 1, 0))
        
        # N adults: Important, on year t-1!
        
        which.alive <- which(ZZad[ite,,t-1]==1) # Which adult FEMALES are alive
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t-1] # Retrieve the activity center for those individuals
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(Xbuf2))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, Xbuf2) # Which inside core buffer
        
        N.ad[ite,t-1] <- sum(which.In,na.rm = T) # Number of (ALL) adults in t-1
        
        # PCR
        pcrmat[ite,t-1] <- R[ite,t-1]/N.ad[ite,t-1]
        
      } # pcr_core
    } # t
  }# ite
  
  results <- list(R,N.ad,pcrmat)
  names(results) <- c("R", "N.ad", "pcrmat")
  
  return(results)
  
}

# Estimate pcr 

pcr_corebuf_fem <- calc.pcr.fem(pcr_core = TRUE)
pcr_all_fem <- calc.pcr.fem(pcr_core = FALSE)

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL")
save(pcr_corebuf_fem, file = "pcr_corebuf_fem.RData")
save(pcr_all_fem, file = "pcr_all_fem.RData")

load("pcr_corebuf_fem.RData")
load("pcr_all_fem.RData")

load("pcr_corebuf.RData")
load("pcr_all.RData")

# Compare pcr females/all individuals

colMeans(pcr_all_fem$pcrmat)
colMeans(pcr_all$pcrmat)

colMeans(pcr_corebuf_fem$pcrmat)
colMeans(pcr_corebuf$pcrmat)

# Conclusion 1: PCR is higher when estimated per female (probably because there are a lower number of females??)

# Now chech the uncertainty in PCR

apply(pcr_all_fem$pcrmat,2,sd)
apply(pcr_all$pcrmat,2,sd)

apply(pcr_corebuf_fem$pcrmat,2,sd)
apply(pcr_corebuf$pcrmat,2,sd)

# Conclusion 2: Uncertainty in PCR is higher when estimated per female. Probably because uncertainty in female abundance is higher? CHECK

## COMPARE UNCERTAINTY IN NUMBER OF FEMALES vS MALES, to prove if there is more uncertainty

# 1. Nº females

ZZadFEM <- myResultsSXYZ$sims.list$z
ZZadFEM[!myResultsSXYZ$sims.list$age.cat %in% c(5) ]  <- 3 # Considers all individuals that are not 5 (ADULTS) as dead
ZZadFEM[!myResultsSXYZ$sims.list$sex %in% c(0) ]  <- 3 # Considers all individuals that are not 0 (females) as dead

adFEM_trapBuf <- matrix(NA,nrow=dim(myResultsSXYZ$sims.list$z)[1],ncol=dim(myResultsSXYZ$sims.list$z)[3]) # nrow = iterations, ncol = year

for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
  for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
    
    which.alive <- which(ZZadFEM[ite,,t]==1) # Select only the individuals alive (z=1)

    which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals

    sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(Xbuf2))) # CONVERT SXY TO SPATIAL POINTS 
    which.In <- over(sp, Xbuf2) # Check which ones are in the buffer
    
    adFEM_trapBuf[ite,t] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
  }}

# 2. Nº males

ZZadMAL <- myResultsSXYZ$sims.list$z
ZZadMAL[!myResultsSXYZ$sims.list$age.cat %in% c(5) ]  <- 3 # Considers all individuals that are not 5 (SADULTS) as dead
ZZadMAL[!myResultsSXYZ$sims.list$sex %in% c(1) ]  <- 3 # Considers all individuals that are not 1 (males) as dead

adMAL_trapBuf <- matrix(NA,nrow=dim(myResultsSXYZ$sims.list$z)[1],ncol=dim(myResultsSXYZ$sims.list$z)[3]) # nrow = iterations, ncol = year

for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
  for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
    
    which.alive <- which(ZZadMAL[ite,,t]==1) # Select only the individuals alive (z=1)
    
    which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
    
    sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(Xbuf2))) # CONVERT SXY TO SPATIAL POINTS 
    which.In <- over(sp, Xbuf2) # Check which ones are in the buffer
    
    adMAL_trapBuf[ite,t] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
  }}

#average number of individuals without the buffer for each year 
colMeans(adFEM_trapBuf)
colMeans(adMAL_trapBuf)

apply(adFEM_trapBuf,2,sd)
apply(adMAL_trapBuf,2,sd)


par(mfrow = c(2,3))
for (t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
  plot(density(adFEM_trapBuf[,t]), main = t)
  abline(v = colMeans(adFEM_trapBuf)[t], col = "blue")
}

par(mfrow = c(2,3))
for (t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
  plot(density(adMAL_trapBuf[,t]), main = t)
  abline(v = colMeans(adMAL_trapBuf)[t], col = "blue")
}


# Sum of individuals alive in total each year (without buffer)
NALL_fem <- apply(ZZadFEM,c(1,3),function(x) sum(x==1))
colMeans(NALL_fem)

NALL_mal <- apply(ZZadMAL,c(1,3),function(x) sum(x==1))
colMeans(NALL_mal)

apply(NALL_fem,2,sd)
apply(NALL_mal,2,sd)

# Coonclusion 3: Uncertainty is higher in the number of estimated females
# That is why uncertainty in PCR is also higher

#Conclusion 4: Sex ratio is skewed to females!! Maybe investigate further?
