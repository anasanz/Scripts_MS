## -------------------------------------------------
##                 Predict future population
##                      Data script 3.1
##                         Model 3.7
## ------------------------------------------------- 

## Our abundance estimates are based in subsetting our population to our trapping area (core buffer)
## To predict abundance non-spatially, we do not have the activity centers in future years, so
## we can't subset it a posteriori and we need use as starting point (year 5) the individuals
## inside the core area.
## To predict, we use the pcr rate. We need to choose if we:
# 1. Estimate pcr for all individuals
# 2. Estimate pcr from individuals in core buffer: R in core/NºAdults in core


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
dim(sampmat1)

which((myResults$sims.list$N == sampmat1[,c(6:10)]) == FALSE) # With this I check that myResults and nimOutput really correspond to each other

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL")
load("myResults_3-3.1_sxy.RData")
sampmat2 <- do.call(rbind, nimOutputSXY)

sampmat <- cbind(sampmat1, sampmat2)

M.aug <- nimConstants$M
Tt <- nimConstants$Nyr

## Unscale sxy coordinates to subset only individuals within core buffer

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")
load("habcoord.RData")

dimnames(myResultsSXYZ$sims.list$sxy)[[3]] <- c('x','y')
myResultsSXYZ$sims.list$sxy <- scaleCoordsToHabitatGrid(coordsData = myResultsSXYZ$sims.list$sxy,## this are your sxy
                                                        coordsHabitatGridCenter = G,# this is your unscaled habitat (as you used when scaling the habitat/detector to the habitat. G?
                                                        scaleToGrid = FALSE)$coordsDataScaled

start.sxy <- myResultsSXYZ$sims.list$sxy[,,,5]

## ---- 1. Set dimensions projection matrix ----

niter <- nrow(sampmat)
t.new <- 6 #6 yrs into future
max.age <- 5 # BEcause we will project into age categories
age.cat <- 3 # age categories to keep track off
Narray <- array(NA, c(age.cat, t.new, niter)) # ASP: For each iteration, project abundance per each category and year

## ---- 2. Get starting values (Values of year 1, so 5th year of our data) ----

# Z

##z, age at T=5 from model
z.which <- grep('z\\[', colnames(sampmat)) # ASP: index columns all z (sampmat matrix)
##ordered: all individuals for yr 1, then all for yr2, etc
z.est <- array(sampmat[,z.which], c(niter,M.aug, Tt))
z.nosuper <- apply(z.est,1:2,sum) == 0 # ASP: Across all years if that individual in that iteration was never detected (sum of z = 0), it was never part of the super pop

# AGE

#age.which<-pmatch(paste('age[', 1:Maug, ', ', Tt,']', sep=''), colnames(sampmat))
age.which <- grep('age.cat\\[', colnames(sampmat))
age.est <- array(sampmat[,age.which], c(niter,M.aug, Tt))
age.est[z.est==0]<-0 # ASP: If is not alive it doesn't have an age (how age model is set up)
age.est.start <- age.est[,,5]

# Count only the individuals that are inside core buffer
# (creating new category -1, which means outside)

for(ite in 1:nrow(sampmat)){
  sp <- SpatialPoints(start.sxy[ite,,],proj4string=CRS(proj4string(Xbuf2))) # Convert to spatial
  which.In <- over(sp, Xbuf2) # Check which ones are in the buffer
  age.est.start[ite,][is.na(which.In)] <- -1 # This -1 includes all age categories, unrecruited or dead individuals that are outside the core area
  }

##get distribution of ages in each iteration at t=5
dim(age.est[,,5])
dim(age.est.start)

ttt <- do.call('rbind', apply(age.est.start, 1, table))

N.age <- matrix(NA, niter, 3) # ASP: Organize into age categories
for (i in 1:niter){
  N.age[i,]<-c(sum(ttt[i,3:4]), sum(ttt[i,5:6]), sum(ttt[i,7:ncol(ttt)]))
} 

Narray[,1,]<-t(N.age) # ASP: Fill N at each age class the first year for each iteration


## ---- 3. Calculate per capita recruitment ----

# Function to estimate pcr in the whole state space (pcr_core = FALSE) 
# or only in core buffer (pcr_core = TRUE)

#calc.pcr <- function(pcr_core = TRUE){
  
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
  
  return(results)
  
}
#
#pcr_corebuf1 <- calc.pcr(pcr_core = TRUE)
#pcr1_corebuf1 <- apply(pcr_corebuf1[[3]],1,mean) # ASP: Mean per capita recruitment per iteration

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL")
load("pcr_corebuf.RData")
load("pcr_all.RData")

pcr1_corebuf <- apply(pcr_corebuf[[3]],1,mean) # ASP: Mean per capita recruitment per iteration
pcr2_all <- apply(pcr_all[[3]],1,mean) # ASP: Mean per capita recruitment per iteration


## ---- 4. Projection, for each iteration ----

itera <- sample(1:nrow(sampmat), 20)

## ---- 4.1. PCR inside core buffer ----

pcr <- pcr1_corebuf

for(ite in 1:length(itera)){
  ageN <- matrix(NA, t.new, max.age) # ASP: Number of individuals of each age per year
  ###make vectors of survival
  phi.vec <- c(rep(sampmat[itera[ite], 'phi.cub'], 2),
             rep(sampmat[itera[ite], 'phi.sub'], 2),
             rep(sampmat[itera[ite], 'phi.ad'], 1)) # ASP: Survival corresponding to each age
  
  ##starting population, ages 
  #ageN[1,] <- c(ttt[itera[ite],2:ncol(ttt)], rep(0, max.age-ncol(ttt)+1)) # ASP: The first row, year 1
  ageN[1,] <- ttt[itera[ite],3:ncol(ttt)] # ASP: The first row, year 1
  
  
  ##move forward
  for (t in 2:t.new){
    ##survivors
    S <- rbinom(length(ageN[t-1,]), ageN[t-1,], phi.vec) #ASP: Simulate 5 values (ages). For the number of individuals of each age, how many survived according to phi vec (their age specific survival)
                                  ### QUESTION, IS THIS RIGHT FOR AGE CATEGORIES OR SHOULD I GET EVERYTHING INTO AGES????? NOT MONITORED
    # if(S[max.age]>0) stop('Up max.age!!') # ASP: Doesn't work with age categories
    
    ##recruits (N at t-1 times pcr) - SHOULD BE ADULTS ONLY!!!
    # B<-rpois(1, pcr[itera[ite]]*sum(Narray[1:3,t-1,itera[ite]]))
    #binomial, assuming there are always 500 available
    
    # ASP: Here modification from simulation to estimate the numer of recruits from the number of adults
    # B <- rbinom(1, 500, (pcr[itera[ite]]*sum(Narray[1:3,t-1,itera[ite]]))/500) # ASP: Number of recruits each year (pcr*nºAdults previous year). 0 age category
    B <- rbinom(1, 500, (pcr[itera[ite]]*Narray[3,t-1,itera[ite]])/500) # ASP: Instead of the sum (all individuals, I take the third value (adults))
    
    aging.ad <- c(S[1:(max.age-2)], sum(S[(max.age-1):max.age]))
    
    ageN[t,]<-c(B, aging.ad) # ASP: The recruits pass to be the cubs of year 1 at time t -> al rest AGING
    
    Narray[,t,itera[ite]]<-c(sum(ageN[t,1:2]),sum(ageN[t,3:4]), ageN[t,5] ) # ASP: Sum the N of the different ages (into age categories)
    
  }
}

subset_N <- Narray[,,itera]
Nall_core <- apply(subset_N, 2:3, sum)

par(mfrow = c(1,2))

plot(1:6, apply(Nall_core[1:6,], 1,mean), type='l', ylim=c(0,400), main = "Core")
for (ite in 1:length(itera)){
  points(1:6, Nall_core[1:6,ite], type='l', col='lightgrey')
}
points(1:6,apply(Nall_core, 1,mean), type='l')

## ---- 4.2. PCR ALL individuals ----

pcr <- pcr2_all

for(ite in 1:length(itera)){
  ageN <- matrix(NA, t.new, max.age) # ASP: Number of individuals of each age per year
  ###make vectors of survival
  phi.vec <- c(rep(sampmat[itera[ite], 'phi.cub'], 2),
               rep(sampmat[itera[ite], 'phi.sub'], 2),
               rep(sampmat[itera[ite], 'phi.ad'], 1)) # ASP: Survival corresponding to each age
  
  ##starting population, ages 
  #ageN[1,] <- c(ttt[itera[ite],2:ncol(ttt)], rep(0, max.age-ncol(ttt)+1)) # ASP: The first row, year 1
  ageN[1,] <- ttt[itera[ite],3:ncol(ttt)] # ASP: The first row, year 1
  
  
  ##move forward
  for (t in 2:t.new){
    ##survivors
    S <- rbinom(length(ageN[t-1,]), ageN[t-1,], phi.vec) #ASP: Simulate 5 values (ages). For the number of individuals of each age, how many survived according to phi vec (their age specific survival)
    ### QUESTION, IS THIS RIGHT FOR AGE CATEGORIES OR SHOULD I GET EVERYTHING INTO AGES????? NOT MONITORED
    # if(S[max.age]>0) stop('Up max.age!!') # ASP: Doesn't work with age categories
    
    ##recruits (N at t-1 times pcr) - SHOULD BE ADULTS ONLY!!!
    # B<-rpois(1, pcr[itera[ite]]*sum(Narray[1:3,t-1,itera[ite]]))
    #binomial, assuming there are always 500 available
    
    # ASP: Here modification from simulation to estimate the numer of recruits from the number of adults
    # B <- rbinom(1, 500, (pcr[itera[ite]]*sum(Narray[1:3,t-1,itera[ite]]))/500) # ASP: Number of recruits each year (pcr*nºAdults previous year). 0 age category
    B <- rbinom(1, 500, (pcr[itera[ite]]*Narray[3,t-1,itera[ite]])/500) # ASP: Instead of the sum (all individuals, I take the third value (adults))
    
    aging.ad <- c(S[1:(max.age-2)], sum(S[(max.age-1):max.age]))
    
    ageN[t,]<-c(B, aging.ad) # ASP: The recruits pass to be the cubs of year 1 at time t -> al rest AGING
    
    Narray[,t,itera[ite]]<-c(sum(ageN[t,1:2]),sum(ageN[t,3:4]), ageN[t,5] ) # ASP: Sum the N of the different ages (into age categories)
    
  }
}

subset_N <- Narray[,,itera]
Nall_all <- apply(subset_N, 2:3, sum)

plot(1:6, apply(Nall_all[1:6,], 1,mean), type='l', ylim=c(0,400), main = "All")
for (ite in 1:length(itera)){
  points(1:6, Nall_all[1:6,ite], type='l', col='lightgrey')
}
points(1:6,apply(Nall_all, 1,mean), type='l')

## Compare:
apply(Nall_all[1:6,], 1,mean)
apply(Nall_core[1:6,], 1,mean)
mean(Nall_core[1,]) # N in year 5 (year 1 of projection) is 66.5 (different from what we obtained) because we are only averaging 20 iterations, not all
# When using only the individuals in the core to estimate pcr(Nº adults core)/Nº recruits in core, it is much higher
# Probably because recruitment stays the same within core and within whole area (all of it happens at the core)
# BUT the number of adults is much higher in whole area, which makes pcr smaller.
# Accordying to this, pcr better at the core? It is about the nº of adults to choose. And we are using the core to give abundance and as starting point

# Explore pcr

plot(1:4, apply(pcr_corebuf[[3]], 2,mean), type='l', ylim=c(0,3), main = "Core")
for (ite in 1:niter){
  points(1:4, pcr_corebuf[[3]][ite,], type='l', col='lightgrey')
}
points(1:4, apply(pcr_corebuf[[3]], 2,mean), type='l')

plot(1:4, apply(pcr_all[[3]], 2,mean), type='l', ylim=c(0,3), main = "All")
for (ite in 1:niter){
  points(1:4, pcr_all[[3]][ite,], type='l', col='lightgrey')
}
points(1:4, apply(pcr_all[[3]], 2,mean), type='l')

apply(pcr_corebuf[[3]], 2,mean)
apply(pcr_all[[3]], 2,mean)

