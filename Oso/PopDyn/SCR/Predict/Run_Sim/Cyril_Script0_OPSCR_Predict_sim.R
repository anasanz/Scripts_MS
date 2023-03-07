
rm(list=ls())
## Load packages
library(nimble)
library(nimbleSCR)
library(basicMCMCplots)
library(R.utils)
library(abind)
library(coda)
# 
# source("C:/Personal_Cloud/OneDrive/Work/PublicCodesGit/Public/SpatialSurvivalOPSCR/sampler_categorical_general1.R")
# source("C:/Personal_Cloud/OneDrive/Work/PublicCodesGit/Public/SpatialSurvivalOPSCR/dcatState2Dead.R")

## ------ SET REQUIRED WORKING DIRECTORIES ------
#source("C:/My_documents/rovquant/analyses/Rgit/RovQuant/Temp/CM/myWorkingDirectories.R")
#source("C:/My_documents/RovQuant/Temp/PD/myWorkingDirectories.R")             
#source("C:/PROJECTS/RovQuant/Temp/RB/myWorkingDirectories.R")   

## ------ SOURCE THE REQUIRED FUNCTIONS ------
#sourceDirectory(dir.function, modifiedOnly = FALSE)

source("D:/MargSalas/Scripts_MS/Functions/ProcessCodaOutput.R")
source("D:/MargSalas/Scripts_MS/Functions/plot.violins2.r")
source("D:/MargSalas/Scripts_MS/Functions/DoScale.r")


## Create habitat grid
coordsHabitatGridCenter <- cbind(rep(seq(23, 8, by = -5), 4),
                                 sort(rep(seq(8, 23, by = 5), 4)))
colnames(coordsHabitatGridCenter) <- c("x","y")

## Create trap grid
coordsObsCenter <- cbind(rep(seq(8, 23, by = 1), 16),
                         sort(rep(seq(23, 8, by = -1), 16)))
colnames(coordsObsCenter) <- c("x","y")

## Plot check
plot(coordsHabitatGridCenter[,"y"] ~ coordsHabitatGridCenter[,"x"],
     xlim = c(5,25), ylim = c(5,25),
     pch = 1, cex = 1.5) 
points(coordsObsCenter[,"y"] ~ coordsObsCenter[,"x"], col="red", pch=16 ) 
par(xpd=TRUE)
legend(x = 7, y = 7,
       legend=c("Habitat window centers", "Observation window centers"),
       pt.cex = c(1.5,1),
       horiz = T,
       pch=c(1,16),
       col=c("black", "red"),
       bty = 'n')

       
habitatMask <- matrix(1, nrow = 4, ncol= 4, byrow = TRUE)

### 1.2 Rescale coordinates

## Rescale coordinates
scaledObjects <- scaleCoordsToHabitatGrid(
  coordsData = coordsObsCenter,
  coordsHabitatGridCenter = coordsHabitatGridCenter)

## Get lower and upper cell coordinates
lowerAndUpperCoords <- getWindowCoords(
  scaledHabGridCenter = scaledObjects$coordsHabitatGridCenterScaled,
  scaledObsGridCenter = scaledObjects$coordsDataScaled,
  plot.check = F)

trapLocal <- getLocalObjects(habitatMask = habitatMask,
                             coords = scaledObjects$coordsDataScaled,
                             dmax = 1,
                             resizeFactor = 1,
                             plot.check = TRUE
)

deadObsLocal <- getLocalObjects(habitatMask = habitatMask,
                             coords = scaledObjects$coordsHabitatGridCenterScaled,
                             dmax = 2,
                             resizeFactor = 1,
                             plot.check = TRUE
)


### 1.3 Define model code
modelCode <- nimbleCode({
  ##--------------------------------------------------------------------------------------------
  ##-----------------------------## 
  ##------ SPATIAL PROCESS ------##  
  ##-----------------------------##
  tau ~ dgamma(0.001, 0.001)
  logHabIntensity[1:numHabWindows] <- mu[1:numHabWindows]
  sumHabInt <- log(sum(mu[1:numHabWindows]))
  ## FIRST YEAR 
  for(i in 1:M){
    s[i, 1:2,1] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],
      upperCoords = upperHabCoords[1:numHabWindows, 1:2],
      logIntensities = logHabIntensity[1:numHabWindows],
      logSumIntensity = sumHabInt,
      habitatGrid = habitatGrid[1:y.max,1:x.max],
      numGridRows =  y.max,
      numGridCols = x.max
    )
  }#i
  
  ## T>1 
  for(t in 2:n.years){
    for(i in 1:M){
      s[i, 1:2, t] ~ dbernppACmovement_normal(lowerCoords = lowerHabCoords[1:numHabWindows, 1:2]
                                                , upperCoords = upperHabCoords[1:numHabWindows, 1:2]
                                                , s = s[i, 1:2, t - 1]
                                                , sd = tau
                                                , baseIntensities = mu[1:numHabWindows]
                                                , habitatGrid =  habitatGrid[1:y.max,1:x.max]
                                                , numGridRows = y.max
                                                , numGridCols = x.max
                                                , numWindows= numHabWindows
      )
    }#i
  }#t
  
  ##--------------------------------------------------------------------------------------------
  ##-------------------------------## 
  ##----- DEMOGRAPHIC PROCESS -----## 
  ##-------------------------------##    
  psi ~ dunif(0,1)   
  omeg1[1:3] <- c(1-psi,psi,0)                                                 
  
  
  
  for(t in 1:n.years1){
    phi[t] ~ dunif(0,1)      
    gamma[t] ~ dunif(0,1)
    
    omega[1,1:3,t] <- c(1-gamma[t],gamma[t],0)
    omega[2,1:3,t] <- c(0,phi[t],1-phi[t])
    omega[3,1:3,t] <- c(0,0,1)
  }#t
  

  for(i in 1:M){
    z[i,1] ~ dcat(omeg1[1:3])
    for(t in 1:n.years1){
      z[i,t+1] ~ dcat(omega[z[i,t],1:3,t])
    }#t
  }#i
  
  
  
  ##---------------------------------------------------------------------------------------------   
  ##-----------------------------##
  ##----- DETECTION PROCESS -----## 
  ##-----------------------------##
  sigma ~ dunif(0,15)
  
  for(t in 1:n.years){
    p0[t] ~ dunif(0,1)
    
    for (i in 1:M){
      
      ## ALIVE DETECTIONS 
       y[i, 1:lengthYCombined,t] ~ dbinomLocal_normal(size = trials[1:n.traps],
                                                 p0 = p0[t],
                                                 s = s[i,1:2,t],
                                                 sigma = sigma,
                                                 trapCoords = trapCoords[1:n.traps,1:2],
                                                 localTrapsIndices = trapIndex[1:n.cells,1:maxNBDets],
                                                 localTrapsNum = nTraps[1:n.cells],
                                                 resizeFactor = resizeFactor,
                                                 habitatGrid = habitatIDDet[1:y.maxDet,1:x.maxDet],
                                                 indicator = isAlive[i,t],
                                                 lengthYCombined = lengthYCombined)
      
      
     
    }#i
  }#t
  
  ##---------------------------------------------------------------------------------------------										
  ##----------------------------------------## 
  ##---------- DERIVED PARAMETERS ----------##
  ##----------------------------------------##
  for(i in 1:M){ 
    isAlive[i,1] <- (z[i,1] == 2) 
    for(t in 1:n.years1){
      isAlive[i,t+1] <- (z[i,t+1] == 2) 
    }
  }
  for(t in 1:n.years){
    N[t] <- sum(isAlive[1:M,t])
  }#t
})

xx=2
#for(xx in 2:20){
### 1.4 Define parameter values to simulate

M <- 500

## ASSIGN SIMULATED VALUES 
p0 <- c( 0.05, 0.04, 0.02, 0.05, 0.045) 
sigma <- 0.4

#expected number of individuals alive (z=2)at t=1
n.individualsT1 <- 100
n.years <- 5
phi <- runif(n.years-1, 0.5,0.8)#rep(0.8, n.years-1)
# mhW <- -2
# mhH <- -2

# betaH <- 1
# betaW <- -1
tau <- 1.2

recruitment <- 0.3  
## calculate the gamma 
gamma <- n.individualsT1/M
Recruit <- n.individualsT1*recruitment ## 40 % of recruitment 
NeverAlive <- n.individualsT1
  
for(t in 2:n.years){
    Navai  <- M- NeverAlive[t-1]# Available to be recruited is Augmented - already alive in t = 1
    gamma[t] <- Recruit/Navai # Individuals recruited in t from all available = RECRUITMENT RATE
    NeverAlive[t] <- NeverAlive[t-1] + Recruit # Sum the recruited individuals (not available anymore to be recruited)
}
psi <- gamma[1]

habCov <- as.numeric(scale(lowerAndUpperCoords$lowerHabCoords[,2])[,1])
habCovImage <- lowerAndUpperCoords$habitatGrid
habCovImage[] <- habCov[lowerAndUpperCoords$habitatGrid]
image(habCovImage)

### plot expected surv/mortality prob as a function of covariate 
# mhH1 <- exp(mhH + betaH*habCov)
# mhW1 <- exp(mhW + betaW*habCov)
#   
# phi <- exp(-(mhH1 + mhW1))
# h <-  (1-phi) * (mhH1/(mhH1 + mhW1))
# w <-  (1-phi) * (mhW1/(mhH1 + mhW1))
#   
# plot(phi ~ habCov, ylim=c(0,1), type="b", ylab="Probability")
# points(h ~ habCov, col="red", type="b")
# points(w ~ habCov, col="blue", type="b")
# legend("topright", legend=c("phi","h","w"), 
#        col=c("black","red","blue"), lty=c(1,1,1))

lengthYCombined <- 1 + trapLocal$numLocalIndicesMax*2


### 1.5 Create data, constants and inits objects 

nimConstants <- list(M = M,
                     n.years = n.years,
                     n.years1 = n.years-1,
                     n.traps = dim(scaledObjects$coordsDataScaled)[1],
                     y.max = dim(habitatMask)[1],
                     x.max = dim(habitatMask)[2],
                     y.maxDet = dim(trapLocal$habitatGrid)[1],
                     x.maxDet = dim(trapLocal$habitatGrid)[2],
                     ResizeFactor = trapLocal$resizeFactor,
                     n.cells = dim(trapLocal$localIndices)[1],
                     maxNBDets = trapLocal$numLocalIndicesMax,
                     trapIndex = trapLocal$localIndices,
                     nTraps = trapLocal$numLocalIndices,
                     habitatIDDet = trapLocal$habitatGrid,
                     lengthYCombined = lengthYCombined,
                     numHabWindows = dim(lowerAndUpperCoords$lowerHabCoords)[1],
                     maxNBDetsReco = deadObsLocal$numLocalIndicesMax)



nimData <- list(trapCoords = scaledObjects$coordsDataScaled,
                trials = rep(1, dim(scaledObjects$coordsDataScaled)[1]),
                lowerHabCoords = lowerAndUpperCoords$lowerHabCoords,
                upperHabCoords = lowerAndUpperCoords$upperHabCoords,
                habitatGrid = lowerAndUpperCoords$habitatGrid,
                #habCov = habCov,
                mu = rep(1, dim(lowerAndUpperCoords$lowerHabCoords)[1]),
                resizeFactor = deadObsLocal$resizeFactor
                )

# We set the parameter values as inits
nimInits <- list(sigma = sigma,
                 p0=p0,
                 tau=tau,
                 gamma = gamma[2:n.years],
                # gamma1 = gamma[1],
                 psi = psi,
                 phi = phi) 




### 1.6 Create NIMBLE model
model <- nimbleModel( code = modelCode,
                      constants = nimConstants,
                      data = nimData,
                      inits = nimInits,
                      check = F,       
                      calculate = T)  

model$calculate()

## Error here doesn't matter. It happens because we don't really provide data?

### 1.7 Simualte SCR data from the NIMBLE model
# FIRST WE GET THE NODES TO SIMULATE
#nodesToSim <- model$getDependencies(c("s", "z"), self=T)
nodesToSim <- model$getDependencies(c("p0","sigma","phi", "psi","tau","mu","gamma"),
                                      self = F,
                                      downstream = T,
                                      returnScalarComponents = TRUE)
# THEN WE SIMULATE THOSE NODES 
#set.seed(1)
model$simulate(nodesToSim, includeData = FALSE)


N <- apply(model$z,2,function(x)sum(x==2))
N
# 
# N.recoveredDead <- apply(model$z,2,function(x)sum(x==3))
# N.recoveredDead

## 2. RUN MCMC WITH NIMBLE
myZ <- model$z
z <- zInits <- model$z
whichDet <- apply(model$y, 3, function(x) x[,1]>0 )
whichNotDet <- apply(model$y, 3, function(x) x[,1]==0 )
apply(whichDet,2,sum)

# give NAS to individuals not detected
z[whichNotDet] <- NA
zInits[whichDet] <- NA

nimData$y <- model$y

nimData$z <- z
nimInits$z <- zInits
nimInits$s <- model$s

# CREATE AND COMPILE THE NIMBLE MODEL
modelfit <- nimbleModel( code = modelCode,
                      constants = nimConstants,
                      data = nimData,
                      inits = nimInits,
                      check = F,
                      calculate = F)
modelfit$calculate()
cmodelfit <- compileNimble(modelfit)
cmodelfit$calculate()
MCMCconf <- configureMCMC(model =  modelfit,
                          monitors = c("p0","sigma","psi","gamma", "N",
                                         "gamma","tau","phi","s","z"),
                          control = list(reflective = TRUE),
                          thin = 1)

MCMC <- buildMCMC(MCMCconf)
cMCMC <- compileNimble(MCMC, project = modelfit, resetFunctions = TRUE)
# RUN THE MCMC 
MCMCRuntime <- system.time(myNimbleOutput <- runMCMC( mcmc = cMCMC,
                                                             nburnin = 100,
                                                             niter = 500,
                                                             nchains = 2,
                                                             samplesAsCodaMCMC = TRUE))


myresults <- ProcessCodaOutput(myNimbleOutput)

## -----------------------------------------------------
##  ASP: This will be the step from which I will start 
##             Full model and results 
## -----------------------------------------------------

#### NOW PREDICT FROM THE MODEL POSTERIORS ####
# WE TAKE THE SAME MOEL, BUT WE ADD MORE YEARS.
# HERE WE PREDICT FOR 2 EXTRA YEARS ##
NExtraYears <- 4
nimConstants$n.years <- nimConstants$n.years+NExtraYears
nimConstants$n.years1 <- nimConstants$n.years1+NExtraYears

nimData1 <- nimData

## WE NEED TO ADD EXTRA TIME DIMENSIONS TO ALL OBJECTS,
# yPredict <- array(0,c(dim(nimData$y)[1:2],dim(nimData$y)[3]+NExtraYears )) # Preguntar a Cyril si esto esta mal? Seria:
yPredict <- array(0,c(dim(nimData$y)[1:2],NExtraYears ))
nimData1$y <- abind(nimData$y,yPredict,along=3) 

##SINCE WE ALSO ADD MORE YEARS, WE INCREASE AUGMENTATION TO MAKE SURE WE DONT RUN OUT OF IDS AVAILABLE. 
NEXTRA <- 400
yAUGM <- array(0,c(NEXTRA,dim(nimData1$y)[2:3]))
nimData1$y <- abind(nimData1$y,yAUGM,along=1)
nimConstants$M <- nimConstants$M+NEXTRA

##NOW WE USE THE POSTERIOR DISTRIBUTION OBTAINED FROM THE FITTED MODEL AS PARAMETER VALUES#
# THEY WILL BE USED FOR SIMULATING THE TRAJECTORY OF THE POPULATION 
#here we take the posterior from the  first iteration 
ite <- 1

# WE SET NA FOR Z FOR YEARS TO PREDICT
zPredict <- array(NA,c(dim(nimData$z)[1], NExtraYears))
nimData1$z <- abind(myresults$sims.list$z[ite,,], zPredict, along=2)
# WE AUGMENT 
zAUGM <- array(1,c(NEXTRA,dim(nimData1$z)[2]))
zAUGM[,(nimConstants$n.years-NExtraYears):nimConstants$n.years] <- NA

nimData1$z <- abind(nimData1$z, zAUGM, along=1)

# WE SET NA FOR s FOR YEARS TO PREDICT
sPredict <- array(NA, c(dim(myresults$sims.list$s)[2:3],NExtraYears ))
nimData1$s <- abind(myresults$sims.list$s[ite,,,], sPredict, along=3)
# WE AUGMENT 
sAUGM <- array(0, c(NEXTRA, dim(nimData1$s)[2:3]))
nimData1$s <- abind(nimData1$s, sAUGM, along=1)


##WE MAKE SOME ASSUMPTIONS FOR THE OTHER PARAMETERS# 
# HERE WE ASSUME THAT PHI IS 0.5 FOR THE YEARS TO PREDICT
nimData1$phi <- c(myresults$sims.list$phi[ite,],rep(0.5,NExtraYears))
# HERE WE ASSUME THAT GAMMA IS 0.2 FOR THE YEARS TO PREDICT
nimData1$gamma <- c(myresults$sims.list$gamma[ite,],rep(0.2,NExtraYears))
nimData1$p0 <- c(myresults$sims.list$p0[ite,], rep(myresults$sims.list$p0[ite,1], NExtraYears))
nimData1$sigma <- c(myresults$sims.list$sigma[ite])
nimData1$tau <- c(myresults$sims.list$tau[ite])
nimData1$psi <- c(myresults$sims.list$psi[ite])

names(nimInits)

##NOW WE BUILD THE MODEL AGAIN WITH THE NEW DATA AND THE YEARS TO PREDICT (EXTRA DIMENSIONS)###
modelSims <- nimbleModel( code = modelCode,
                      constants = nimConstants,
                      data = nimData1,
                      #inits = nimInits,
                      check = F,
                      calculate = F)
nimData1$trapCoords
names(nimData1)

# paste("z[",c(1:dim(model$z)[1]), " ,", c(1:dim(model$z)[2]),"]")
nodesToSim <- modelSims$getDependencies(c("p0","sigma",
                                      "gamma","tau","psi","phi"),
                                    self = F,
                                    downstream = T,
                                    returnScalarComponents = TRUE)

modelSims$simulate(nodes = nodesToSim,#c(zNodes,sNodes,yNodes),
                   includeData = F)#---if TRUE: want to simulate new values also for nodes considered as data

# WE CHECK, HERE WE SEE THAT Z FOR THE FIRST 5 FIRST YEARS IS EQUAL TO THE ONES SET AS DATA, BUT WE HAVE VALUES FOR Z 
modelSims$z == nimData1$z
# FOR THE LAST TWO YEARS, THE MODEL PREDICTED Z VALUES! 
modelSims$z
# P
modelSims$phi

# WE CAN SIMULATE SEVERAL TIMES 
N <-list()
for(i in 1:3){
  modelSims$simulate(nodes = nodesToSim,#c(zNodes,sNodes,yNodes),
                     includeData = F)#---if TRUE: want to simulate new values also for nodes considered as data

  N[[i]] <- apply(modelSims$z,2,function(x)sum(x==2))
}
N# only the last 2 years of N changes, it means that it works n is predicted for the last 2 years everytime


##  COMPILE THE MODEL FOR FASTER SIMULATIONS ##
cmodelSims <- compileNimble(modelSims)



# WE FIND THE S AND Z NODES THAT SHOULDNT BE PREDICTED FOR. 
samplerConfList <- modelSims$getNodeNames()
zNodes <- samplerConfList[grep("z",samplerConfList)]
## FIND Z NODES FROM THE FIRST YEAR (SHOULD NOT BE REPLACED)
zNodes <- zNodes[-grep( ", 1]",zNodes)]
zNodes <- zNodes[-grep( ", 2]",zNodes)]
zNodes <- zNodes[-grep( ", 3]",zNodes)]
zNodes <- zNodes[-grep( ", 4]",zNodes)]
zNodes <- zNodes[-grep( ", 5]",zNodes)]

sNodes <- samplerConfList[grep("s\\[",samplerConfList)]
## FIND Z NODES FROM THE FIRST YEAR (SHOULD NOT BE REPLACED)
sNodes <- sNodes[-grep( ", 1]",sNodes)]
sNodes <- sNodes[-grep( ", 2]",sNodes)]
sNodes <- sNodes[-grep( ", 3]",sNodes)]
sNodes <- sNodes[-grep( ", 4]",sNodes)]
sNodes <- sNodes[-grep( ", 5]",sNodes)]



### now we do it for our simulated scenarios ##
# we can picK 3 values of phi# 
phiVal <- c(0.3,0.5,0.8)

NN <- list()
# we take random posteriors to use for simulations as parameter values for the first 5 years.  
itera <- sample(1:dim(myresults$sims.list$z)[1], 15)

for(j in 1:length(phiVal)){
  N <- list()
  for(ite in 1:length(itera)){
    # WE UPDATE THE Z VALUES USING THE POSTERIORS PREDICTED Z FOR THE FIVE FIRST YEARS AND THEN 
    # USE NA FOR THE years to predict#
    
    # WE SET NA FOR Z FOR YEARS TO PREDICT
    zPredict <- array(NA,c(dim(nimData$z)[1],NExtraYears ))
    nimData1$z <- abind(myresults$sims.list$z[ite,,], zPredict, along=2)
    # WE AUGMENT 
    zAUGM <- array(1,c(NEXTRA,dim(nimData1$z)[2])) # The extra individuals (400)
    zAUGM[,(nimConstants$n.years-NExtraYears):nimConstants$n.years] <- NA # They were alive before
    nimData1$z <- abind(nimData1$z, zAUGM, along=1)
    ###
    zNodes <- samplerConfList[grep("z",samplerConfList)]
    # we set the values in the model
    values(cmodelSims, zNodes) <- nimData1$z # Fill the values predicted by the model by new values where the extra years are NA
    
    # WE SET NA FOR s FOR YEARS TO PREDICT
    sPredict <- array(NA, c(dim(myresults$sims.list$s)[2:3],NExtraYears ))
    nimData1$s <- abind(myresults$sims.list$s[ite,,,], sPredict, along=3)
    # WE AUGMENT 
    sAUGM <- array(0, c(NEXTRA, dim(nimData1$s)[2:3]))
    nimData1$s <- abind(nimData1$s, sAUGM, along=1)
    
    sNodes <- samplerConfList[grep("s\\[",samplerConfList)]
    # we set the values in the model
    values(cmodelSims, sNodes) <- nimData1$s
    
    ##phi
    values(cmodelSims, c("phi[5]","phi[6]","phi[7]","phi[8]")) <- rep(phiVal[j], 4)
    ##gamma
    values(cmodelSims, c("gamma[5]","gamma[6]","gamma[7]","gamma[8]")) <- c(0.2,0.25,0.30,0.35)
    
    
    #now we simulate 
    cmodelSims$simulate(nodes = nodesToSim,#c(zNodes,sNodes,yNodes),
                       includeData = F)#---if TRUE: want to simulate new values also for nodes considered as data
    
    N[[ite]] <- apply(cmodelSims$z,2,function(x)sum(x==2))
  }
NN[[j]] <- N
}


### plot  

plot(-1000,xlim=c(0,nimConstants$n.years+1),ylim=c(0,500),ylab="N",xlab="Years")
years <- c(1:nimConstants$n.years)
offset <- c(-0.2,0,0.2)
col <-  gray.colors(3, start =0.1, end = 0.7)
Ntmp <- do.call(rbind,NN[[j]])

for(t in 1:5){
plot.violins2(dat.list =list(Ntmp[,t]),
             at = t+offset[j],
             x=t,
             border=col[j], 
             col=col[j],
             alpha = 0.5,
             add=T)
}

for(j in 1:length(phiVal)){
  Ntmp <- do.call(rbind,NN[[j]])
  #for(i in 1:nrow(Ntmp)){
    for(t in 6:nimConstants$n.years){
  plot.violins2(dat.list =list(Ntmp[,t]),
               at = t+offset[j],
               x=t,
               border=col[j], 
               col=col[j],
              alpha = 0.5,add=T)
    }
#}
}

legend("topleft",legend = phiVal,title="phi",fill=col)


