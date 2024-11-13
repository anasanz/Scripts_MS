
SCRhab_covsp0_sigSex <- nimbleCode({
  
  ###data augmentation parameter
  psi ~ dunif(0,1)
  
  ##movement parameters: detection model and between-year AC movement model
  sigma[1]~dunif(0,5) # Sex-specific sigma (1 = Females; 2 = Males)
  sigma[2]~dunif(0,5)
  
  ##detection parameters - lam0 (baseline detection probability)
  p0 ~ dunif(0,1)
  
  b.effort1~dnorm(0, 0.01)
  b.effort2~dnorm(0, 0.01)
  b.trap~dnorm(0, 0.01)
  b.bh~dnorm(0, 0.01)
  
  omega~dunif(0,1) # Prior for sex latent variable
  
  ##effect of habitat cov on density
  beta.dens ~ dnorm(0, 0.01)
  
  ##model for density surface - vectorized, plus some other stuff for the function
  mu1[1:numHabWindows] <- exp(beta.dens * habDens[1:numHabWindows]) # Expected dens in each cell
  sumHabIntensity <- sum(mu1[1:numHabWindows]) # Sum of expected densities
  logHabIntensity[1:numHabWindows] <- log(mu1[1:numHabWindows]) # Backtransformed expected density
  logSumHabIntensity <- log(sumHabIntensity) # Back-transformed sum of expected densities?
  
  ####observation model
  for (i in 1:M){
    #activity centers according to density surface
    sxy[i, 1:2] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],#not used; from getWindowCoords()
      upperCoords = upperHabCoords[1:numHabWindows, 1:2],#not used; getWindowCoords()
      logIntensities = logHabIntensity[1:numHabWindows], #from model
      logSumIntensity = logSumHabIntensity, #from model
      habitatGrid = habitatGrid[1:numGridRows,1:numGridCols],#from getWindowCoords()
      numGridRows =  numGridRows, #calculated from habitatGrid (data)
      numGridCols = numGridCols
    )
    
    #alive state
    z[i]~dbern(psi)
    
    
    #detection model
    
    sex[i]~dbern(omega)     # Sex is a latent variable
    
    for (k in 1:K){
      
    logit(p.eff[i,1:J,k]) <- p0 + b.effort1*effort[1:J,k,1] + b.effort2*effort[1:J,k,2] + b.trap*trap[1:J] +
      b.bh * prevcap[i,1:J,k]

    y[i,1:maxDetNums,k] ~ dbinomLocal_normal(detNums = detNums[i,k],#getSparseY()$detNums
                                             detIndices=detIndices[i,1:maxDetNums,k],#getSparseY()$detIndices
                                             size = ones[1:J], ##NOW: always 1, because we model each occasion separately
                                             p0Traps = p.eff[i,1:J,k], #model parameter
                                             sigma = sigma[sex[i]+1], #model parameter
                                             s=sxy[i,1:2], #model parameter
                                             trapCoords=X.sc[1:J,1:2], #trap coordinates (data)
                                             localTrapsIndices=localTrapsIndex[1:numHabWindows,1:MaxLocalTraps], #from getLocalTraps()
                                             localTrapsNum=localTrapsNum[1:numHabWindows], #from getLocalTraps()
                                             resizeFactor=1, #no resizing
                                             habitatGrid=habitatGridDet[1:numGridRows,1:numGridCols],#from getLocalTraps()
                                             indicator=z[i]) #model parameter
    }
  }
  #total abundance in state-space
  N<-sum(z[1:M])
})


