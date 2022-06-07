SCR0<-nimbleCode({
  
    ###data augmentation parameter
    psi~dunif(0,1)
    ##detection parameters - sigma (movement)
    sigma~dunif(0,5) #adjust to units of trap array and space use of species
    ##detection parameters - lam0 (baseline detection probability)
    p0~dunif(0,1)
    
  ####observation model
      for (i in 1:M){
        #activity centers (min and max coordinates for state space)
        S[i,1]~dunif(xmin, xmax)
        S[i,2]~dunif(ymin, ymax)
        #alive state
        z[i]~dbern(psi)
        
        for (j in 1:J){
          #distance, detection function, effective detection probability
          D[i,j]<-sqrt((X[j,1]-S[i,1])^2 + (X[j,2]-S[i,2])^2)
          p[i,j]<-p0 * exp(-D[i,j]^2/(2*sigma^2))
          p.eff[i,j]<-z[i]*p[i,j]
          #detection model
          y[i,j]~dbinom(p.eff[i,j], K)#K=number of occasions
        }}
    #total abundance in state-space
    N<-sum(z[1:M])
})


SCRhab<-nimbleCode({
  
  ###data augmentation parameter
  psi ~ dunif(0,1)
  ##detection parameters - sigma (movement)
  sigma ~ dunif(0,5) #adjust to units of trap array and space use of species
  ##detection parameters - lam0 (baseline detection probability)
  p0 ~ dunif(0,1)
  
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
      
      #detection model; data are from getSparseY()$y
      #maxDetNums = getSparseY()$maxDetNums
      #Note that in this example, the detection and the habitat grid have the same dimensions
      #otherwise numHabWindows would need to be provided separately for the detection grid
      y[i,1:maxDetNums,1] ~ dbinomLocal_normal(detNums = detNums[i,1],#getSparseY()$detNums
                                    detIndices=detIndices[i,1:maxDetNums,1],#getSparseY()$detIndices
                                    size=K[1:J], ##number of trials per trap (data)
                                    p0 = p0, #model parameter
                                    sigma= sigma, #model parameter
                                    s=sxy[i,1:2], #model parameter
                                    trapCoords=X.sc[1:J,1:2], #trap coordinates (data)
                                    localTrapsIndices=localTrapsIndex[1:numHabWindows,1:MaxLocalTraps], #from getLocalTraps()
                                    localTrapsNum=localTrapsNum[1:numHabWindows], #from getLocalTraps()
                                    resizeFactor=1, #no resizing
                                    habitatGrid=habitatGridDet[1:numGridRows,1:numGridCols],#from getLocalTraps()
                                    indicator=z[i]) #model parameter
    }
  #total abundance in state-space
  N<-sum(z[1:M])
})


