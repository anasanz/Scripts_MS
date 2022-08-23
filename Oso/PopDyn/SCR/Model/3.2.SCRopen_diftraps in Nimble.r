
SCRhab.Open.diftraps<-nimbleCode({
  
  ### PRIORS ###
  
  ###recruitment probabilities
  for (t in 1:Nyr){
    gamma[t]~dunif(0,1)  
  }
  
  ##movement parameters: detection model and between-year AC movement model
  sigma~dunif(0,5) #adjust to units of trap array and space use of species
  sigD~dunif(0,5)  #dispersal Kernel SD, adjust to units of trap array
  
  ##detection parameter - p0 (baseline detection probability)
  p0~dunif(0,1)
  
  ##survival probability
  phi~dunif(0,1)
  
  ##effect of habitat cov on density
  beta.dens ~ dnorm(0, 0.01)
  
  ##model for density surface - vectorized, plus some other stuff for the function
  ##for initial AC 
  mu1[1:numHabWindows] <- exp(beta.dens * habDens[1:numHabWindows])
  sumHabIntensity <- sum(mu1[1:numHabWindows])
  logHabIntensity[1:numHabWindows] <- log(mu1[1:numHabWindows])
  logSumHabIntensity <- log(sumHabIntensity)
  
  ##  FIRST YEAR
  ####activity centers, alive state yr 1
  for (i in 1:M){
    
    #alive state yr 1 (gamma[1] = recruitment prob yr 1)
    z[i,1]~dbern(gamma[1])
    avail[i,1]<-1-z[i,1] #available for recruitment next yr?
    #ever alive?
    alive[i]<-sum(z[i,1:Nyr])>0
    
    
    #activity centers according to density surface
    sxy[i, 1:2,1] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],#not used; from getWindowCoords()
      upperCoords = upperHabCoords[1:numHabWindows, 1:2],#not used; getWindowCoords()
      logIntensities = logHabIntensity[1:numHabWindows], #from model
      logSumIntensity = logSumHabIntensity, #from model
      habitatGrid = habitatGrid[1:numGridRows,1:numGridCols],#from getWindowCoords()
      numGridRows =  numGridRows, #calculated from habitatGrid (data)
      numGridCols = numGridCols
    )
  } # end M loop for yr 1 ACs
  
  ## FOLLOWING YEARS
  ###activity centers, demographic model, t>1
  for (t in 2:Nyr){
    ##for observed inds, model movement of ACs between years
    for (i in 1:nobs){
      sxy[i, 1:2, t] ~ dbernppACmovement_normal(
        lowerCoords            = lowerHabCoords[1:numHabWindows, 1:2],#data getWindowCoords()
        upperCoords            = upperHabCoords[1:numHabWindows, 1:2],#data getWindowCoords()
        s                      = sxy[i, 1:2, t-1], #model parameter
        sd                     = sigD, #model parameter
        baseIntensities        = mu1[1:numHabWindows], 
        habitatGrid            = habitatGrid[1:numGridRows,1:numGridCols],
        numGridRows            = numGridRows,
        numGridCols            = numGridCols,
        numWindows             = numHabWindows
      )
    }
    
    #for never observed, always generate random ACs each year
    for (i in (nobs+1):M){
      sxy[i, 1:2, t] ~ dbernppAC(
        lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],#not used; from getWindowCoords()
        upperCoords = upperHabCoords[1:numHabWindows, 1:2],#not used; getWindowCoords()
        logIntensities = logHabIntensity[1:numHabWindows], #from model
        logSumIntensity = logSumHabIntensity, #from model
        habitatGrid = habitatGrid[1:numGridRows,1:numGridCols],#from getWindowCoords()
        numGridRows =  numGridRows, #calculated from habitatGrid (data)
        numGridCols = numGridCols
      )
    }
    
    ###demographic model
    for (i in 1:M){
      z[i,t]~dbern(phi*z[i,t-1] + gamma[t]*avail[i,t-1])
      avl[i,t]<-sum(z[i,1:t])>0 #if true, individual was alive, no longer avail NEXT year
      avail[i,t]<-1-avl[i,t] #available for recruitment in following year
    }
    
    #number of recruits and per capita recruitment rate
    R[t]<-sum(avail[1:M,t-1] * z[1:M,t])
    pc.gam[t]<-R[t]/N[t-1]
    
  }#end yr loop for ACs, demographic model
  
  ##detection model
  for (t in 1:Nyr){
    for(i in 1:M){
      #detection model; data are from getSparseY()$y
      #maxDetNums = getSparseY()$maxDetNums
      #Note that in this example, the detection and the habitat grid have the same dimensions
      #otherwise numHabWindows would need to be provided separately for the detection grid
      y[i,1:maxDetNums,t]~dbinomLocal_normal(detNums = detNums[i,t],#getSparseY()$detNums
                                             detIndices = detIndices[i,1:maxDetNums,t],#getSparseY()$detIndices; ASP: Links with trapID
                                             size = K[1:J[t]], ##number of trials per trap (data); ASP: different each year
                                             p0 = p0, #model parameter
                                             sigma = sigma, #model parameter
                                             s = sxy[i,1:2,t], #model parameter
                                             trapCoords = X.sc[1:J[t],1:2,t], #trap coordinates (data); ASP: Year specific trap array
                                             localTrapsIndices = localTrapsIndex[1:numHabWindows,1:MaxLocalTraps[t],t], #from getLocalTraps()
                                             localTrapsNum = localTrapsNum[1:numHabWindows,t], #from getLocalTraps()
                                             resizeFactor = 1, #no resizing
                                             habitatGrid = habitatGridDet[1:numGridRows,1:numGridCols],#from getLocalTraps()
                                             indicator = z[i,t]) #model parameter
    }#end ind loop
    
    #total abundance in state-space per year
    N[t]<-sum(z[1:M,t])
    
  }
  
  Nsuper<-sum(alive[1:M])
})


