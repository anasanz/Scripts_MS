JS.SCRhab.Open.diftraps.sex.effortTrapBhCov.fast<-nimbleCode({
  
  ### PRIORS ###
  ###recruitment probabilities
  for (t in 1:Nyr){
    gamma[t]~dunif(0,1)  
  }
  
  ##movement parameters: detection model and between-year AC movement model
  sigma[1]~dunif(0,5) #adjust to units of trap array and space use of species
  sigma[2]~dunif(0,5) #adjust to units of trap array and space use of species
  sigD~dunif(0,5)  #dispersal Kernel SD, adjust to units of trap array
  
  ##sex ratio
  omega~dunif(0,1)
  
  ##detection parameter - p0 (baseline detection probability), per age category
  p0 ~ dbeta(1,1)            # detection per age class
  # Covariate effects on p (b.effort1, b.effort2, b.trap)
  for(c in 1:nTrapCovs){ 
    trapBetas[c] ~ dnorm(0, 0.01)
  }
  b.bh~dnorm(0, 0.01)
  
  
  ##survival probability, spatial cov and additive age effect
  mu.phi~ dnorm(0,0.01)  #baseline, adult survival
  beta.cov~dnorm(0,0.01) #effect of spatial covariate
  
  ##effect of habitat cov on density
  beta.dens ~ dnorm(0, 0.01)
  
  ##model for density surface - vectorized, plus some other stuff for the function
  ##for initial AC
  mu1[1:numHabWindows] <- exp(beta.dens * habDens[1:numHabWindows]) 
  sumHabIntensity <- sum(mu1[1:numHabWindows])
  logHabIntensity[1:numHabWindows] <- log(mu1[1:numHabWindows])
  logSumHabIntensity <- log(sumHabIntensity)
  
  ##  FIRST YEAR
  ####activity centers, alive state, age yr 1
  for (i in 1:M){ #
    sex[i]~dbern(omega)
    z[i,1] ~ dbern(gamma[1])  #
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
  } #end M loop yr 1
  
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
      
      # get habitat covariate value at previous ac location
      COV[i,t]<-getcovAC(habitatGrid= habitatGrid[1:numGridRows,1:numGridCols],
                         habcov = habSurv[1:numHabWindows], #habSurv is data
                         AC = sxy[i, 1:2, t-1])
      
      # State process, only survival
      logit(phi.eff[i,t])<-mu.phi + beta.cov*COV[i,t]   # effect of covariate 
      
      z[i,t]~dbern(phi.eff[i,t]*z[i,t-1] + gamma[t]*avail[i,t-1])
      avl[i,t]<-sum(z[i,1:t])>0 #if true, individual was alive, no longer avail NEXT year
      avail[i,t]<-1-avl[i,t] #available for recruitment in following year
    }
    #number of recruits and per capita recruitment rate
    R[t]<-sum(avail[1:M,t-1] * z[1:M,t])
    pc.gam[t]<-R[t]/N[t-1]
  }#end yr loop for ACs, demographic model
  
  
  ##################################################################################################################################
  ##detection model
  for (t in 1:Nyr){
    
    for (k in 1:K){
      
      for(i in 1:M){
        
        
        #calculate baseline detection probability as function of effort, 
        # which varies by year, trap and occasion
        #logit(p.eff[i,1:J[t],k,t]) <- p.ad + b.effort1*effort[1:J[t],k,t,1] + b.effort2*effort[1:J[t],k,t,2] + b.trap*trap[1:J[t],t] +
        #  b.bh * prevcap[i,1:J[t],k,t]
        
        #detection model; data are from getSparseY()$y
        #maxDetNums = getSparseY()$maxDetNums
        #Note that in this example, the detection and the habitat grid have the same dimensions
        #otherwise numHabWindows would need to be provided separately for the detection grid
        y[i,1:lengthYCombined[t],k,t]~dbinomLocal_normalBear(#detNums = detNums[i,k,t],#getSparseY()$detNums
                                        #detIndices = detIndices[i,1:maxDetNums[t],k,t],#getSparseY()$detIndices; ASP: Links with trapID
                                        size = ones[1:J[t]], ##NOW: always 1, because we model each occasion separately
                                        #p0Traps = p.eff[i,1:J[t],k,t], #model parameter
                                        p0 = p0,
                                        sigma = sigma[sex[i]+1], #model parameter
                                        s = sxy[i,1:2,t], #model parameter
                                        trapCoords = X.sc[1:J[t],1:2,t], #trap coordinates (data); ASP: Year specific trap array
                                        localTrapsIndices = localTrapsIndex[1:numHabWindows,1:MaxLocalTraps[t],t], #from getLocalTraps()
                                        localTrapsNum = localTrapsNum[1:numHabWindows,t], #from getLocalTraps()
                                        #resizeFactor = 1, #no resizing
                                        lengthYCombined = lengthYCombined[t],
                                        habitatGrid = habitatGridDet[1:numGridRows,1:numGridCols],#from getLocalTraps()
                                        trapCovs = effort[1:J[t],k,t,1:nTrapCovs],
                                        trapBetas = trapBetas[1:nTrapCovs],
                                        indTrapCov = prevcap[i,1:J[t],k,t],
                                        indTrapBeta = b.bh,
                                        indicator = z[i,t]) #model parameter
      }#end individual loop
    }#end occasion loop
    #total abundance in state-space per year
    N[t]<-sum(z[1:M,t])
  } #end year loop for observation model
  
  Nsuper<-sum(alive[1:M])
  
})


getcovAC <- nimbleFunction(
  run=function(
    habitatGrid     = double(2),
    habcov          = double(1),
    AC              = double(1) #current activity center
  ){
    returnType(double(0))
    ## Find which windowcurrent AC falls within
    windowInd <- habitatGrid[trunc(AC[2])+1, trunc(AC[1])+1]
    ## return covariate value in that window
    return(habcov[windowInd])
  })