
# This model is a combination of the age-structured JS model (code_JS_AGEcatV3.1)
# and the open SCR model with different trap arrays/year (SCR in Nimble_diftraps)

SCRhab.Open.diftraps.age.effortTrapCov<-nimbleCode({
  
  ### PRIORS ###
  psi~ dbeta(1,1)           # data augmentation
  
  # recruitment prob at k 
  beta[1:Nyr] ~ ddirch(b[1:Nyr])
  
  eta[1] <- beta[1]
  for(k in 2:Nyr){
    eta[k] <- beta[k]/(1-sum(beta[1:(k-1)]))
  }
  
  # starting age distribution is not yet recruited (age = 0), or age if recruited
  piAGE[1:max.age] ~ ddirch(a[1:max.age])
  piAGEuncond[1:(max.age+1)] <- c( (1-eta[1] ), eta[1]*piAGE[1:max.age] )  
  
  ##movement parameters: detection model and between-year AC movement model
  sigma~dunif(0,5) #adjust to units of trap array and space use of species
  sigD~dunif(0,5)  #dispersal Kernel SD, adjust to units of trap array
  
  ##detection parameter - p0 (baseline detection probability), per age category
  p.ad ~ dnorm(0, 0.01)            # detection per age class. Now in logit scale
  p.sub ~ dnorm(0, 0.01) 
  p.cub ~ dnorm(0, 0.01) 
  
  p0[1]<-0  #not recruited yet, placeholder to make indexing work
  p0[2]<-p.cub
  p0[3]<-p.cub
  p0[4]<-p.sub
  p0[5]<-p.sub
  p0[6]<-p.ad
  
  b.effort1~dnorm(0, 0.01)
  b.effort2~dnorm(0, 0.01)
  b.trap~dnorm(0, 0.01)
  
  ##survival probability, per age category
  phi.ad ~ dbeta(1,1)          # survival adults
  phi.sub ~ dbeta(1,1)          # survival subadults
  phi.cub ~ dbeta(1,1)          # survival cubs
  
  phi[1]<-0  #not recruited yet, placeholder to make indexing work
  phi[2]<-phi.cub
  phi[3]<-phi.cub
  phi[4]<-phi.sub
  phi[5]<-phi.sub
  phi[6]<-phi.ad
  
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
  for (i in 1:M){
    
    w[i] ~ dbern(psi) # part of superpopulation?
    
    # Age process
    agePlusOne[i] ~ dcat(piAGEuncond[1:(max.age+1)]) 
    age[i,1] <- agePlusOne[i]-1 #age 0 = not yet entered, 5 = adult
    age.cat[i,1]<-age[i,1] #age category, everything above 5 == 5
    
    # State process
    u[i,1] <- step(agePlusOne[i]-1.1)  # alive if age[i,1] >0
    z[i,1] <- u[i,1]*w[i]  
    
    # derived stuff
    avail[i,1] <- 1- u[i,1]            # still available for recruitment. 
    recruit[i,1] <- z[i,1]             # recruited at k
    
    
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
      # State process
      u[i,t] ~ dbern( u[i,t-1]*phi[ (age.cat[i,t-1]+1) ] + avail[i,t-1]*eta[t] )   #
      z[i,t] <- u[i,t]*w[i]  
      
      # Age process
      age[i,t] <- age[i,t-1] + max(u[i,1:t]) # ages by one year after recruitment
      ##make sure age>5 get converted to age class 5 (adult)
      age.cat[i,t]<-min(age[i,t], max.age)
      
      # derived stuff
      avail[i,t] <- 1- max(u[i,1:t])       
      recruit[i,t] <- equals(z[i,t]-z[i,t-1],1) # recruited at k
      
    }
    
  }#end yr loop for ACs, demographic model
  
  #derived population level stuff
  for (t in 1:Nyr){
    N[t] <- sum(z[1:M,t])               # Annual abundance
    B[t] <- sum(recruit[1:M,t])         # Number of entries
    
    for(c in 1:max.age){ # Abundance per age class.
      N.age[c,t] <- sum(age.cat[1:M,t]==c & z[1:M,t]== 1 )
    } # c
    
    N.cub[t] <- sum(N.age[1:2,t])
    N.sub[t] <- sum(N.age[3:4,t])
    N.ad[t] <- N.age[5,t]
    
  } #t
  
  # In case it doesn't work the sum
  #for(c in 1:max.age){
  #  for (i in 1:M){
  #  isAgeAlive[i,t,c] <- (age.cat[i,t]==c) & (z[i,t]==1)  
  #  }
  #  N.age[c,t] <- sum(isAgeAlive[1:M,t,c])
  #}
  
  Nsuper <- sum(w[1:M])            # Superpopulation size
  
  ##################################################################################################################################
  ##detection model
  for (t in 1:Nyr){
    
    for (k in 1:K){
      
    for(i in 1:M){
      
      logit(p.eff[i,1:J[t],k,t]) <- p0[age.cat[i,t]+1] + b.effort1*effort[1:J[t],k,t,1] + b.effort2*effort[1:J[t],k,t,2] + b.trap*trap[1:J[t],t]
      
      y[i,1:maxDetNums[t],k,t]~dbinomLocal_normal(detNums = detNums[i,k,t],#getSparseY()$detNums
                                             detIndices = detIndices[i,1:maxDetNums[t],k,t],#getSparseY()$detIndices; ASP: Links with trapID
                                             size = ones[1:J[t]], ##NOW: always 1, because we model each occasion separately
                                             p0Traps = p.eff[i,1:J[t],k,t], #model parameter
                                             sigma = sigma, #model parameter
                                             s = sxy[i,1:2,t], #model parameter
                                             trapCoords = X.sc[1:J[t],1:2,t], #trap coordinates (data); ASP: Year specific trap array
                                             localTrapsIndices = localTrapsIndex[1:numHabWindows,1:MaxLocalTraps[t],t], #from getLocalTraps()
                                             localTrapsNum = localTrapsNum[1:numHabWindows,t], #from getLocalTraps()
                                             resizeFactor = 1, #no resizing
                                             habitatGrid = habitatGridDet[1:numGridRows,1:numGridCols],#from getLocalTraps()
                                             indicator = z[i,t]) #model parameter
    }#end occasion loop
    }#end ind loop
  }
  })