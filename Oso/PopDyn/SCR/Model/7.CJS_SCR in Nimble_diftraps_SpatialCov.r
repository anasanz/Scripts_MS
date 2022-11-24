
CJS.SCRhab.Open.diftraps<-nimbleCode({
  
  ### PRIORS ###
  
  ##movement parameters: detection model and between-year AC movement model
  sigma~dunif(0,5) #adjust to units of trap array and space use of species
  sigD~dunif(0,5)  #dispersal Kernel SD, adjust to units of trap array
  
  ##detection parameter - p0 (baseline detection probability), per age category
  p.ad ~ dbeta(1,1)            # detection per age class
  #p.sub ~ dbeta(1,1) 
  #p.cub ~ dbeta(1,1) 
  #clunky, but easier to identify monitored nodes
  #not in current model
  #p0[1]<-p.cub
  #p0[2]<-p.sub
  #p0[3]<-p.ad 
  
  ##survival probability, spatial cov and additive age effect
  mu.phi~ dnorm(0,0.01)  #baseline, adult survival
  #beta.sub~dnorm(0,0.01) #effect of being subadult, not in current model
  #beta.cub~dnorm(0,0.01) #effect of being cub, not in current model
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
  for (i in 1:nobs){ #loop only over observed
    
    ##state stuff for first observation year
    
    z[i,first[i]] <- 1  #always alive when first observed 
    
    #activity centers according to density surface
    sxy[i, 1:2,first[i]] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],#not used; from getWindowCoords()
      upperCoords = upperHabCoords[1:numHabWindows, 1:2],#not used; getWindowCoords()
      logIntensities = logHabIntensity[1:numHabWindows], #from model
      logSumIntensity = logSumHabIntensity, #from model
      habitatGrid = habitatGrid[1:numGridRows,1:numGridCols],#from getWindowCoords()
      numGridRows =  numGridRows, #calculated from habitatGrid (data)
      numGridCols = numGridCols
    )
    
    ## FOLLOWING YEARS
    ###activity centers, demographic model, t>1
    for (t in first.po[i]:Nyr){
      #movement of activity centers
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
      
      ###demographic model
      
      # get habitat covariate value at previous ac location
      COV[i,t]<-getcovAC(habitatGrid= habitatGrid[1:numGridRows,1:numGridCols],
                         habcov = habSurv[1:numHabWindows], #habSurv is data
                         AC = sxy[i, 1:2, t-1])
      
      # State process, only survival
      logit(phi.eff[i,t])<-mu.phi + #beta.cub*is.cub[i,f] + #is.cub binary data
        #beta.sub*is.sub[i,t] +               #is.sub binary data
        beta.cov*COV[i,t]                    # effect of covariate 
      
      z[i,t] ~ dbern(z[i,t-1]*phi.eff[i,t])
      
    }#end yr loop for ACs, demographic model
    
    
    ##################################################################################################################################
    ##detection model
    for (t in first[i]:Nyr){
      
      #detection model; data are from getSparseY()$y
      #maxDetNums = getSparseY()$maxDetNums
      #Note that in this example, the detection and the habitat grid have the same dimensions
      #otherwise numHabWindows would need to be provided separately for the detection grid
      y[i,1:maxDetNums,t]~dbinomLocal_normal(detNums = detNums[i,t],#getSparseY()$detNums
                                             detIndices = detIndices[i,1:maxDetNums,t],#getSparseY()$detIndices; ASP: Links with trapID
                                             size = K[1:J[t]], ##number of trials per trap (data); ASP: different each year
                                             p0 = p.ad, #p0[age.cat[i,t]], #age.cat is data, 1=cub, 2=sub, 3=adult
                                             sigma = sigma, #model parameter
                                             s = sxy[i,1:2,t], #model parameter
                                             trapCoords = X.sc[1:J[t],1:2,t], #trap coordinates (data); ASP: Year specific trap array
                                             localTrapsIndices = localTrapsIndex[1:numHabWindows,1:MaxLocalTraps[t],t], #from getLocalTraps()
                                             localTrapsNum = localTrapsNum[1:numHabWindows,t], #from getLocalTraps()
                                             resizeFactor = 1, #no resizing
                                             habitatGrid = habitatGridDet[1:numGridRows,1:numGridCols],#from getLocalTraps()
                                             indicator = z[i,t]) #model parameter
    }#end year loop for observation model
    
  } #end individual loop
  
  
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