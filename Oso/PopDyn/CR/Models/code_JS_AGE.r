
## DIFFERENCE BETWEEN NORMAL JS AND JS-AGE-STRUCTURED MODEL
# In normal JS the state process depends on whether the individivual is recruited or survives.
# In Age-structure it depends on whether the individivual is recruited or survives AND HAS AN AGE
# BUT because we are modeling an augmented population, we don't know the age of the augmented individuals
# That is why we need the distribution of age classes, to assign augmented individuals to an age class.

code_JS_AGE <- nimbleCode({
  
  # Priors and constraints
  psi~ dbeta(1,1)           # data augmentation
  p ~ dbeta(1,1)            # detection
  phi ~ dbeta(1,1)          # survival
  
  # recruitment prob at k 
  beta[1:K] ~ ddirch(b[1:K])
  
  eta[1] <- beta[1]
  for(k in 2:K){
    eta[k] <- beta[k]/(1-sum(beta[1:(k-1)]))
  }
  
  # starting age distribution is not yet recruited (age = 0), or age if recruited
  piAGE[1:max.age] ~ ddirch(a[1:max.age])
  piAGEuncond[1:(max.age+1)] <- c( (1-eta[1] ), eta[1]*piAGE[1:max.age] )  
  
  ## ASP: piAGE is the age distribution 
  # piAGEuncond adds the year when the individual was still not recruited to the 
  # age distribution, because we are modeled augmented individuals as well, that may or
  # may not be recruited. This is actually why we provide to the model agePlusOne, 
  # to index it properly. Then an agePlusOne = 2 is an age = 1, but is indexed as the second category
  # in piAGEuncond (as the first category is UNRECRUITED)
  # RS: add 0 in first place for category 'dead' - no-one starts dead
  
  ## eta[1]*piAGE[1:max.age]
  # The probability of being in an age class in year 1 is conditional on the probability
  # of being recruited in year 1
  
  # Likelihoods 
  for (i in 1:M){
    w[i] ~ dbern(psi)
    
    # Define latent state at first occasion
    # Age process
    # Because you are modeling both real and augmented individuals, you need to place
    # the augmented individuals on an age class depending on the distribution of ages.
    agePlusOne[i] ~ dcat(piAGEuncond[1:(max.age+1)]) 
    age[i,1] <- agePlusOne[i]-1
    trueAge[i,1] <- age[i,1]*z[i,1]    # age if alive (zero if not yet entered or dead)
    # This trueAge is not super useful, but they probably do it because they need it as a covariate
    # in the model including the quadratic effect of age on survival.
    
    # State process
    # The FIRST year, it depends on the age structure of the population. So you construct an age
    # structure (piAGE), only for the first year.
    u[i,1] <- step(agePlusOne[i]-1.1)  # alive if age[i,1] >0 (STEP = EITHER 0-1)
    z[i,1] <- u[i,1]*w[i]  
    
    # Observation process
    y[i,1] ~ dbern(z[i,1]*p)
    
    # derived stuff
    avail[i,1] <- 1- u[i,1]            # still available for recruitment. 
    recruit[i,1] <- z[i,1]             # recruited at k
    
    # for occasions >1     
    for (t in 2:K){
      
      # State process
      u[i,t] ~ dbern(u[i,t-1]*phi + avail[i,t-1]*eta[t])   
      z[i,t] <- u[i,t]*w[i]  
      trueAge[i,t] <- age[i,t]*z[i,t]    
      
      # Age process
      age[i,t] <- age[i,t-1] + max(u[i,1:t]) # ages by one year after recruitment
      
      # Observation process
      y[i,t] ~ dbern(z[i,t]*p)
      
      # derived stuff
      avail[i,t] <- 1- max(u[i,1:t])       
      recruit[i,t] <- equals(z[i,t]-z[i,t-1],1) # recruited at k; check how many individuals are recruited each year
    } #t
  } #i
  
  #derived population level stuff
  for (t in 1:K){
    N[t] <- sum(z[1:M,t])               # Annual abundance
    B[t] <- sum(recruit[1:M,t])         # Number of entries
  } #t
  
  Nsuper <- sum(w[1:M])            # Superpopulation size
  
  # Optional Freeman-Tukey GOF evaluating observed and expected counts.
  for(t in 1:K){
    expectedCount[t] <- N[t]*p
    Count.new[t] ~ dbin(p, N[t])
    
    ## discrepancy from original
    E.dat[t] <- (sqrt(Counts[t])-sqrt(expectedCount[t]) )^2
    
    ## discrepancy from new
    E.new[t] <- (sqrt(Count.new[t])-sqrt(expectedCount[t]) )^2
  }
  fit.dat <- sum(E.dat[1:K])
  fit.new <- sum(E.new[1:K])
  
})#END model