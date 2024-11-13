code_JS <- nimbleCode({
  
  # Priors and constraints
  psi ~ dbeta(1,1)          # data augmentation
  p ~ dbeta(1,1)            # detection probability
  phi ~ dbeta(1,1)          # survival
  
  # recruitment prob at k 
  beta[1:K] ~ ddirch(b[1:K])
  
  eta[1] <- beta[1]
  for(k in 2:K){
    eta[k] <- beta[k]/(1-sum(beta[1:(k-1)]))
  }
  
  
  # Likelihoods 
  for (i in 1:M){
    w[i] ~ dbern(psi)
    
    # At first occasion
    # State process
    u[i,1] ~ dbern(eta[1])   
    z[i,1] <- u[i,1]*w[i]
    
    # Observation process
    y[i,1] ~ dbern(z[i,1]*p)
    
    # derived stuff
    avail[i,1] <- 1- u[i,1]            # still available for 'recruitment'. 
    recruit[i,1] <- z[i,1]             # 'recruited' at k
    
    # for occasions >1     
    for (t in 2:K){
      
      # State process
      u[i,t] ~ dbern(u[i,t-1]*phi + avail[i,t-1]*eta[t])   
      z[i,t] <-u[i,t]*w[i]   
      
      # Observation process
      y[i,t] ~ dbern(z[i,t]*p)
      
      # derived stuff
      avail[i,t] <- 1- max(u[i,1:t])            # still available for 'recruitment'. 
      recruit[i,t] <- equals(z[i,t]-z[i,t-1],1) # 'recruited' at k
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