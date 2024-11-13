code_JS_AGE <- nimbleCode({
  
  # Priors and constraints
  psi~ dbeta(1,1)           # data augmentation
  
  p.ad ~ dbeta(1,1)            # detection per age class
  p.sub ~ dbeta(1,1) 
  p.cub ~ dbeta(1,1) 
  
  p[1]<-0  #not recruited yet, placeholder to make indexing work
  p[2]<-p.cub
  p[3]<-p.cub
  p[4]<-p.sub
  p[5]<-p.sub
  p[6]<-p.ad
  
  #phi ~ dbeta(1,1) 
  phi.ad ~ dbeta(1,1)          # survival adults
  phi.sub ~ dbeta(1,1)          # survival subadults
  phi.cub ~ dbeta(1,1)          # survival cubs

  phi[1]<-0  #not recruited yet, placeholder to make indexing work
  phi[2]<-phi.cub
  phi[3]<-phi.cub
  phi[4]<-phi.sub
  phi[5]<-phi.sub
  phi[6]<-phi.ad
  
  # recruitment prob at k 
  beta[1:K] ~ ddirch(b[1:K])
  
  eta[1] <- beta[1]
  for(k in 2:K){
    eta[k] <- beta[k]/(1-sum(beta[1:(k-1)]))
  }
  
  # starting age distribution is not yet recruited (age = 0), or age if recruited
  piAGE[1:max.age] ~ ddirch(a[1:max.age])
  piAGEuncond[1:(max.age+1)] <- c( (1-eta[1] ), eta[1]*piAGE[1:max.age] )  
  
  
  # Likelihoods 
  for (i in 1:M){
    w[i] ~ dbern(psi)
    
    # Define latent state at first occasion
    # Age process
    agePlusOne[i] ~ dcat(piAGEuncond[1:(max.age+1)]) 
    age[i,1] <- agePlusOne[i]-1 #age 0 = not yet entered, 5 = adult
    age.cat[i,1]<-age[i,1] #age category, everything above 5 == 5
    #trueAge[i,1] <- age[i,1]*z[i,1]    # age if alive (zero if not yet entered or dead)
    
    # State process
    u[i,1] <- step(agePlusOne[i]-1.1)  # alive if age[i,1] > 0
    z[i,1] <- u[i,1]*w[i]  
    
    # Observation process
    y[i,1] ~ dbern(z[i,1]*p[ (age.cat[i,1]+1) ])
    
    # derived stuff
    avail[i,1] <- 1- u[i,1]            # still available for recruitment. 
    recruit[i,1] <- z[i,1]             # recruited at k
    
    # for occasions >1     
    for (t in 2:K){
      
      # State process
      
      ## This is the main difference with the V2 of the model: The state process depends
      # on survival from the previous year, and this survival is specific of each age 
      # category (indesec properly into the different years). We sum 1 because otherwise
      # in age category 0 the model would crash when multiplying by phi. Like this, the
      # second category (phi[2]) corresponds to the first year (phi.cub1)
      
      u[i,t] ~ dbern( u[i,t-1]*phi[ (age.cat[i,t-1]+1) ] + avail[i,t-1]*eta[t] )   #
      z[i,t] <- u[i,t]*w[i]  
      
      #trueAge[i,t] <- age.cat[i,t]*z[i,t]    
      
      # Age process
      age[i,t] <- age[i,t-1] + max(u[i,1:t]) # ages by one year after recruitment
      ##make sure age>5 get converted to age class 5 (adult)
      age.cat[i,t]<-min(age[i,t], max.age) ## ASP: If age is bigger than max.age, it
                                          # takes the max.age instead (to ensure that
                                          # age categories are maximun 5)
      
      # Observation process
      y[i,t] ~ dbern(z[i,t]*p[(age.cat[i,t]+1)])
      
      # derived stuff
      avail[i,t] <- 1- max(u[i,1:t])       
      recruit[i,t] <- equals(z[i,t]-z[i,t-1],1) # recruited at k
    } #t
  } #i
  
  #derived population level stuff
  for (t in 1:K){
    N[t] <- sum(z[1:M,t])               # Annual abundance
    B[t] <- sum(recruit[1:M,t])         # Number of entries
  } #t
  
  Nsuper <- sum(w[1:M])            # Superpopulation size
  
  # # Optional Freeman-Tukey GOF evaluating observed and expected counts.
  # for(t in 1:K){
  #   expectedCount[t] <- N[t]*p
  #   Count.new[t] ~ dbin(p, N[t])
  #   
  #   ## discrepancy from original
  #   E.dat[t] <- (sqrt(Counts[t])-sqrt(expectedCount[t]) )^2
  #   
  #   ## discrepancy from new
  #   E.new[t] <- (sqrt(Count.new[t])-sqrt(expectedCount[t]) )^2
  # }
  # fit.dat <- sum(E.dat[1:K])
  # fit.new <- sum(E.new[1:K])
  
})#END model