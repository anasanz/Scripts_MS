code_JS_AGE <- nimbleCode({
  
  # Priors and constraints
  psi~ dbeta(1,1)           # data augmentation
  
  p.ad ~ dbeta(1,1)            # detection per age class
  p.sub ~ dbeta(1,1) 
  p.cub ~ dbeta(1,1) 
  
  p[1]<-0
  p[2]<-0
  p[3]<-p.cub
  p[4]<-p.cub
  p[5]<-p.sub
  p[6]<-p.sub
  p[7]<-p.ad

  phi.ad ~ dbeta(1,1)          # survival adults
  phi.sub ~ dbeta(1,1)          # survival subadults
  phi.cub ~ dbeta(1,1)          # survival cubs
  
  # recruitment prob at k 
  beta[1:K] ~ ddirch(b[1:K])
  
  eta[1] <- beta[1]
  for(k in 2:K){
    eta[k] <- beta[k]/(1-sum(beta[1:(k-1)]))
  }
  
  ##state transition probabilities
  ##make 1 = dead, 2 = unrecruited, 3-7 = age classes, for ease of tracking who is alive
  ##start out with fixed transition probs (not individual specific)
  ##consider adding sex in the future
    for (t in 2:K){
    ##from death
    ps[1,t,1]<-1        #remain dead
    ps[1,t,2]<-0        #go back to unrecruited
    ps[1,t,3]<-0        #become cub1
    ps[1,t,4]<-0        #become cub2 
    ps[1,t,5]<-0        #become subadult1
    ps[1,t,6]<-0        #become subadult2
    ps[1,t,7]<-0        #become adult
    
    ##from unrecruited
  ps[2,t,1]<-0 #recruit into death
  ps[2,t,2]<-1-eta[t]   #remain unrecruited
  ps[2,t,3]<-eta[t]   #recruit into 'cub1' age class
  ps[2,t,4]<-0   #recruit into cub2 age class
  ps[2,t,5]<-0   #recruit into subadult1 age class
  ps[2,t,6]<-0   #recruit into subadult2 age class
  ps[2,t,7]<-0   #recruit into adult age class
  
    ##from cub1
  ps[3,t,1]<-1-phi.cub   #die
  ps[3,t,2]<-0           # go back to unrecruited
  ps[3,t,3]<-0           #remain cub1 
  ps[3,t,4]<-phi.cub     #become cub2
  ps[3,t,5]<-0   #become subdult1 
  ps[3,t,6]<-0   #become subdult2 
  ps[3,t,7]<-0   #become adult1 
  
    ##from cub2
  ps[4,t,1]<-1-phi.cub   #die
  ps[4,t,2]<-0           #go back to unrecruited 
  ps[4,t,3]<-0           #go back to cub1
  ps[4,t,4]<-0           #remain cub2
  ps[4,t,5]<-phi.cub     #become subadult1 
  ps[4,t,6]<-0           #become subadult2
  ps[4,t,7]<-0           #become adult
  
    ##from subadult1
  ps[5,t,1]<-1-phi.sub   #die
  ps[5,t,2]<-0        #go back to unrecruited 
  ps[5,t,3]<-0        #become cub1
  ps[5,t,4]<-0        #become cub2
  ps[5,t,5]<-0        #remain subadult1
  ps[5,t,6]<-phi.sub  #become subadult2
  ps[5,t,7]<-0        #become adult 
  
  ##from subadult2
  ps[6,t,1]<-1-phi.sub   #die
  ps[6,t,2]<-0        #go back to unrecruited 
  ps[6,t,3]<-0        #become cub1
  ps[6,t,4]<-0        #become cub2
  ps[6,t,5]<-0        #become subadult1
  ps[6,t,6]<-0        #remain subadult2
  ps[6,t,7]<-phi.sub  #become adult 
  
  ##from adult
  ps[7,t,1]<-1-phi.ad   #die
  ps[7,t,2]<-0        #go back to unrecruited 
  ps[7,t,3]<-0        #become cub1
  ps[7,t,4]<-0        #become cub2
  ps[7,t,5]<-0        #become subadult1
  ps[7,t,6]<-0        #become subadult2
  ps[7,t,7]<-phi.ad   #remain adult 
  
    }
  
  
  # starting age class distribution 
  # age here refers to age category, max.age=5 (adult)
  # uncond account for not yet recruited's and deads
  # where 1=dead, 2= not recruited and 3-7 = alive age classes
  # no dead animals in yr 1
  piAGE[1:max.age] ~ ddirch(a[1:max.age])
  
  #add 0 in first place for category 'dead' - no-one starts dead
  piAGEuncond[1:(max.age+2)] <- c( 0, (1-eta[1] ), eta[1]*piAGE[1:max.age] )  
  
  # Likelihoods 
  for (i in 1:M){
    w[i] ~ dbern(psi) #part of superpopulation (ie, ever alive)?
    
    # Define latent state at first occasion
    # Age process (no-one starts as 'dead')
    age.cat[i] ~ dcat(piAGEuncond[1:(max.age+2)]) 
    #1=dead (no-one), 2=not recruited, 3-7=alive age classes
    
    age[i,1] <- age.cat[i]              #true age class (ie, 2=unrecr, 3 = cub1 to 7 = adult)
    #trueAge[i,1] <- age[i,1]*z[i,1]    # age if alive at t=1 (zero if not yet entered)
    
    # State process
    u[i,1] <- step(age[i,1]-2.1)  # alive if age[i,1] >2 
    z[i,1] <- u[i,1]*w[i]  #z = alive state in each year, alive if in alive age class and part of super-population
    
    # Observation process
    y[i,1] ~ dbern(z[i,1]*p[age[i,1]])
    
    # derived stuff
    #avail[i,1] <- 1- u[i,1]            # still available for recruitment. 
    recruit[i,1] <- z[i,1]             # recruited at k=1
    #dead[i,1]<-0                       #dead, yes or no; no-one starts out dead
    
    # for occasions >1     
    for (t in 2:K){
      
      # State process: aging and survival (via state transitions)
      age[i,t] ~ dcat(ps[age[i,t-1],t, 1:7 ]) #1 means dead, 2 means not yet recruited 
      #age[i,t] <- age.cat[i,t]-1  
      #dead[i,t]<-step(age.cat[i,t]-1.1)            #dead
      
      ## ASP: So in this model, aging depends on the probability of transitioning 
      # from the age class you are in, to any other category (vector of prob, although
      # many are 0). This transition probabilities are linked to survival and recruitment.
      
      u[i,t] <- step(age[i,t]-2.1)             #alive (in cat 3 or higher)
      
      z[i,t] <- u[i,t] *w[i]                      #alive, part of superpop
      #trueAge[i,t] <- age[i,t]*z[i,t]    
      
      # Observation process
      y[i,t] ~ dbern(z[i,t]*p[age[i,t]])
      
      # derived stuff
      #avail[i,t] <- 1- max(u[i,1:t])       
      recruit[i,t] <- equals(z[i,t]-z[i,t-1],1) # recruited at k; ASP: check how many individuals are recruited each year
    } #t
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