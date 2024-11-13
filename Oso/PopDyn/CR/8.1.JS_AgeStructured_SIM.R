##
##
## Simulation and models described in Hostetter NJ, NJ Lunn, ES Richardson, EV Regehr, and SJ Converse: 
## Age-structured Jolly-Seber model expands inference and improves parameter estimation from capture-recapture data.
##
## see also Royle and Dorazio (2008) and Kery and Schaub (2012) for examples of 
## Jolly-Seber models using data augmentation
##
## models include
## JS without age data ("NoAge")
## JS with age data and constant survival ("ConstantSurvival")
## JS with age data and survival is a quadratic function of age ("AgeSpecificSurvival")
##
## 12-Jan-2020
## Although these data have been processed successfully on a computer system at the U.S. Geological Survey (USGS), no warranty expressed or implied is made regarding the display or utility of the data for other purposes, nor on all computer systems, nor shall the act of distribution constitute any such warranty. The USGS or the U.S. Government shall not be held liable for improper or incorrect use of the data described and/or contained herein.


## Packages
library(nimble)
library(coda)
set.seed(1234)

## SELECT MODEL
run.model <- "NoAge"   ## "NoAge", "ConstantSurvival", "AgeSpecificSurvival",
				
## MCMC settings
nc <- 3; AdaptInterval <- 500; nb <- 20000; ni <- 100000+nb; nthin <- 10


## Define stuff
K <- 7                          # years/primary periods
Nsuper <- 400                   # superpopulation
M <- 700                        # Augmented population
psi <- Nsuper/M                 # augmentation parameter (for reference)
beta <- c(0.40, rep(0.10,K-1))  # recruitment
phi0 <- 0.85                    # Survival (if constant) or intercept at centered age (centered_age = 5) if a function of age
centered_age <- 5               # centered age for phi0
alpha0 <- qlogis(phi0)          # Survival intercept
if(run.model == "AgeSpecificSurvival"){
 alpha1 <- -0.50             # Survival parameter * centered_age[i,t] 
 alpha2 <- -0.20             # Survival parameter * centered_age[i,t]^2
 }else{ # not quadratic models
 alpha1 <- 0.00              # Survival parameter * centered_age[i,t] 
 alpha2 <- 0.00              # Survival parameter * centered_age[i,t]^2
}
p <- 0.25                       # detection probability if alive 
piAGE <- c(0.27, 0.17, 0.14, 0.12, 0.11, 0.09, 0.06, 0.03, 0.01)  # p of AGE CLASSES; reasonable values derived from quadratic simulation
max.age <- 10                   # (J) maximum allowed age at time of first sampling
if(length(piAGE) < max.age ) piAGE <- c(piAGE, rep(0, max.age-length(piAGE))  ) 

#age-specific survival (extended to > max.age just in case an individual survives to a higher age)
phi <-  plogis(alpha0 + alpha1*(1:(max.age+K)-centered_age) + alpha2 * (1:(max.age+K)-centered_age)^2 )  

## placeholders
zTRUE <- yTRUE <- ageTRUE <- alive <- avail <- matrix(0L, Nsuper, K)

## generate realizations 

## recruitment
B <- rmultinom(1, Nsuper, beta) # Generate no. of entering ind. per occasion
ent.occ <- rep(1:K, B) ## ASP: Generate the entering occasion for each individual

## Age at occasion 1
ageTRUE[which(ent.occ==1), 1] <- rep(1:max.age, rmultinom(1, B[1], piAGE)) # Years of id entering the population
# The max age is 10, but there are 9 age classes = 9 probabilities (0 prob of being in last). Last age class not generated
## ASP: For example for the first occasion, thera are 152 individuals entering the population (generated in B). These indiv.
## are distributed into the different age classes according to piAGE

## Survival and aging
for(i in 1:Nsuper){
 if(ent.occ[i] > 1)  ageTRUE[i, ent.occ[i]] <- 1 # Second occasion, all individuals have 1 
 zTRUE[i,ent.occ[i]] <- 1

 if(ent.occ[i] < K){ # All except individuals entering the last year (because it is to calculate survival)
  for(k in (ent.occ[i]+1):K){ # Calculate survival every year from entering occasion to last occasion
   phi.it <- plogis(alpha0 + alpha1*(ageTRUE[i, k-1]-centered_age) + alpha2 * (ageTRUE[i, k-1]-centered_age)^2 ) # Survival prob
   zTRUE[i,k] <- rbinom(1, 1, zTRUE[i,k-1]*phi.it) # survival (the individual survives to the next year? 1/0)
   ageTRUE[i, k] <- ageTRUE[i, k-1] +1             # aging (we add a year because it survives)
  }
 }
}

## Detection given alive
for(k in 1:K){
 yTRUE[,k] <- rbinom(Nsuper,1,p*zTRUE[,k])           # detection 1/0
}

N <- apply(zTRUE, 2, sum)   # realized annual abundances
lambda <- N[2:K]/N[1:(K-1)] # realized annual growth rates


### DATA
# total n ever captured
captured <- which(apply(yTRUE,1,sum) > 0) #which individuals were captured
n <- length(captured)

# captures per year
nt <- apply(yTRUE,2,sum)

# capture histories
y <- yTRUE[captured, ]
first <- apply(y,1,function(x) min(which(x>0))) # Vector with year of first cap
last <- apply(y,1,function(x) max(which(x>0)))  # last cap

# In this example, we assume ages are known for captured individuals. 
age <- matrix(NA, n, K)
for(i in 1:n){
 age[i, ]<-ageTRUE[captured[i], ]
}

#year of recruitment for observed individuals (individual added to the population)
r <- apply(age,1,function(x) min(which(x!=0))) 

#Augment datasets
nz <- M-n
yAug <- rbind(y, matrix(0, nz,K)) # Capture stories augmented with undetected augmented ind
ageAug <- rbind(age, matrix(NA, nz,K))


# Known z: state process is informed by both y and age data 
zdatAGE <- zdatNoAGE <- matrix(NA, M, K)
for(i in 1:n){
 zdatNoAGE[i, first[i]:last[i]] <- 1         # when ignoring age data, we known an individuals is alive between first and last observation
 #if age is used to inform z
 zdatAGE[i, r[i]:last[i]] <- 1               # alive between known recruitment year and last observation
 if(r[i]>1)  zdatAGE[i, 1:(r[i]-1)] <- 0     # Not entered prior to age==1
}

## initial values for z
#without age data
 zstNoAGE <- zdatNoAGE
  zstNoAGE[is.na(zdatNoAGE)] <- 0
  zstNoAGE[!is.na(zdatNoAGE)] <- NA

#with age data
 zstAGE <- zdatAGE
 for(i in 1:n){
  if(last[i]<K) zstAGE[i,(last[i]+1):K] <- 0
 }
 zstAGE[(n+1):M,] <- 0  # start augmented individuals as not entered to prevent error messages.
 zstAGE[!is.na(zdatAGE)] <- NA





#############################################
## MODELS IN NIMBLE


## JS without Age data (See Royle and Dorazio 2008 and Kery and Schaub 2012)
  code_JS_SuperPop <- nimbleCode({
    
    # Priors and constraints
    psi ~ dbeta(1,1)             # probability in super population 
    p ~ dbeta(1,1)               # Prior for mean capture
    phi0 ~ dbeta(1,1)            # Prior for survival at theoretical age 0
    alpha0 <- logit(phi0)  

   # recruitment rate from 'not entered population' at t 
    beta[1:K] ~ ddirch(b[1:K]) # beta = entry probability at each occasion

    eta[1] <- beta[1] # eta = conditional probabilities (probability of entry at k conditional on having not yet entered)
    for(k in 2:K){
     eta[k] <- beta[k]/(1-sum(beta[1:(k-1)]))
    }
    
    ## ASP: eta is basically the probability of being recruited

    # Likelihoods 
    for (i in 1:M){
      w[i] ~ dbern(psi)
      
      ## For occasion 1:

      # state process
      u[i,1] ~ dbern(eta[1])## ASP: Recruited or not recruited? Generated from prob. recruitment
      z[i,1] <- u[i,1]*w[i] # Whether it has entered * whether it exists
      
      # Observation process
      y[i,1] ~ dbern(z[i,1]*p) # Whether its alive * whether is detected
      
      # derived stuff
      avail[i,1] <- 1- u[i,1]            # still available for 'recruitment'. ie not yet recruited
                                        # If recruited (u = 1), is not available (avail = 0)
 
      ## for occasions >1     
      for (t in 2:K){

        # State process
        u[i,t] ~ dbern(u[i,t-1]*phi[i,t] + avail[i,t-1]*eta[t]) # an individual may either enter if not previously entered or survive if previously present   
        logit(phi[i,t]) <- alpha0 
        z[i,t] <- u[i,t]*w[i]

        # Observation process
        y[i,t] ~ dbern(z[i,t]*p)

        # derived stuff
        avail[i,t] <- 1- max(u[i,1:t])       # still available for 'recruitment'. ie not yet recruited
      } #t
    } #i
    
    ## Derived population level stuff
    # Annual abundance
    for (t in 1:K){
      N[t] <- sum(z[1:M,t])               
    } #t

    # Annual growth rate
    for (t in 1:(K-1)){
      lambda[t] <- N[t+1]/N[t]               
    } #t

    # Superpopulation size
    Nsuper <- sum(w[1:M])            
	
  })#END model




## JS with Age data and Constant survival
  code_JS_SuperPop_AgeData <- nimbleCode({
    
    # Priors and constraints
    psi ~ dbeta(1,1)             # probability in super population 
    p ~ dbeta(1,1)               # Prior for mean capture
    phi0 ~ dbeta(1,1)            # Prior for survival at theoretical age 0
    alpha0 <- logit(phi0)  

    # CONSTRUCT THE PROBABILITIES
    # starting age distribution is not yet recruited (age = 0), or age if recruited
    piAGE[1:max.age] ~ ddirch(a[1:max.age]) # Probability of having an age (age distribution)
    piAGEuncond[1:(max.age+1)] <- c( (1-eta[1]), eta[1]*piAGE[1:max.age] )  # age distribution conditional on alive at t=1
        # 1-eta is the conditional probability of non having entered
        # Probability of having entered (recruited) and have an age, conditional on the probability of non having entered (1-eta) (??? Dont think so now)
   
    # recruitment rate from 'not entered population' at t  
    beta[1:K] ~ ddirch(b[1:K])

    eta[1] <- beta[1]
    for(k in 2:K){
     eta[k] <- beta[k]/(1-sum(beta[1:(k-1)]))
    }

    # Likelihoods 
    for (i in 1:M){
      w[i] ~ dbern(psi)

      # initial ages: come from the distribution of ages
      agePlusOne[i] ~ dcat(piAGEuncond[1:(max.age+1)]) # where agePlusOne are data (construct the probability distrib of ages)
      age[i,1] <- (agePlusOne[i]-1) 

      # state process
      ## ASP: The difference between an age-structured model and a simple model, is 
      ## 1) in the state process.
      #### In the simple model, the state process in the first occasion is informed by whether the individual was present or not
      #### In the age model, it is informed by the age directly(?)
      ## 2) The age (known, data) paralelly serves to construct the probabilities of each age class, and know the age structure of the population
      ## 3) If you know the year of recruitment from the age data (z initial values) yo can better informe the model on whether the individual
      #### was already on the population or not
      
      u[i,1] <- step(age[i,1]-.1)
      z[i,1] <- u[i,1]*w[i]
      
      # Observation process
      y[i,1] ~ dbern(z[i,1]*p)
      
      # derived stuff
      avail[i,1] <- 1- u[i,1]            # still available for 'recruitment'. ie not yet recruited
 
      # for occasions >1     
      for (t in 2:K){

        # State process
        u[i,t] ~ dbern(u[i,t-1]*phi[i,t] + avail[i,t-1]*eta[t])   
        logit(phi[i,t]) <- alpha0 
        z[i,t] <- u[i,t]*w[i]

        # Age process
        age[i,t] <- age[i,t-1] + max(u[i,1:t]) # ages by one year after recruitment (NIMBLE allows this syntax)

        # Observation process
        y[i,t] ~ dbern(z[i,t]*p)

        # derived stuff
        avail[i,t] <- 1- max(u[i,1:t])       # still available for 'recruitment'. ie not yet recruited
      } #t
    } #i
    
    ## Derived population level stuff
    # Annual abundance
    for (t in 1:K){
      N[t] <- sum(z[1:M,t])               
    } #t

    # Annual growth rate
    for (t in 1:(K-1)){
      lambda[t] <- N[t+1]/N[t]               
    } #t

    # Superpopulation size
    Nsuper <- sum(w[1:M])            
	

  })#END model





## Quadratic age survival
  code_JS_SuperPop_AgeSpecificSurvival <- nimbleCode({
    
    # Priors and constraints
    psi ~ dbeta(1,1)             # probability in super population 
    p ~ dbeta(1,1)               # Prior for mean capture
    phi0 ~ dbeta(1,1)            # Prior for survival at theoretical age 0
    alpha0 <- logit(phi0)
    alpha1 ~ dnorm(0, sd=10)  # relationship with age[i,t]
    alpha2 ~ dnorm(0, sd=10)  # relationship with age[i,t]^2

    # starting age distribution is not yet recruited (age = 0), or age if recruited
    piAGE[1:max.age] ~ ddirch(a[1:max.age])
    piAGEuncond[1:(max.age+1)] <- c( (1-eta[1]), eta[1]*piAGE[1:max.age] )  # age distribution conditional on alive at t=1

   # recruitment rate from 'not entered population' at t 
    beta[1:K] ~ ddirch(b[1:K])

    eta[1] <- beta[1]
    for(k in 2:K){
     eta[k] <- beta[k]/(1-sum(beta[1:(k-1)]))
    }


    # Likelihoods 
    for (i in 1:M){
      w[i] ~ dbern(psi)

      # initial ages
      agePlusOne[i] ~ dcat(piAGEuncond[1:(max.age+1)]) # where agePlusOne are data
      age[i,1] <- (agePlusOne[i]-1) 

      # state process
      u[i,1] <- step(age[i,1]-.1)
      z[i,1] <- u[i,1]*w[i]
      
      # Observation process
      y[i,1] ~ dbern(z[i,1]*p)
      
      # derived stuff
      avail[i,1] <- 1- u[i,1]            # still available for 'recruitment'. ie not yet recruited
 
      # for occasions >1     
      for (t in 2:K){

        # State process
        u[i,t] ~ dbern(u[i,t-1]*phi[i,t] + avail[i,t-1]*eta[t])   
        logit(phi[i,t]) <- alpha0 + alpha1*(age[i,t-1]-centered_age) + alpha2*(age[i,t-1]-centered_age)^2
        z[i,t] <- u[i,t]*w[i]

        # Age process
        age[i,t] <- age[i,t-1] + max(u[i,1:t]) # ages by one year after recruitment (NIMBLE allows this syntax)

        # Observation process
        y[i,t] ~ dbern(z[i,t]*p)

        # derived stuff
        avail[i,t] <- 1- max(u[i,1:t])       # still available for 'recruitment'. ie not yet recruited
      } #t
    } #i
    
    ## Derived population level stuff
    # Annual abundance
    for (t in 1:K){
      N[t] <- sum(z[1:M,t])               
    } #t

    # Annual growth rate
    for (t in 1:(K-1)){
      lambda[t] <- N[t+1]/N[t]               
    } #t

    # Superpopulation size
    Nsuper <- sum(w[1:M])  

  })#END model
  




################ END MODELS

 
################ Run analyses 
  if(run.model == "AgeSpecificSurvival"){
   nim.data<- list(u=zdatAGE,y = yAug, agePlusOne=ageAug[,1]+1, age=ageAug, w=c(rep(1, n), rep(NA, nz)) , centered_age=centered_age, a=rep(1, max.age), b=rep(1, K) ) 
   nim.constants<- list(K = K, M = M, max.age=max.age)
   params <- c("beta","p", "phi0", "alpha0", "alpha1", "alpha2", "piAGE","psi", "Nsuper", "N", "lambda")
   inits <- list(phi0 = .8, p = .3, u = zstAGE, w=c(rep(NA, n), rep(0, nz)), alpha1=0, alpha2=0, beta=rep(1/K, K), piAGE=rep(1/max.age, max.age))  
  }

  if(run.model == "ConstantSurvival"){
   nim.data<- list(u=zdatAGE,y = yAug, agePlusOne=ageAug[,1]+1, age=ageAug, w=c(rep(1, n), rep(NA, nz)), a=rep(1, max.age), b=rep(1, K) ) 
   nim.constants<- list(K = K, M = M, max.age=max.age)
   params <- c("beta","p", "phi0", "piAGE","psi", "Nsuper", "N", "lambda")
   inits <- list(phi0 = .8, p = .3, u = zstAGE, w=c(rep(NA, n), rep(0, nz)), beta=rep(1/K, K), piAGE=rep(1/max.age, max.age) )  
  }

  if(run.model == "NoAge"){
   nim.data<- list(u=zdatNoAGE,y = yAug, w=c(rep(1, n), rep(NA, nz)), b=rep(1, K) ) 
   nim.constants<- list(K = K, M = M)
   params <- c("beta","p", "phi0", "psi", "Nsuper", "N", "lambda")
   inits <- list(phi0 = .8, p = .3, u = zstNoAGE, w=c(rep(NA, n), rep(0, nz)), beta=rep(1/K, K))  
  }

  ## Nimble steps
  rm(list=c("Rmodel","conf","Rmcmc","Cmodel","Cmcmc"))  # just in case previous versions are stored.
  
  if(run.model == "NoAge")               Rmodel <- nimbleModel(code=code_JS_SuperPop, constants=nim.constants, data=nim.data, calculate=F, check=F, inits=inits)
  if(run.model == "ConstantSurvival")    Rmodel <- nimbleModel(code=code_JS_SuperPop_AgeData, constants=nim.constants, data=nim.data, calculate=F, check=F, inits=inits)
  if(run.model == "AgeSpecificSurvival") Rmodel <- nimbleModel(code=code_JS_SuperPop_AgeSpecificSurvival, constants=nim.constants, data=nim.data, calculate=F, check=F, inits=inits)

  conf  <- configureMCMC(Rmodel,monitors=params,control = list(adaptInterval = AdaptInterval), thin=nthin, useConjugacy = TRUE)
  Rmcmc <- buildMCMC(conf)  #produces an uncompiled R mcmc function
  Cmodel<- compileNimble(Rmodel)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
  out   <- runMCMC(Cmcmc, niter = ni, nburnin = nb, nchains = nc, inits=inits,
                   setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)

## Summarize (e.g., ...)
round(summary(out)$q,2)
traceplot(out[,"N[1]"])


## END