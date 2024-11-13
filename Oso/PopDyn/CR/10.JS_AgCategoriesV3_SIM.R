
## Simulate and analyze JS data with age categories
## mirroring bear analysis
## Question: can we estimate cub survival if we don't observe most cubs that
## die young

## Packages
library(nimble)
library(MCMCvis)
library(coda)
set.seed(1234)

				
## MCMC settings
nc <- 3; nAdapt <- 100; nb <- 5000; ni <- 10000+nb; nthin <- 5 # demonstrate model

## model code

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/CR/Models")
source('code_JS_AGEcatV3.r')

## Define stuff
K <- 10                          # years/primary periods
Nsuper <- 100                   # superpopulation
M <- 200                        # Augmented population
psi <- Nsuper/M                 # augmentation parameter (for reference)
beta <- c(0.25, rep(0.75/(K-1),K-1))  # recruitment
phi0 <- c(0.5,0.5, 0.85, 0.85,0.9)       # Survival cubs, subadults, adults
alpha0 <- qlogis(phi0)          # Survival intercept

p <- c(0.25,0.25, 0.5, 0.5,0.6)                       # detection probability (cub, subadult, adult)

piAGE <- c(0.23, 0.09, 0.04, 0.14, 0.5)  # estimates from bear model
max.age <- 5                   # maximum age category, adult

#age-specific survival (extended to >max.age just in case an individual survives to a higher age)
phi <-  phi0 


#niter<-10

#for (ii in 2:niter){
## placeholders
zTRUE <- yTRUE <- ageTRUE <- alive <- avail <- matrix(0L, Nsuper, K)

## generate realizations 

## recruitment
B <- rmultinom(1, Nsuper, beta) # Generate no. of entering ind. per occasion
ent.occ <- rep(1:K, B)

## Age at occasion 1
ageTRUE[which(ent.occ==1), 1] <- rep(1:max.age, rmultinom(1, B[1], piAGE)) 
## ASP: Generate age structure for the first 20 individuals that are recruited in (year 1)

## Survival and aging
for(i in 1:Nsuper){
 if(ent.occ[i] > 1)  ageTRUE[i, ent.occ[i]] <- 1
 zTRUE[i,ent.occ[i]] <- 1

 if(ent.occ[i] < K){ 
   
  for(k in (ent.occ[i]+1) : K){
    
   phi.it <- phi[ageTRUE[i, k-1]] ## ASP: Index survival of the age category the previous year
   zTRUE[i,k] <- rbinom(1, 1, zTRUE[i,k-1]*phi.it) ## ASP: Whether it survives or not depends on
                                                  # being alive last year(zTRUE) and the probability 
                                                  # of surviving for its age class
   newage <- ageTRUE[i, k-1] + 1 ## ASP: Sum a year to the previous age
   ageTRUE[i, k] <- ifelse(newage > 5, 5, newage) # aging. ASP: Attribute age to current year, the 
                                                  # Maximum age is the max age category (5)
  }
 }
}

## Detection given alive
for(i in 1:Nsuper){
  for(k in 1:K){
    yTRUE[i,k] <- rbinom(1,1,p[ageTRUE[i, k]]*zTRUE[i,k])           # detection 1/0
  }
}
##NA warnings come from 0 age, all of these are 0 detections since animal is not
##yet alive

yTRUE[is.na(yTRUE)] <- 0

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
first <- apply(y,1,function(x) min(which(x>0))) # first cap
last <- apply(y,1,function(x) max(which(x>0)))  # last cap

# In this example, we assume ages are known for captured individuals. 
age <- matrix(NA, n, K)
for(i in 1:n){
 age[i, ]<-ageTRUE[captured[i], ]
}

#year of recruitment for observed individuals
r <- apply(age,1,function(x) min(which(x!=0)))

##turn into category, with 1=not recruited
ageTRUE.cat<- ageTRUE + 1

##ignore death state for this version, is separate from aging data
#ageTRUE.cat[zTRUE==0 & ageTRUE>0]<-1 #dead

age.cat <- age + 1

##to use as data in model, set to NA after last observation
##though in this version we only use yr 1 as data
for (i in 1:n){
  if (last[i] < K) age.cat[i,(last[i]+1):K]<-NA
}

###Q1: How many cubs do we observe versus how many we know were there 
###    from back-filling age data

nobs <- table(age.cat[y==1])/table(age.cat[age.cat>1])

# table(age.cat[y==1]) ## ASP: Number of individuals from each age category we observe (throughout the study)
# table(age.cat[age.cat>2])

##about 50% for all age classes, that does not match data where it's about 30%
##for youngest cubs
##maybe rectified with Ana's data fix to extract cub detections from mom's detection

## ASP: ?? Question for RS: How do you know in data that 30% of detected cubs?
## Because we only detect individuals when they are older? Check script 9 (data)

#Augment datasets
nz <- M-n
yAug <- rbind(y, matrix(0, nz,K))
ageMatAug <- rbind(age.cat, matrix(NA, nz, K))


# Known z: state process is informed by both y and age data 
zdatAGE <- matrix(NA, M, K)
for(i in 1:n){
 #if age is used to inform z
 zdatAGE[i, r[i]:last[i]] <- 1               # alive between known recruitment year and last observation
 if(r[i] > 1)  zdatAGE[i, 1:(r[i]-1)] <- 0     # Not entered prior to age==1
}

###############################################################################
## inits for age
## initial values for z

## ---- ASP: u (alive state = zstAGE) ----

#with age data
 zstAGE <- zdatAGE
 for(i in 1:n){
  if(last[i] < K) zstAGE[i,(last[i]+1):K] <- 0 ## ASP: Start as dead after the last observation
 }
 zstAGE[(n+1):M,] <- 0  # start augmented individuals as not entered to prevent error messages.
 zstAGE[!is.na(zdatAGE)] <- NA ## ASP: Inverse matrix: we put NA in the individuals for which we have data
                              ## BEcause it feeds from the data already?
 
 
 #randomly generate augmented individuals as part of superpop
 w.in <- c(rep(NA, n), rbinom(nz, 1, 0.5))
 
 ##generate entrance occasion for all individuals
 ent.occ.aug <- sample(1:K, nz, replace = TRUE)
 
 ## ASP: GEnerate age category of entry individuals
 age.cat.in <- c(rep(NA, n), rep(1, nz)) #1=not yet recruited (ASP, 1 is the default category)
 
 #randomly assign starting age category in year 1 to all individuals 
 #alive in year 1
 ## ASP: For the augmented individuals that entered in year 1, assign an age cat
 age.cat.in[(n+1):M][ent.occ.aug == 1] <- sample(2:(max.age+1),sum(ent.occ.aug==1) , 
                                             prob = piAGE,
                                             replace = TRUE) ## ASP: Samples from a category with probabilities
 
 ## ASP: Only assign age category for year 1, because it is what is used in the model to construct the 
 ## probabilities
 
 ##adjust starting values for u accordingly
 
 for(i in (n+1):M){
   zstAGE[i,ent.occ.aug[i-n] : K] <- 1
 }  ## ASP: From the entered occasion to K, fill with one

 ## ---- ASP: agePlusOne (Age structure the year 1 = age.in)  ----
 ## INITIAL VALUES FOR AGE CATEGORIES
 
 age.in <- ageMatAug
 
 ## ASP: OBSERVED
 ##for observed, continue series or stick with adult
 ##not really used in this version, only 1st year needs inits
 
 for (i in 1:n){
   
   if (last[i] < K){
     
     #continue series, then set all > 6 to 6
     nyr <-length((last[i]+1):K) # Number of years to be filled
     age.in[i, (last[i]+1):K] <- (ageMatAug[i,(last[i])]+1) : (ageMatAug[i,(last[i])]+nyr) ## ASP: Fill up sequence from last age category to last occasion
     age.in[i,][age.in[i,]>6]<-6} ## ASP: All ages > 6 correspond to the age category 6
 }
 
 ## ASP: AUGMENTED
 for (i in (n+1):M){
   if (ent.occ.aug[i-n] == K){ ## ASP: If the entered occasion is the last occasion (K) 
     
     age.in[i,K] <- 2 ## ASP: fill this occasion with age category 2 (year 1)
     age.in[i,1:(K-1)] <- 1 ## ASP: and before set as age category 1 (non-recruited)
     
   } else { ## ASP: If the entered occasion is NOT the last occasion (K)
     
     nyr <- length((ent.occ.aug[i-n]+1):K) ## ASP: Number of years to fill (from entered year to K)
     
     #start at age.cat.in on entry occasion, then continue series; set values >7 to 7
     age.in[i,ent.occ.aug[i-n]:K] <- age.cat.in[i] : (age.cat.in[i]+nyr) ## ASP: Fill up
     age.in[i,][age.in[i,] > 6] <- 6
     
     ##fill values before entry occasion with 2's
     if (ent.occ.aug[i-n] > 1){
       age.in[i,1:(ent.occ.aug[i-n]-1)]<-1
     }
   }}
 ##set observed data to NA
 age.in[!is.na(zdatAGE)] <- NA ## ASP: Inverse matrix: we put NA in the individuals for which we have data
                                ## BEcause it feeds from the data already?
 
 
 ## Format for NIMBLE with age
 nim.data<- list(y = yAug, agePlusOne = ageMatAug[,1], #age=ageMatAug,
                 w = c(rep(1, n), rep(NA, nz)), u = zdatAGE,
                 b = rep(1,K), a = rep(1, max.age) ) #age=ageMatAugMinusOne, 
 
 nim.constants<- list(K = K, M = M, max.age = max.age) #, Counts=apply(ydatAug, 2,sum)
 
 params <- c("p.ad", "p.sub","p.cub","phi.ad","phi.sub","phi.cub", 
             "beta", "psi", "piAGE", "Nsuper", "N") #,"fit.dat", "fit.new"
 
 #inits <- list(phi = .9, p = .5, psi=.6, u = zstAGE, w=c(rep(NA, n), rep(0, nz)), 
 #              beta=rep(1/K, K),agePlusOne=c(rep(NA, n), rep(1, nz)) )  
 
 inits <- list(phi.ad = 0.9, phi.sub = 0.8, phi.cub = 0.8, p.ad = .5, p.sub = .5, p.cub = .25, 
               psi = .6, w = w.in, agePlusOne  = age.in[,1],#age.in[,1],
               beta = c(0.15,rep(0.85/(K-1), K-1)),#age.cat=age.cat.in,age=age.in,
               piAGE = piAGE, u = zstAGE) #c(rep(NA, n), rep(0, nz)) c(rep(NA, n), rep(1, nz))
#phi.ad = 0.9,phi.sub=0.8, phi.cub=0.8,
  ###save data, raw and input, and inits
 
 #dat<-list(zTRUE=zTRUE, ageTRUE.cat=ageTRUE.cat, yTRUE=yTRUE, nim.data=nim.data,
#           nim.inits=inits)
 
 #saveRDS(dat, paste('Sim/Data_', ii, '.rds', sep=''))
 
  ## Nimble steps
 
  Rmodel <- nimbleModel(code=code_JS_AGE, constants=nim.constants, 
                        data=nim.data, calculate=TRUE, check=TRUE, inits=inits)

  conf  <- configureMCMC(Rmodel,monitors=params,control = list(adaptInterval = nAdapt), 
                         thin=nthin, useConjugacy = TRUE)
  
  Rmcmc <- buildMCMC(conf)  #produces an uncompiled R mcmc function
  Cmodel<- compileNimble(Rmodel)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
  out   <- runMCMC(Cmcmc, niter = ni, nburnin = nb, nchains = nc, inits=inits,
                   setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)

  summ <- MCMCsummary(out)
  
  #saveRDS(list(out=out, summ=summ), paste('Sim/Results_', ii, '.rds', sep=''))
#} #end iteration loop


## Summarize (e.g., ...)
round(summary(out)$q,2)
traceplot(out[,"Nsuper"])


## END