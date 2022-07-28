
## -------------------------------------------------
##          Age structured model Hostetter
## ------------------------------------------------- 


library(nimble)
library(MCMCvis)

#setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Functions/bpa")
# source("known.state.cjs.r")
# source("simul.cjs.r")
# source("cjs.init.z.r")

setwd("D:/MargSalas/Scripts_MS/Functions")
source("capt_hist_bear.r")


## ---- Monitoring data and capture history ----

setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2020")
os <- read.csv("Data_os_96_20.csv", header = TRUE, row.names = NULL)
os <- os[,-1] 

os_id <- os[which(os$Confirmed_Individual != "Indetermined"), ] # Only identified
os_id <- os_id[-which(os_id$Confirmed_Individual == "Camille_AspeOuest"), ] # Remove because I don't know age

# Capture history
HairStCH <- capt_hist_bear(data = os_id, 
                           method = "Sampling_station", 
                           obs_type = c("Hair"))

chb <- as.data.frame(HairStCH$capt.hist)
rownames(chb) <- chb$Confirmed_Individual
chb <- as.matrix(chb[,c(3:13)]) # Keep onlu from 2010 - 2020 (systematic sampling)
#chb <- as.matrix(chb[,c(3:12)])
###PROBLEM WHEN INDIVIDUALS ARE FIRST DETECTED THE LAST OCCASION, REMOVE THOSE ONES!
# The individuals to remove are: Blizzard,
#rem <- c("Blizzard", "New 17-02", "New 18-03", "New 20-03", "New 20-10", "New 20-11", "Sardo")
#chb <- chb[-which(rownames(chb) %in% rem), ]
#chb <- chb[which(apply(chb,1,function(x) sum(x)) == 0), ]


n <- dim(chb)[1]
K <- dim(chb)[2]
#nz <- 135 # I augment the population three times the number of detected individuals
nz <- 114 # I augment the population three times the number of detected individuals
#M <- 180
M <- n+nz
ydatAug <- rbind(chb, matrix(0, nz,K))

## ---- Age structure ----

setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2020")
info <- read.csv("Info_individuals.csv", header = TRUE, row.names = NULL, sep = ";") # Load info for year of birth
info <- info[,c(4,8)]

# Matrix with exact ages
x <- matrix(NA, nrow = nrow(chb), ncol = ncol(chb)) # Matrix to fill
colnames(x) <- colnames(chb)
rownames(x) <- rownames(chb)

for (i in 1:nrow(x)){
  birth <- as.numeric(info$Year_birth[which(info$ID %in% rownames(x)[i])])
  if (birth >= 2010){
    x[i,which(colnames(x) %in% birth):K] <- 0:(K-which(colnames(x) %in% birth))
  } else { x[i,] <- (2010-birth):((2010-birth)+(K-1)) }
}

# In polar bear case study, age 1 is corresponding to 2 year-old bears (independent) and so on
# But in our case it corresponds to 1st year individuals

# Sum 1 so that the first year of life is age = 1
age <- x+1
#set known unborn individuals to 0
age[is.na(age)]<-0

###convert to 'raw' age categories: 
### age = 0: not yet recruited (category 2)
### age = 1: cub1 ( category 3), 
### age = 2: cub2 (category 4), 
### age = 3: subadult 1 (category 5)
### age = 4: subadult 2 (category 6)
### age >=5: adult (category 7)
### category 1: dead

age.cat<-age
age.cat[age==0]<-2
age.cat[age==1]<-3
age.cat[age==2]<-4
age.cat[age==3]<-5
age.cat[age==4]<-6
age.cat[age>=5]<-7

ageMatAug <- rbind(age.cat, matrix(NA, nz,K)) #Matrix with age categories augmented

##in ageMatAug, set age to NA after last observation of individual - below!
##then provide as data

# here, this is highest age category excluding dead and unrecruited, so 5
max.age <- 5

## ---- Known z ----

# Some index
first <- apply(chb,1,function(x) min(which(x>0))) # Vector with year of first cap
last <- apply(chb,1,function(x) max(which(x>0)))  # last cap
r <- apply(age.cat,1,function(x) min(which(x>2))) # year of recruitment for observed individuals (individual added to the population)

for (i in 1:n){
  if (last[i]<K)
ageMatAug[i, (last[i]+1):K]<-NA
}

# State process is informed by both y and age data 

zdatAGE <- matrix(NA, M, K)
for(i in 1:n){
  #zdatNoAGE[i, first[i]:last[i]] <- 1         # when ignoring age data, we known an individuals is alive between first and last observation
  #if age is used to inform z
  zdatAGE[i, r[i]:last[i]] <- 1               # alive between known recruitment year and last observation
  if(r[i]>1)  zdatAGE[i, 1:(r[i]-1)] <- 0     # Not entered prior to age==1
}

## ASP: Both ageMatAug and zdatAGE have NA in the same observations: The ones where we don't know if its alive

## ---- Initial values for z ----

##these aren't random nodes, I wonder if these starting values matter...
##generate but also try without

#with age data
zstAGE <- zdatAGE
for(i in 1:n){
  if(last[i] < K) zstAGE[i,(last[i]+1):K] <- 1 ## ASP: if the last capture is less than the number of occasions, fill as alive until the last oc
}
zstAGE[(n+1):M,] <- 0  # start augmented individuals as not entered to prevent error messages.
zstAGE[!is.na(zdatAGE)] <- NA ## ASP: Put as NA the ones that we want to use as data as initial values (I dont know why)

#randomly generate augmented individuals as part of superpop
w.in <- c(rep(NA, n), rbinom(nz, 1, 0.5))

##generate entrance occasion for all individuals
ent.occ <- sample(1:K, M-n, replace=TRUE)
age.cat.in <- c(rep(NA, n), rep(2, nz)) #2=not yet recruited
#randomly assign starting age category in year 1 to all individuals 
#alive in year 1
age.cat.in[(n+1):M][ent.occ==1]<-sample(3:(max.age+2),sum(ent.occ==1) , replace=TRUE)

##adjust starting values for u accordingly
for(i in (n+1):M){
  zstAGE[i,ent.occ[i-n]:K]<-1
}

##initial values for starting age distribution
is.alive.1 <- ageMatAug[,1][which(ageMatAug[,1]>2)] ## ASP: Those that have entered in the population
##no-one in cat 3 so create apprx. proportions (age cat 5)
piAGE.in<-c(0.2, 0.1,0.05 , 0.1, 0.55) 

###possibly provide inits for age to avoid problems at initialization

age.in <- ageMatAug

##for observed, continue series or stick with adult, bc z values are initiated at 1

for (i in 1:n){
  if (last[i] < K){
    #continue series, then set all >7 to 7
    nyr <- length((last[i]+1):K)
    age.in[i, (last[i]+1):K]<-(ageMatAug[i,(last[i])]+1) : (ageMatAug[i,(last[i])]+nyr)
    age.in[i,][age.in[i,]>7]<-7}
}
## ASP: No entiendo muy bien, me parece dejarlo como está

for (i in (n+1):M){
  if (ent.occ[i-n] == K){
    age.in[i,K]<-3
    age.in[i,1:(K-1)]<-2
  } else {
  nyr <- length((ent.occ[i-n]+1):K)
  #start in cat 3 on entry occasion, then continue series; set values >7 to 7
  age.in[i,ent.occ[i-n]:K] <- age.cat.in[i] : (age.cat.in[i]+nyr)
  age.in[i,][age.in[i,]>7]<-7
  
  ##fill values before entry occasion with 2's
  if (ent.occ[i-n]>1){
    age.in[i,1:(ent.occ[i-n]-1)]<-2
  }
  }}

age.in[!is.na(zdatAGE)]<-NA

## ---- MCMC settings ----
#nc <- 6; nAdapt <- 500; nb <- 20000; ni <- 200000+nb; nthin <- 10 ##  (warning, hours of run time with full settings)
nc <- 3; nAdapt <- 100; nb <- 5000; ni <- 10000+nb; nthin <- 5 # demonstrate model

nc <- 1; nAdapt <- 50; nb <- 100; ni <- 200+nb; nthin <- 1 #test


## ---- Run models ----

# Model 2: Age

source('code_JS_AGEcatV2.r')

## Format for NIMBLE with age
nim.data<- list(u=zdatAGE,y = ydatAug, age.cat=ageMatAug[,1], age=ageMatAug,
                w=c(rep(1, n), rep(NA, nz)),
                b=rep(1,K), a=rep(1, max.age) ) #age=ageMatAugMinusOne, 
nim.constants<- list(K=K, M=M, max.age=max.age) #, Counts=apply(ydatAug, 2,sum)
params <- c("p.ad", "p.sub","p.cub","phi.ad","phi.sub","phi.cub", 
            "beta", "psi", "piAGE", "Nsuper", "N") #,"fit.dat", "fit.new"

#inits <- list(phi = .9, p = .5, psi=.6, u = zstAGE, w=c(rep(NA, n), rep(0, nz)), 
#              beta=rep(1/K, K),agePlusOne=c(rep(NA, n), rep(1, nz)) )  

inits <- list(phi.ad = 0.9,phi.sub=0.8, phi.cub=0.8, p.ad = .5, p.sub = .5, p.cub = .25, 
              psi=.6, u = zstAGE, w=w.in, 
              beta=c(rep(0.09,10), .1) ,age.cat=age.cat.in,age=age.in,
              piAGE=piAGE.in) 

## RUN!
## Nimble steps
Rmodel <- nimbleModel(code=code_JS_AGE, constants=nim.constants, data=nim.data, calculate=F, check=F, inits=inits)
conf <- configureMCMC(Rmodel,monitors=params,control = list(adaptInterval = nAdapt), thin=nthin, useConjugacy = TRUE)
Rmcmc <- buildMCMC(conf)  
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

out <- runMCMC(Cmcmc, niter = ni, nburnin = nb, nchains = nc, inits=inits, thin=nthin,
               setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)


##results make sense, similar to previous model

summary(out)
str(out)
head(out$chain1)
summ<-MCMCsummary(out)
par(mfrow = c(1,1))
MCMCtrace(out,
          pdf = FALSE,
          ind = TRUE,
          Rhat = TRUE,
          n.eff = TRUE)
saveRDS(out, 'BearsAgeCat.rds')

out2<-readRDS('BearsAgeCat.rds')
MCMCsummary(out2)
################################################################################
nim.data<- list(u=zdatAGE,y = ydatAug, age=ageMatAugMinusOne, agePlusOne=ageMatAugMinusOne[,1]+1, w=c(rep(1, n), rep(NA, nz)),
                b=rep(1,K), a=rep(1, max.age)) 
nim.constants<- list(K=K, M=M, max.age=max.age, centered_age=centered_age, Counts=apply(ydatAug, 2,sum))
params <- c("p", "phi0","alpha0", "alpha1", "alpha2", #"phi.age", 
            "beta", "piAGE","psi", "Nsuper", "N","fit.dat", "fit.new") # save trueAge to derive estimated age structure our side the model.
inits <- list(alpha1 = 0, alpha2 = 0, p = .5, psi=.6, u = zstAGE, w=w.in, 
              beta=c(rep(0.09,10), .1) ,agePlusOne=agePlusOne.in ) 


source('code_JS_AGE_survival.r')
## RUN!
## Nimble steps
Rmodel <- nimbleModel(code=code_JS_AgeSurvival, constants=nim.constants, data=nim.data, calculate=F, check=F, inits=inits)
conf <- configureMCMC(Rmodel,monitors=params,control = list(adaptInterval = nAdapt), thin=nthin, useConjugacy = TRUE)
Rmcmc <- buildMCMC(conf)  
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

out  <- runMCMC(Cmcmc, niter = ni, nburnin = nb, nchains = nc, inits=inits, thin=nthin,
                setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)

summ<-MCMCsummary(out)
##this model does not seem to be estimable
saveRDS(out, 'BearsAgeOnSurvival.rds')
