
## -------------------------------------------------
##          Age structured model Hostetter
## ------------------------------------------------- 


library(nimble)
library(MCMCvis)

setwd("D:/MargSalas/Scripts_MS/Functions/bpa")
source("known.state.cjs.r")
source("simul.cjs.r")
source("cjs.init.z.r")

setwd("D:/MargSalas/Scripts_MS/Functions")
source("capt_hist_bear.r")


## ---- Monitoring data and capture history ----

setwd("D:/MargSalas/Oso/Datos/Tablas_finales")
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
#chb <- as.matrix(chb[,c(3:13)]) # Keep onlu from 2010 - 2020 (systematic sampling)
chb <- as.matrix(chb[,c(3:12)])
###PROBLEM WHEN INDIVIDUALS ARE FIRST DETECTED THE LAST OCCASION, REMOVE THOSE ONES!
# The individuals to remove are: Blizzard,
rem <- c("Blizzard", "New 17-02", "New 18-03", "New 20-03", "New 20-10", "New 20-11", "Sardo")
chb <- chb[-which(rownames(chb) %in% rem), ]
#chb <- chb[which(apply(chb,1,function(x) sum(x)) == 0), ]


n <- dim(chb)[1]
K <- dim(chb)[2]
#nz <- 135 # I augment the population three times the number of detected individuals
nz <- 114 # I augment the population three times the number of detected individuals
#M <- 180
M <- 152
ydatAug <- rbind(chb, matrix(0, nz,K))

## ---- Age structure ----

setwd("D:/MargSalas/Oso/Datos")
info <- read.csv("Info_individuals.csv", header = TRUE, row.names = NULL, sep = ";") # Load info for year of birth
info <- info[,c(4,8)]

# Matrix with exact ages
x <- matrix(NA, nrow = nrow(chb), ncol = ncol(chb)) # Matrix to fill
colnames(x) <- colnames(chb)
rownames(x) <- rownames(chb)

for (i in 1:nrow(x)){
  birth <- as.numeric(info$Year_birth[which(info$ID %in% rownames(x)[i])])
  if (birth > 2010){
    x[i,which(colnames(x) %in% birth):K] <- 0:(K-which(colnames(x) %in% birth))
  } else { x[i,] <- (2010-birth):((2010-birth)+(K-1)) }
}

# In polar bear case study, age 1 is corresponding to 2 year-old bears (independent) and so on
# But in our case it corresponds to 1st year individuals

# Sum 1 so that the first year of life is age = 1
age <- x+1
ageMatAugMinusOne <- rbind(age, matrix(NA, nz,K)) # It would make more sense that is called ageAug (fill in nimble object = age)

# max.age: maximum age in occasion 1 (minus 1) --> BUT I GUESS IT SHOULD BE +1?
#max.age <- max(age[,1], na.rm = TRUE)
max.age <- 23

# centered_age: median observed age (minus 1)
centered_age <- median(age, na.rm = TRUE) 


## ---- Known z ----

# Some index
first <- apply(chb,1,function(x) min(which(x>0))) # Vector with year of first cap
last <- apply(chb,1,function(x) max(which(x>0)))  # last cap
r <- apply(age,1,function(x) min(which(x!=0))) # year of recruitment for observed individuals (individual added to the population)

# State process is informed by both y and age data 

zdatAGE <- zdatNoAGE <- matrix(NA, M, K)
for(i in 1:n){
  zdatNoAGE[i, first[i]:last[i]] <- 1         # when ignoring age data, we known an individuals is alive between first and last observation
  #if age is used to inform z
  zdatAGE[i, r[i]:last[i]] <- 1               # alive between known recruitment year and last observation
  if(r[i]>1)  zdatAGE[i, 1:(r[i]-1)] <- 0     # Not entered prior to age==1
}

## ---- Initial values for z ----

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

zstAGE[38,3]

## ---- MCMC settings ----
nc <- 6; nAdapt <- 500; nb <- 20000; ni <- 200000+nb; nthin <- 10 ##  (warning, hours of run time with full settings)
nc <- 2; nAdapt <- 500; nb <- 100; ni <- 500+nb; nthin <- 1 # demonstrate model



## ---- Run models ----

# Model 1: No age

setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Oso/PopDyn/CR/Models")
source('code_JS.r')

## Format for NIMBLE without age
nim.data<- list(u=zdatNoAGE,y = ydatAug, w=c(rep(1, n), rep(NA, nz)), b=rep(1,K)  ) 
nim.constants<- list(K=K, M=M, Counts=apply(ydatAug, 2,sum))
params <- c("p", "phi", "beta", "psi", "Nsuper", "N","fit.dat", "fit.new")
inits <- list(phi = .80, p = .2, u = zstNoAGE , w=c(rep(NA, n), rep(0, nz)), psi = .5, beta=rep(1/K, K) )
  
## RUN!
## Nimble steps
Rmodel <- nimbleModel(code=code_JS, constants=nim.constants, data=nim.data, calculate=F, check=F, inits=inits)
conf <- configureMCMC(Rmodel,monitors=params, control = list(adaptInterval = nAdapt), thin=nthin, useConjugacy = TRUE)
Rmcmc <- buildMCMC(conf)  #produces an uncompiled R mcmc function
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
# Start: 19:20
out <- runMCMC(Cmcmc, niter = ni, nburnin = nb, nchains = nc, inits=inits, thin=nthin,
                 setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)
  
#setwd("~/Scripts_MS/Oso/PopDyn/CR/Results")
setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Oso/PopDyn/CR/Results")
save(out, file = "out_noAge.RData")



# Model 2: Age

setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Oso/PopDyn/CR/Models")
source('code_JS_AGE.r')

## Format for NIMBLE with age
nim.data<- list(u=zdatAGE,y = ydatAug, age=ageMatAugMinusOne, agePlusOne=ageMatAugMinusOne[,1]+1, w=c(rep(1, n), rep(NA, nz)),
                b=rep(1,K), a=rep(1, max.age) ) 
nim.constants<- list(K=K, M=M, max.age=max.age, Counts=apply(ydatAug, 2,sum))
params <- c("p", "phi", "beta", "psi", "piAGE", "Nsuper", "N","fit.dat", "fit.new")   
inits <- list(phi = .9, p = .2, psi=.6, u = zstAGE, w=c(rep(NA, n), rep(0, nz)), 
              beta=rep(1/K, K),agePlusOne=c(rep(NA, n), rep(1, nz)) )  

ageMatAugMinusOne


## RUN!
## Nimble steps
Rmodel <- nimbleModel(code=code_JS_AGE, constants=nim.constants, data=nim.data, calculate=F, check=F, inits=inits)
conf <- configureMCMC(Rmodel,monitors=params,control = list(adaptInterval = nAdapt), thin=nthin, useConjugacy = TRUE)
Rmcmc <- buildMCMC(conf)  
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

out <- runMCMC(Cmcmc, niter = ni, nburnin = nb, nchains = nc, inits=inits, thin=nthin,
               setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)


#Running chain 1 ...
#warning: logProb of data node agePlusOne[33]: logProb is -Inf.
#warning: logProb of data node u[23, 3]: logProb is -Inf.
#warning: logProb of data node u[37, 4]: logProb is -Inf.
#warning: logProb of data node u[6, 5]: logProb is -Inf.
#warning: logProb of data node u[36, 5]: logProb is -Inf.
#warning: logProb of data node u[18, 6]: logProb is -Inf.
#|-------------|-------------|-------------|-------------|
#  |-------------------------------------------------------|
#  Running chain 2 ...
#warning: logProb of data node agePlusOne[33]: logProb is -Inf.
#warning: logProb of data node u[23, 3]: logProb is -Inf.
#warning: logProb of data node u[36, 5]: logProb is -Inf.
#warning: problem initializing stochastic node Count.new[10]: logProb is -Inf.
#|-------------|-------------|-------------|-------------|
#  |-------------------------------------------------------|
#  
setwd("~/Scripts_MS/Oso/PopDyn/CR/Results")
save(out, file = "out_Age.RData")





summary(out)
str(out)
head(out$chain1)
MCMCsummary(out)
par(mfrow = c(1,1))
MCMCtrace(out,
          pdf = FALSE,
          ind = TRUE,
          Rhat = TRUE,
          n.eff = TRUE)

