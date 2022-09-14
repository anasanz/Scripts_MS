## -------------------------------------------------
##          Age structured model Hostetter
##    Age categories modeled similar to continuous age model (indexing)
## ------------------------------------------------- 

#### 10.V3.1. Survival and detection specific to originally established age categories 
## Cub = 1+2 
## Subadult = 3+4
## Adult = 5+ 

## The data used:
##  - 2010-2021
##  - ID of cubs detected with the mother ARE included
##  - Included information of death individuals

rm(list = ls())

library(nimble)
library(MCMCvis)

#setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Functions/bpa")
# source("known.state.cjs.r")
# source("simul.cjs.r")
# source("cjs.init.z.r")

setwd("D:/MargSalas/Scripts_MS/Functions")
source("capt_hist_bear.r")


## ---- Monitoring data and capture history ----

setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022/Copia_versionPrev")
os <- read.csv("Data_os_96_21_cubLocations.csv", header = TRUE, row.names = NULL)
os <- os[,-1] 

os_id <- os[which(os$Confirmed_Individual != "Indetermined"), ] # Only identified
os_id <- os_id[-which(os_id$Confirmed_Individual == "Camille_AspeOuest"), ] # Remove because I don't know age

# Capture history
HairStCH <- capt_hist_bear(data = os_id, 
                           method = "Sampling_station", 
                           obs_type = c("Hair"))

chb <- as.data.frame(HairStCH$capt.hist)
rownames(chb) <- chb$Confirmed_Individual
chb <- as.matrix(chb[,c(4:15)]) 

n <- dim(chb)[1]
K <- dim(chb)[2]

nz <- 192 # I augment the population three times the number of detected individuals
M <- n+nz
ydatAug <- rbind(chb, matrix(0, nz,K))

## ---- Age structure ----

setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
info <- readxl::read_xlsx("D:/MargSalas/Oso/Datos/Tablas_finales/2022/Info_individuals_2021.xlsx", sheet = 1)

info <- info[,c(4,8,10)] # I keep column 10 (confirmed_death) because it is the least restrictive

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

# Sum 1 so that the first year of life is age = 1
age <- x+1
#set known unborn individuals to 0
age[is.na(age)]<-0

###convert to 'raw' age categories: 
### age = 0: not yet recruited (category 1)
### age = 1: cub1 ( category 2), 
### age = 2: cub2 (category 3), 
### age = 3: subadult 1 (category 4)
### age = 4: subadult 2 (category 5)
### age >=5: adult (category 6)
### category 1: dead

age.cat<-age
age.cat[age==0]<-1
age.cat[age==1]<-2
age.cat[age==2]<-3
age.cat[age==3]<-4
age.cat[age==4]<-5
age.cat[age>=5]<-6


ageMatAug <- rbind(age.cat, matrix(NA, nz,K)) #Matrix with age categories augmented

# here, this is highest age category excluding dead and unrecruited, so 5
max.age <- 5

## ---- Known z ----

# To include the death information into the z, I create an age matrix including the deaths (value 0)
# The function of this matrix is ONLY to then create the last_alive vector
age.cat.z <- age.cat
for (i in 1:nrow(age.cat)){
  death <- info$Confirmed_death[which(info$ID %in% rownames(age.cat.z)[i])]
  if (death == "Alive") next
  if (death == 2021 | death == 2022) next
  death <- as.numeric(death) + 1 # We don't know when it died exactly, so we set it as death the year after
  # I could dig into this and get the specific date for some individuals
  age.cat.z[i,which(colnames(age.cat.z) %in% death):K] <- 0 
}

# Index
first <- apply(chb,1,function(x) min(which(x>0))) # Vector with year of first cap
r <- apply(age.cat,1,function(x) min(which(x>1))) # year of recruitment for observed individuals (individual added to the population)

# Combined vector to inform of the death or last capture (if death not available) of individuals
last_cap <- apply(chb,1,function(x) max(which(x>0)))  # last capture
last_alive <- apply(age.cat.z,1,function(x) min(which(x == 0))-1) # last occasion alive

last <- last_cap
last[which(last_alive != Inf)] <- last_alive[which(last_alive != Inf)] 

## I THINK THIS STEP IS NOT NEEDED IN THIS MODEL AS IT ONLY MATTERS THE AGE STRUCTURE AT YEAR 1 RIGHT?
##in ageMatAug, set age to NA after last observation of individual 
## ** This would be ONLY for the ones not death no??
##then provide as data
#for (i in 1:n){
#  if (info$Confirmed_death[which(info$ID %in% rownames(age.cat)[i])] != "Alive") next # If it's death I don't introduce NA
#  if (last[i]<K)
#    ageMatAug[i, (last[i]+1):K]<-NA
#}

# State process is informed by y, age data and death recoveries

zdatAGE <- matrix(NA, M, K)
for(i in 1:n){
  #zdatNoAGE[i, first[i]:last[i]] <- 1         # when ignoring age data, we known an individuals is alive between first and last observation
  #if age is used to inform z
  zdatAGE[i, r[i]:last[i]] <- 1               # alive between known recruitment year and last capture or last occasion alive when known
  if(r[i]>1)  zdatAGE[i, 1:(r[i]-1)] <- 0     # Not entered prior to age==1
  if(last_alive[i] != Inf) zdatAGE[i, (last_alive[i]+1):K] <- 0 # Death after the last occasion alive (when available, so when last_alive !=Inf)
}

## ---- Initial values ----
## ----- 1. Initial values for z -----

#with age data
zstAGE <- zdatAGE
for(i in 1:n){
  if(last[i] < K) zstAGE[i,(last[i]+1):K] <- 0 ## ASP: Start as dead after the last observation
} 

zstAGE[(n+1):M,] <- 0  # start augmented individuals as not entered to prevent error messages.
zstAGE[!is.na(zdatAGE)] <- NA ## ASP: Put as NA the ones that we want to use as data (model reasons)

#randomly generate augmented individuals as part of superpop
w.in <- c(rep(NA, n), rbinom(nz, 1, 0.5))

##generate entrance occasion for all individuals
ent.occ <- sample(1:K, M-n, replace=TRUE)

## ASP: GEnerate age category of entry individuals
age.cat.in <- c(rep(NA, n), rep(1, nz)) #2=not yet recruited

#randomly assign starting age category in year 1 to all individuals 
#alive in year 1
age.cat.in[(n+1):M][ent.occ==1]<-sample(2:(max.age+1),sum(ent.occ==1) , replace=TRUE)

##adjust starting values for u accordingly
for(i in (n+1):M){
  zstAGE[i,ent.occ[i-n]:K] <- 1
}

## ----- 2. Initial values for starting age distribution -----

## ASP: PROB. OF AGE CATEGORIES: check observed age categories of first year
ageMatAug[,1][which(ageMatAug[,1]>1)] ## ASP: Those that have entered in the population
##no-one in cat 3 so create apprx. proportions (age cat 4)
piAGE.in<-c(0.2, 0.1, 0.05, 0.1, 0.55)

###possibly provide inits for age to avoid problems at initialization
age.in <- ageMatAug

## ASP: OBSERVED
##for observed, continue series or stick with adult, bc z values are initiated at 1
##not really used in this version, only 1st year needs inits 
##(that is why ageMatAug doesn't include information on following years like deaths or detection info)

for (i in 1:n){
  if (last[i] < K & last_alive[i] == Inf){ # Only for the individuals without death info
    #continue series, then set all >6 to 6
    nyr <- length((last[i]+1):K)
    age.in[i, (last[i]+1):K]<-(ageMatAug[i,(last[i])]+1) : (ageMatAug[i,(last[i])]+nyr)
    age.in[i,][age.in[i,]>6]<-6}
}

## ASP: AUGMENTED
for (i in (n+1):M){
  if (ent.occ[i-n] == K){ # ASP: If the entered occasion is the last occasion (K)  
    
    age.in[i,K] <- 2 ## ASP: fill this occasion with age category 2 (year 1)
    age.in[i,1:(K-1)] <- 1 ## ASP: and before set as age category 1 (non-recruited)
  
    } else {
      
    nyr <- length((ent.occ[i-n]+1):K) # ASP: Number of years to fill (from entered year to K)
    
    #start in cat 3 on entry occasion, then continue series; set values >6 to 6
    age.in[i,ent.occ[i-n]:K] <- age.cat.in[i] : (age.cat.in[i]+nyr)
    age.in[i,][age.in[i,] > 6] <- 6
    
    ##fill values before entry occasion with 1's
    if (ent.occ[i-n] > 1){
      age.in[i,1:(ent.occ[i-n]-1)] <- 1
    }
  }}

age.in[!is.na(zdatAGE)]<-NA ## ASP: Inverse matrix: we put NA in the individuals for which we have data

## ---- MCMC settings ----
#nc <- 6; nAdapt <- 500; nb <- 20000; ni <- 200000+nb; nthin <- 10 ##  (warning, hours of run time with full settings)
nc <- 3; nAdapt <- 100; nb <- 5000; ni <- 10000+nb; nthin <- 5 # demonstrate model

#nc <- 1; nAdapt <- 50; nb <- 100; ni <- 200+nb; nthin <- 1 #test


## ---- Run models ----


setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/CR/Models")
source('code_JS_AGEcatV3.1.R')

## Format for NIMBLE with age

nim.data<- list(y = ydatAug, agePlusOne = ageMatAug[,1], 
                w = c(rep(1, n), rep(NA, nz)), u = zdatAGE,
                b = rep(1,K), a = rep(1, max.age) ) 

nim.constants<- list(K=K, M=M, max.age=max.age) #, Counts=apply(ydatAug, 2,sum)

params <- c("p.ad", "p.sub","p.cub","phi.ad","phi.sub","phi.cub", 
            "beta", "psi", "piAGE", "Nsuper", "N", "B") #,"fit.dat", "fit.new"

inits <- list(phi.ad = 0.9, phi.sub = 0.8, phi.cub = 0.8, p.ad = .5, p.sub = .5, p.cub = .25, 
              psi = .6, w = w.in, agePlusOne  = age.in[,1],
              beta = c(0.15,rep(0.85/(K-1), K-1)),
              piAGE = piAGE.in, u = zstAGE)

Rmodel <- nimbleModel(code=code_JS_AGE, constants=nim.constants, 
                      data=nim.data, calculate=TRUE, check=TRUE, inits=inits)

conf <- configureMCMC(Rmodel,monitors=params,control = list(adaptInterval = nAdapt), 
                       thin=nthin, useConjugacy = TRUE)

Rmcmc <- buildMCMC(conf)  #produces an uncompiled R mcmc function
Cmodel<- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
out   <- runMCMC(Cmcmc, niter = ni, nburnin = nb, nchains = nc, inits=inits,
                 setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)

summ <- MCMCsummary(out,
                    round = 3)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/CR/Results/10.JS_V3.1")
saveRDS(list(out=out, summ=summ),"10.JS_V3.1.rds")

## ---- Plot results ----

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/CR/Results/10.JS_V3.1")
out <- readRDS(file = "10.JS_V3.1.rds")
out_sum2 <- out$summ

library(MCMCvis)

MCMCtrace(out$out,
          ind = TRUE,
          pdf = TRUE, 
          open_pdf = FALSE, 
          filename = '10.JS_V3.1_Traceplot', 
          wd = 'D:/MargSalas/Scripts_MS/Oso/PopDyn/CR/Results/10.JS_V3.1')

source("D:/PhD/MyScripts_PhD/Ch. 2-3/Ch. 3/Results/Functions/ProcessCodaOutput.R")
source("D:/PhD/MyScripts_PhD/Ch. 2-3/Ch. 3/Results/Functions/plot.violins3.r")
source("D:/PhD/MyScripts_PhD/Ch. 2-3/Ch. 3/Results/Functions/DoScale.r")
library(coda)

out.list<- list()
out.list[[1]] <- as.mcmc(out$out$chain1)
out.list[[2]] <- as.mcmc(out$out$chain2)
out.list[[3]] <- as.mcmc(out$out$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)
out$sims.list$N

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/CR/Results/10.JS_V3.1")
pdf("10.JS_V3.1Param.pdf", 7, 7)

par(mfrow = c(2,2),
    oma = c(2,4,2,1),
    mar = c(3,3,2,3))

# Detection probability

params_p <- c("p.cub", "p.sub", 'p.ad')

plot(1, ylim = c(0.5, length(params_p)+0.5), 
     xlim = c(0.2,0.8), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Detection probability",
     cex.axis = 0.8)

axis(2, c(1:length(params_p)), labels = c("Cub(1+2)", "Subad(3+4)", "Adult(>4)"), las = 2, cex.axis = 1)

for (i in 1:length(params_p)){
  plot.violins3(list(out$sims.list[names(out$sims.list) %in% params_p[i]][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("purple"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "black",
              horizontal = TRUE)}

# Survival

params_phi <- c("phi.cub", "phi.sub", 'phi.ad')

plot(1, ylim = c(0.5, length(params_phi)+0.5), 
     xlim = c(0.5,1), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Survival probability",
     cex.axis = 0.8)

axis(2, c(1:length(params_phi)), labels = c("Cub(1+2)", "Subad(3+4)", "Adult(>4)"), las = 2, cex.axis = 1)

for (i in 1:length(params_phi)){
  plot.violins3(list(out$sims.list[names(out$sims.list) %in% params_phi[i]][[1]]),
                x = i,
                at = i,
                violin.width = 0.2,
                plot.ci = 0.95,
                col = c("purple"),
                add = T,
                alpha = 0.3,
                scale.width = FALSE,
                border.col = "black",
                horizontal = TRUE)}

# Age structure

plot(1, ylim = c(0.5, ncol(out$sims.list$piAGE)+0.5), 
     xlim = c(0,1), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Age classes probabilities",
     cex.axis = 0.8)

axis(2, c(1:ncol(out$sims.list$piAGE)), labels = c("1y", "2y", "3y", "4y", ">5y"), las = 2, cex.axis = 1)

ncol(out$sims.list$piAGE)

for (i in 1:ncol(out$sims.list$piAGE)){
  plot.violins3(list(out$sims.list$piAGE[ ,i]),
                x = i,
                at = i,
                violin.width = 0.2,
                plot.ci = 0.95,
                col = c("purple"),
                add = T,
                alpha = 0.3,
                scale.width = FALSE,
                border.col = "black",
                horizontal = TRUE)}

# Abundance

plot(1, ylim = c(0.5, ncol(out$sims.list$N)+0.5), 
     xlim = c(0,max(out$sims.list$N)), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Abundance",
     cex.axis = 0.8)

axis(2, c(1:ncol(out$sims.list$N)), labels = c(2010:2021), las = 2, cex.axis = 1)


for (i in 1:ncol(out$sims.list$N)){
  plot.violins3(list(out$sims.list$N[ ,i]),
                x = i,
                at = i,
                violin.width = 0.4,
                plot.ci = 0.95,
                col = c("purple"),
                add = T,
                alpha = 0.3,
                scale.width = FALSE,
                border.col = "black",
                horizontal = TRUE)}
dev.off()












