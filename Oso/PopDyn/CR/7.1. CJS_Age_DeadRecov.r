
## -------------------------------------------------
##    MODEL 7.7 BPA: CJS with AGE EFFECTS
##        + Dead recovery in latent state
## ------------------------------------------------- 

# I want to test if including the deaths the model works:
# 1. I test it in the data already simulated under the model Phi[ageclass]_p0
#   -> It works!
# 2. I analyze my data with the model Phi[ageclass]_p[ageclass] 

rm(list = ls())

library(dplyr)
library(tidyr)
library(rjags)
library(jagsUI)

setwd("D:/MargSalas/Scripts_MS/Functions/bpa")
source("known.state.cjs.r")
source("simul.cjs.r")
source("cjs.init.z.r")
source("cjs.init.zdead.r")

setwd("D:/MargSalas/Scripts_MS/Functions")
source("capt_hist_bear.r")

## ---- Simulation ----

# Define parameter values
n.occasions <- 10                   # Number of capture occasions
marked.j <- rep(200, n.occasions-1) # Annual number of newly marked juveniles
marked.a <- rep(30, n.occasions-1)  # Annual number of newly marked adults
phi.juv <- 0.3                      # Juvenile annual survival
phi.ad <- 0.65                      # Adult annual survival
p <- rep(0.5, n.occasions-1)        # Recapture
phi.j <- c(phi.juv, rep(phi.ad, n.occasions-2))
phi.a <- rep(phi.ad, n.occasions-1)

# Define matrices with survival and recapture probabilities
PHI.J <- matrix(0, ncol = n.occasions-1, nrow = sum(marked.j))

for (i in 1:length(marked.j)){
  PHI.J[(sum(marked.j[1:i])-marked.j[i]+1):sum(marked.j[1:i]),i:(n.occasions-1)] <- matrix(rep(phi.j[1:(n.occasions-i)],marked.j[i]), ncol = n.occasions-i, byrow = TRUE)
}
P.J <- matrix(rep(p, sum(marked.j)), ncol = n.occasions-1, nrow = sum(marked.j), byrow = TRUE)
PHI.A <- matrix(rep(phi.a, sum(marked.a)), ncol = n.occasions-1, nrow = sum(marked.a), byrow = TRUE)
P.A <- matrix(rep(p, sum(marked.a)), ncol = n.occasions-1, nrow = sum(marked.a), byrow = TRUE)

# Apply simulation function
CH.J <- simul.cjs(PHI.J, P.J, marked.j)
CH.A <- simul.cjs(PHI.A, P.A, marked.a) 

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f.j <- apply(CH.J, 1, get.first)
f.a <- apply(CH.A, 1, get.first)

# Create matrices X indicating age classes
x.j <- matrix(NA, ncol = dim(CH.J)[2]-1, nrow = dim(CH.J)[1])
for (i in 1:dim(CH.J)[1]){
  for (t in f.j[i]:(dim(CH.J)[2]-1)){
    x.j[i,t] <- 2
    x.j[i,f.j[i]] <- 1   
  } #t
} #i

x.a <- matrix(NA, ncol = dim(CH.A)[2]-1, nrow = dim(CH.A)[1])
for (i in 1:dim(CH.A)[1]){
  for (t in f.a[i]:(dim(CH.A)[2]-1)){
    x.a[i,t] <- 2
  } #t
} #i

CH <- rbind(CH.J, CH.A)
f <- c(f.j, f.a)
x <- rbind(x.j, x.a)

# Set known state for z, including dead recovery
get.last <- function(x) max(which(!is.na(x)))

zknown <- known.state.cjs(CH)
dead <- sample(dim(zknown)[1],1000) # Dead individuals
last_alive <- apply(zknown[dead,], 1, get.last) # Last year they were alive

z <- zknown
for (i in 1:length(dead)){
  if (last_alive[i] == -Inf | last_alive[i] == 10 ) next # If its only detected the first or last occasion, problems
    z[dead[i], (last_alive[i]+1):n.occasions] <- 0
  }

# Bundle data
bugs.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = z, x = x)

# Specify model in BUGS language
setwd("D:/MargSalas/Oso/SCR/Model")

cat("
model {
# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      phi[i,t] <- beta[x[i,t]]
      p[i,t] <- mean.p
      } #t
   } #i
for (u in 1:2){
   beta[u] ~ dunif(0, 1)              # Priors for age-specific survival
   }
mean.p ~ dunif(0, 1)                  # Prior for mean recapture
# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- 1
   for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
      } #t
   } #i
}",fill = TRUE, file = "cjs-age.txt")


# Initial values

cjs.init.zdead <- function(ch,f){
  for (i in 1:dim(ch)[1]){
    if (sum(ch[i,])==1) next
    n2 <- max(which(ch[i,] == 1 | z[i,] == 0))
    ch[i,f[i]:n2] <- NA
  }
  for (i in 1:dim(ch)[1]){
    ch[i,1:f[i]] <- NA
  }
  return(ch)
}# Latent state (include dead individuals)

inits <- function(){list(z = cjs.init.zdead(CH, f), beta = runif(2, 0, 1), mean.p = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("beta", "mean.p")

# MCMC settings
ni <- 2000
nt <- 3
nb <- 1000
nc <- 3

# Call WinBUGS from R (BRT 3 min)
cjs.age <- jags(bugs.data, inits, parameters, "cjs-age.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(cjs.age, digits = 3)

## ---- Data ----

# Load monitoring data and capture history
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
chb <- as.matrix(chb[,c(3:13)]) # Keep onlu from 2010 - 2020 (systematic sampling)

nind <- dim(chb)[1]
n.occasions <- dim(chb)[2]

# First observation
get.first <- function(x) min(which(x!=0)) # Important: Only initial values after first capture 
f_bear <- apply(chb, 1, get.first)

# Age class

setwd("D:/MargSalas/Oso/Datos")
info <- read.csv("Info_individuals.csv", header = TRUE, row.names = NULL, sep = ";") # Load info for year of birth
info <- info[,c(4,8,9)]

# Matrix with exact ages
x <- matrix(NA, nrow = nrow(chb), ncol = ncol(chb)) # Matrix to fill
colnames(x) <- colnames(chb)
rownames(x) <- rownames(chb)

for (i in 1:nrow(x)){
  birth <- as.numeric(info$Year_birth[which(info$ID %in% rownames(x)[i])])
  if (birth > 2010){
    x[i,which(colnames(x) %in% birth):n.occasions] <- 0:(n.occasions-which(colnames(x) %in% birth))
  } else { x[i,] <- (2010-birth):((2010-birth)+(n.occasions-1)) }
}

# Matrix with reclassified ages
x_cat <- matrix(NA, nrow = nrow(chb), ncol = ncol(chb)) # Matrix to fill
colnames(x_cat) <- colnames(chb)
rownames(x_cat) <- rownames(chb)

for (i in 1:length(x)){
  if(is.na(x[i])) next
  if(x[i] < 2){ x_cat[i] <- 1 }
  else if (x[i] > 1 & x[i] < 5){ x_cat[i] <- 2
  } else {x_cat[i] <- 3}
}
# Select from 2010 to 2019 because it is a co-variate on survival, 
# which happens between years
x_cat <- x_cat[,c(1:10)]

# Latent state: 
z <- known.state.cjs(chb) # Known state without including deaths

x_death <- matrix(NA, nrow = nrow(chb), ncol = ncol(chb)) # Create matrix to include year of death
colnames(x_death) <- colnames(chb)
rownames(x_death) <- rownames(chb)

for (i in 1:nrow(x_death)){
  death <- as.numeric(info$Year_death[which(info$ID %in% rownames(x)[i])])
  if (is.na(death) | death > 2019) next
  x_death[i,(which(colnames(x_death) %in% death)+1):n.occasions] <- 0}

z[which(x_death == 0)] <- 0 # Known state including deaths

# I don't know why, the model doesn't work if I include the deaths of the individuals
# detected once and death afterwards. 
det1 <- which(apply(chb,1,sum)==1) # Only 4 (Fifonet, Patoune, S28-SLO3, Soulane) from which I know the date of death
z[det1,] <- NA                    # So I don't introduce latent state on those for now


#### Model 2: Phi[age_class] & p[age_class] ####

# Bundle data
data_bear <- list(y = chb, f = f_bear, nind = dim(chb)[1], n.occasions = dim(chb)[2], z = z, x = x_cat)

# Model
setwd("D:/MargSalas/Oso/SCR/Model")

cat("
model {
# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      phi[i,t] <- beta.phi[x[i,t]]
      p[i,t] <- beta.p[x[i,t]]
      } #t
   } #i
for (u in 1:3){
   beta.phi[u] ~ dunif(0, 1)}              # Priors for age-specific survival
for (u in 1:3){
   beta.p[u] ~ dunif(0, 1)}              # Priors for age-specific p recapture
   
# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- 1
   for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
      } #t
   } #i
}",fill = TRUE, file = "cjs_phi[ageclass]_p[ageclass].txt")

# Initial values
inits <- function(){list(z = cjs.init.zdead(chb, f_bear), beta.phi = runif(3, 0, 1), beta.p = runif(3, 0, 1))}  

# Parameters monitored
parameters <- c("beta.phi", "beta.p")

# MCMC settings
ni <- 2000
nt <- 10
nb <- 1000
nc <- 3

# Call WinBUGS from R (BRT 3 min)
setwd("D:/MargSalas/Oso/SCR/Model")
outHairSt.ageclass <- jags(data_bear, inits, parameters, "cjs_phi[ageclass]_p[ageclass].txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(outHairSt.ageclass, digits = 3)

# Summarize posteriors 
par(mfrow = c(1,1)) #Phi
hist(outHairSt.ageclass$sims.list$beta.phi[,1], nclass = 50, col = adjustcolor("purple", alpha.f = 0.2), border =  adjustcolor("purple", alpha.f = 0.2), main = "", xlab = "Phi (Cub = purple; Subadult = green; Adult = orange)", las = 1, xlim = c(0.5,1))
abline(v = mean(outHairSt.ageclass$sims.list$beta.phi[,1]), col = "purple", lwd = 3)

hist(outHairSt.ageclass$sims.list$beta.phi[,2], nclass = 50, col = adjustcolor("green", alpha.f = 0.2), border = adjustcolor("green", alpha.f = 0.2),  add = TRUE)
abline(v = mean(outHairSt.ageclass$sims.list$beta.phi[,2]), col = "darkgreen", lwd = 3)

hist(outHairSt.ageclass$sims.list$beta.phi[,3], nclass = 50, col = adjustcolor("orange", alpha.f = 0.2), border = adjustcolor("orange", alpha.f = 0.2),  add = TRUE)
abline(v = mean(outHairSt.ageclass$sims.list$beta.phi[,3]), col = "darkorange", lwd = 3)

par(mfrow = c(1,1)) # P detection
hist(outHairSt.ageclass$sims.list$beta.p[,1], nclass = 50, col = adjustcolor("purple", alpha.f = 0.2), border =  adjustcolor("purple", alpha.f = 0.2), main = "", xlab = "P (Cub = purple; Subadult = green; Adult = orange)", las = 1, xlim = c(0.2,1))
abline(v = mean(outHairSt.ageclass$sims.list$beta.p[,1]), col = "purple", lwd = 3)

hist(outHairSt.ageclass$sims.list$beta.p[,2], nclass = 50, col = adjustcolor("green", alpha.f = 0.2), border = adjustcolor("green", alpha.f = 0.2),  add = TRUE)
abline(v = mean(outHairSt.ageclass$sims.list$beta.p[,2]), col = "darkgreen", lwd = 3)

hist(outHairSt.ageclass$sims.list$beta.p[,3], nclass = 50, col = adjustcolor("orange", alpha.f = 0.2), border = adjustcolor("orange", alpha.f = 0.2),  add = TRUE)
abline(v = mean(outHairSt.ageclass$sims.list$beta.p[,3]), col = "darkorange", lwd = 3)


# Save outputs
setwd("D:/MargSalas/Oso/SCR/Exp_analysis")
write.csv(outHairSt.ageclass$summary, file = "outHairSt_phi[ageclass]_p[ageclass]_DeathRecov.csv")

