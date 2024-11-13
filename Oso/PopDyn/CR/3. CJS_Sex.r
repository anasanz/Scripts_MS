## -------------------------------------------------
##    MODEL 7.5.1 BPA: CJS with fixed group effects 
##                  GROUP = SEX
## ------------------------------------------------- 


rm(list = ls())

library(dplyr)
library(tidyr)
library(rgdal)
library(raster)
library(rjags)
library(jagsUI)

setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Functions")
source("known.state.cjs.r")
source("simul.cjs.r")
source("cjs.init.z.r")
source("capt_hist_bear.r")


## ---- Simulation ----

# Define parameter values
n.occasions <- 12                  # Number of capture occasions
marked <- rep(30, n.occasions-1)   # Annual number of newly marked individuals
phi.f <- rep(0.65, n.occasions-1)  # Survival of females
p.f <- rep(0.6, n.occasions-1)     # Recapture of females
phi.m <- rep(0.8, n.occasions-1)   # Survival of males
p.m <- rep(0.3, n.occasions-1)     # Reacpture of males

# Define matrices with survival and recapture probabilities
PHI.F <- matrix(phi.f, ncol = n.occasions-1, nrow = sum(marked))
P.F <- matrix(p.f, ncol = n.occasions-1, nrow = sum(marked))
PHI.M <- matrix(phi.m, ncol = n.occasions-1, nrow = sum(marked))
P.M <- matrix(p.m, ncol = n.occasions-1, nrow = sum(marked))

# Simulate capture-histories

CH.F <- simul.cjs(PHI.F, P.F, marked)
CH.M <- simul.cjs(PHI.M, P.M, marked)

# Merge capture-histories by row
CH <- rbind(CH.F, CH.M)

# Create group variable to index phi and p into group
group <- c(rep(1, dim(CH.F)[1]), rep(2, dim(CH.M)[1]))

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)


# Bundle data
bugs.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH), g = length(unique(group)), group = group)


# Specify model in BUGS language
setwd("D:/MargSalas/Oso/SCR/Model")
cat("
model {

# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      phi[i,t] <- phi.g[group[i]]
      p[i,t] <- p.g[group[i]]
      } #t
   } #i
for (u in 1:g){
   phi.g[u] ~ dunif(0, 1)              # Priors for group-specific survival
   p.g[u] ~ dunif(0, 1)                # Priors for group-specific recapture
   }

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
}",fill = TRUE, file = "cjs-group.txt")


# Initial values

inits <- function(){list(z = cjs.init.z(CH, f), phi.g = runif(length(unique(group)), 0, 1), p.g = runif(length(unique(group)), 0, 1))}  

# Parameters monitored
parameters <- c("phi.g", "p.g")

# MCMC settings
ni <- 5000
nt <- 3
nb <- 2000
nc <- 3

# Call WinBUGS from R (BRT 2 min)
setwd("D:/MargSalas/Oso/SCR/Model")
cjs.group <- jags(bugs.data, inits, parameters, "cjs-group.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

# Summarize posteriors
print(cjs.group, digits = 3)


## ---- Data ----

# Load monitoring data and capture history
setwd("D:/MargSalas/Oso/Datos/Tablas_finales")
os <- read.csv("Data_os_96_20.csv", header = TRUE, row.names = NULL)
os <- os[,-1] 

os_id <- os[which(os$Confirmed_Individual != "Indetermined"), ] # Only identified

HairStCH <- capt_hist_bear(data = os_id, 
                           method = "Sampling_station", 
                           obs_type = c("Hair"))

# Create sex vector to index
setwd("D:/MargSalas/Oso/Datos")
info <- read.csv("Info_individuals.csv", header = TRUE, row.names = NULL, sep = ";") # Load info to join sex
info <- info[,c(4,5)]
colnames(info)[1] <- "Confirmed_Individual"

info2 <- left_join(HairStCH$capt.hist,info, by = "Confirmed_Individual")
info2 <- info2[-46,-2]
info2$Sex <- as.numeric(as.factor(info2$Sex))

sex <- NULL
for (i in 1:nrow(info2)){
  sex <- c(sex, rep(info2$Sex[i],ncol(info2)-2))
}


# Arrange data for JAGS
# Capture history
chb <- as.matrix(info2[,c(2:12)]) # Keep onlu from 2010 - 2020 (systematic sampling)

# First observation
get.first <- function(x) min(which(x!=0)) # Important: Only initial values after first capture 
f_bear <- apply(chb, 1, get.first)

# JAGS data
setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Functions")
source("known.state.cjs.r")

data_bear <- list(y = chb, f = f_bear, nind = dim(chb)[1], n.occasions = dim(chb)[2], z = known.state.cjs(chb), g = length(unique(sex)), group = sex)

# Initial values

setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Functions")
source("cjs.init.z.r")

inits <- function(){list(z = cjs.init.z(chb, f_bear), phi.g = runif(length(unique(sex)), 0, 1), p.g = runif(length(unique(sex)), 0, 1))}  

# Parameters monitored
parameters <- c("phi.g", "p.g")

# MCMC settings
ni <- 10000
nt <- 10
nb <- 5000
nc <- 3

# Call JAGS from R 
setwd("D:/MargSalas/Oso/SCR/Model")
outHairSt.sex <- jags(data_bear, inits, parameters, "cjs-group.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(outHairSt.sex, digits = 3)


# Summarize posteriors 
par(mfrow = c(1,1)) #Phi
hist(outHairSt.sex$sims.list$phi.g[,1], nclass = 50, col = adjustcolor("purple", alpha.f = 0.2), border =  adjustcolor("purple", alpha.f = 0.2), main = "", xlab = "Phi (F = purple; M = green)", las = 1, xlim = c(0.5,1))
abline(v = mean(outHairSt.sex$sims.list$phi.g[,1]), col = "purple", lwd = 3)

hist(outHairSt.sex$sims.list$phi.g[,2], nclass = 50, col = adjustcolor("green", alpha.f = 0.2), border = adjustcolor("green", alpha.f = 0.2),  add = TRUE)
abline(v = mean(outHairSt.sex$sims.list$phi.g[,2]), col = "darkgreen", lwd = 3)

par(mfrow = c(1,1)) # P detection
hist(outHairSt.sex$sims.list$p.g[,1], nclass = 50, col = adjustcolor("purple", alpha.f = 0.2), border =  adjustcolor("purple", alpha.f = 0.2), main = "", xlab = "P (F = purple; M = green)", las = 1, xlim = c(0.5,1))
abline(v = mean(outHairSt.sex$sims.list$p.g[,1]), col = "purple", lwd = 3)

hist(outHairSt.sex$sims.list$p.g[,2], nclass = 50, col = adjustcolor("green", alpha.f = 0.2), border = adjustcolor("green", alpha.f = 0.2),  add = TRUE)
abline(v = mean(outHairSt.sex$sims.list$p.g[,2]), col = "darkgreen", lwd = 3)


# Save outputs
setwd("D:/MargSalas/Oso/SCR/Exp_analysis")
write.csv(outHairSt.sex$summary, file = "outHairSt_sex.csv")
