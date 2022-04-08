
## -------------------------------------------------
##      CJS model on bear data (different methods)
##      7.3. Models with constant parameters (BPA)
## ------------------------------------------------- 

library(rjags)
library(jagsUI)


## ---- Simulation ----

# Define parameter values
n.occasions <- 6                   # Number of capture occasions
marked <- rep(50, n.occasions-1)   # Annual number of newly marked individuals
phi <- rep(0.65, n.occasions-1)
p <- rep(0.4, n.occasions-1)

# Define matrices with survival and recapture probabilities
PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked))
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked))

# Define function to simulate a capture-history (CH) matrix
simul.cjs <- function(PHI, P, marked){
  n.occasions <- dim(PHI)[2] + 1
  CH <- matrix(0, ncol = n.occasions, nrow = sum(marked))
  # Define a vector with the occasion of marking
  mark.occ <- rep(1:length(marked), marked[1:length(marked)])
  # Fill the CH matrix
  for (i in 1:sum(marked)){
    CH[i, mark.occ[i]] <- 1       # Write an 1 at the release occasion
    if (mark.occ[i]==n.occasions) next
    for (t in (mark.occ[i]+1):n.occasions){
      # Bernoulli trial: does individual survive occasion?
      sur <- rbinom(1, 1, PHI[i,t-1])
      if (sur==0) break		# If dead, move to next individual 
      # Bernoulli trial: is individual recaptured? 
      rp <- rbinom(1, 1, P[i,t-1])
      if (rp==1) CH[i,t] <- 1
    } #t
  } #i
  return(CH)
}

# Execute function
CH <- simul.cjs(PHI, P, marked)

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0)) # Important: Only initial values after first capture 
f <- apply(CH, 1, get.first) #             -> The initial capture process is not modeled at CJS

# Specify model in BUGS language
setwd("D:/MargSalas/Oso/Model")
cat("
model {

# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      phi[i,t] <- mean.phi
      p[i,t] <- mean.p
      } #t
   } #i

mean.phi ~ dunif(0, 1)         # Prior for mean survival
mean.p ~ dunif(0, 1)           # Prior for mean recapture

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
}",fill = TRUE, file = "cjs-c-c.txt")

# Function to create a matrix with information about known latent state z
# So that z doesn't need to be calculated at each iteration, we give the information that we have to the model
known.state.cjs <- function(ch){
  state <- ch
  for (i in 1:dim(ch)[1]){
    n1 <- min(which(ch[i,]==1))
    n2 <- max(which(ch[i,]==1))
    state[i,n1:n2] <- 1
    state[i,n1] <- NA
  }
  state[state==0] <- NA
  return(state)
}

# Bundle data
bugs.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH))

# Function to create a matrix of initial values for latent state z
# Important: - Only initial values after first capture(The initial capture process is not modeled at CJS)
#             - Because we gave information to z: We should not give initial values for those elements of z whose value was specified in the data

cjs.init.z <- function(ch,f){
  for (i in 1:dim(ch)[1]){
    if (sum(ch[i,])==1) next
    n2 <- max(which(ch[i,]==1))
    ch[i,f[i]:n2] <- NA
  }
  for (i in 1:dim(ch)[1]){
    ch[i,1:f[i]] <- NA
  }
  return(ch)
}

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("mean.phi", "mean.p")

# MCMC settings
ni <- 10000
nt <- 6
nb <- 5000
nc <- 3

# Call JAGS from R 
cjs.c.c <- jags(bugs.data, inits, parameters, "cjs-c-c.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

# Summarize posteriors
print(cjs.c.c, digits = 3)

## ---- Data ----

setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Functions")
source("capt_hist_bear.r")

# Load monitoring data and capture history
setwd("D:/MargSalas/Oso/Datos/Tablas_finales")
os <- read.csv("Data_os_96_20.csv", header = TRUE, row.names = NULL)
os <- os[,-1] 

os_id <- os[which(os$Confirmed_Individual != "Indetermined"), ] # Only identified

## ---- 1. Camera trap Sampling Station----

CamSampStCH <- capt_hist_bear(data = os_id, 
                              method = "Sampling_station", 
                              obs_type = c("Photo", "Photo/Video", "Video"))

# Arrange data for JAGS
# Capture history
chb <- as.data.frame(CamSampStCH$capt.hist)
chb <- as.matrix(chb[,c(4:14)]) # Keep onlu from 2010 - 2020 (systematic sampling)
chb <- chb[-which(apply(chb,1,sum) == 0), ] # Delete individuals only detected before 2010

# First observation
get.first <- function(x) min(which(x!=0)) # Important: Only initial values after first capture 
f_bear <- apply(chb, 1, get.first)

# JAGS data
data_bear <- list(y = chb, f = f_bear, nind = dim(chb)[1], n.occasions = dim(chb)[2], z = known.state.cjs(chb))

# Initial values
inits <- function(){list(z = cjs.init.z(chb, f_bear), mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("mean.phi", "mean.p")

# MCMC settings
ni <- 10000
nt <- 10
nb <- 5000
nc <- 3

# Call JAGS from R 
setwd("D:/MargSalas/Oso/Model")
outCamera <- jags(data_bear, inits, parameters, "cjs-c-c.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(outCamera, digits = 3)


# Save outputs
setwd("D:/MargSalas/Oso/Exp_analysis")
write.csv(CamSampStCH$capt.hist, file = "CamSampStCH.csv")
write.csv(outCamera$summary, file = "outCamera.csv")

CamSampStCH$Nident


## ---- 2. Hair Sampling Station----

HairStCH <- capt_hist_bear(data = os_id, 
                              method = "Sampling_station", 
                              obs_type = c("Hair"))

# Arrange data for JAGS
# Capture history
chb <- as.data.frame(HairStCH$capt.hist)
chb <- as.matrix(chb[,c(4:14)]) # Keep onlu from 2010 - 2020 (systematic sampling)

# First observation
get.first <- function(x) min(which(x!=0)) # Important: Only initial values after first capture 
f_bear <- apply(chb, 1, get.first)

# JAGS data
data_bear <- list(y = chb, f = f_bear, nind = dim(chb)[1], n.occasions = dim(chb)[2], z = known.state.cjs(chb))

# Initial values
inits <- function(){list(z = cjs.init.z(chb, f_bear), mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("mean.phi", "mean.p")

# MCMC settings
ni <- 10000
nt <- 10
nb <- 5000
nc <- 3

# Call JAGS from R 
setwd("D:/MargSalas/Oso/Model")
outHairSt <- jags(data_bear, inits, parameters, "cjs-c-c.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(outHairSt, digits = 3)


# Save outputs
setwd("D:/MargSalas/Oso/Exp_analysis")
write.csv(HairStCH$capt.hist, file = "HairStCH.csv")
write.csv(outHairSt$summary, file = "outHairSt.csv")

HairStCH$Nident
HairStCH$Nsurv
HairStCH$capt.hist

## ---- 3. Hair Transect ----

HairTransCH <- capt_hist_bear(data = os_id, 
                              method = "Transect", 
                              obs_type = c("Hair"))

# Arrange data for JAGS
# Capture history
chb <- as.data.frame(HairTransCH$capt.hist)
chb <- as.matrix(chb[,c(3:13)]) # Keep onlu from 2010 - 2020 (systematic sampling)

# First observation
get.first <- function(x) min(which(x!=0)) # Important: Only initial values after first capture 
f_bear <- apply(chb, 1, get.first)

# JAGS data
data_bear <- list(y = chb, f = f_bear, nind = dim(chb)[1], n.occasions = dim(chb)[2], z = known.state.cjs(chb))

# Initial values
inits <- function(){list(z = cjs.init.z(chb, f_bear), mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("mean.phi", "mean.p")

# MCMC settings
ni <- 10000
nt <- 10
nb <- 5000
nc <- 3

# Call JAGS from R 
setwd("D:/MargSalas/Oso/Model")
outHairTr <- jags(data_bear, inits, parameters, "cjs-c-c.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
print(outHairTr, digits = 3)

# Save outputs
setwd("D:/MargSalas/Oso/Exp_analysis")
write.csv(HairTransCH$capt.hist, file = "HairTransCH.csv")
write.csv(outHairTr$summary, file = "outHairTr.csv")
HairTransCH$Nident

