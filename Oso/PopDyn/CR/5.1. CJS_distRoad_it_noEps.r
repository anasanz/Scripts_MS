## -------------------------------------------------
##    MODEL 7.4.3 BPA: CJS with temporal covariates 
##          Variation: Temporal + ID co-variate 
##        Covariate: Distance to roads --> REMOVE EPSILON (the model can't estimate it anyway)
## -------------------------------------------------


rm(list = ls())

library(dplyr)
library(tidyr)
library(rgdal)
library(raster)
library(rjags)
library(jagsUI)

## ---- SIMULATION ----

# Define parameter values
n.occasions <- 20                  # Number of capture occasions
marked <- rep(15, n.occasions-1)   # Annual number of newly marked individuals
mean.phi <- 0.65
p <- rep(0.4, n.occasions-1)
beta <- -0.3                       # Slope of survival-winter relationship	
#r.var <- 0.2                       # Residual temporal variance

# Draw annual survival probabilities

#eps <- matrix(rnorm((n.occasions-1), 0, 1^0.5), nrow = sum(marked), ncol = n.occasions-1, byrow = TRUE)
winter <- matrix(rnorm((n.occasions-1)*sum(marked), 0, 1^0.5), nrow = sum(marked), ncol = n.occasions-1)
logit.phi <- qlogis(mean.phi) + beta*winter 
              #+ eps
phi <- plogis(logit.phi)

# Define matrices with survival and recapture probabilities
PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked), byrow = TRUE)
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked))

# Simulate capture-histories

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

CH <- simul.cjs(PHI, P, marked)

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)


# Specify model in BUGS language

setwd("D:/MargSalas/Oso/SCR/Model")
cat("
model {

# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      logit(phi[i,t]) <- mu + beta*x[i,t] 
      # + epsilon[t] 
      p[i,t] <- mean.p
      } #t
   } #i
   
  #for (t in 1:(n.occasions-1)){
     #epsilon[t] ~ dnorm(0, tau)
     #phi.est[t] <- 1 / (1+exp(-mu-beta*x[t]-epsilon[t])) # Yearly survival
     #} #t

mu ~ dnorm(0, 0.001)                     # Prior for logit of mean survival
mean.phi <- 1 / (1+exp(-mu))             # Logit transformation
beta ~ dnorm(0, 0.001)I(-10, 10)         # Prior for slope parameter
# sigma ~ dunif(0, 10)                     # Prior on standard deviation
# tau <- pow(sigma, -2)
# sigma2 <- pow(sigma, 2)                  # Residual temporal variance
mean.p ~ dunif(0, 1)                     # Prior for mean recapture

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
}",fill = TRUE, file = "cjs-cov-raneff[it]_noEps.txt")

# Bundle data

setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Functions")
source("known.state.cjs.r")

bugs.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH), x = winter)

# Initial values
setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Functions")
source("cjs.init.z.r")


inits <- function(){list(z = cjs.init.z(CH, f), mu = rnorm(1), beta = runif(1, -5, 5), mean.p = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("mean.phi", "mean.p", "beta")

# MCMC settings
ni <- 50000
nt <- 10
nb <- 10000
nc <- 3

# Call WinBUGS from R (BRT 12 min)

setwd("D:/MargSalas/Oso/SCR/Model")
cjs.cov.it.noEps <- jags(bugs.data, inits, parameters, "cjs-cov-raneff[it]_noEps.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

# Summarize posteriors+
print(cjs.cov.it.noEps, digits = 3)

# Produce graph
par(mfrow = c(1, 2), las = 1)
hist(cjs.cov.it.noEps$sims.list$beta, nclass = 25, col = "gray", main = "", xlab = expression(beta), ylab = "Frequency")
abline(v = -0.3, col = "red", lwd = 2)
hist(cjs.cov.it.noEps$sims.list$sigma2, nclass = 50, col = "gray", main = "", xlab = expression(sigma^2), ylab = "Frequency", xlim=c(0, 3))
abline(v = 0.2, col = "red", lwd = 2)


## ---- DATA ----

# Load monitoring data

setwd("D:/MargSalas/Oso/Datos/Tablas_finales")
os <- read.csv("Data_os_96_20.csv", header = TRUE, row.names = NULL)
os <- os[,-1] 

os_id <- os[which(os$Confirmed_Individual != "Indetermined"), ] # Only identified
os_id$Method[which(os_id$Method == "sampling_station")] <- "Sampling_station"
os_id$Obs_type[which(os_id$Obs_type == "hair")] <- "Hair"

os_id <- os_id[which(os_id$Region == "Catalunya"), ]

coordinates(os_id) <- os_id[,c("x_long","y_lat")] # Spatial object
os_id@proj4string <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# Load and calculate distance to roads (1º and 2º)

# road1 <- raster("D:/MargSalas/Oso/Datos/GIS/Variables/Catalunya/Autopistes_ln.tif")
# road2 <- raster("D:/MargSalas/Oso/Datos/GIS/Variables/Catalunya/ViaPref_ln.tif")
# road3 <- raster("D:/MargSalas/Oso/Datos/GIS/Variables/Catalunya/ViaConven_ln.tif")
# 
# road1_4326 <- projectRaster(road1, crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
# road2_4326 <- projectRaster(road2, crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
# road3_4326 <- projectRaster(road3, crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
# 
# 
# asph <- brick(road1_4326, road2_4326, road3_4326)
# dist_asphalted <- calc(asph, function(x){min(x)})
# 
# setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Catalunya")
# writeRaster(dist_asphalted, filename = 'dist_roads_4326', format = 'GTiff')

# Extract coordinates
setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Catalunya")
road <- raster("D:/MargSalas/Oso/Datos/GIS/Variables/Catalunya/dist_roads_4326.tif")

coord <- os_id[ ,c("x_long","y_lat")] 
cells <- cellFromXY(road, coord) # 1. Tells the number of the cells where the coord. fall
dist_road <- road[cells] # 2. Gets the raster values of those cells

# Capture history

setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Functions")
source("capt_hist_bear.r")

d <- os_id@data

HairStCH <- capt_hist_bear(data = d, 
                           method = "Sampling_station", 
                           obs_type = c("Hair"))

ch <- data.frame(HairStCH$capt.hist[,c(1,3:13)])
rownames(ch) <- ch$Confirmed_Individual
ch <- ch[,-c(1)]
colnames(ch) <- c(2010:2020)

# --- Mean distance to paths per year and individual ----

# Data frame to fill

dist <- matrix(NA, nrow = nrow(ch), ncol = ncol(ch))
rownames(dist) <- rownames(ch)
colnames(dist) <- colnames(ch)

# Calculate distances

d <- cbind(d, dist_road)

dist_sum <- d %>%
  group_by(Confirmed_Individual,Year) %>%
  summarise(
    mean = mean(dist_road, na.rm = TRUE))

for (i in 1:nrow(dist_sum)){
  dist[which(rownames(dist) %in% dist_sum$Confirmed_Individual[i]),which(colnames(dist) %in% dist_sum$Year[i])] <- dist_sum$mean[i] 
}

# Select from 2010 to 2019 because it is a co-variate on survival, 
# which happens between years

dist <- dist[,c(1:10)]

# ---- Arrange data for JAGS ----

# Capture history
chb <- as.data.frame(HairStCH$capt.hist)
chb <- as.matrix(chb[,c(3:13)]) # Keep onlu from 2010 - 2020 (systematic sampling)

# Covariate distance to path

# There can't be NA in co-variate, so simulate values
mean_years <- apply(dist,2,mean, na.rm = TRUE)
sd_years <- apply(dist,2,sd, na.rm = TRUE)
number <- colSums(is.na(dist)) # How many values generate

for (i in 1:length(number)){
  col <- dist[,i] 
  col[is.na(col)] <- rnorm(number[i],mean_years[i],sd_years[i])
  dist[,i] <- col
}

# Generate randomly (not from the distribution)
random_dist <- unique(matrix(dist))[-1]
dist[which(is.na(dist))] <- sample(random_dist, length(which(is.na(dist))), replace = TRUE) # No NA in covariate


# Standardize
dist_mean <- mean(dist)
dist_sd <- sd(dist)
dist_sc <- (dist - dist_mean) / dist_sd


# First observation
get.first <- function(x) min(which(x!=0)) # Important: Only initial values after first capture 
f_bear <- apply(chb, 1, get.first)

# JAGS data
setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Functions")
source("known.state.cjs.r")

data_bear <- list(y = chb, f = f_bear, nind = dim(chb)[1], n.occasions = dim(chb)[2], z = known.state.cjs(chb), x = dist_sc)

# Specify model in BUGS language

# Specify model in BUGS language

setwd("D:/MargSalas/Oso/SCR/Model")
cat("
model {

# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      logit(phi[i,t]) <- mu + beta*x[i,t] 
      # + epsilon[t] 
      p[i,t] <- mean.p
      } #t
   } #i
   
  #for (t in 1:(n.occasions-1)){
     #epsilon[t] ~ dnorm(0, tau)
     #phi.est[t] <- 1 / (1+exp(-mu-beta*x[t]-epsilon[t])) # Yearly survival
     #} #t

mu ~ dnorm(0, 0.001)                     # Prior for logit of mean survival
mean.phi <- 1 / (1+exp(-mu))             # Logit transformation
beta ~ dnorm(0, 0.001)I(-10, 10)         # Prior for slope parameter
# sigma ~ dunif(0, 10)                     # Prior on standard deviation
# tau <- pow(sigma, -2)
# sigma2 <- pow(sigma, 2)                  # Residual temporal variance
mean.p ~ dunif(0, 1)                     # Prior for mean recapture

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
}",fill = TRUE, file = "cjs-cov-raneff[it]_noEps.txt")


# Initial values

setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Functions")
source("cjs.init.z.r")

inits <- function(){list(z = cjs.init.z(chb, f_bear), mu = rnorm(1), sigma = runif(1, 0, 5), beta = runif(1, -5, 5), mean.p = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("mean.phi", "mean.p", "beta")

# MCMC settings
ni <- 20000
nt <- 10
nb <- 10000
nc <- 3

# Call JAGS from R 
setwd("D:/MargSalas/Oso/SCR/Model")

outHairSt_distRoad_it_noEps<- jags(data_bear, inits, parameters, "cjs-cov-raneff[it]_noEps.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(outHairSt_distRoad_it_noEps, digits = 3)

# Summarize posteriors 
par(mfrow = c(1,1))
hist(outHairSt_distRoad_it_noEps$sims.list$mean.phi, nclass = 50, col = "gray", main = "", xlab = "Phi", las = 1, xlim = c(0,1))
abline(v = mean(outHairSt_distRoad_it_noEps$sims.list$mean.phi), col = "blue", lwd = 3)

hist(outHairSt_distRoad_it_noEps$sims.list$beta, nclass = 50, col = "gray", main = "", xlab = "Dist_road", las = 1, xlim = c(-2,10))
abline(v = mean(outHairSt_distRoad_it_noEps$sims.list$beta), col = "blue", lwd = 3)

# Save outputs
setwd("D:/MargSalas/Oso/SCR/Exp_analysis")
write.csv(outHairSt_distRoad_it_noEps$summary, file = "outHairSt_Cat_DistRoad[it]_noEps.csv")

