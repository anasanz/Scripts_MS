## -------------------------------------------------
##      Size estimation PTS using Bayesian
##              Close population model
## ------------------------------------------------- 
rm(list=ls())
setwd("D:/MargSalas/Ganga/Data/CMR")

## Function to estimate the mode
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
## -------------------------------------------------
##                1. M0: Constant p
## ------------------------------------------------- 


library(jagsUI)
library(rjags)

## ---- Simulation ----

data.fn <- function(N = 100, p = 0.5, T = 3){
  yfull <- yobs <- array(NA, dim = c(N, T))
  for (j in 1:T){
    yfull[,j] <- rbinom(n = N, size = 1, prob = p)
  }
  ever.detected <- apply(yfull, 1, max)
  C <- sum(ever.detected)
  yobs <- yfull[ever.detected == 1,]
  cat(C, "out of", N, "animals present were detected.\n")
  return(list(N = N, p = p, C = C, T = T, yfull = yfull, yobs = yobs))
}

data <- data.fn()

str(data)

# Augment data set by 150 potential individuals
nz <- 150
yaug <- rbind(data$yobs, array(0, dim = c(nz, data$T)))

# Specify model in BUGS language
sink("model_m0.txt")
cat("
model {

# Priors
omega ~ dunif(0, 1)
p ~ dunif(0, 1)

# Likelihood
for (i in 1:M){
   z[i] ~ dbern(omega)			# Inclusion indicators (probability that exists)
   for (j in 1:T){
      yaug[i,j] ~ dbern(p.eff[i,j])
      p.eff[i,j] <- z[i] * p		# Can only be detected if z=1
      } #j
   } #i

# Derived quantities
N <- sum(z[])
}",fill = TRUE)
sink()

# Bundle data
win.data <- list(yaug = yaug, M = nrow(yaug), T = ncol(yaug))

# Initial values
inits <- function() list(z = rep(1, nrow(yaug)), p = runif(1, 0, 1))

# Parameters monitored
params <- c("N", "p", "omega")

# MCMC settings
ni <- 2500
nt <- 2
nb <- 500
nc <- 3

# Call WinBUGS from R (BRT <1 min)
out_sim_M0 <- jags(win.data, inits, params, "model_m0.txt", n.chains = nc, 
                   n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

# Summarize posteriors
print(out_sim_M0, dig = 3)
hist(out_sim_M0$sims.list$N, nclass = 50, col = "gray", main = "", xlab = "Population size", las = 1, xlim = c(80, 150))
abline(v = data$C, lwd = 3)
abline(v = mean(out_sim_M0$sims.list$N), col = "blue", lwd = 3)

## ---- Data ----

setwd("D:/MargSalas/Ganga/Data/CMR")

load("cr_sandgrouse_2022_new.RData")
capt.hist <- as.data.frame(capt.hist$ch)
colnames(capt.hist)[1] <- "ch"
n_oc <- 5

data_ganga <- matrix(data = as.numeric(do.call("rbind", strsplit(as.character(capt.hist$ch), "", fixed = TRUE))), nrow = nrow(capt.hist), ncol = n_oc)
T = ncol(data_ganga)

# Augment data set by 150 potential individuals
nz <- 150
yaug <- rbind(data_ganga, array(0, dim = c(nz, T)))

# Bundle data
win.data <- list(yaug = yaug, M = nrow(yaug), T = ncol(yaug))

# Call WinBUGS from R (BRT <1 min)
out_dat_m0 <- jags(win.data, inits, params, "model_m0.txt", n.chains = nc, 
                   n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

setwd("D:/MargSalas/Ganga/Results/CR/Preliminary_2022/new_samples")
write.csv(out_dat_m0$summary, file = "out_dat_m0.csv")

# Summarize posteriors 
print(out_dat_m0, dig = 3) # p very low, 0.13
hist(out_dat_m0$sims.list$N, nclass = 50, col = "gray", main = "m0", xlab = "Population size", las = 1, xlim = c(35, 210))
abline(v = nrow(data_ganga), lwd = 3)
abline(v = mean(out_dat_m0$sims.list$N), col = "blue", lwd = 3)
abline(v = getmode(out_dat_m0$sims.list$N), col = "red", lwd = 3)



## -------------------------------------------------
##            2. Mt: p varies each occasion
## ------------------------------------------------- 

## ---- Simulation ----

data.fn <- function(N = 100, mean.p = 0.5, T = 3, time.eff = runif(T, -2, 2)){
  yfull <- yobs <- array(NA, dim = c(N, T) )
  p.vec <- array(NA, dim = T)
  for (j in 1:T){
    p <- plogis(log(mean.p / (1-mean.p)) + time.eff[j])
    yfull[,j] <- rbinom(n = N, size = 1, prob = p)
    p.vec[j] <- p
  }
  ever.detected <- apply(yfull, 1, max)
  C <- sum(ever.detected)
  yobs <- yfull[ever.detected == 1,]
  cat(C, "out of", N, "animals present were detected.\n")
  return(list(N = N, p.vec = p.vec, C = C, T = T, yfull = yfull, yobs = yobs))
}

data <- data.fn()

# Augment data set
nz <- 150
yaug <- rbind(data$yobs, array(0, dim = c(nz, data$T)))

# Specify model in BUGS language
sink("model_mt.txt")
cat("
model {
# Priors
omega ~ dunif(0, 1)
for (i in 1:T){
   p[i] ~ dunif(0, 1)
   }

# Likelihood
for (i in 1:M){
   z[i] ~ dbern(omega)
   for (j in 1:T){
      yaug[i,j] ~ dbern(p.eff[i,j])
      p.eff[i,j] <- z[i] * p[j]
      } #j
   } #i

# Derived quantities
N <- sum(z[])
} # end model
",fill = TRUE)
sink()

# Bundle data
win.data <- list(yaug = yaug, M = nrow(yaug), T = ncol(yaug))

# Initial values
inits <- function() list(z = rep(1, nrow(yaug)), p = runif(data$T, 0, 1))

# Parameters monitored
params <- c("N", "p", "omega")

# MCMC settings
ni <- 2500
nt <- 2
nb <- 500
nc <- 3

# Call WinBUGS from R (BRT <1 min)
out_sim_mt <- jags(win.data, inits, params, "model_mt.txt", n.chains = nc, 
                   n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

# Summarize posteriors
print(out_sim_mt, dig = 3)
par(mfrow = c(1,1))
hist(out_sim_mt$sims.list$N, nclass = 40, col = "gray", main = "", xlab = "Population size", las = 1, xlim = c(70, 150))
abline(v = data$C, col = "black", lwd = 3)
abline(v = mean(out_sim_mt$sims.list$N), col = "blue", lwd = 3)


## ---- Data ----

setwd("D:/MargSalas/Ganga/Data/CMR")

load("cr_sandgrouse_2022_new.RData")
capt.hist <- as.data.frame(capt.hist$ch)
colnames(capt.hist)[1] <- "ch"
n_oc <- 5

data_ganga <- matrix(data = as.numeric(do.call("rbind", strsplit(as.character(capt.hist$ch), "", fixed = TRUE))), nrow = nrow(capt.hist), ncol = n_oc)
T = ncol(data_ganga)

# Augment data set by 150 potential individuals
nz <- 150
yaug <- rbind(data_ganga, array(0, dim = c(nz, T)))

# Bundle data
win.data <- list(yaug = yaug, M = nrow(yaug), T = ncol(yaug))

inits <- function() list(z = rep(1, nrow(yaug)), p = runif(T, 0, 1))

# Call WinBUGS from R (BRT <1 min)
out_dat_mt <- jags(win.data, inits, params, "model_mt.txt", n.chains = nc, 
                   n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

setwd("D:/MargSalas/Ganga/Results/CR/Preliminary_2022/new_samples")
write.csv(out_dat_mt$summary, file = "out_dat_mt.csv")

# Summarize posteriors 
print(out_dat_mt, dig = 3) # p very low, 0.13
hist(out_dat_mt$sims.list$N, nclass = 50, col = "gray", main = "MT", xlab = "Population size", las = 1, xlim = c(35, 200))
abline(v = nrow(data_ganga), lwd = 3)
abline(v = mean(out_dat_mt$sims.list$N), col = "blue", lwd = 3)
abline(v = getmode(out_dat_mt$sims.list$N), col = "red", lwd = 3)


#######  PLOT M0 AND MT  #######

setwd("D:/MargSalas/Ganga/Results/CR/Preliminary_2022/new_samples")

pdf("Posterior distribution new samples.pdf", 7,5)
par(mfrow = c(1,1))
hist(out_dat_m0$sims.list$N, nclass = 50, col = "gray", main = "m0", xlab = "Population size", las = 1, xlim = c(35, 210))
abline(v = nrow(data_ganga), lwd = 3)
abline(v = mean(out_dat_m0$sims.list$N), col = "blue", lwd = 3)
abline(v = getmode(out_dat_m0$sims.list$N), col = "red", lwd = 3)

hist(out_dat_mt$sims.list$N, nclass = 50, col = "gray", main = "Mt", xlab = "Population size", las = 1, xlim = c(35, 200))
abline(v = nrow(data_ganga), lwd = 3)
abline(v = mean(out_dat_mt$sims.list$N), col = "blue", lwd = 3)
abline(v = getmode(out_dat_mt$sims.list$N), col = "red", lwd = 3)

dev.off()

## -------------------------------------------------
##            3. Mb: trap response 
##  (p varies with number of times/if an individual was captured before)
##      This is not our case, our sampling occasions are independent
## ------------------------------------------------- 

## -------------------------------------------------
##            4. Mh: p varies per individual
##   Not enough N, and not enough recaptures of id
## ------------------------------------------------- 
## ---- Simulation ----

# Define function to simulate data under Mh
data.fn <- function(N = 100, mean.p = 0.5, T = 3, sd = 1){
  yfull <- yobs <- array(NA, dim = c(N, T) )
  mean.lp <- log(mean.p / (1-mean.p))
  p.vec <- plogis(mean.lp+ rnorm(N, 0, sd))
  
  for (i in 1:N){
    yfull[i,] <- rbinom(n = T, size = 1, prob = p.vec[i])
  }
  
  ever.detected <- apply(yfull, 1, max)
  C <- sum(ever.detected)
  yobs <- yfull[ever.detected == 1,]
  cat(C, "out of", N, "animals present were detected.\n")
  hist(p.vec, xlim = c(0,1), nclass = 20, col = "gray", main = "", xlab = "Detection probability", las = 1)
  return(list(N = N, p.vec = p.vec, mean.lp = mean.lp, C = C, T = T, yfull = yfull, yobs = yobs))
}

data <- data.fn()

# Aggregate capture-histories and augment data set
y <- sort(apply(data$yobs, 1, sum), decreasing = TRUE)
nz <- 300
yaug <- c(y, rep(0, nz))
yaug

# Specify model in BUGS language
sink("model_mh.txt")
cat("
model {

# Priors
omega ~ dunif(0, 1)
mean.lp <- logit(mean.p)
mean.p ~ dunif(0, 1)
tau <- 1 / (sd * sd)
sd ~ dunif(0, 5)

# Likelihood
for (i in 1:M){
   z[i] ~ dbern(omega)
   logit(p[i]) <- eps[i]
   eps[i] ~ dnorm(mean.lp, tau)	
   p.eff[i] <- z[i] * p[i]
   y[i] ~ dbin(p.eff[i], T)
   }

# Derived quantities
N <- sum(z[])}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(y = yaug, M = length(yaug), T = ncol(data$yobs))

# Initial values
inits <- function() list(z = rep(1, length(yaug)), sd = runif(1, 0.1, 0.9))

# Parameters monitored
params <- c("N", "mean.p", "sd", "omega")

# MCMC settings
ni <- 25000
nt <- 2
nb <- 5000
nc <- 3

# Call WinBUGS from R (BRT 6 min)
out_sim_mh <- jags(win.data, inits, params, "model_mh.txt", n.chains = nc, 
                   n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

# Summarize posteriors
print(out_sim_mh, dig = 3)
hist(out_sim_mh$sims.list$N, nclass = 50, col = "gray", main = "", xlab = "Population size", las = 1, xlim = c(80, 250))
abline(v = data$C, col = "black", lwd = 3)
abline(v = mean(out_sim_mh$sims.list$N), col = "blue", lwd = 3)

## ---- Data ----

setwd("D:/MargSalas/Ganga/Data/CMR")

load("cr_sandgrouse_new.RData")
capt.hist <- as.data.frame(capt.hist$ch)
colnames(capt.hist)[1] <- "ch"

data_ganga <- matrix(data = as.numeric(do.call("rbind", strsplit(as.character(capt.hist$ch), "", fixed = TRUE))), nrow = nrow(capt.hist), ncol = 3)
T = ncol(data_ganga)

# Augment data set by 150 potential individuals
nz <- 150
yaug <- rbind(data_ganga, array(0, dim = c(nz, T)))

# Bundle data
win.data <- list(yaug = yaug, M = nrow(yaug), T = ncol(yaug))

# Initial values
inits <- function() list(z = rep(1, nrow(yaug)), sd = runif(1, 0.1, 0.9))


# Call WinBUGS from R (BRT <1 min)
out_dat_mh <- jags(win.data, inits, params, "model_mh.txt", n.chains = nc, 
                   n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

# Summarize posteriors 
print(out_dat_mh, dig = 3) # Super super impreciso, lo he probado pero no tiene sentido meterlo para cada individuo
hist(out_dat_mh$sims.list$N, nclass = 50, col = "gray", main = "", xlab = "Population size", las = 1, xlim = c(0, 200))
abline(v = nrow(data_ganga), lwd = 3)
abline(v = mean(out_dat_mh$sims.list$N), col = "blue", lwd = 3) 
