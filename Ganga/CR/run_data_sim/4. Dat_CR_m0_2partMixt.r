

## -------------------------------------------------
##                1. M0: Constant p in time
##                    Sex-specific p
## ------------------------------------------------- 

rm(list=ls())

library(jagsUI)
library(rjags)

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


# Initial values
inits <- function() list(z = rep(1, nrow(yaug)), p = runif(2, 0, 1),
                         group = c(rbinom(nrow(data_ganga) + nz,1,0.5)))


# Parameters monitored
params <- c("N", "p", "omega", "psi")

# MCMC settings
ni <- 20000
nt <- 2
nb <- 1000
nc <- 3

# Call WinBUGS from R (BRT <1 min)
out_sim_M0_2partMixt <- jags(win.data, inits, params, "model_m0_2partMixt.txt", n.chains = nc, 
                        n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

# There is not a group with higher p than the other, and estimates are not that different.
# So better to use the sex model (even if I have to double check with RS the sim script)
