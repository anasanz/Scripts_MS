

## -------------------------------------------------
##                1. M0: Constant p in time
##                    Sex-specific p
## ------------------------------------------------- 

rm(list=ls())

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

# Load sex

setwd("D:/MargSalas/Ganga/Data/CMR")
load("id_sex_sandgrouse_2022.RData") # It is in the same order than capt.hist
id_sex$sex[which(id_sex$sex == "F")] <- 0
id_sex$sex[which(id_sex$sex == "M")] <- 1
id_sex$sex[which(id_sex$sex == "X")] <- NA

sex <- as.numeric(id_sex$sex)
sexAug <- c(sex, rep(NA, nz))

# Bundle data
win.data <- list(yaug = yaug, M = nrow(yaug), T = ncol(yaug), sex = sexAug)


# Initial values
inits <- function() list(z = rep(1, nrow(yaug)), p = runif(2, 0, 1),
                         sex = c(rep(NA,nrow(data_ganga)), rbinom(nz,1,0.5)))

# Parameters monitored
params <- c("N", "p", "omega", "psi")

# MCMC settings
ni <- 2500
nt <- 2
nb <- 500
nc <- 3

# Call WinBUGS from R (BRT <1 min)
out_sim_M0_psex <- jags(win.data, inits, params, "model_m0_pSex.txt", n.chains = nc, 
                        n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

