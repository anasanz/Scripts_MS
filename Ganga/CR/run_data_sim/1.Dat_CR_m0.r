

## -------------------------------------------------
##                1. M0: Constant p
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
inits <- function() list(z = rep(1, nrow(yaug)), p = runif(1, 0, 1))

# Parameters monitored
params <- c("N", "p", "omega")

# MCMC settings
ni <- 2500
nt <- 2
nb <- 500
nc <- 3

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

