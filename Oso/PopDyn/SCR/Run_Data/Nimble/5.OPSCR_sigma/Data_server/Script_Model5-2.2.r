rm(list = ls())

# Load packages

library(nimble)
library(nimbleSCR)
library(parallel)

setwd("D:/MargSalas/Scripts_MS/Stats/Nimble")
#source('dbinomLocal_normalBear.R')
source('dbinomLocal_normalBear_rbinom2.R')


# Load data
setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/5.OPSCR_sigma/Data_server")
#setwd("~/data/data/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/5.OPSCR_sigma/Data_server")

load("Data_Model5-2.2.RData")

#### RUN IN PARALLEL ####
detectCores()

##start cluster w/ 3 cores (for 3 chains)
this_cluster <- makeCluster(3)

old <- Sys.time()

chain_output <- parLapply(cl = this_cluster, X = 1:3, 
                          fun = run_MCMC_allcode,      ##function in "Parallel Nimble function.R"
                          data = nimData,              ##your data list
                          code = modelcode,   ##your model code
                          inits = inits,                 ##your inits function
                          constants = nimConstants,      ##your list of constants
                          params = params,               ##your vector with params to monitor
                          niter = 150000,                  ##iterations per chain
                          nburnin = 100000,                ##burn-in
                          nthin = 10,                  ##thinning, main parameters
                          Tt = Tt,                     ##additional objects needed within inits
                          z.in = z.in,
                          S.in.sc_coords = S.in.sc_coords,
                          sex.in = sex.in
                          
)
new <- Sys.time() - old

## ALWAYS close cluster when model is done
stopCluster(this_cluster)

# Save

#setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/4.openSCRdenscov_Effort")
setwd("~/data/data/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/5.OPSCR_sigma")
save(chain_output, file = "Results_Model5-2.2.RData")