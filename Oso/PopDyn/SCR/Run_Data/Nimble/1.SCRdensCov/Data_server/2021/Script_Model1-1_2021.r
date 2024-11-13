rm(list = ls())

# Load packages

library(nimble)
library(nimbleSCR)
library(parallel)

# Load data
#setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/1.SCRdensCov/Data_server/2021")
setwd("~/data/data/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/1.SCRdensCov/Data_server/2021")

covs <- c("forest", "slope", "logDistcore", "roads1")

for (xxx in 1:length(covs)) {
  
  load(paste("Data_Model1-1_", covs[xxx], ".RData", sep = ""))
  
  #### RUN IN PARALLEL ####
  detectCores()
  
  ##start cluster w/ 3 cores (for 3 chains)
  this_cluster <- makeCluster(3)
  
  old <- Sys.time()
  
  chain_output <- parLapply(cl = this_cluster, X = 1:3, 
                            fun = run_MCMC_allcode,      ##function in "Parallel Nimble function.R"
                            data = nimData,              ##your data list
                            code = model,   ##your model code
                            inits = inits,                 ##your inits function
                            constants = nimConstants,      ##your list of constants
                            params = params,               ##your vector with params to monitor
                            niter = 150000,                  ##iterations per chain
                            nburnin = 100000,                ##burn-in
                            nthin = 10,                  ##thinning, main parameters
                            z.in = z.in,
                            S.in = S.in,
                            sex.in = sex.in
  )
  new <- Sys.time() - old
  
  ## ALWAYS close cluster when model is done
  stopCluster(this_cluster)
  
  # Save
  
  #setwd(paste("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_covsp0/", covs[xxx], sep = ""))
  setwd(paste("~/data/data/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_covsp0/", covs[xxx], sep = ""))
  
  save(chain_output, file = paste("Results2021_Model1-1_", covs[xxx], ".RData", sep = ""))
}