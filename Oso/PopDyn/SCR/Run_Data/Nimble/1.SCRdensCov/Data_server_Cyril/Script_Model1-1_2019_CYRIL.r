rm(list = ls())

# Load packages

library(nimble)
library(nimbleSCR)
library(parallel)

# Load data
setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/1.SCRdensCov/Data_server_Cyril")


covs <- c("forest", "slope", "logDistcore", "roads1")

for (xxx in 1:length(covs)) {
  
  load(paste("Data_Model1-1_2019_", covs[xxx], ".RData", sep = ""))
  
  model <- nimbleModel(model, constants = nimConstants, 
                       data=nimData, inits=inits(), check = FALSE)
  
  cmodel <- compileNimble(model)       
  
  conf.mcmc <- configureMCMC(model, monitors = params, thin=10)
  
  mcmc <- buildMCMC(conf.mcmc)
  
  cmcmc <- compileNimble(mcmc, project = cmodel, resetFunctions = TRUE)
  
  samp <- runMCMC(cmcmc, niter = 150000, nburnin = 100000, nchains=3, inits = inits) 
  
  # Save
  
  setwd(paste("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_covsp0/", covs[xxx], sep = ""))

  save(chain_output, file = paste("Results2019_Model1-1_", covs[xxx], ".RData", sep = ""))
}