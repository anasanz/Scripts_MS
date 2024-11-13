

rm(list = ls())

# Load packages

library(nimble)
library(nimbleSCR)

# Load data
#setwd("/home/csic/dde/jol/Documentos/OPSCR/")
setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/3.openSCR_Age/Data_server")

load("Data_Model3-3.1_CYRIL_allparams.RData")

setwd("D:/MargSalas/Scripts_MS/Functions/Nimble")
source('dbinomLocal_normalBear_rbinom2.R')
source('getSexSigma.R')


model <- nimbleModel(modelcode, constants = nimConstants,
                     data=nimData, inits=inits(), check = FALSE)

cmodel <- compileNimble(model)      

conf.mcmc <- configureMCMC(model, monitors = params, thin=10, monitors2 = params2, thin2 = 20)

mcmc <- buildMCMC(conf.mcmc)

cmcmc <- compileNimble(mcmc, project = cmodel, resetFunctions = TRUE)

samp <- runMCMC(cmcmc, niter = 150000, nburnin = 100000, nchains=3, inits = inits) 
