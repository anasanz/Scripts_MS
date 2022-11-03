
# Model for Cyril

rm(list = ls())

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_prefinal_1721/Cyril")

load("data.RData")


#(1) set up model

model <- nimbleModel(SCRhab.Open.diftraps.age, constants = nimConstants, 
                     data=nimData, inits=inits(), check = FALSE)

cmodel <- compileNimble(model)       

conf.mcmc <- configureMCMC(model, monitors = params, thin=10)

mcmc <- buildMCMC(conf.mcmc)

cmcmc <- compileNimble(mcmc, project = cmodel, resetFunctions = TRUE)

samp <- runMCMC(cmcmc, niter = 150000, nburnin = 100000, nchains=3, inits = inits) 


setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_prefinal_1721/Cyril")

save(samp, file = "sampOpenSCR_diftraps_age2.RData")


