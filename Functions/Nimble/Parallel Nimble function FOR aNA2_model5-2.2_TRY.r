###########################################################################################################
#### compile all model setup steps into a single function that is then executed in parallel ###############

##seed: random number generator seed for each core
##data: list with data (ie, observations)
##code: Nimble model code
##inits: initial values function
##constants: list with constants for the model (indices, covariates, etc)
##params: character vector with parameters to monitor
##niter: number of iterations for each chain
## nburnin: how many iterations to discard as burnin
##nthin: thinning 


##NOTE: custom functions and distributions have to be defined inside the function and passed to the
##      global environment; these are the same functions as in "Nimble custom function.R"

#################################################################################################################
## I have added Tt as an argument to the function and am passing it to the global environment further down in the 
## function
## You will have to do the same for all R objects that you use inside the inits() function
## and then provide these to the function when you execute the parallel computing
#################################################################################################################

run_MCMC_allcode <- function(seed, data, code, inits, constants, params, niter, nburnin, nthin, nthin2, Tt, z.in, S.in.sc_coords, sex.in, dbinomLocal_normalBear) {
  
  library(nimble)
  library(nimbleSCR)
  #### add any other libraries we use, like nimbleSCR
  
  ## pass objects for inits to global environment
  assign('Tt', Tt, envir = .GlobalEnv)
  #assign('piAGE.in', piAGE.in, envir = .GlobalEnv)
  #assign('zstAGE', zstAGE, envir = .GlobalEnv)
  #assign('w.in', w.in, envir = .GlobalEnv)
  #assign('age.cat.in', age.cat.in, envir = .GlobalEnv)
  assign('S.in.sc_coords', S.in.sc_coords, envir = .GlobalEnv)
  assign('z.in', z.in, envir = .GlobalEnv)
  assign('sex.in', z.in, envir = .GlobalEnv)
  
  ### Add function that speeds up by including p calculation
  
  setwd("D:/MargSalas/Scripts_MS/Stats/Nimble")
  #setwd("~/data/data/Scripts_MS/Stats/Nimble")
  source('dbinomLocal_normalBear.R')
  
 
  ### pass functions to global environment
  assign('dbinomLocal_normalBear', dbinomLocal_normalBear, envir = .GlobalEnv)
  
  ###################################################################################
  ### implement model 
  ### objects need to refer to function definition
  
  #(1) set up model
  model <- nimbleModel(code = code, constants = constants, data=data, check = FALSE)
  
  #(2) Compile model in c++
  cmodel <- compileNimble(model)       
  
  # (3) Configure MCMC - on an uncompiled model
  conf.mcmc<-configureMCMC(model, monitors = params, thin=nthin)
  
  # (4) Build the MCMC sampler based on configurations
  mcmc <- buildMCMC(conf.mcmc)
  
  # (5) Compile sampler in c++ together with compiled model
  cmcmc <- compileNimble(mcmc, project = cmodel, resetFunctions = TRUE)
  
  # (6) Run - single chain!!!
  samp <- runMCMC(cmcmc, niter = niter, nburnin = nburnin, nchains=1, inits = inits) 
  
  ### define what function returns (MCMC samples)
  return(samp)
  
} #end wrapper function


