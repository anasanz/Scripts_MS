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

run_MCMC_allcode <- function(seed, data, code, inits, constants, params, niter, nburnin, nthin, nthin2) {
  
  library(nimble)
  library(nimbleSCR)
  #inits <- inits
  Tt<- Tt
  
  
  #### add any other libraries we use, like nibleSCR
  
  ###################################################################################
  ## include custom functions, distributions - same as in "Nimble model and function"
  
  ## If we have no custom functions (ie, we only use the functions from the packages above), delete
  ## this all the way to the next row of hashtags
  
  
  #dBernM2 <- nimbleFunction( run = function(x = double(1), ##observation
  #                                          detectionProb=double(1),   ##p
  #                                          occProb = double(0), ##psi
  #                                          log = double(0)) {
  #  
  #  returnType(double())
  #  
  #  nocc<-length(x) ##number of occasions
  #  
  #  probDetectionHistoryGivenOccupied <- 1
  #  probDetectionHistoryGivenUnoccupied <- 1
  #  
  #  for (k in 1:nocc){
  #    
  #    if(x[k] == 1) {
  #      probDetectionHistoryGivenOccupied <-
  #        probDetectionHistoryGivenOccupied * detectionProb[k]
  #      probDetectionHistoryGivenUnoccupied <- 0
  #    } else {
  #      probDetectionHistoryGivenOccupied <-
  #        probDetectionHistoryGivenOccupied * (1-detectionProb[k])
  #    }
  #    
  #  }
  #  ans <- log(occProb *
  #               probDetectionHistoryGivenOccupied +
  #               (1-occProb) *
  #               probDetectionHistoryGivenUnoccupied)
  #  
  #  if(log) return(ans)
  #  return(exp(ans))
  #  
  #})
  #
  #
  #rBernM2 <- nimbleFunction( run = function(n = integer(0), ##default
  #                                          detectionProb=double(1),   ##p
  #                                          occProb = double(0)) {
  #  
  #  returnType(double(1))
  #  
  #  nocc<-length(detectionProb) ##number of occasions
  #  z<-rbinom(1, 1, occProb)
  #  x<-rbinom(nocc, 1, detectionProb*z)
  #  return(x)
  #  
  #})
  #
  #
  #
  #registerDistributions(list(dBernM2 = list(
  #  BUGSdist = "dBernM2(detectionProb, occProb)",
  #  Rdist = "dBernM2(detectionProb,occProb)",
  #  types = c('value = double(1)',
  #            'detectionProb = double(1)',
  #            'occProb = double(0)'
  #  )
  #))
  #)
  #
  ### pass functions to global environment
  #assign('dBernM2', dBernM2, envir = .GlobalEnv)
  #assign('rBernM2', rBernM2, envir = .GlobalEnv)
  #
  ###################################################################################
  ### implement model 
  ### objects need to refer to function definition
  
  #(1) set up model
  model <- nimbleModel(code = code, constants = constants, data=data, check = FALSE)
  
  #(2) Compile model in c++
  cmodel <- compileNimble(model)       
  
  # (3) Configure MCMC - on an uncompiled model
  conf.mcmc<-configureMCMC(model, monitors = params, thin = nthin)
  
  # (4) Build the MCMC sampler based on configurations
  mcmc <- buildMCMC(conf.mcmc)
  
  # (5) Compile sampler in c++ together with compiled model
  cmcmc <- compileNimble(mcmc, project = cmodel, resetFunctions = TRUE)
  
  # (6) Run - single chain!!!
  samp <- runMCMC(cmcmc, niter = niter, nburnin = nburnin, nchains=1, inits = inits) 
  
  ### define what function returns (MCMC samples)
  return(samp)
  
} #end wrapper function



