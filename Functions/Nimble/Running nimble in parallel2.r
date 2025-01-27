
############################################################################################
### running a model in parallel ############################################################

##source code to run model in parallel 
source("Parallel Nimble function FOR aNA.R") ##sorry, caps lock...

### set up data, constants, parameters, inits as for the simple model run

##NOTE: because all model setup steps have to be part of the parallel run, this isn't 3x faster than
##      the regular run and probably most beneficial for long runs, not for test runs


##start cluster w/ 3 cores (for 3 chains)
this_cluster <- makeCluster(3)


##run wrapper function in parallel - cl and X need to be this way
## cl defines which cluster to use, X provides random number seeds to each core

##added argument Tt at the end; you'll need to add all other objects used within inits

chain_output <- parLapply(cl = this_cluster, X = 1:3, 
                          fun = run_MCMC_allcode,      ##function in "Parallel Nimble function.R"
                          data = nimData,              ##your data list
                          code=SCRhab.Open.diftraps,   ##your model code
                          inits=inits,                 ##your inits function
                          constants=nimConstants,      ##your list of constants
                          params=params,               ##your vector with params to monitor
                          niter=100000,                  ##iterations per chain
                          nburnin=50000,                ##burn-in
                          nthin=10,                     ##thinning, main parameters
                          Tt=Tt				##additional objects needed within inits
                          
)

## ALWAYS close cluster when model is done
stopCluster(this_cluster)

### output is a list 





