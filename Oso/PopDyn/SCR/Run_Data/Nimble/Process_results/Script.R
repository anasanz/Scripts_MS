##########################################
###--- LOAD LIBRARIES AND FUNCTION 
##########################################
rm(list=ls())
library("nimble")
library(nimbleSCR)

setwd("/mnt/users/cymi/Ana/")
source("SaveStateModel.R")
path <- "/mnt/users/cymi/Ana/NimbleInFiles/"
path2 <- "/mnt/users/cymi/Ana/NimbleProcFiles/"
path3 <- "/mnt/users/cymi/Ana/NimbleOutFiles/"



##########################################
###--- FITTING THE MODELS
##########################################
pausesecs <- sample(seq(5,10, by=0.1), 1)
Sys.sleep(pausesecs)


## ----- PICK A FILE AT RANDOM -----

file.list <- list.files(path)
#save(file.list,file="file.list.RData")
#while(length(file.list)>0  ){
set <- sample(file.list, 1)
print(set)

## ----- LOAD THE FILE -----
load(paste(path, set, sep=""))
file.rename(from = paste(path, set, sep=""), to = paste(path2, set, sep=""))
file.remove(paste(path, set, sep=""))                                    #---ACTIVATE!


## ----- CREATE NIMBLE OBJETS -----
#nimbleOptions(clearNimbleFunctionsAfterCompiling = TRUE)
model <- nimbleModel( code = SCRhab.Open.diftraps.age
                      , constants = nimConstants
                      , data = nimData
                      , inits = inits()
                      , check = FALSE       
                      , calculate = FALSE)

cmodel <- compileNimble(model,showCompilerOutput = T)
calc <- cmodel$calculate()  
print(calc)

if(calc =="-Inf"){
   file.remove(paste(path2, set, sep=""))                                     
   stop("There is an INF")
}

conf <- configureMCMC(model,control = list(reflective = TRUE, adaptScaleOnly = TRUE),
                      useConjugacy = FALSE)
conf$addMonitors(params)   ##


Rmcmc <- buildMCMC(conf)

compiledList <- compileNimble(list(model=model, mcmc=Rmcmc))#,resetFunctions = TRUE
Cmodel <- compiledList$model
Cmcmc <- compiledList$mcmc

## ----- GENERATE MCMC SAMPLES -----
bite.size <- 1000# number of iterations in each bite to be exported and cleared
#set.seed(0)


for(nb in 1:150){
   print(nb)
   ####################################
   ## save samples1000, or otherwise ##
   ####################################
   
   ## now try to free up memory:
   ## remove the samples1000 variable,
   ## reduce the internal mvSamples object to 0 rows,
   ## and run R's garbage collector
   
   ## continue same run of the MCMC for another 1000 iterations, using reset = FALSE
   # print memory info
   memInfo <- gc()
   print(memInfo) 
   
   if(nb==1){
      stm <- proc.time()
      Cmcmc$run(bite.size)
      RunTime <- proc.time()-stm 
   }
   
   if(nb>1){
      stm <- proc.time()
      Cmcmc$run(bite.size, reset = FALSE)
      RunTime <- proc.time()-stm 
      
      #remove the last .rds file 
      pathrds <- file.path(path3, dirName, fileNameRds)
      file.remove(pathrds)                                     #---ACTIVATE!
      
   }
   this.sample <- as.matrix(Cmcmc$mvSamples)

   
   
   ## EXPORT NIMBLE OUTPUT WITH NAMES 
   dirName <- paste("NimbleOutFOR", set, sep="")
   fileName <- paste("bite_", nb, ".RData", sep="")
   fileNameRds <- paste("bite_", nb, ".rds", sep="")
   
   
   path <- file.path(path3, dirName)
   dir.create(path)
   

   path <- file.path(path3, dirName, fileName)
   save(this.sample, RunTime, file=path)
   
   # SAVE STATE OF THE MODEL 

   stateList <- list(modelState = getModelState(Cmodel),
                     mcmcState = getMCMCstate(conf, Cmcmc))
   path1 <- file.path(path3, dirName, fileNameRds)
   saveRDS(stateList, file = path1 )
   
   file.remove(paste(path2, set, sep=""))                                     #---ACTIVATE!
   
   rm("this.sample")
   Cmcmc$mvSamples$resize(0)
   gc()
   
}


rm(list=ls())

