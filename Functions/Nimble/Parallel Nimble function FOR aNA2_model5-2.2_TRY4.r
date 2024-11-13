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

run_MCMC_allcode <- function(seed, data, code, inits, constants, params, niter, nburnin, nthin, nthin2, Tt, z.in, S.in.sc_coords, sex.in) {
  
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
  
  dbinomLocal_normalBear <- nimbleFunction(
    
    run = function( x = double(1),
                    
                    # detNums = double(0, default = -999),
                    
                    # detIndices = double(1),
                    
                    size = double(1),
                    
                    p0 = double(0, default = -999),
                    
                    #p0Traps = double(1),
                    
                    sigma = double(0),
                    
                    s = double(1),
                    
                    trapCoords = double(2),
                    
                    localTrapsIndices = double(2),
                    
                    localTrapsNum = double(1),
                    
                    #resizeFactor = double(0, default = 1),
                    
                    habitatGrid = double(2),
                    
                    indicator = double(0),
                    
                    lengthYCombined = double(0, default = 0),
                    
                    trapCovs = double(2),
                    
                    trapBetas = double(1),
                    
                    indTrapCovs = double(2),
                    
                    indTrapBetas = double(1),
                    
                    log = integer(0, default = 0)
                    
    ) {
      
      ## Specify return type
      
      returnType(double(0))
      
      
      
      ## Deal with cases where detection info is combined in one vector
      
      # if(detNums==-999){
      
      detNums <- x[1]
      
      nMaxDetectors <- (lengthYCombined-1)/2
      
      detIndices1 <- x[(nMaxDetectors+2):lengthYCombined]
      
      x1 <- x[2:(nMaxDetectors+1)]
      
      # }else{
      
      #   x1 <- x
      
      #   detIndices1 <- detIndices
      
      # }
      
      
      
      ## Shortcut if the current individual is not available for detection
      
      if(indicator == 0){
        
        if(detNums == 0){
          
          if(log == 0) return(1.0)
          
          else return(0.0)
          
        } else {
          
          if(log == 0) return(0.0)
          
          else return(-Inf)
          
        }
        
      }
      
      resizeFactor <- 1
      
      ## Retrieve the index of the habitat cell where the current AC is
      
      sID <- habitatGrid[trunc(s[2]/resizeFactor)+1, trunc(s[1]/resizeFactor)+1]
      
      
      
      ## Retrieve the indices of the local traps surrounding the selected habita grid cell
      
      theseLocalTraps <- localTrapsIndices[sID,1:localTrapsNum[sID]]
      
      
      
      ## CHECK IF DETECTIONS ARE WITHIN THE LIST OF LOCAL TRAPS
      
      if(detNums > 0){
        
        for(r in 1:detNums){
          
          if(sum(detIndices1[r] == theseLocalTraps) == 0){
            
            if(log == 0) return(0.0)
            
            else return(-Inf)
            
          }
          
        }
        
      }
      
      
      
      ## Calculate the log-probability of the vector of detections
      
      alpha <- -1.0 / (2.0 * sigma * sigma)
      
      logProb <- 0.0
      
      detIndices1 <- c(detIndices1, 0)
      
      count <- 1
      
      
      
      
      
      #if(p0==-999){# when p0 is provide through p0Traps
      
      for(r in 1:localTrapsNum[sID]){
        
        if(theseLocalTraps[r] == detIndices1[count]){
          
          d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
          
          #p <- p0Traps[theseLocalTraps[r]] * exp(alpha * d2)
          
          pZero <- ilogit(logit(p0) +
                            
                            #indTrapBeta*indTrapCov[theseLocalTraps[r]] +
                            
                            inprod(indTrapBetas, indTrapCovs[theseLocalTraps[r],]) +
                            
                            inprod(trapBetas, trapCovs[theseLocalTraps[r],]))
          
          p <- pZero * exp(alpha * d2)
          
          
          
          logProb <-  logProb + dbinom(x1[count], prob = p, size = size[theseLocalTraps[r]], log = TRUE)
          
          count <- count + 1
          
        }else{
          
          d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
          
          #p <- p0Traps[theseLocalTraps[r]] * exp(alpha * d2)
          
          pZero <- ilogit(logit(p0) +
                            
                            #indTrapBeta*indTrapCov[theseLocalTraps[r]] +
                            
                            inprod(indTrapBetas, indTrapCovs[theseLocalTraps[r],]) +
                            
                            inprod(trapBetas, trapCovs[theseLocalTraps[r],]))
          
          p <- pZero * exp(alpha * d2)        
          
          logProb <- logProb + dbinom(0, prob = p, size = size[theseLocalTraps[r]], log = TRUE)
          
          
          
        }
        
      }
      
      # #}else{# when p0 is provide through p0
      
      #   for(r in 1:localTrapsNum[sID]){
      
      #     if(theseLocalTraps[r] == detIndices1[count]){
      
      #       d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
      
      #       #p <- p0Traps[theseLocalTraps[r]] * exp(alpha * d2)
      
      #       pZero <- ilogit(logit(p0) +
      
      #                         indTrapBeta*indTrapCovs[theseLocalTraps[r]] +
      
      #                         inprod(trapBetas, trapCovs[theseLocalTraps[r],]))
      
      #       p <- pZero * exp(alpha * d2)
      
      #       logProb <-  logProb + dbinom(x1[count], prob = p, size = size[theseLocalTraps[r]], log = TRUE)
      
      #       count <- count + 1
      
      #     }else{
      
      #       d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
      
      #       #p <- p0Traps[theseLocalTraps[r]] * exp(alpha * d2)
      
      #       pZero <- ilogit(logit(p0) +
      
      #                         indTrapBeta*indTrapCovs[theseLocalTraps[r]] +
      
      #                         inprod(trapBetas, trapCovs[theseLocalTraps[r],]))
      
      #       p <- pZero * exp(alpha * d2)
      
      #       logProb <- logProb + dbinom(0, prob = p, size = size[theseLocalTraps[r]], log = TRUE)
      
      #     }
      
      #   }
      
      # }
      
      
      
      
      
      ## Return the probability of the vector of detections (or log-probability if required)
      
      if(log)return(logProb)
      
      return(exp(logProb))
      
    })
  
  
  
  
  
  #' @rdname dbinomLocal_normalBear
  
  #' @export
  
  rbinomLocal_normalBear <- nimbleFunction(
    
    run = function( n = double(0, default = 1),
                    
                    # detNums = double(0, default = -999),
                    
                    # detIndices = double(1),
                    
                    size = double(1),
                    
                    p0 = double(0, default = -999),
                    
                    # p0Traps = double(1),
                    
                    sigma = double(0),
                    
                    s = double(1),
                    
                    trapCoords = double(2),
                    
                    localTrapsIndices = double(2),
                    
                    localTrapsNum = double(1),
                    
                    #resizeFactor = double(0, default = 1),
                    
                    habitatGrid = double(2),
                    
                    indicator = double(0),
                    
                    lengthYCombined = double(0, default = 0),
                    
                    trapCovs = double(2),
                    
                    trapBetas = double(1),
                    
                    indTrapCovs = double(2),
                    
                    indTrapBetas = double(1)
                    
    ) {
      
      ## Specify return type
      
      returnType(double(1))
      
      #  if(detNums >= 0) stop("Random generation for the rbinomLocal_normal distribution is not currently supported without combining all individual detections information in one vector. See 'getSparseY()'")
      
      
      
      #========================================================
      
      # RETURN TYPE DECLARATION
      
      if(n!=1){print("rbinomLocal_normal only allows n = 1; using n = 1")}
      
      alpha <- -1.0 / (2.0 * sigma * sigma)
      
      nMAxDetections <- (lengthYCombined-1)/2
      
      ## SHORTCUT IF INDIVIDUAL IS NOT AVAILABLE FOR DETECTION
      
      if(indicator == 0){return(rep(0.0, lengthYCombined))}
      
      
      
      ## RETRIEVE THE ID OF THE HABITAT WINDOW THE CURRENT sxy FALLS IN FROM THE HABITAT_ID MATRIX
      
      resizeFactor <- 1
      
      sID <- habitatGrid[trunc(s[2]/resizeFactor)+1, trunc(s[1]/resizeFactor)+1]
      
      
      
      ## RETRIEVE THE IDs OF THE RELEVANT DETECTORS
      
      theseLocalTraps <- localTrapsIndices[sID, 1:localTrapsNum[sID]]
      
      
      
      ## INITIALIZE THE OUTPUT VECTOR OF DETECTIONS
      
      detectOut <- rep(0, localTrapsNum[sID])
      
      ys <- rep(-1, nMAxDetections)
      
      dets <- rep(-1, nMAxDetections)
      
      count <- 1
      
      
      
      ## SAMPLE THE DETECTION HISTORY (FOR RELEVANT DETECTORS ONLY)
      
      #if(p0==-999){## when p0 is provided through p0Traps
      
      for(r in 1:localTrapsNum[sID]){
        
        d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
        
        
        
        pZero <- ilogit(logit(p0) +
                          
                          #indTrapBeta*indTrapCov[theseLocalTraps[r]] +
                          
                          inprod(indTrapBetas, indTrapCovs[theseLocalTraps[r],]) +
                          
                          inprod(trapBetas, trapCovs[theseLocalTraps[r],]))
        
        
        p <- pZero * exp(alpha * d2)
        
        
        
        # Draw the observation at detector j from a binomial distribution with probability p
        
        detectOut[r] <- rbinom(1, size[theseLocalTraps[r]], p)
        
        if(detectOut[r] >0){
          
          if(nMAxDetections<count){stop("Simulated individual detections occur at more traps than what can be stored within x.\n

                                           You may need to augment the size of the x object with the argument 'nMaxTraps' from the getSparseY() function")}
          
          ys[count] <- detectOut[r]
          
          dets[count] <- theseLocalTraps[r]
          
          count <- count + 1
          
        }#if
        
      }#r
      
      
      
      count <- count - 1
      
      
      
      
      
      out <- rep(-1, lengthYCombined)
      
      
      
      out[1] <- count
      
      if(count >= 1){
        
        out[2:(count+1)] <- ys[1:count]
        
        out[(nMAxDetections+2):(nMAxDetections+count+1)] <- dets[1:count]
        
      }
      
      ## OUTPUT
      
      out <- rep(-1, lengthYCombined)
      
      
      return(out)
      
      
      #return()
      
    })
  
  
  registerDistributions(
    
    list(
      
      dbinomLocal_normalBear = list(
        
        BUGSdist ='dbinomLocal_normalBear(size, p0       , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, habitatGrid, indicator, lengthYCombined, trapCovs, trapBetas, indTrapCovs, indTrapBetas)',
        
        # Rdist = c('dbinomLocal_normalBear(detNums       , detIndices    , size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined, trapCovs, trapBetas, indTrapCovs, indTrapBetas)',
        
        #           'dbinomLocal_normalBear(detNums       , detIndices    , size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined = 0, trapCovs, trapBetas, indTrapCovs, indTrapBetas)',
        
        #           'dbinomLocal_normalBear(detNums = -999, detIndices    , size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined, trapCovs, trapBetas, indTrapCovs, indTrapBetas)',
        
        #           'dbinomLocal_normalBear(detNums = -999, detIndices = s, size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined, trapCovs, trapBetas, indTrapCovs, indTrapBetas)',
        
        #           'dbinomLocal_normalBear(detNums = -999, detIndices = s, size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined, trapCovs, trapBetas, indTrapCovs, indTrapBetas)',
        
        #           'dbinomLocal_normalBear(detNums = -999, detIndices    , size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined = 0, trapCovs, trapBetas, indTrapCovs, indTrapBetas)',
        
        #           'dbinomLocal_normalBear(detNums = -999, detIndices = s, size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0, trapCovs, trapBetas, indTrapCovs, indTrapBetas)',
        
        #           'dbinomLocal_normalBear(detNums       , detIndices    , size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined, trapCovs, trapBetas, indTrapCovs, indTrapBetas)',
        
        #           'dbinomLocal_normalBear(detNums = -999, detIndices    , size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined, trapCovs, trapBetas, indTrapCovs, indTrapBetas)',
        
        #           'dbinomLocal_normalBear(detNums = -999, detIndices    , size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0, trapCovs, trapBetas, indTrapCovs, indTrapBetas)',
        
        #           'dbinomLocal_normalBear(detNums       , detIndices    , size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0, trapCovs, trapBetas, indTrapCovs, indTrapBetas)',
        
        #          
        
        #           'dbinomLocal_normalBear(detNums       , detIndices    , size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined = 0, trapCovs, trapBetas, indTrapCovs, indTrapBetas)',
        
        #           'dbinomLocal_normalBear(detNums       , detIndices    , size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined, trapCovs, trapBetas, indTrapCovs, indTrapBetas)',
        
        #           'dbinomLocal_normalBear(detNums = -999, detIndices    , size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined, trapCovs, trapBetas, indTrapCovs, indTrapBetas)',
        
        #           'dbinomLocal_normalBear(detNums = -999, detIndices = s, size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined, trapCovs, trapBetas, indTrapCovs, indTrapBetas)',
        
        #           'dbinomLocal_normalBear(detNums = -999, detIndices = s, size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined, trapCovs, trapBetas, indTrapCovs, indTrapBetas)',
        
        #           'dbinomLocal_normalBear(detNums = -999, detIndices    , size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined = 0, trapCovs, trapBetas, indTrapCovs, indTrapBetas)',
        
        #           'dbinomLocal_normalBear(detNums = -999, detIndices = s, size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0, trapCovs, trapBetas, indTrapCovs, indTrapBetas)',
        
        #           'dbinomLocal_normalBear(detNums       , detIndices    , size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined, trapCovs, trapBetas, indTrapCovs, indTrapBetas)',
        
        #           'dbinomLocal_normalBear(detNums = -999, detIndices    , size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined, trapCovs, trapBetas, indTrapCovs, indTrapBetas)',
        
        #           'dbinomLocal_normalBear(detNums = -999, detIndices    , size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0, trapCovs, trapBetas, indTrapCovs, indTrapBetas)',
        
        #           'dbinomLocal_normalBear(detNums       , detIndices    , size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0, trapCovs, trapBetas, indTrapCovs, indTrapBetas)'
        
        # ),
        
        types = c('value = double(1)', 'size = double(1)', 's = double(1)', 'trapCoords = double(2)', 'localTrapsIndices = double(2)', 'localTrapsNum = double(1)', 'habitatGrid = double(2)',"trapCovs = double(2)",
                  
                  'trapBetas = double(1)', 'indTrapCovs = double(2)', 'indTrapBetas = double(1)'),
        
        discrete = TRUE,
        
        mixedSizes = TRUE,
        
        pqAvail = FALSE
        
      )
      
    ),
    
    verbose = T)
  
  ### pass functions to global environment
  assign('dbinomLocal_normalBear', dbinomLocal_normalBear, envir = .GlobalEnv)
  assign('rbinomLocal_normalBear', rbinomLocal_normalBear, envir = .GlobalEnv)
  
  
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


