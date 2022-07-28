
## -------------------------------------------------
##        Understand nimbleSCR functions
## ------------------------------------------------- 

library(nimbleSCR)

## ---- getSparseY to prepare dataset to speed up Nimble ----

## EXAMPLE DATA

#' # I. DATA SET UP 
coordsHabitatGridCenter <- matrix(c(0.5, 3.5,
                                    1.5, 3.5,
                                    2.5, 3.5,
                                    3.5, 3.5,
                                    0.5, 2.5,
                                    1.5, 2.5,
                                    2.5, 2.5,
                                    3.5, 2.5,
                                    0.5, 1.5,
                                    1.5, 1.5,
                                    2.5, 1.5,
                                    3.5, 1.5,
                                    0.5, 0.5,
                                    1.5, 0.5,
                                    2.5, 0.5,
                                    3.5, 0.5), ncol=2,byrow = TRUE)
colnames(coordsHabitatGridCenter) <- c("x","y")
# CREATE OBSERVATION WINDOWS
trapCoords <- matrix(c(1.5, 1.5, 2.5, 1.5, 1.5, 2.5, 2.5, 2.5), nrow = 4, byrow = TRUE)
colnames(trapCoords) <- c("x","y")
# PLOT CHECK
plot(coordsHabitatGridCenter[,"y"]~coordsHabitatGridCenter[,"x"],pch=16) 
points(trapCoords[,"y"]~trapCoords[,"x"],col="red",pch=16) 

# PARAMETERS
p0 <- 0.2
sigma <- 2
indicator <- 1 
# WE CONSIDER 2 INDIVIDUALS
y <- matrix(c(0, 1, 1, 0,
              0, 1, 0, 1),ncol=4,nrow=2) ## ASP: Capture history: Number of id x number of traps
s <- matrix(c(0.5, 1,
              1.6, 2.3),ncol=2,nrow=2)

# RESCALE COORDINATES 
ScaledtrapCoords <- scaleCoordsToHabitatGrid(coordsData =  trapCoords,
                                             coordsHabitatGridCenter = coordsHabitatGridCenter)
ScaledtrapCoords<- ScaledtrapCoords$coordsDataScaled
habitatMask <- matrix(1, nrow = 4, ncol=4, byrow = TRUE)


# CREATE LOCAL OBJECTS 
TrapLocal <- getLocalObjects(habitatMask = habitatMask,
                             coords = ScaledtrapCoords,
                             dmax=2.5,
                             resizeFactor = 1,
                             plot.check = TRUE)

# GET SPARSE MATRIX 
SparseY <- getSparseY(y)

getSparseY <- function( 
                        x = y
                        ,
                        noDetections = -1
                        ,
                        nMaxTraps = NULL 
                        
){
  # IF 2D OBSERVATION MATRIX, CONVERT TO 3D ARRAY
  if(length(dim(x))==2){
    Y <- array(x, c(dim(x),1))
  }else{Y <- x}
  
  # RETRIEVE THE NUMBER OF DETECTIONS FOR EACH ID
  detNums <- apply(Y, c(1,3), function(x) length(which(x>0)))
  
  
  
  # INITIALIAZE EMPTY ARRAYS FOR THE DETECTIONS AND DETECTOR INDICES
  detIndices <- array(-1, c(dim(Y)[1],max(detNums), dim(Y)[3]))
  ySparse <- array(-1, c(dim(Y)[1],max(detNums), dim(Y)[3]))
  
  # FILL IN THE ARRAYS
  for(t in 1:dim(Y)[3]){
    for(i in 1:dim(Y)[1]){
      if(detNums[i,t] > 0){
        # GET WHERE (DETECTOR ID) DETECTIONS OCCUR
        detIndices[i,1:detNums[i,t],t] <- which(Y[i, ,t] > 0)
        # GET NUMBER OF DETECTIONS 
        ySparse[i,1:detNums[i,t],t] <- Y[i,which(Y[i, ,t] > 0),t] 
      }
    }
  }   
  
  # IF DETECTIONS SHOULD BE BINDED TO ALLOW NUMBER LARGE NUMBER OF DETECTORS AT WHICH DETECTIONS CAN OCCUR
  if(!is.null(nMaxTraps)){
    if(max(detNums) >= nMaxTraps){
      
      nMaxTraps <-max(detNums)*2
      print("Warnings! nMaxTraps was less than or equal to maxDetNums (the maximum number of spatial recaptures). We have set nMaxTraps to 2*maxDetNums")
    }
  }else{# if nMaxTraps is not provided, the default value is max(detNums)*2
    nMaxTraps <- max(detNums)*2
  }
  
  increaseSize <- nMaxTraps 
  yCombined <- array(NA, c(dim(ySparse)[1], 1 + (increaseSize)*2, dim(ySparse)[3])) 
  for(t in 1:dim(yCombined)[3]){
    yCombined[,,t] <- cbind(detNums[,t],
                            ySparse[,,t],
                            matrix(-1, nrow=dim(ySparse[,,t])[1], ncol=increaseSize-max(detNums)),
                            detIndices[,,t],
                            matrix(-1, nrow=dim(ySparse[,,t])[1], ncol=increaseSize-max(detNums))
    )
  }
  
  
  return(list( y = ySparse, ## ASP: Fills in the 1 in each detection per individual
               detIndices = detIndices, ## ASP: Number of rows = max number of detections. Which trap?
               detNums = detNums, ## ASP: Number of detections (summing 1 ch) per individual
               maxDetNums = max(detNums),
               yCombined = yCombined,
               lengthYCombined = dim(yCombined)[2]))
}

## ---- getLocalObjects to perform local evaluation and assign closest traps to activity centers ----

## Run function with example data

getLocalObjects <- function( habitatMask = habitatMask
                             ,
                             coords = ScaledtrapCoords
                             ,
                             dmax=2.5
                             ,
                             resizeFactor = 1
                             ,
                             plot.check = TRUE
){
  ## STORE THE COORDINATES OF THE ORIGINAL HABITAT CELLS
  oldCoords <- which(habitatMask == 1, arr.ind = T) - 0.5
  oldCoords <- cbind(oldCoords[,2], oldCoords[,1])
  names(oldCoords) <- c("x", "y")
  
  ## GET ORIGIN FOR THE NEW (RESIZED) HABITAT MATRIX 
  origin <- (resizeFactor/2) 
  
  ## GET DIMENSIONS FOR THE NEW HABITAT MATRIX
  xMax <- ceiling(dim(habitatMask)[2]/resizeFactor) * resizeFactor
  yMax <- ceiling(dim(habitatMask)[1]/resizeFactor) * resizeFactor
  
  ## GET COORDINATES FOR THE NEW HABITAT CELLS
  xCoords <- seq(origin, xMax, by = resizeFactor)
  yCoords <- seq(origin, yMax, by = resizeFactor)
  habitatCoords <- expand.grid(list(x = xCoords, y = yCoords))
  habitatCoords <- habitatCoords[order(habitatCoords[ ,1], habitatCoords[ ,2]), ]
  
  ## GET UPPER AND LOWER COORDINATES FOR THE NEW HABITAT CELLS
  habUpCoords <- habitatCoords + resizeFactor/2
  habLoCoords <- habitatCoords - resizeFactor/2
  
  ## CHECK WHICH "NEW CELLS" CONTAIN AT LEAST ONE "OLD CELL"
  isIn <- unlist(lapply(1:dim(habUpCoords)[1], function(c){
    sum((habLoCoords[c,1] <= oldCoords[ ,1]) *
          (habLoCoords[c,2] <= oldCoords[ ,2]) *
          (habUpCoords[c, 1] > oldCoords[ ,1]) *
          (habUpCoords[c,2] > oldCoords[ ,2])) > 0
  })) 
  
  ## REMOVE NEW HABITAT CELL COORDINATES THAT ARE NOT HABITAT  
  habitatCoords <- habitatCoords[isIn, ]
  
  ## CREATE AN EMPTY MATRIX OF NEW HABITAT CELL IDs
  habitatID <- matrix(0, nrow = length(yCoords), ncol = length(xCoords))
  for(c in 1:dim(habitatCoords)[1]){
    habitatID[trunc(habitatCoords[c,2]/resizeFactor)+1, trunc(habitatCoords[c,1]/resizeFactor)+1] <- c
  }
  
  ## DETERMINE WHICH POINTS ARE WITHIN dmax OF THE CENTER OF EACH NEW HABITAT CELL
  
  x = habitatCoords[1,]
  
  localIndices  <- apply(habitatCoords, 1, function(x){ ## ASP: Apply to each row (coords of habitat each cell)
    D <- sqrt((x[1] - coords[,1])^2 + (x[2] - coords[,2])^2) # Distance from the cell to each trap
    which(D < dmax) ## ASP: These traps are below the dmax from the habitat cell
  })
  
  ## ASP: Gives the traps ID that are below the dmax from each cell
 
  ## MAKE SURE IT ALWAYS RETURN A LIST
  if(class(localIndices) == "matrix"){
    localIndices <- lapply(1:dim(localIndices)[2], function(x) localIndices[ ,x])
  }
  
  ## GET NUMBER OF DETECTORS WITHIN dmax OF EACH NEW HABITAT CELL
  numLocalIndices <- unlist(lapply(localIndices, function(x) length(x)))
    ## ASP: Number of traps that are near the cell
  maxLocalIndices <- max(numLocalIndices)
  
  ## FOR ALL HABITAT GRIDS, THE LOCAL EVALUATION SHOULD BE LARGE ENOUGH TO OVERLAP WITH >0 TRAP
  if(any(numLocalIndices %in% 0 )){
    stop("dmax value too small. All habitat grid cells should have at least one local object within a radius of dmax.")
  }
  ## STORE LOCAL  INDICES IN A MATRIX 
  Index <- matrix(0, nrow = length(localIndices), ncol = maxLocalIndices)
  for(j in 1:length(localIndices)){
    if(length(localIndices[[j]])!=0){
      Index[j, 1:numLocalIndices[j]] <- localIndices[[j]]
    }
  }
  
  
  ## PLOT CHECK 
  if(plot.check){
    SXY <- as.numeric(habitatCoords[sample(1:dim(habitatCoords)[1], size = 1), ])
    sxyID <- habitatID[trunc(SXY[2]/resizeFactor)+1, trunc(SXY[1]/resizeFactor)+1]
    index <- Index[sxyID, 1:numLocalIndices[sxyID]]
    
    plot(habitatCoords[ ,2] ~ habitatCoords[ ,1], pch = 16, cex = 0.1)
    points(habitatCoords[sxyID,2] ~ habitatCoords[sxyID,1], pch = 16, cex = 0.4, col = "orange")
    points(coords[ ,2] ~ coords[ ,1], pch = 16, cex = 0.2, col = "red")
    points(coords[index,2] ~ coords[index,1], pch = 16, cex = 0.4, col = "blue")
    points(SXY[2] ~ SXY[1], bg = "red", pch = 21, cex = 1.2)
  }
  
  
  
  ## OUTPUT LIST
  output <- list( habitatGrid = habitatID,
                  localIndices = Index,
                  numLocalIndices = numLocalIndices,
                  numLocalIndicesMax = maxLocalIndices,
                  resizeFactor = resizeFactor)
  
  return(output)
}

## ---- Bernoulli point process: distribution of activity centers ----

## RANDOM GENERATION (rbernppAC): Generates 1 AC

# logIntensities (model) from simulated beta (beta.dens = 1)
logIntensities <- log(exp(beta.dens * habDens[1:numHabWindows])) 
logSumHabIntensity <- log(sum(exp(beta.dens * habDens[1:numHabWindows])))
# I guess this needs to be converted to probabilities like in the model and some to one??
# So I use p.cell:
p.cell <- exp(beta.d*X.d)/sum(exp(beta.d*X.d))

# I need to introduce it as the log of pcell, because then it exp within the function to obtain probabilities
logpcell <- log(p.cell)

n = 1 # Only able to generate one activity center
lowerCoords     = lowerHabCoords
upperCoords     = upperHabCoords
logIntensities  = logpcell # Needs beta.dens from model, introduced with simulated
logSumIntensity = logSumHabIntensity # Needs beta.dens from model
habitatGrid     = habitatGrid[1:numGridRows,1:numGridCols]
numGridRows     = numGridRows
numGridCols     = numGridCols

rbernppAC <- nimbleFunction(
  run = function(
    n               = integer(0),
    lowerCoords     = double(2),
    upperCoords     = double(2),
    logIntensities  = double(1),
    logSumIntensity = double(0),
    habitatGrid     = double(2),
    numGridRows     = integer(0),
    numGridCols     = integer(0)
  ) {
    if(n <= 0) stop("The number of requested samples must be above zero")
    else if(n > 1) print("rbernppAC only allows n = 1; using n = 1")
    ## Simulate window index
    windowInd <- rcat(1, exp(logIntensities)) # This choses a cell from 144 categories?
    numDims <- 2 
    ## A uniform distribution is used within the target window (so that is not right in center??)
    outCoordinates <- lowerCoords[windowInd,] + 
      runif(numDims, 0.0, 1.0) * (upperCoords[windowInd,] - lowerCoords[windowInd,])
    return(outCoordinates)
    returnType(double(1)) 
  }
)

## DENSITY FUNCTION (dbernppAC) of the Bernoulli point process for the distribution of AC

x               = outCoordinates # I guess????
lowerCoords     = lowerHabCoords
upperCoords     = upperHabCoords
logIntensities  = logIntensities
logSumIntensity = logSumHabIntensity
habitatGrid     = habitatGrid[1:numGridRows,1:numGridCols]
numGridRows     = numGridRows
numGridCols     = numGridCols
#log             = integer(0, default = 0)

dbernppAC <- nimbleFunction(
  run = function(
    x               = double(1),
    lowerCoords     = double(2),
    upperCoords     = double(2),
    logIntensities  = double(1),
    logSumIntensity = double(0),
    habitatGrid     = double(2),
    numGridRows     = integer(0),
    numGridCols     = integer(0),
    log             = integer(0, default = 0)
  ) {
    ## Check if the point falls within the habitat: 
    ## Getting numGridRows and numGridCols using the following code takes some time
    ## and may cause inefficiency if the function is called repeatedly in a loop.
    ## numGridRows <- dim(habitatGrid)[1]
    ## numGridCols <- dim(habitatGrid)[2]
    ## In addition, when the true habitat grid has one row/column we need to inflate it for use in NIMBLE model code. 
    ## In this case, we have problems using the code above.
    ## Note that we need to rescale the habitat gird to ensure x and y coordinates start from 0
    ## and each window is of size 1x1. So the following code works correctly. 
    if(min(x) < 0 | x[2] >= numGridRows | x[1] >= numGridCols) {
      if(log) return(-Inf) 
      else return(0.0)
    }
    ## Find which window point x falls within
    windowInd <- habitatGrid[trunc(x[2])+1, trunc(x[1])+1] # Good! Its cell 15
    ## windowInd == 0 means this window is not defined as habitat
    if(windowInd == 0) {
      if(log) return(-Inf)
      else return(0.0)
    }
    ## Log probability density 
    logProb <- logIntensities[windowInd] - logSumIntensity
    
    if(log) return(logProb)
    else return(exp(logProb))
    returnType(double(0))
  }

)

## ---- Binomial SCR detection process ----

## model and simulate binomial observations (x) of a single individual over a set of detectors defined by
## their coordinates (trapCoords). The distribution assumes that an individual’s detection probability
## at any detector follows a half-normal function of the distance between the individual’s activity
## center (s) and the detector location

## ASP: From AC location, calculate the probability of detection 
## -> Which depends on the distance between the AC and X
## -> Follows a binomial distribution

## EXAMPLE DATA

#' # I. DATA SET UP 
coordsHabitatGridCenter <- matrix(c(0.5, 3.5,
                                   1.5, 3.5,
                                   2.5, 3.5,
                                   3.5, 3.5,
                                   0.5, 2.5,
                                   1.5, 2.5,
                                   2.5, 2.5,
                                   3.5, 2.5,
                                   0.5, 1.5,
                                   1.5, 1.5,
                                   2.5, 1.5,
                                   3.5, 1.5,
                                   0.5, 0.5,
                                   1.5, 0.5,
                                   2.5, 0.5,
                                   3.5, 0.5), ncol=2,byrow = TRUE)
colnames(coordsHabitatGridCenter) <- c("x","y")
# CREATE OBSERVATION WINDOWS
trapCoords <- matrix(c(1.5, 1.5, 2.5, 1.5, 1.5, 2.5, 2.5, 2.5), nrow = 4, byrow = TRUE)
colnames(trapCoords) <- c("x","y")
# PLOT CHECK
plot(coordsHabitatGridCenter[,"y"]~coordsHabitatGridCenter[,"x"],pch=16) 
points(trapCoords[,"y"]~trapCoords[,"x"],col="red",pch=16) 

# PARAMETERS
p0 <- 0.2
sigma <- 2
indicator <- 1 
# WE CONSIDER 2 INDIVIDUALS
y <- matrix(c(0, 1, 1, 0,
             0, 1, 0, 1),ncol=4,nrow=2)
s <- matrix(c(0.5, 1,
             1.6, 2.3),ncol=2,nrow=2)

# RESCALE COORDINATES 
ScaledtrapCoords <- scaleCoordsToHabitatGrid(coordsData =  trapCoords,
                                            coordsHabitatGridCenter = coordsHabitatGridCenter)
ScaledtrapCoords<- ScaledtrapCoords$coordsDataScaled
habitatMask <- matrix(1, nrow = 4, ncol=4, byrow = TRUE)


# CREATE LOCAL OBJECTS 
TrapLocal <- getLocalObjects(habitatMask = habitatMask,
                                  coords = ScaledtrapCoords,
                                  dmax=2.5,
                                  resizeFactor = 1,
                                  plot.check = TRUE)

# GET SPARSE MATRIX 
SparseY <- getSparseY(y)

#########################################################################
## RANDOM GENERATION

# Take the first individual 
i = 1
# This generates the detection data for an individual: From the location of the AC
# locates the closer traps, calculates distance and the p of detection
# With this p generates the observation process through the binomial model

rbinomLocal_normal <- nimbleFunction(
  run = function( 
                  #n = double(0, default = 1),
                  #detNums = double(0, default = -999),
                  #detIndices = double(1),
                  #size = double(1),
                  #p0 = double(0, default = -999),
                  #p0Traps = double(1),
                  #sigma = double(0),
                  #s = double(1),
                  #trapCoords = double(2),
                  #localTrapsIndices = double(2),
                  #localTrapsNum = double(1),
                  #resizeFactor = double(0, default = 1),
                  #habitatGrid = double(2),
                  #indicator = double(0),
                  #lengthYCombined = double(0, default = 0)
                  
                  n=1
                  ,
                  size=rep(1,4)
                  ,
                  p0 = p0
                  ,
                  sigma= sigma
                  , 
                  s=s[i,1:2]   ## ASP: This comes from within model (point process)
                  ,
                  trapCoords=ScaledtrapCoords
                  ,
                  localTrapsIndices=TrapLocal$localIndices
                  ,
                  localTrapsNum=TrapLocal$numLocalIndices
                  ,
                  resizeFactor=TrapLocal$resizeFactor
                  ,
                  habitatGrid=TrapLocal$habitatGrid
                  ,
                  indicator=indicator
                  ,
                  lengthYCombined = SparseY$lengthYCombined
  ) {
    ## Specify return type
    returnType(double(1))
    if(detNums >= 0) stop("Random generation for the rbinomLocal_normal distribution is not currently supported without combining all individual detections information in one vector. See 'getSparseY()'")
    
    #========================================================
    # RETURN TYPE DECLARATION
    if(n!=1){print("rbinomLocal_normal only allows n = 1; using n = 1")}
    # returnType(double(3))
    # len <- 2*MAX + 1
    ## GET NECESSARY INFO
    alpha <- -1.0 / (2.0 * sigma * sigma)
    # n.detectors <- dim(detector.xy)[1]
    # nMAxDetections <- length(detIndices)
    nMAxDetections <- (lengthYCombined-1)/2
    ## SHORTCUT IF INDIVIDUAL IS NOT AVAILABLE FOR DETECTION
    #if(indicator == 0){return(rep(0.0, 2*nMAxDetections + 1))}
    if(indicator == 0){return(rep(0.0, lengthYCombined))}
    
    ## RETRIEVE THE ID OF THE HABITAT WINDOW THE CURRENT sxy FALLS IN FROM THE HABITAT_ID MATRIX
    sID <- habitatGrid[trunc(s[2]/resizeFactor)+1, trunc(s[1]/resizeFactor)+1]
    
    ## RETRIEVE THE IDs OF THE RELEVANT DETECTORS
    theseLocalTraps <- localTrapsIndices[sID, 1:localTrapsNum[sID]]
    
    ## INITIALIZE THE OUTPUT VECTOR OF DETECTIONS
    detectOut <- rep(0, localTrapsNum[sID])
    ys <- rep(-1, nMAxDetections)
    dets <- rep(-1, nMAxDetections)
    count <- 1
    
    ## SAMPLE THE DETECTION HISTORY (FOR RELEVANT DETECTORS ONLY)
    if(p0==-999){## when p0 is provided through p0Traps
        ## ASP: If there are covariates that are trap specific (different p0 per trap)
      
      for(r in 1:localTrapsNum[sID]){
        d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
        p <- p0Traps[theseLocalTraps[r]] * exp(alpha * d2)
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
    }else{## when p0 is provided through p0
      ## ASP: No covariates, the case of the example
      
      for(r in 1:localTrapsNum[sID]){   
        ## ASP: ex, 4 traps available at 2.5 dmax from cell 2
        d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
          ## ASP: Distance from the activity center to the available trap r 
        p <- p0 * exp(alpha * d2)
          ## ASP: Probability of detection at that trap
        
        # Draw the observation at detector j from a binomial distribution with probability p
        detectOut[r] <- rbinom(1, size[theseLocalTraps[r]], p) 
          ## ASP: size[theseLocalTraps[r]] are the number of trials, for example the occasions
        
        if(detectOut[r] >0){
          if(nMAxDetections<count){stop("Simulated individual detections occur at more traps than what can be stored within x.\n
                                          You may need to augment the size of the x object with the argument 'nMaxTraps' from the getSparseY() function")}
          ys[count] <- detectOut[r]
          dets[count] <- theseLocalTraps[r]
          count <- count + 1
        }#if
      }#r 
      
    }
    count <- count - 1
    
    
    # out <- rep(-1, 2*nMAxDetections + 1)
    out <- rep(-1, lengthYCombined)
    
    out[1] <- count
    if(count >= 1){
      out[2:(count+1)] <- ys[1:count]
      out[(nMAxDetections+2):(nMAxDetections+count+1)] <- dets[1:count]
    }
    ## OUTPUT
    return(out)
    
    ## ASP: output of the function for each individual (activity center):
      # out[1]: Total number of spatial recaptures
      # out[2:...x(maximun number of detections, lets say 5)]: How many detections in each trap
      # out[5:10]: ids of the traps
    # So in one vector per individual you get all the capture information. There have to be always -1
    # because it rarely fills up (otherwise there wouldn't be space to store, it's a Nimble thing)
    
    # If is not alive, it only gives -1. Indicator?
    
  })


########################################################################################

dbinomLocal_normal <- nimbleFunction(
  run = function( 
                  #x = double(1),
                  #detNums = double(0, default = -999),
                  #detIndices = double(1),
                  #size = double(1),
                  #p0 = double(0, default = -999),
                  #p0Traps = double(1),
                  #sigma = double(0),
                  #s = double(1),
                  #trapCoords = double(2),
                  #localTrapsIndices = double(2),
                  #localTrapsNum = double(1),
                  #resizeFactor = double(0, default = 1),
                  #habitatGrid = double(2),
                  #indicator = double(0),
                  #lengthYCombined = double(0, default = 0),
                  #log = integer(0, default = 0)
                  
                  x=SparseY$y[i,,1]
                  ,
                  detNums=SparseY$detNums[i] ## ASP: Number of total detections of the individual
                  ,
                  detIndices=SparseY$detIndices[i,,1] ## ASP: In which traps are the detections
                  ,
                  size=rep(1,4) ## ASP: Four traps
                  ,
                  p0 = p0
                  ,
                  sigma= sigma
                  , 
                  s=s[i,1:2]
                  ,
                  trapCoords=ScaledtrapCoords
                  ,
                  localTrapsIndices=TrapLocal$localIndices
                  ,
                  localTrapsNum=TrapLocal$numLocalIndices
                  ,
                  resizeFactor=TrapLocal$resizeFactor
                  ,
                  habitatGrid=TrapLocal$habitatGrid
                  ,
                  indicator=indicator
  ) {
    ## Specify return type
    returnType(double(0))
    
    ## Deal with cases where detection info is combined in one vector 
    if(detNums==-999){
      detNums <- x[1]
      nMaxDetectors <- (lengthYCombined-1)/2
      detIndices1 <- x[(nMaxDetectors+2):lengthYCombined]
      x1 <- x[2:(nMaxDetectors+1)]
    }else{
      x1 <- x
      detIndices1 <- detIndices
    }
    
    ## Shortcut if the current individual is not available for detection
      ## ASP: Because z = 0??
    
    if(indicator == 0){
      if(detNums == 0){
        if(log == 0) return(1.0)
        else return(0.0)
      } else {
        if(log == 0) return(0.0)
        else return(-Inf)
      }
    }
    
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
    }    ## ASP: Detections are supossed to be always in the list of local traps right?
        # to detect errors?
      
    ## Calculate the log-probability of the vector of detections
    alpha <- -1.0 / (2.0 * sigma * sigma)
    logProb <- 0.0 
    detIndices1 <- c(detIndices1, 0)
    count <- 1 
    
    
    if(p0==-999){# when p0 is provided through p0Traps
      for(r in 1:localTrapsNum[sID]){
        if(theseLocalTraps[r] == detIndices1[count]){ 
          d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
          p <- p0Traps[theseLocalTraps[r]] * exp(alpha * d2)
          logProb <-  logProb + dbinom(x1[count], prob = p, size = size[theseLocalTraps[r]], log = TRUE)
          count <- count + 1
        }else{
          d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
          p <- p0Traps[theseLocalTraps[r]] * exp(alpha * d2)
          logProb <- logProb + dbinom(0, prob = p, size = size[theseLocalTraps[r]], log = TRUE)
          
        }
      }
    }else{# when p0 is provide through p0
      for(r in 1:localTrapsNum[sID]){ # For each trap
        if(theseLocalTraps[r] == detIndices1[count]){ 
            ## ASP: For the first one?
          d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
            ## ASP: Distance from the trap to the activity center
          p <- p0 * exp(alpha * d2)
            ## ASP: Probability of detection at that trap (according to how is from the activity center)
          logProb <-  logProb + dbinom(x1[count], prob = p, size = size[theseLocalTraps[r]], log = TRUE)
            ## ASP: x1 is the detections vector (nº detections in a vector)
            ## dbinom (x1 = 1: what is the probability of getting 1 detection/how likely 
            ##        prob = when the prob. of success is p,
            ##        with a number of trials = size (occasions))
            
          count <- count + 1
          ## ASP: This is to index the detections so that it fits with the localtraps
          ## the first fits but then it doesn't fit anymore?
        }else{
          d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
            ## ASP: Distance from the activity center to the trap
          p <- p0 * exp(alpha * d2)
          logProb <- logProb + dbinom(0, prob = p, size = size[theseLocalTraps[r]], log = TRUE)
          ## ASP: dbinom (x1 = 1: what is the probability of getting 0 detection/how likely..why now a 0?
        }
      }
    }
    
    ## ASP: This function sums one logprob per trap.
    ## ??DOUBT When you do dbinom con 0 or 1 (when theseLocalTraps[r] == detIndices1[count]).
    ## BUT I can live without this
    
    ## Return the probability of the vector of detections (or log-probability if required)
    if(log)return(logProb)
    return(exp(logProb))
  })

