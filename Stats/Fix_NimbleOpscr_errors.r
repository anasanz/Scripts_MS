
## -------------------------------------------------
##               Fix errors nimble opscr
## ------------------------------------------------- 


#(1) set up model

model <- nimbleModel(SCRhab.Open.diftraps.age.effortTrapBhCov.sigsexage, constants = nimConstants, 
                     data=nimData, inits=inits(), check = FALSE)
##ignore error message, only due to missing initial values at this stage

model$initializeInfo()

model$calculate()

# IF Likelihood is -Inf
# Normally the error is in the following:

which(is.infinite(model$logProb_sex),arr.ind = T)
which(is.infinite(model$logProb_sxy),arr.ind = T)
which(is.infinite(model$logProb_u),arr.ind = T)
which(is.infinite(model$logProb_w),arr.ind = T)
which(is.infinite(model$logProb_y),arr.ind = T)

# Sometimes when the error is in y, the problem is in the local evaluation (increase dmax)
# In this case it is:
#      dim1 dim2 dim3 dim4
#[1,]   65    1    5    7
# It doesnt get solved when increasing dmax, so you can check what it is by:
# Runing the Nimble function with the x = what is in the left side of the tilde
# (y[i,1:lengthYCombined[t],k,t]~dbinomLocal_normalBear)

i=65
k=5
t=7

dbinomLocal_normalBear(
  x=model$y[i,1:nimConstants$lengthYCombined[t],k,t]
  ,
  size = model$ones[1:nimConstants$J[t]]
  , ##NOW: always 1, because we model each occasion separately
  p0 = model$p0[model$age.cat[i,t]+1]
  ,
  sigma = model$sigma[model$sex.age[i,t]+1]
  , #model parameter
  s = model$sxy[i,1:2,t]
  , #model parameter
  trapCoords = nimData$X.sc[1:nimConstants$J[t],1:2,t]
  , #trap coordinates (data); ASP: Year specific trap array
  localTrapsIndices = nimData$localTrapsIndex[1:nimConstants$numHabWindows,1:nimConstants$MaxLocalTraps[t],t]
  , #from getLocalTraps()
  localTrapsNum = nimData$localTrapsNum[1:nimConstants$numHabWindows,t]
  , #from getLocalTraps()
  lengthYCombined = nimConstants$lengthYCombined[t]
  ,
  habitatGrid = nimData$habitatGridDet[1:numGridRows,1:numGridCols]
  ,#from getLocalTraps()
  trapCovs = nimData$effort[1:nimConstants$J[t],k,t,1:nimConstants$nTrapCovs]
  ,
  trapBetas = model$trapBetas[1:nimConstants$nTrapCovs]
  
  ,
  indTrapCov = model$prevcap[i,1:nimConstants$J[t],k,t]
  ,
  indTrapBeta = model$b.bh
  ,
  indicator = model$z[i,t]
) #

# In this example, we see by running: model$y[i,1:nimConstants$lengthYCombined[t],k,t] you see that it has been detected
# But that year, it was not alive according to: nimData$u[i,]

# You can fix it, or set up the problematic thing (ex. the problematic year to 1) and see if it calculates it
model$u[i,t] <- 1
model$calculate()

