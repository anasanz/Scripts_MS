

## -------------------------------------------------
##      Predict from simulation script RS (run1)
##              SPATIAL PROJECTION
## ------------------------------------------------- 

#rm(list = ls())

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Predict")
samp2<-readRDS('TestOutputAGE.rds')
sampmat<-do.call(rbind, samp2)
dim(sampmat)

##create starting z's and sxy's based on model output and new M and time frame
t.new<-5 #sim 5 addl years
M.new<-400 #new augmentation limit

##use estimates of z, age etc from last year of estimation model only
z.start<-matrix(NA, M.new, 1+t.new)
age.start<-matrix(NA, M.new, 1+t.new)#

##fill z.in with estimates of z form fitted model; set addl individuals in years
##covered by model to 0 - known NOT to be alive
#first, make matrix of z estimates
z.which1<-pmatch(paste('z[', 1:Maug, ', ', Tt,']', sep=''), colnames(sampmat)) # ASP: index columns z for LAST YEAR (sampmat matrix)
z.which<-grep('z\\[', colnames(sampmat)) # ASP: index columns all z (sampmat matrix)
##ordered: all individuals for yr 1, then all for yr2, etc
z.est<-matrix(sampmat[1,z.which], Maug, Tt) # ASP: All z
z.start[1:Maug, 1]<-sampmat[1,z.which1] # ASP: Starting z for projection: z of year 5 
z.start[(Maug+1):M.new, 1]<-0 # ASP: Augmented as not yet recruited

##fill age.start with estimates of age form fitted model; set addl individuals in years
##covered by model to 0 - known NOT to be alive = age 0
#first, make matrix of age estimates
age.which<-pmatch(paste('age[', 1:Maug, ', ', Tt,']', sep=''), colnames(sampmat)) # ASP: index columns AGE for LAST YEAR
age.est<-sampmat[1,age.which] # ASP: First iteration, Ages last year
##set everyone not alive at all to all-0, ie, can be recruited in the NEW recruitment 
##model (they were never part of the superpopulation -> ASP: ASP: because of how the model is formulated z=u*w
z.nosuper<-which(apply(z.est,1,sum)==0)
age.est[z.nosuper]<-0
age.est[age.est>5]<-5 #set to max age category
age.start[1:Maug, 1]<-age.est # ASP: Starting AGE for projection: AGE of year 5 
age.start[(Maug+1):M.new, 1]<-0 # ASP: Augmented as not yet recruited


sxy.start<-array(NA, c(M.new, 2, 1+t.new))
##make array of activity centers from model output
s.which<-grep('sxy', colnames(sampmat)) # ASP: index columns all sxy (sampmat matrix)
##all individuals X Yr1, then all inds Y Yr1; all inds X yr 2, then y yr 2 etc ## ASP: This is the order the sxy vector is placed (s.which)

indx<-(Tt*Maug*2-(Maug*2-1)):(Tt*Maug*2) # ASP: Index the LAST year of sxy. *2 because there are the double of columns: x and y
sxy.start[1:Maug,,1]<-matrix(sampmat[ 1, s.which[indx] ], Maug, 2) # ASP: Starting SXY for projection: SXY of year 5 

##calculate per capita recruitment
#R.which<-grep('R', colnames(sampmat))
##didn't monitor R... calculate from z's
R<-matrix(NA, nrow(sampmat), Tt-1)
for (t in 2:Tt){
  zwt<-paste('z[', 1:Maug,', ', t, ']', sep='' ) # ASP: names z all individuals on year t
  zwtm<-paste('z[', 1:Maug,', ', t-1, ']', sep='' ) # ASP: names z all individuals on year t-1
  
  for (nn in 1:nrow(sampmat)){ # ASP: Each iteration
  R[nn,t-1] <- sum(ifelse((sampmat[nn,zwt]-sampmat[nn,zwtm]) == 1, 1, 0)) # ASP: Identified recruited ind (1 at t - 0 at t - 1 = recruited 1); sum all to get total numer of recruited from year t-1 to t
  }
}

##ASP: PCR -> Average number of individuals recruited per Adult from t-1 to t. 
N.which<-grep('N', colnames(sampmat))[1:(Tt-1)] # ASP: Index columns N the four first years (to calculate pcr)
pcr<-mean(R[1,]/sampmat[1,N.which[1:4]]) # ASP: Mean per capita recruitment (for ONE iteration). Example?

pcrmat<-matrix(NA, nrow(sampmat), 4)
for (ite in 1:nrow(sampmat)){
  pcrmat[ite,]<-R[ite,]/sampmat[ite,N.which[1:4]]
}

phi.ad<-sampmat[1,'phi.ad']
phi.cub<-sampmat[1,'phi.cub']
phi.sub<-sampmat[1,'phi.sub']
beta.dens<-sampmat[1,'beta.dens']
sigD<-sampmat[1,'sigD']

##calculate psi and unconditional age structure at t=1, which don't really matter 
## as they  don't apply to the projection years
##but need to be in the model for proper structure
##won't be changed in simulations
psi<-mean(sampmat[,'N[5]'])/M.new # ASP: This is the data augmentation parameter (N last year/AUG), not really understand but it doesn't matter?
##unconditional age structure at t=1 (old t=5, from single iteration)
pi.uncond<-table(age.start[,1])/M.new # ASP: Reminder -> it is unconditional because it adds the prob. of not being recruited for augmented individuals
  
##set up constants
nimConstants <- list(
  M=(M.new),numHabWindows=numHabWindows, 
  numGridRows=numGridRows, numGridCols=numGridCols, Nyr=(1+t.new),
  max.age=max.age
)

##set up data
nimData<-list(habDens=X.d, 
              lowerHabCoords=lowerHabCoords, upperHabCoords=upperHabCoords,
              habitatGrid=habitatGrid,
              z=z.start,
              sxy=sxy.start,
              age=age.start+1,
              age.cat=age.start,
              u=z.start,
              pcr=pcr, ##transformed per capita recruitment
              phi.ad = phi.ad,
              phi.cub = phi.cub,
              phi.sub = phi.sub,
              beta.dens = beta.dens,
              sigD=sigD,
              psi=psi,
              pi.uncond=pi.uncond
)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Predict")
source('SCR in Nimble_diftraps_Age_simRS(model1&2).r')

# Create model using projection code (per capita recruitment, no observations)
model <- nimbleModel( code = SCRhab.Open.diftraps.age.PR,
                      constants = nimConstants,
                      data = nimData,
                      #inits = nimInits,
                      check = F,       
                      calculate = T)  

nodesToSim <- model$getDependencies(c("sigD","phi.ad","phi.cub","phi.sub",
                                      "pcr", "beta.dens", 'pi.uncond', 'psi'),
                                    self = F,
                                    downstream = T,
                                    returnScalarComponents = TRUE)
model$simulate(nodesToSim, includeData = FALSE)
N <- apply(model$z,2,sum)


##compile for faster simulation
cmodelSims <- compileNimble(model)

# WE FIND THE S AND Z NODES THAT SHOULDNT BE PREDICTED FOR. 
samplerConfList <- model$getNodeNames()
zNodes <- samplerConfList[grep("z\\[",samplerConfList)]
sNodes <- samplerConfList[grep("sxy",samplerConfList)]
uNodes <- samplerConfList[grep("u\\[",samplerConfList)]
ageNodes <- samplerConfList[grep("age\\[",samplerConfList)]

##get some random iterations from posterior
itera <- sample(1:nrow(sampmat), 20)
nimData1 <- nimData
Nmat<-Rmat<-matrix(NA, 20, 1+t.new)

for(ite in 1:length(itera)){
  # WE UPDATE THE Z VALUES USING THE POSTERIORS PREDICTED Z FOR THE FIVE FIRST YEARS AND THEN 
  # USE NA FOR THE years to predict#
  
  # WE SET NA FOR Z FOR YEARS TO PREDICT - from posterior
  z.est<-matrix(sampmat[itera[ite],z.which], Maug, Tt)
  z.start[1:Maug, 1]<-sampmat[itera[ite],z.which1]
  z.start[(Maug+1):M.new, 1]<-0
  
  nimData1$z <- z.start
  nimData1$u <- z.start
  
  # we set the values in the model - yr 1 aren't nodes in model (only data)
  values(cmodelSims, zNodes) <- nimData1$z#[,2:10] # Fill the values predicted by the model by new values where the extra years are NA
  values(cmodelSims, uNodes) <- nimData1$z#[,2:10] # Fill the values predicted by the model by new values where the extra years are NA
  
  # WE SET SXY
  sxy.start[1:Maug,,1]<-matrix(sampmat[ itera[ite], s.which[indx] ], Maug, 2)
  nimData1$sxy <- sxy.start
  
  # we set the values in the model
  values(cmodelSims, sNodes) <- nimData1$sxy
  
  
  ###set age
  age.est<-sampmat[itera[ite],age.which]
  ##set everyone not alive at all to all-0, ie, can be recruited in the NEW recruitment 
  ##model (they were never part of the superpopulation)
  z.nosuper<-which(apply(z.est,1,sum)==0)
  age.est[z.nosuper]<-0
  age.est[age.est>5]<-5 #set to max age category
  age.start[1:Maug, 1]<-age.est
  age.start[(Maug+1):M.new, 1]<-0
  
  nimData1$age <- age.start
  values(cmodelSims, ageNodes) <- age.start#[,2:(Tt+t.new)]
  
  ##set pcr, phi, sigD, beta.dens
  ##not sure if needed to set in data and cmodelSims...
  nimData1$pcr<-mean(R[itera[ite],]/sampmat[itera[ite],N.which[1:4]])
  values(cmodelSims,"pcr") <- nimData1$pcr
  
  nimData1$phi.ad<-sampmat[itera[ite],'phi.ad']
  values(cmodelSims,"phi.ad") <- nimData1$phi.ad
  nimData1$phi.cub<-sampmat[itera[ite],'phi.cub']
  values(cmodelSims,"phi.cub") <- nimData1$phi.cub  
  nimData1$phi.sub<-sampmat[itera[ite],'phi.sub']
  values(cmodelSims,"phi.sub") <- nimData1$phi.sub
  
  nimData1$sigD<-sampmat[itera[ite],'sigD']
  values(cmodelSims,"sigD") <- nimData1$sigD
  
  nimData1$beta.dens<-sampmat[itera[ite],'beta.dens']
  values(cmodelSims,"beta.dens") <- nimData1$beta.dens
  
  #now we simulate 
  cmodelSims$simulate(nodes = nodesToSim,#c(zNodes,sNodes,yNodes),
                      includeData = F)#---if TRUE: want to simulate new values also for nodes considered as data
  
  Nmat[ite,] <- apply(cmodelSims$z,2,sum)
  Rmat[ite,] <- cmodelSims$R
}

##plot trajectory
plot(1:(1+t.new), apply(Nmat,2,mean), type='l', ylim=range(Nmat, na.rm=TRUE))
for (ite in 1:20){
  points(1:(1+t.new), Nmat[ite,], type='l', col='lightgrey')
}
points(1:(1+t.new), apply(Nmat,2,mean), type='l')
