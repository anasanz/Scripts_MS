## -------------------------------------------------
##                openSCR + denscov + Age structure
##                    Our data 2017-2021
##                        ONLY Syst
##                  Different trap arrays per year
##                  Effort covariates + sigma[sex]
## ------------------------------------------------- 
rm(list = ls())

library(nimble)
library(MCMCvis)
library(nimbleSCR)
library(raster)
library(rgeos)
library(oSCR)
library(terra)
library(sp)
library(dplyr)
library(parallel)

setwd("D:/MargSalas/Scripts_MS/Stats/Nimble")
#source('dbinomLocal_normalBear.R')
source('dbinomLocal_normalBear_rbinom2.R')

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")

#---- 1. LOAD THE DETECTION DATA ---- 

load("edf1721.RData")

# As the model is set, there can be only one capture per trap and occasion
# --> The number of trials is 7 (From may to November): So max number of captures per trap = 7
# To fix it in this edf, remove duplicates
edf <- edf[-which(duplicated(edf)), ]

load("tdf2017_effort.RData")
load("tdf2018_effort.RData")
load("tdf2019_effort.RData")
load("tdf2020_effort.RData")
load("tdf2021_effort.RData")

tdf_all <- rbind(tdf2017[,1:3], tdf2018[,1:3], tdf2019[,1:3],
                 tdf2020[,1:3], tdf2021[,1:3]) # Join to define state space
rownames(tdf_all) <- 1:nrow(tdf_all)


# We remove Nere and Goiat, which are two individuals moving a lot that don't represent
# the rest of the population.

edf <- edf[-which(edf$ind %in% c("Nere", "Goiat")), ]

# Remove cubs from dataset

setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
#setwd("~/Data_server")

info <- readxl::read_excel("info_individuals_2021.xlsx")
#info <- readxl::read_xlsx("D:/MargSalas/Oso/Datos/Tablas_finales/2022/Info_individuals_2021.xlsx", sheet = 1)

info <- info[,c(4,8,10)]

# Create matrix with age data

# Matrix with exact ages
Tt <- 5

x <- matrix(NA, nrow = length(unique(edf$ind)), ncol = Tt) # Matrix to fill
colnames(x) <- c(2017, 2018, 2019, 2020, 2021)
rownames(x) <- unique(edf$ind)

# While filling, 0 counts as an age (first year)

for (i in 1:nrow(x)){
  birth <- as.numeric(info$Year_birth[which(info$ID %in% rownames(x)[i])])
  if (birth >= 2017){
    x[i,which(colnames(x) %in% birth):Tt] <- 0:(Tt-which(colnames(x) %in% birth))
  } else { x[i,] <- (2017-birth):((2017-birth)+(Tt-1)) }
}

### USE THIS TO REMOVE INDIVIDUALS THAT ARE ONLY CUBS IN THE DATASET
cubs <- rownames(x)[apply(x,1,sum, na.rm  = TRUE) < 2]
edf <- edf[-which(edf$ind %in% cubs), ]

# There is an error in the model: All u the first year infinite likelihood.
# Remove all detections of cubs, see if that is causing the problem

delete <- NULL

for (i in 1:nrow(edf)){
  tmp1 <- edf[i,]
  tmpage <- x[which(rownames(x) %in% tmp1$ind), tmp1$session]
  if (tmpage < 2 | is.na(tmpage)){
    delete <- c(delete,i)
  }
  
}

edf <- edf[-delete,]

#---- 2. DEFINE THE TRAP AND THE HABITAT  ---- 
#----   2.1 GET TRAPS---- 
X <- tdf_all[,c(2,3)]
colnames(X) <- c('x', 'y')
J <- dim(X)[1]

#----   2.2 DEFINE STATE SPACE EXTENT ---- 
# State space coordinates
# Buffer: 25000 (used by Maelis, also ~3*sigma in pre-analysis where sig = 6640)
xmin <- min(X[,1]) - 25000
ymin <- min(X[,2]) - 25000
xmax <- max(X[,1]) + 25000
ymax <- max(X[,2]) + 25000
e <- as(raster::extent(xmin, xmax, ymin, ymax), "SpatialPolygons") # Extent of state space


#----   2.3 GET A RASTER FOR THE HABITAT ---- 
# USE A FOREST RASTER TO GET A BASIS FOR THE HABITAT RASTER
# Set up a raster file 
setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
#setwd("~/Data_server/Variables_hrscale")
distcore <- raster("logDistcore_hrbear.tif")

# Crop it to extent of state-space
habitat.r <- crop(distcore, e) 

#----   2.4 DEFINE THE BUFFER AREA AND CUT WHAT IS NOT HABITAT---- 
# Buffer around traps (5*sigma = 33200)
Xpoints <- X
coordinates(Xpoints) <- Xpoints[,c(1,2)]
Xbuf <- gBuffer(Xpoints, width = 25000)

distcoreMask <- rasterize(Xbuf, habitat.r, mask = TRUE)

#RETAIN HABITAT COORDINATES THAT ARE WITHIN THE HABITAT
G <- coordinates(distcoreMask)[!is.na(distcoreMask[]),]
colnames(G) <- c("x","y")
#PLOT CHECK 
#plot(distcoreMask)
#points(Xpoints)
#points(G[,2]~G[,1],col="red",cex=0.5)

#----   2.5 RESCALE HABITAT AND TRAP COORDINATES ---- 
###scale X and G so that bottom left corner of state space is origin (0,0)
sc.coord <- scaleCoordsToHabitatGrid(coordsData = X,
                                     coordsHabitatGridCenter = G)

##this returns S also in row order (all x for a given y) but starting top left 
## corner 
G.sc <- sc.coord$coordsHabitatGridCenterScaled
X.sc <- sc.coord$coordsDataScaled

###get cell coordinates for G.sc
windowCoords <- getWindowCoords(G.sc, plot.check = FALSE)
habitatGrid <- windowCoords$habitatGrid

#----   2.6 HABITAT COVARIATES ---- 
#[CM] THIS MIGHT BE WHERE THE ISSUE WAS
# AS YOU DID FOR THE COORDINATES THE HABITAT, YOU NEED TO SELECT HAB COORDS THAT ARE CONSIDERED AS HABITAT
X.d <- values(distcoreMask)[!is.na(distcoreMask[])]

# Scale
X.d_mean <- mean(X.d)
X.d_sd <- sd(X.d)
X.d_sc <- (X.d - X.d_mean) / X.d_sd

#----   2.7 GET THE LOCAL DETECTORS OBJECTS  ---- 
# USE THE HABITAT GRID PROVIDED BY GETWINDOWCOORDS
habitatMask <- habitatGrid
habitatMask[habitatMask>0] <- 1 #TURN CELL ID TO 1 TO DEFINE THE HABITAT

# Format trap array for nimble: ARRAY WITH TRAP MATRIX PER YEAR

Jyear <- c(nrow(tdf2017), nrow(tdf2018), nrow(tdf2019), nrow(tdf2020), nrow(tdf2021)) # Traps indexes
idx <- 1:nrow(tdf_all)
J.year <- list(idx[1]:Jyear[1], 
               (idx[Jyear[1]]+1):(idx[Jyear[1]] + idx[Jyear[2]]),
               (idx[Jyear[1]] + idx[Jyear[2]]+1): (idx[Jyear[1]] + idx[Jyear[2]] + idx[Jyear[3]]),
               (idx[Jyear[1]] + idx[Jyear[2]] + idx[Jyear[3]] + 1) : (idx[Jyear[1]] + idx[Jyear[2]] + idx[Jyear[3]] + idx[Jyear[4]]),
               (idx[Jyear[1]] + idx[Jyear[2]] + idx[Jyear[3]] + idx[Jyear[4]] + 1) : (idx[Jyear[1]] + idx[Jyear[2]] + idx[Jyear[3]] + idx[Jyear[4]] + idx[Jyear[5]]))


Tt <- 5
Yrs <- seq(1:Tt)

Xt <- Xt.sc <- list() 
for (t in 1:Tt){
  Xt.sc[[t]] <- X.sc[J.year[[t]],] # For the getLocalTraps function
  Xt[[t]] <- X[J.year[[t]],] # To simulate yearly detection data later on
}

Xt.sc.array <- array(NA, c(max(Jyear), 2, Tt))
for (t in 1:Tt){
  Xt.sc.array[1:Jyear[t],,t] <- as.matrix(Xt.sc[[t]])
}

Xt.sc.array[1:Jyear[1], 1:2, 1] # Example: to get traps from year 1

# Visualize traps per year before local evaluation
#for (t in 1:Tt){
#  plot(distcoreMask, main = t)
#  points(Xpoints[J.year[[t]], ])
#  points(G[,2]~G[,1], col="red", cex=0.5)
#} # We will need to use a different dmax per year (in year 1 there are no traps in a region)


# Get one localtraps per year
##determine which traps are within some threshold distance of each habitat grid 


# 5 times log sigma
(5*6640)/5000

dmax <- c(20,10,10,10,9) ## PROBLEM: First year the local evaluation is almost useless

localTraps <- localTrapsNum.l <- MaxLocalTraps.l <- list()
for (t in 1:Tt){
  localTraps[[t]] <- getLocalObjects(habitatMask, Xt.sc[[t]], resizeFactor = 1, dmax = dmax[[t]], plot.check = FALSE)
  localTrapsNum.l[[t]]  <- localTraps[[t]]$numLocalIndices
  MaxLocalTraps.l[[t]]  <- localTraps[[t]]$numLocalIndicesMax
}

##check that there are always local traps available
lapply(localTrapsNum.l, function(x){min(x)})
##yes, every cell has at least some local traps

# Store objects of localtraps as numeric objects because nimble doesn't support lists

localTrapsNum <- do.call(cbind, localTrapsNum.l) # Matrix with cols = years
habitatGridDet <- localTraps[[1]]$habitatGrid # The habitat grid is the same for all
MaxLocalTraps <- unlist(MaxLocalTraps.l) # This will also act as an index to map the dimensions of localTrapsIndex in the model

localTrapsIndex <- array(NA, c(max(habitatGridDet), max(MaxLocalTraps), Tt)) # Array with dimensions of maximun maxlocaltraps, will need to subset in the model
for (t in 1:Tt){
  localTrapsIndex[,1:MaxLocalTraps[t],t] <- localTraps[[t]]$localIndices
}

##some characteristics of S, not affected by changing trap array
numHabWindows <- sum(habitatGrid !=0) #number of cells in S
numGridRows <- nrow(localTraps[[1]]$habitatGrid) # I take the first year of local traps but it doesn't matter, all the same
numGridCols <- ncol(localTraps[[1]]$habitatGrid)


#the following two are still passed to function but no longer used by it
lowerHabCoords <- windowCoords$lowerHabCoords
upperHabCoords <- windowCoords$upperHabCoords

#----   2.8. EFFORT COVARIATE ---- 

K <- 7 # 7 occasions

effort<-array(NA, c(max(Jyear), K, Tt))
tdf_list <- list(tdf2017, tdf2018, tdf2019, tdf2020, tdf2021)

for (t in 1:Tt){
  effort[1:Jyear[t],,t] <- as.numeric(as.matrix(tdf_list[[t]][1:Jyear[t], 6:12]))
}

# Create dummy variable

effort.dummy <- array(1, c(max(Jyear), K, Tt, 2)) # 4th dimension includes 2 arrays: one per level (intercept doesnt count)

# effort.dummy[,,,1] =0  effort.dummy[,,,2] =0 ==> 1 visit in France (cat 1): Intercept, no need to add
# effort.dummy[,,,1] =1  effort.dummy[,,,2] =0 ==> 2 visit in France (cat 2): Multiply b.effort1*array #1 
# effort.dummy[,,,1] =0  effort.dummy[,,,2] =1 ==> Spain (cat 3): Multiply b.effort2*array #2

for (t in 1:Tt){
  tmp <-tmp2 <- tmp3 <- effort[,,t]
  
  # Dummy variable 2 visits in France (only de 2 appear as 1)
  tmp2[tmp2[] %in% c(1,3)] <- 0
  tmp2[tmp2[] %in% c(2)] <- 1
  effort.dummy[,,t,1] <- tmp2
  
  # Dummy variable trap in Spain (only de 3 appear as 1)
  tmp3[tmp3[] %in% c(1,2)] <- 0
  tmp3[tmp3[] %in% c(3)] <- 1
  effort.dummy[,,t,2] <- tmp3
}

#----   2.9. TRAP COVARIATE ---- 

# 3D array to store it together with effort (fast run)
trap <- array(NA, c(max(Jyear), K, Tt)) 

for (t in 1:Tt){
  for(k in 1:K){
    trap[1:Jyear[t],k,t] <- as.numeric(as.factor(as.matrix(tdf_list[[t]][1:Jyear[t], 4])))
  }
}
trap[trap[] %in% c(2)] <- 0

#---- TRAP + EFFORT AS ARRAY FOR FAST FUNCTION ----#

effortarray <- array(1, c(max(Jyear), K, Tt, 3))
effortarray[,,,1] <- effort.dummy[,,,1]
effortarray[,,,2] <- effort.dummy[,,,2]
effortarray[,,,3] <- trap


#----   2.10. SEX COVARIATE ---- 

setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
info <- readxl::read_xlsx("D:/MargSalas/Oso/Datos/Tablas_finales/2022/info_individuals_2021.xlsx", sheet = 1)

#setwd("~/data/data/Data_server")
#info <- readxl::read_xlsx("~/data/data/Data_server/info_individuals_2021.xlsx", sheet = 1)

info <- info[,c(4,5)]
colnames(info)[1] <- "ind"

# Arrange in the same order than detections data frame below
indsex <- data.frame(ind = unique(edf$ind))
indsex <- left_join(indsex,info, by = "ind")

indsex$Sex[indsex$Sex == "F"] <- 0
indsex$Sex[indsex$Sex == "M"] <- 1

sex <- as.numeric(indsex$Sex)


#---- 3. DETECTION DATA   ---- 


#----   3.1 MAKE Y   ---- 

n <- length(unique(edf$ind))

Y <- array(0, c(n, max(Jyear), K, Tt))
rownames(Y) <- unique(edf$ind)
xx <- edf[,c(1,2,3,4)]

for (t in 1:Tt){
  xxt <- xx[which(xx$session == t), ]
  
  for (k in 1:K){
    xxtk <- xxt[which(xxt$occ == k), ]
    
    if(nrow(xxtk) == 0) next # If there were no detections in that occasion
    
    for (obs in 1:nrow(xxtk)) {
      Y[xxtk[obs, 2], xxtk[obs, 4], xxtk[obs, 3], t] <- 1 # ASP: 1 because it can only be detected once per trap and occasion
    }
  }
}

max(Y) # Check the number max number of detections (it can't be higher than K)

#----   3.1.1 BEHAVIOURAL RESPONSE COVARIATE FROM Y   ---- 

prevcap <- list()
for (s in 1:dim(Y)[4]) {
  Ys <- Y[,,,s]
  prevcap[[s]] <- array(0, dim = c(dim(Ys)[1], dim(Ys)[2], 
                                   dim(Ys)[3]))
  first <- matrix(0, dim(Ys)[1], dim(Ys)[2])
  for (i in 1:dim(Ys)[1]) {
    for (j in 1:dim(Ys)[2]) {
      if (sum(Ys[i, j, ]) > 0) {
        first[i, j] <- min((1:(dim(Ys)[3]))[Ys[i, j, 
        ] > 0])
        prevcap[[s]][i, j, 1:first[i, j]] <- 0
        if (first[i, j] < dim(Ys)[3]) 
          prevcap[[s]][i, j, (first[i, j] + 1):(dim(Ys)[3])] <- 1
      }
    }
  }
  # zeros <- array(0, c(1, dim(prevcap[[s]])[2], dim(prevcap[[s]])[3]))
  # prevcap[[s]] <- abind(prevcap[[s]], zeros, along = 1)
}

##
prevcapArray <- array(0,dim(Y))
for(t in 1:dim(Y)[4]){
  prevcapArray[,,,t] <- prevcap[[t]]
}


#----   3.2 AUGMENT Y   ---- 
##augment observed data to size M
M <- 400
nz <- M-n

y.in <- array(0, c(M, max(Jyear), K, Tt))
y.in[1:n,,,] <- Y

prevcapArray.in <- array(0, c(M, max(Jyear), K, Tt))
prevcapArray.in[1:n,,,] <- prevcapArray

##augment observed sex variable to size M (it becomes latent)
sex <- c(sex,rep(NA,length((n+1):M)))

#----   3.3 USE SPARSE FORMAT FOR Y   ---- 
##change to 'sparse' format - speeds up computation by reducing file size
## getSparseY cannot handle 4d arrays, so loop over years to get a 3d array per year
y.sparse <- list()
for (t in 1:Tt){
  y.sparse[[t]] <- getSparseY(y.in[,,,t]) 
}

##extract pieces to be passed to Nimble 

max.max <- max(sapply(y.sparse, function(x)x$lengthYCombined)) # ASP: Maximun of the maxDetNums

y.sp  <- array(NA, c(M, max.max, K, Tt))
for (t in 1:Tt){
  y.sp[,1:y.sparse[[t]]$lengthYCombined,,t] <- y.sparse[[t]]$yCombined
}

lengthYCombined <- sapply(y.sparse, function(x)x$lengthYCombined)

##number of trials per trap - now always 1 but needs to be passed to Nimble
# always 1, because we model each occasion separately
ones <- rep(1, max(Jyear))

## ---- 3.4 INCLUDE AGE DATA ----

# Ro-construct age matrix, now without cubs

setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
#setwd("~/Data_server")

info <- readxl::read_excel("info_individuals_2021.xlsx")
#info <- readxl::read_xlsx("D:/MargSalas/Oso/Datos/Tablas_finales/2022/Info_individuals_2021.xlsx", sheet = 1)

info <- info[,c(4,8,10)]

# Create matrix with age data

# Matrix with exact ages
x <- matrix(NA, nrow = nrow(Y), ncol = Tt) # Matrix to fill
colnames(x) <- c(2017, 2018, 2019, 2020, 2021)
rownames(x) <- rownames(Y)

# While filling, 0 counts as an age (first year)

for (i in 1:nrow(x)){
  birth <- as.numeric(info$Year_birth[which(info$ID %in% rownames(x)[i])])
  if (birth >= 2017){
    x[i,which(colnames(x) %in% birth):Tt] <- 0:(Tt-which(colnames(x) %in% birth))
  } else { x[i,] <- (2017-birth):((2017-birth)+(Tt-1)) }
}

# Sum 1 so that the first year of life is age = 1
age <- x+1
#set known unborn individuals to 0
age[is.na(age)]<-0

###convert to 'raw' age categories -> NEW FORMULATION WHEN DELETING cubs
### age = 0: not yet recruited (category 1)
### age = 1: not yet recruited (category 1) 
### age = 2: not yet recruited (category 1) 
### age = 3: subadult 1 (category 2)
### age = 4: subadult 2 (category 3)
### age >=5: adult (category 4)
### category 1: dead

age.cat<-age
age.cat[age==0]<-1
age.cat[age==1]<-1
age.cat[age==2]<-1
age.cat[age==3]<-2
age.cat[age==4]<-3
age.cat[age>=5]<-4

ageMatAug <- matrix(NA,M,Tt)
ageMatAug[1:n,] <- age.cat #Matrix with age categories augmented

# here, this is highest age category excluding dead and unrecruited, so 3 WITHOUT CUBS
max.age <- 3

## ---- Known z ----

# To include the death information into the z, I create an age matrix including the deaths (value 0)
# The function of this matrix is ONLY to then create the last_alive vector
age.cat.z <- age.cat
for (i in 1:nrow(age.cat)){
  death <- info$Confirmed_death[which(info$ID %in% rownames(age.cat.z)[i])]
  if (death == "Alive") next
  if (death == 2021 | death == 2022) next #This has to be modified depending on the years of study
  death <- as.numeric(death) + 1 # We don't know when it died exactly, so we set it as death the year after
  # I could dig into this and get the specific date for some individuals
  age.cat.z[i,which(colnames(age.cat.z) %in% death):Tt] <- 0 
}

# Check month of death of: Hvala, Nougat, Fifonet, Cachou, Plume, Gribouille                               

##some stuff for known z state below
#year of recruitment for observed individuals
#r <- apply(age,1,function(x) min(which(x!=0)))
r <- apply(age,1,function(x) min(which(x>2))) # Now the recruitment year is that when the age is > 2 (come directly as subadults, age 3)
Y2d <- apply(Y, c(1,4), sum, na.rm=T)

# Combined vector to inform of the death or last capture (if death not available) of individuals
last_cap <- apply(Y2d,1,function(x) max(which(x>0)))  # last year with detection
last_alive <- apply(age.cat.z,1,function(x) min(which(x == 0))-1) # last occasion alive

last <- last_cap
last[which(last_alive != Inf)] <- last_alive[which(last_alive != Inf)] 

# Known z: state process is informed by both y, age data AND death recoveries

zdatAGE <- matrix(NA, M, Tt)
for(i in 1:n){
  #zdatNoAGE[i, first[i]:last[i]] <- 1         # when ignoring age data, we known an individuals is alive between first and last observation
  #if age is used to inform z
  zdatAGE[i, r[i]:last[i]] <- 1               # alive between known recruitment year and last capture or last occasion alive when known
  if(r[i]>1)  zdatAGE[i, 1:(r[i]-1)] <- 0     # Not entered prior to age==1
  if(last_alive[i] != Inf) zdatAGE[i, (last_alive[i]+1):Tt] <- 0 # Death after the last occasion alive (when available, so when last_alive !=Inf)
}


#---- 4. FIT NIMBLE MODEL    ---- 

############################################################################################
### running a model in parallel ############################################################

##source code to run model in parallel 
setwd("D:/MargSalas/Scripts_MS/Stats/Nimble")
source("Parallel Nimble function FOR aNA2_model3-2.r") ##sorry, caps lock...

#----   4.1 CONSTANT AND DATA    ---- 

##compile constants
nimConstants <- list(
  M = M,
  J = Jyear,
  numHabWindows = numHabWindows, 
  numGridRows = numGridRows,
  numGridCols = numGridCols,
  lengthYCombined = lengthYCombined,
  MaxLocalTraps = MaxLocalTraps,
  nobs = n, 
  Nyr = Tt,
  max.age = max.age,
  K = K,
  effort = effortarray,
  nTrapCovs = dim(effortarray)[4],
  prevcap = prevcapArray.in
)

##compile data
nimData <- list(habDens = X.d_sc,
                y = y.sp,
                lowerHabCoords = lowerHabCoords,
                upperHabCoords = upperHabCoords,
                habitatGrid = habitatGrid,
                ones=ones,           # ASP: ones substitutes K=rep(K, n.max.traps)
                X.sc = Xt.sc.array,
                habitatGridDet = habitatGridDet,
                localTrapsIndex = localTrapsIndex, 
                localTrapsNum = localTrapsNum,
                sex = sex,
                agePlusOne = ageMatAug[,1],       #known age in yr 1
                w = c(rep(1, n), rep(NA, nz)),    #membership in superpopulation
                u = zdatAGE,                      #known alive states
                b = rep(1,K), a = rep(1, max.age) #prior params for recruitment, age distribution
)


#----   4.2 INITIAL VALUES  ---- 
#----     4.2.1 Z AND AGE COMPONENT   ---- 

piAGE.in <- c(0.25, 0.25, 0.50)

zstAGE <- zdatAGE # ASP: u (alive state = zstAGE)
for(i in 1:n){
  if(last[i]<Tt) zstAGE[i,(last[i]+1):Tt] <- 0 ## ASP: Start as dead after the last observation
}
zstAGE[(n+1):M,] <- 0  # start augmented individuals as not entered to prevent error messages.
zstAGE[!is.na(zdatAGE)] <- NA #set known z values to NA

#because setting all augmented guys to 0 did not work:
#randomly generate augmented individuals as part of superpop
w.in <- c(rep(NA, n), rbinom(nz, 1, 0.5))

##generate entrance occasion for all individuals
ent.occ.aug <- sample(1:Tt, nz, replace=TRUE)

## ASP: GEnerate age category of entry individuals
age.cat.in <- c(rep(NA, n), rep(1, nz)) #1=not yet recruited

#randomly assign starting age category in year 1 to all individuals 
#alive in year 1
age.cat.in[(n+1):M][ent.occ.aug==1] <- sample(2:(max.age+1),sum(ent.occ.aug==1) , 
                                              prob=piAGE.in,
                                              replace=TRUE)
## ASP: Only assign age category for year 1, because it is what is used in the model to construct the 
## probabilities

##adjust starting values for u accordingly
for(i in (n+1):M){
  zstAGE[i,ent.occ.aug[i-n]:Tt]<-1
}


#----     4.2.2 SXY  ---- 

##because of local evaluation of possible detectors, activity center initial 
##values have to be specified, eg average capture location in 'model' space

S.in <- array(NA, c(M, 2, Tt))

for ( i in 1:n){
  for (t in 1:Tt){
    caps <- which(apply(Y[i,,,t],1,sum) > 0) ## ASP: Get in which traps the ind was captured at year t
    # With 3d data Need to sum the rows (all occasions) to get idx right
    
    if (length(caps)==0) next #fill in missing ACs with reasonable values later
    if (length(caps)==1){ ## ASP: If its only in 1 trap, a put the trap location as AC
      S.in[i,,t]<- as.numeric(Xt[[t]][caps,])
    }else{ 
      S.in[i,,t]<- as.numeric(apply(Xt[[t]][caps,],2,mean))}  ## ASP: If its > 1 trap average location
  }
}

##fill in missing ACs as average of 'observed' ACs in nearby time step


for (i in 1:n){
  #which ACs unobserved
  nac <- which(is.na(S.in[i,1,]))
  wac <- (1:Tt)[-nac] ## ASP: Only the years that are observed
  #for those, use mean observed
  
  if(length(wac)==1) {
    S.in[i,1,nac]<-S.in[i,1,wac] ## ASP: Use the mean of the observed to fill unobserved
    S.in[i,2,nac]<-S.in[i,2,wac]   
  } else {
    ran.ac<-sample(wac,length(nac), replace=TRUE)
    S.in[i,1,nac]<-S.in[i,1,ran.ac] ## ASP: Use the mean of the observed to fill unobserved
    S.in[i,2,nac]<-S.in[i,2,ran.ac] 
  }
  
}

##random ACs for individuals never observed
for(i in (n+1) : M){
  for (t in 1:Tt){
    ssg<-sample(1:length(X.d), 1)
    S.in[i,,t]<-G[ssg,]
  }
}

colnames(S.in) <- c('x', 'y')
S.in.sc <- scaleCoordsToHabitatGrid(S.in, G)

#----     4.2.3 SEX  ---- 

sex.in <- c(rep(NA,n), rep(0,length((n+1):M))) # NA for known values and 
sex.in[sample((n+1):M, 30, replace = FALSE)] <- 1 # random sex for augmented individuals


#----     4.2.2 COMPILE INITIAL VALUES   ---- 

S.in.sc_coords <- S.in.sc$coordsDataScaled

inits<-function(){list(beta=c(0.15,rep(0.85/(Tt-1), Tt-1)), 
                       sigma=runif(2,0.5, 1.5),
                       p.ad=runif(1,0,0.5),
                       p.sub=runif(1,0,0.5),
                       #p.cub=runif(1,0,0.5),
                       trapBetas = runif(3, 0.5,1),
                       b.bh = runif(1, 0.5,1),
                       omega = runif(1, 0.5,1),
                       sex = sex.in,
                       phi.ad=runif(1,0.5,1),
                       phi.sub=runif(1,0.5,1),
                       #phi.cub=runif(1,0.5,1),
                       beta.dens = runif(1,-0.1, 0.1), 
                       piAGE=piAGE.in,
                       u = zstAGE,
                       w=w.in,
                       agePlusOne=age.cat.in,
                       psi=0.7,
                       sxy=S.in.sc_coords,
                       sigD=runif(1, 1.5, 2.5))}


##source model code
##I prefer working on code in a separate script but you can also have everything in
##one script and just execute the code

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Model")
source('3.6.1.SCRopen_diftraps_AgeNOcubs_trapsBhCov_sigsex in Nimble.r')

##determine which parameters to monitor
params<-c("p.ad", "p.sub","phi.ad","phi.sub", 
          "beta", "psi", "piAGE", "Nsuper", "N", "B", "N.sub", "N.ad",
          'sigma', 'beta.dens', 'sigD',
          'trapBetas', 'b.bh', 'omega')

params2 <- c("sxy", "z", "age.cat") 


###### SAVE FOR RUNNING #####

modelcode = SCRhab.Open.diftraps.ageNOcubs.effortTrapBhCov.sigsex

# For Cyril
setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/3.openSCR_Age/Data_server")
save(nimData, nimConstants, 
     inits, Tt, sex.in, piAGE.in, zstAGE, w.in, age.cat.in ,S.in.sc_coords, 
     params, params2, modelcode, file = "Data_Model3-4_CYRIL.RData")

#### OPTION 1: PARALLEL ####
detectCores()

##start cluster w/ 3 cores (for 3 chains)
this_cluster <- makeCluster(3)


##run wrapper function in parallel - cl and X need to be this way
## cl defines which cluster to use, X provides random number seeds to each core

old <- Sys.time()

chain_output <- parLapply(cl = this_cluster, X = 1:3, 
                          fun = run_MCMC_allcode,      ##function in "Parallel Nimble function.R"
                          data = nimData,              ##your data list
                          code = SCRhab.Open.diftraps.age,   ##your model code
                          inits = inits,                 ##your inits function
                          constants = nimConstants,      ##your list of constants
                          params = params,               ##your vector with params to monitor
                          niter = 150000,                  ##iterations per chain
                          nburnin = 100000,                ##burn-in
                          nthin = 10,                  ##thinning, main parameters
                          Tt = Tt,                     ##additional objects needed within inits
                          piAGE.in = piAGE.in,
                          zstAGE = zstAGE,
                          w.in = w.in,
                          age.cat.in = age.cat.in,
                          S.in.sc_coords = S.in.sc_coords 
)
new <- Sys.time() - old

## ALWAYS close cluster when model is done
stopCluster(this_cluster)

### output is a list 

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/3.openSCRdenscov_Age")
#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/3.openSCRdenscov_Age")
save(chain_output, file = "sampOpenSCR_diftraps_age_2021.RData")




#### OPTION 2: NO PARALLEL (TO TRY INITIAL VALUES AND SEE IF MODEL WORKS) ####

#(1) set up model

model <- nimbleModel(SCRhab.Open.diftraps.ageNOcubs.effortTrapBhCov.sigsex, constants = nimConstants, 
                     data=nimData, inits=inits(), check = FALSE)
##ignore error message, only due to missing initial values at this stage

model$initializeInfo()
#THEN YOU SHOULD ALWAYS CHECK THAT THE MODEL IS ABLE TO CALCULATE A LIKELIHOOD (RETURN A VALUE) GIVEN THE INITIAL VALUES PROVIDED
model$calculate()#
# IF A -INF OF NA IS RETURNED YOU CAN CHECK WHERE THE PROBLEM COMES FROM WITH (AND THEN TRY TO FIX IT UNTIL THE -INF DISAPEARS):

model$logProb_u # ALL U in the first year have infinite likelihood (relation with age structure the first year?) Or how recruitment is set? I have set recruitment so that they enter as subadults.
model$logProb_sxy[is.na(model$logProb_sxy)]
model$logProb_y[is.na(model$logProb_y)]

#(2) Compile model in c++
#     In complex models, this step can take a while (as well as step 5)
#     Much longer than in JAGS, but the model typically runs much faster
cmodel <- compileNimble(model)       

# (3) Configure MCMC - on an uncompiled model - this step allows setting which quantities to monitor
#     Also, nimble allows two sets of monitors, these can be thinned at different rates
#     all of which is more important in complex models but not to start with
conf.mcmc <- configureMCMC(model, monitors = params, thin=10)

# (4) Build the MCMC sampler based on configurations
mcmc <- buildMCMC(conf.mcmc)

# (5) Compile sampler in c++ together with compiled model
cmcmc <- compileNimble(mcmc, project = cmodel, resetFunctions = TRUE)

# (6) Run (monitor time just for fun) [takes 20 seconds on my computer]
system.time(
  (samp <- runMCMC(cmcmc, niter = 10, nburnin = 5, nchains=3, inits = inits) )
)

## -------------------------------------------------
##                   PROCESS RESULTS
## ------------------------------------------------- 

library(MCMCvis)
library(rgdal)
library(coda)
library(nimbleSCR)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/3.openSCRdenscov_Age/2021/Cyril/3-4")
load("myResults_3-4_param.RData")
summary(nimOutput)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")
load("habcoord.RData")


#----  1. WHICH STUDY AREA WILL I KEEP TO ESTIMATE ABUNDANCE? ---- 

setwd("D:/MargSalas/Oso/Datos/GIS/Countries")
Xbuf <- readOGR("Buffer_statespace.shp")
osbuf <- readOGR("Buffer_8500_allobs.shp")
Xbuf2 <- readOGR("Buffer_8500_traps.shp")

# Plots to decide

#plot(Xbuf, col = "lightblue")
#plot(Xbuf2, col = adjustcolor("pink", alpha = 0.5), add = TRUE)
#plot(osbuf, col = adjustcolor("green", alpha = 0.5), add = TRUE)

#---- 2.  ESTIMATE ABUNDANCE IN BUFFER ---- 

# Load results
setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/3.openSCRdenscov_Age/2021/Cyril/3-4")
load("myResults_3-4_sxy.RData")

dim(nimOutputSXY[[1]]) #  iterations * 6000 elements (e.g., z[1,5]) ???
5*300 + 5*300 + 5*300*2 # 6000 elements

dim(myResultsSXYZ$sims.list$sxy) 
dim(myResultsSXYZ$sims.list$z)

# Unscale the sxy coordinates

dimnames(myResultsSXYZ$sims.list$sxy)[[3]] <- c('x','y')
myResultsSXYZ$sims.list$sxy <- scaleCoordsToHabitatGrid(coordsData = myResultsSXYZ$sims.list$sxy,## this are your sxy
                                                        coordsHabitatGridCenter = G,# this is your unscaled habitat (as you used when scaling the habitat/detector to the habitat. G?
                                                        scaleToGrid = FALSE)$coordsDataScaled

# Each iteration is a reality of the model, so for each iteration I take the activity center of alive individuals, test
# if its in the buffer and sum them to get abundance.

# Matrix to store abundance in the buffer each iteration and year
NIn <- matrix(NA,nrow=dim(myResultsSXYZ$sims.list$z)[1],ncol=dim(myResultsSXYZ$sims.list$z)[3]) # nrow = iterations, ncol = year

ite=1
t=1

for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
  for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
    
    which.alive <- which(myResultsSXYZ$sims.list$z[ite,,t]==1) # Select only the individuals alive (z=1)
    
    which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
    
    sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(osbuf))) # CONVERT SXY TO SPATIAL POINTS 
    
    which.In <- over(sp, osbuf) # Check which ones are in the buffer
    
    NIn[ite,t] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
  }
}

#average number of individuals without the buffer for each year 
colMeans(NIn)
NIn[,1]#posterior distrib for the first year

par(mfrow = c(2,3))
for (t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
  plot(density(NIn[,t]), main = t)
  abline(v = colMeans(NIn)[t], col = "blue")
}

# Sum of individuals alive in total each year (without buffer)
NALL <- apply(myResultsSXYZ$sims.list$z,c(1,3),function(x) sum(x==1))
colMeans(NALL)


