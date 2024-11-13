
## -------------------------------------------------
##               Map density from model 3.1 
##          Explore abundances per age class and sex
## ------------------------------------------------- 

rm(list = ls())

library(Rcpp)
library(RcppArmadillo)
library(RcppProgress)
library(rgdal)
library(raster)

setwd("D:/MargSalas/Scripts_MS/Functions/Nimble")

# Load functions
#sourceCpp("GetSpaceUse_PD.cpp")
sourceCpp("GetDensity_PD.cpp")
source("getDensityInput.R")

# Load buffers
setwd("D:/MargSalas/Oso/Datos/GIS/Countries")
Xbuf <- readOGR("Buffer_statespace.shp")
#Xbuf2 <- readOGR("Buffer_8500_traps.shp")
Xbuf2 <- readOGR("Buffer_8500_traps_sxyObs.shp") # This sampling buffer includes AC of observed individuals a bit outside the trapping array

# Conver to rasters

setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
r <- raster("logDistcore_hrbear.tif")
values(r) <- NA
r <- crop(r,Xbuf) 

Xbuf_raster <- rasterize(Xbuf, r)
values(Xbuf_raster)[which(is.na(values(Xbuf_raster)))] <- 0 # Raster of 0 and 1

Xbuf2_raster <- rasterize(Xbuf2, r)
values(Xbuf2_raster)[which(values(Xbuf2_raster) == 1)] <- 2
values(Xbuf2_raster)[which(is.na(values(Xbuf2_raster)))] <- 0 # Raster of 0 and 2

ras <- overlay(Xbuf_raster, Xbuf2_raster, fun = max)
f <- rasterize(Xbuf, ras, mask = TRUE)
#values(f)[which(values(f) == 1)] <- 0
#values(f)[which(values(f) == 2)] <- 1

# Load posterior distribution
library(nimbleSCR) # Load nimbleSCR here, otherwise it gets in conflict with raster package, weird

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL")
load("myResults_3-3.1_sxy.RData")

# Load original habitat coordinates
setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")
load("habcoord.RData")

#give dim names to your posteriors sxy
dim(myResultsSXYZ$sims.list$sxy) # YOURPOSTERIORSXY 
dimnames(myResultsSXYZ$sims.list$sxy)[[3]] <- c("x","y")

## first rescale the coordinates to the original scale 
myResultsSXYZ$sims.list$sxy <- scaleCoordsToHabitatGrid(coordsData = myResultsSXYZ$sims.list$sxy,
  coordsHabitatGridCenter = G,
  scaleToGrid = FALSE)$coordsDataScaled

f1 <- f
f1[f1%in%c(1,2)] <-  1
f[] <- as.factor(as.character(f[]))

##GET OBJECTS IN SHAPE
densityInputRegions <- getDensityInput( 
  regions = f,## THIS  A RASTER FILE WITH 0/1 HABITAT VS BUFFER, OR 1/2 fRANCE SPAIN... WHATEVER YOU WANT. 
  habitat = f1,## here put the same than regions argument. 
  s = myResultsSXYZ$sims.list$sxy,
  plot.check = TRUE)

## extract density
yearnames <- c("2017", "2018", "2019", "2020", "2021")

## ---- Explore how to get density and store abundance in sampling buffer and ss per age class and sex ----

# Storage:
summary_class <- array(NA, dim = c(3,5,5)) # To store summary of each age /sex class and year
colnames(summary_class) <- c("mean","median", "mode", "95%CILow", "95%CIHigh")
rownames(summary_class) <- c("StateSpace", "SamplingBuffer", "Total")

all_class <- list() # Store all summary_class

######  ALL INDIVIDUALS  #####

DensityCountriesRegions <- list()
n.years = 5
for(t in 1:n.years){
  DensityCountriesRegions[[t]] <- GetDensity_PD(
    sx = densityInputRegions$sx[,,t],# X COORDINATES
    sy =  densityInputRegions$sy[,,t],# Y COORDINATES
    z = myResultsSXYZ$sims.list$z[,,t],# Z 
    IDmx = densityInputRegions$habitat.id,
    aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
    regionID = densityInputRegions$regions.rgmx,
    returnPosteriorCells = F)
}#t

#mean density in each cell 
DensityCountriesRegions[[1]]$MeanCell

#plot()
ACDens<-list()
for(t in 1:n.years){
  ACDens[[t]] <- densityInputRegions$regions.r
  ACDens[[t]][] <- NA
  ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  plot(ACDens[[t]], main = paste(yearnames[t], "ALL"))
}

for(t in 1:n.years){
summary_class[,,t] <- DensityCountriesRegions[[t]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
}

all_class$all <- summary_class

######  CUBS  #####

ZZcubs <- myResultsSXYZ$sims.list$z
ZZcubs[!myResultsSXYZ$sims.list$age.cat %in% c(1,2) ]  <- 3 # Considers all individuals that are not 1 and 2 as dead

DensityCountriesRegions <- list()
n.years = 5
for(t in 1:n.years){
  DensityCountriesRegions[[t]] <- GetDensity_PD(
    sx = densityInputRegions$sx[,,t],# X COORDINATES
    sy =  densityInputRegions$sy[,,t],# Y COORDINATES
    z = ZZcubs[,,t],# Z 
    IDmx = densityInputRegions$habitat.id,
    aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
    regionID = densityInputRegions$regions.rgmx,
    returnPosteriorCells = F)
}#t

#mean density in each cell 
DensityCountriesRegions[[1]]$MeanCell

#plot()
ACDens<-list()
for(t in 1:n.years){
  ACDens[[t]] <- densityInputRegions$regions.r
  ACDens[[t]][] <- NA
  ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  plot(ACDens[[t]], main = paste(yearnames[t], "CUBS"))
}

for(t in 1:n.years){
  summary_class[,,t] <- DensityCountriesRegions[[t]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
}

all_class$cubs <- summary_class

######  SUBADULTS  #####

ZZsub <- myResultsSXYZ$sims.list$z
ZZsub[!myResultsSXYZ$sims.list$age.cat %in% c(3,4) ]  <- 3 # Considers all individuals that are not 3 and 4 (SUBADULTS) as dead

DensityCountriesRegions <- list()
n.years = 5
for(t in 1:n.years){
  DensityCountriesRegions[[t]] <- GetDensity_PD(
    sx = densityInputRegions$sx[,,t],# X COORDINATES
    sy =  densityInputRegions$sy[,,t],# Y COORDINATES
    z = ZZsub[,,t],# Z 
    IDmx = densityInputRegions$habitat.id,
    aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
    regionID = densityInputRegions$regions.rgmx,
    returnPosteriorCells = F)
}#t

#mean density in each cell 
DensityCountriesRegions[[1]]$MeanCell

#plot()
ACDens<-list()
for(t in 1:n.years){
  ACDens[[t]] <- densityInputRegions$regions.r
  ACDens[[t]][] <- NA
  ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  plot(ACDens[[t]], main = paste(yearnames[t], "SUBADULTS"))
}

for(t in 1:n.years){
  summary_class[,,t] <- DensityCountriesRegions[[t]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
}

all_class$sub <- summary_class


######  ADULTS  #####

ZZad <- myResultsSXYZ$sims.list$z
ZZad[!myResultsSXYZ$sims.list$age.cat %in% c(5) ]  <- 3 # Considers all individuals that are not 5 (adults) as dead

DensityCountriesRegions <- list()
n.years = 5
for(t in 1:n.years){
  DensityCountriesRegions[[t]] <- GetDensity_PD(
    sx = densityInputRegions$sx[,,t],# X COORDINATES
    sy =  densityInputRegions$sy[,,t],# Y COORDINATES
    z = ZZad[,,t],# Z 
    IDmx = densityInputRegions$habitat.id,
    aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
    regionID = densityInputRegions$regions.rgmx,
    returnPosteriorCells = F)
}#t

#mean density in each cell 
DensityCountriesRegions[[1]]$MeanCell

#plot()
ACDens<-list()
for(t in 1:n.years){
  ACDens[[t]] <- densityInputRegions$regions.r
  ACDens[[t]][] <- NA
  ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  plot(ACDens[[t]], main = paste(yearnames[t], "ADULTS"))
}

for(t in 1:n.years){
  summary_class[,,t] <- DensityCountriesRegions[[t]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
}

all_class$ad <- summary_class

######  FEMALES  #####
# Sex: 0 Females; 1 Males
ZZfem <- myResultsSXYZ$sims.list$z
ZZfem[!myResultsSXYZ$sims.list$sex %in% c(0) ]  <- 3 # Considers all individuals that are not 0 (females) as dead

DensityCountriesRegions <- list()
n.years = 5
for(t in 1:n.years){
  DensityCountriesRegions[[t]] <- GetDensity_PD(
    sx = densityInputRegions$sx[,,t],# X COORDINATES
    sy =  densityInputRegions$sy[,,t],# Y COORDINATES
    z = ZZfem[,,t],# Z 
    IDmx = densityInputRegions$habitat.id,
    aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
    regionID = densityInputRegions$regions.rgmx,
    returnPosteriorCells = F)
}#t

#mean density in each cell 
DensityCountriesRegions[[1]]$MeanCell

#plot()
ACDens<-list()
for(t in 1:n.years){
  ACDens[[t]] <- densityInputRegions$regions.r
  ACDens[[t]][] <- NA
  ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  plot(ACDens[[t]], main = paste(yearnames[t], "FEMALES"))
}

for(t in 1:n.years){
  summary_class[,,t] <- DensityCountriesRegions[[t]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
}

all_class$females <- summary_class

######  MALES  #####
# Sex: 0 Females; 1 Males
ZZmal <- myResultsSXYZ$sims.list$z
ZZmal[!myResultsSXYZ$sims.list$sex %in% c(1) ]  <- 3 # Considers all individuals that are not 0 (females) as dead

DensityCountriesRegions <- list()
n.years = 5
for(t in 1:n.years){
  DensityCountriesRegions[[t]] <- GetDensity_PD(
    sx = densityInputRegions$sx[,,t],# X COORDINATES
    sy =  densityInputRegions$sy[,,t],# Y COORDINATES
    z = ZZmal[,,t],# Z 
    IDmx = densityInputRegions$habitat.id,
    aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
    regionID = densityInputRegions$regions.rgmx,
    returnPosteriorCells = F)
}#t

#mean density in each cell 
DensityCountriesRegions[[1]]$MeanCell

#plot()
ACDens<-list()
for(t in 1:n.years){
  ACDens[[t]] <- densityInputRegions$regions.r
  ACDens[[t]][] <- NA
  ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  plot(ACDens[[t]], main = paste(yearnames[t], "MALES"))
}

for(t in 1:n.years){
  summary_class[,,t] <- DensityCountriesRegions[[t]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
}

all_class$males <- summary_class

######  ADULT FEMALES  #####
# Sex: 0 Females; 1 Males

ZZadFEM <- myResultsSXYZ$sims.list$z
ZZadFEM[!myResultsSXYZ$sims.list$age.cat %in% c(5) ]  <- 3 # Considers all individuals that are not 5 (ADULTS) as dead
ZZadFEM[!myResultsSXYZ$sims.list$sex %in% c(0) ]  <- 3 # Considers all individuals that dont have a sex 0 (females) as dead

DensityCountriesRegions <- list()
n.years = 5
for(t in 1:n.years){
  DensityCountriesRegions[[t]] <- GetDensity_PD(
    sx = densityInputRegions$sx[,,t],# X COORDINATES
    sy =  densityInputRegions$sy[,,t],# Y COORDINATES
    z = ZZadFEM[,,t],# Z 
    IDmx = densityInputRegions$habitat.id,
    aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
    regionID = densityInputRegions$regions.rgmx,
    returnPosteriorCells = F)
}#t

#mean density in each cell 
DensityCountriesRegions[[1]]$MeanCell

#plot()
ACDens<-list()
for(t in 1:n.years){
  ACDens[[t]] <- densityInputRegions$regions.r
  ACDens[[t]][] <- NA
  ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  plot(ACDens[[t]], main = paste(yearnames[t], "FEMALES"))
}

for(t in 1:n.years){
  summary_class[,,t] <- DensityCountriesRegions[[t]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
}

all_class$adFemales <- summary_class

######  ADULT MALES  #####
# Sex: 0 Females; 1 Males

ZZadMAL <- myResultsSXYZ$sims.list$z
ZZadMAL[!myResultsSXYZ$sims.list$age.cat %in% c(5) ]  <- 3 # Considers all individuals that are not 5 (ADULTS) as dead
ZZadMAL[!myResultsSXYZ$sims.list$sex %in% c(1) ]  <- 3 # Considers all individuals that dont have a sex 1 (males) as dead

DensityCountriesRegions <- list()
n.years = 5
for(t in 1:n.years){
  DensityCountriesRegions[[t]] <- GetDensity_PD(
    sx = densityInputRegions$sx[,,t],# X COORDINATES
    sy =  densityInputRegions$sy[,,t],# Y COORDINATES
    z = ZZadMAL[,,t],# Z 
    IDmx = densityInputRegions$habitat.id,
    aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
    regionID = densityInputRegions$regions.rgmx,
    returnPosteriorCells = F)
}#t

#mean density in each cell 
DensityCountriesRegions[[1]]$MeanCell

#plot()
ACDens<-list()
for(t in 1:n.years){
  ACDens[[t]] <- densityInputRegions$regions.r
  ACDens[[t]][] <- NA
  ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  plot(ACDens[[t]], main = paste(yearnames[t], "MALES"))
}

for(t in 1:n.years){
  summary_class[,,t] <- DensityCountriesRegions[[t]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
}

all_class$adMales <- summary_class

## ---- Methods ----
###### Abundance estimates and BCI ######
all_class$all

###### Sex ratio (proportion of males) ######
apply(all_class$males, c(1,2), sum)[2,1]/apply(all_class$all, c(1,2), sum)[2,1] # Mean
apply(all_class$males, c(1,2), sum)[2,4]/apply(all_class$all, c(1,2), sum)[2,4] # Low CI
apply(all_class$males, c(1,2), sum)[2,5]/apply(all_class$all, c(1,2), sum)[2,5] # Up CI


## ---- Descriptive results (Methods): Explore sex ratio ----

######  ALL INDIVIDUALS  #####

## Check if males and females sum the total number of individuals in sampling buffer
all_class$males[2,1,5] + all_class$females[2,1,5] == all_class$all[2,1,5] # YES

## Sex ratio per year

sex.ratio <- matrix(NA, nrow = n.years, ncol = 2)
rownames(sex.ratio) <- yearnames
colnames(sex.ratio) <- c("females", "males")
for(t in 1:n.years){
  sex.ratio[t,1] <- all_class$females[2,1,t]/all_class$all[2,1,t]
  sex.ratio[t,2] <- all_class$males[2,1,t]/all_class$all[2,1,t]
} # Sex structure is shifting along years

## Total sex ratio only in buffer: Very skewed to females! 70/30

# This led us to include the ACs of observed males into the buffer:
# (new Xbuf2: "Buffer_8500_traps_sxyObs.shp")
# Now it is a bit less skewed as we include more males: 66/34

apply(all_class$males, c(1,2), sum)[2,1]/apply(all_class$all, c(1,2), sum)[2,1]
apply(all_class$females, c(1,2), sum)[2,1]/apply(all_class$all, c(1,2), sum)[2,1] # Mean


# To write in METHODS:
## Proportion males
apply(all_class$males, c(1,2), sum)[2,1]/apply(all_class$all, c(1,2), sum)[2,1] # Mean
apply(all_class$males, c(1,2), sum)[2,4]/apply(all_class$all, c(1,2), sum)[2,4] # Low CI
apply(all_class$males, c(1,2), sum)[2,5]/apply(all_class$all, c(1,2), sum)[2,5] # Up CI


# Get 
## Total sex ratio in all ss: Also quite skewed to females: 68/32

apply(all_class$males, c(1,2), sum)[3,1]/apply(all_class$all, c(1,2), sum)[3,1]
apply(all_class$females, c(1,2), sum)[3,1]/apply(all_class$all, c(1,2), sum)[3,1]

## Check if we leave outside a higher proportion of males than females (and that is why sex ratio is skewed)

apply(all_class$males, c(1,2), sum)[1,1]/apply(all_class$males, c(1,2), sum)[3,1] # Proportion of males outside
apply(all_class$females, c(1,2), sum)[1,1]/apply(all_class$females, c(1,2), sum)[3,1] # Proportion of males outside

prop.outside <- sex.ratio
for(t in 1:n.years){
  prop.outside[t,1] <- all_class$females[1,1,t]/all_class$females[3,1,t]
  prop.outside[t,2] <- all_class$males[1,1,t]/all_class$males[3,1,t]
} # We leave a slight highest proportion of males outside, but not so many really (5% more)


######  ADULTs  #####

## Check if adult females and males sum the total number of adults in sampling buffer
all_class$adMales[2,1,5] + all_class$adFemales[2,1,5] == all_class$ad[2,1,5] # YES

## Sex ratio adults per year

sex.ratio.adults <- matrix(NA, nrow = n.years, ncol = 2)
rownames(sex.ratio.adults) <- yearnames
colnames(sex.ratio.adults) <- c("females", "males")
for(t in 1:n.years){
  sex.ratio.adults[t,1] <- all_class$adFemales[2,1,t]/all_class$ad[2,1,t]
  sex.ratio.adults[t,2] <- all_class$adMales[2,1,t]/all_class$ad[2,1,t]
} # Sex structure is shifting along years

## Total sex ratio adults: Very skewed to females, even more than in all individuals! 75/25
apply(all_class$adFemales, c(1,2), sum)[2,1]/apply(all_class$ad, c(1,2), sum)[2,1]
apply(all_class$adMales, c(1,2), sum)[2,1]/apply(all_class$ad, c(1,2), sum)[2,1]

## Total sex ratio adults in all ss: Also quite skewed to females: 72/27

apply(all_class$adMales, c(1,2), sum)[3,1]/apply(all_class$ad, c(1,2), sum)[3,1]
apply(all_class$adFemales, c(1,2), sum)[3,1]/apply(all_class$ad, c(1,2), sum)[3,1]

## Check if we leave outside a higher proportion of males than females (and that is why sex ratio is skewed)

apply(all_class$adMales, c(1,2), sum)[1,1]/apply(all_class$adMales, c(1,2), sum)[3,1] # Proportion of males outside
apply(all_class$adFemales, c(1,2), sum)[1,1]/apply(all_class$adFemales, c(1,2), sum)[3,1] # Proportion of females outside

prop.outside.adults <- sex.ratio
for(t in 1:n.years){
  prop.outside.adults[t,1] <- all_class$adFemales[1,1,t]/all_class$adFemales[3,1,t]
  prop.outside.adults[t,2] <- all_class$adMales[1,1,t]/all_class$adMales[3,1,t]
} # We leave a slight highest proportion of males outside, but not so many really (5% more)


## ---- Explore ONLY OBSERVED INDIVIDUALS  ----
# SUBSET: ONLY SEEN INDIVIDUALS (FIRST 61)

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL")
load("myResults_3-3.1_sxy.RData")

dim(myResultsSXYZ$sims.list) # YOURPOSTERIORSXY 

myResultsSubset<- myResultsSXYZ

dim(myResultsSubset$sims.list$age.cat)
myResultsSubset$sims.list$age.cat <- myResultsSubset$sims.list$age.cat[ ,1:61,]

dim(myResultsSubset$sims.list$sex)
myResultsSubset$sims.list$sex <- myResultsSubset$sims.list$sex[ ,1:61]

dim(myResultsSubset$sims.list$sxy)
myResultsSubset$sims.list$sxy <- myResultsSubset$sims.list$sxy[ ,1:61,,]

dim(myResultsSubset$sims.list$z)
myResultsSubset$sims.list$z <- myResultsSubset$sims.list$z[ ,1:61,]

#give dim names to your posteriors sxy
dim(myResultsSubset$sims.list$sxy) # YOURPOSTERIORSXY 
dimnames(myResultsSubset$sims.list$sxy)[[3]] <- c("x","y")

## first rescale the coordinates to the original scale 
myResultsSubset$sims.list$sxy <- scaleCoordsToHabitatGrid(coordsData = myResultsSubset$sims.list$sxy,
                                                          coordsHabitatGridCenter = G,
                                                          scaleToGrid = FALSE)$coordsDataScaled

##GET OBJECTS IN SHAPE
densityInputRegions <- getDensityInput( 
  regions = f,## THIS  A RASTER FILE WITH 0/1 HABITAT VS BUFFER, OR 1/2 fRANCE SPAIN... WHATEVER YOU WANT. 
  habitat = f1,## here put the same than regions argument. 
  s = myResultsSubset$sims.list$sxy,
  plot.check = TRUE)

# Storage:
summary_class_observed <- array(NA, dim = c(3,5,5)) # To store summary of each age /sex class and year
colnames(summary_class_observed) <- c("mean","median", "mode", "95%CILow", "95%CIHigh")
rownames(summary_class_observed) <- c("StateSpace", "SamplingBuffer", "Total")

all_class_observed <- list() # Store all summary_class

######  ALL INDIVIDUALS  #####

DensityCountriesRegions <- list()
n.years = 5
for(t in 1:n.years){
  DensityCountriesRegions[[t]] <- GetDensity_PD(
    sx = densityInputRegions$sx[,,t],# X COORDINATES
    sy =  densityInputRegions$sy[,,t],# Y COORDINATES
    z = myResultsSubset$sims.list$z[,,t],# Z 
    IDmx = densityInputRegions$habitat.id,
    aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
    regionID = densityInputRegions$regions.rgmx,
    returnPosteriorCells = F)
}#t


#plot()
ACDens<-list()
for(t in 1:n.years){
  ACDens[[t]] <- densityInputRegions$regions.r
  ACDens[[t]][] <- NA
  ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  plot(ACDens[[t]], xaxt = "n", yaxt = "n", bty="n")
}


for(t in 1:n.years){
  summary_class_observed[,,t] <- DensityCountriesRegions[[t]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
}

all_class_observed$all <- summary_class_observed

######  FEMALES  ###### 
# Sex: 0 Females; 1 Males
ZZfem <- myResultsSubset$sims.list$z
ZZfem[!myResultsSubset$sims.list$sex %in% c(0) ]  <- 3 # Considers all individuals that are not 0 (females) as dead

DensityCountriesRegions <- list()
n.years = 5
for(t in 1:n.years){
  DensityCountriesRegions[[t]] <- GetDensity_PD(
    sx = densityInputRegions$sx[,,t],# X COORDINATES
    sy =  densityInputRegions$sy[,,t],# Y COORDINATES
    z = ZZfem[,,t],# Z 
    IDmx = densityInputRegions$habitat.id,
    aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
    regionID = densityInputRegions$regions.rgmx,
    returnPosteriorCells = F)
}#t


#plot()
ACDens<-list()
for(t in 1:n.years){
  ACDens[[t]] <- densityInputRegions$regions.r
  ACDens[[t]][] <- NA
  ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  plot(ACDens[[t]], main = paste(yearnames[t], "FEMALES"))
}

for(t in 1:n.years){
  summary_class_observed[,,t] <- DensityCountriesRegions[[t]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
}

all_class_observed$females <- summary_class_observed

######  MALES  ###### 
# Sex: 0 Females; 1 Males
ZZmal <- myResultsSubset$sims.list$z
ZZmal[!myResultsSubset$sims.list$sex %in% c(1) ]  <- 3 # Considers all individuals that are not 0 (females) as dead

DensityCountriesRegions <- list()
n.years = 5
for(t in 1:n.years){
  DensityCountriesRegions[[t]] <- GetDensity_PD(
    sx = densityInputRegions$sx[,,t],# X COORDINATES
    sy =  densityInputRegions$sy[,,t],# Y COORDINATES
    z = ZZmal[,,t],# Z 
    IDmx = densityInputRegions$habitat.id,
    aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
    regionID = densityInputRegions$regions.rgmx,
    returnPosteriorCells = F)
}#t


#plot()
ACDens<-list()
for(t in 1:n.years){
  ACDens[[t]] <- densityInputRegions$regions.r
  ACDens[[t]][] <- NA
  ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  plot(ACDens[[t]], main = paste(yearnames[t], "MALES"))
}

for(t in 1:n.years){
  summary_class_observed[,,t] <- DensityCountriesRegions[[t]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
}

all_class_observed$males <- summary_class_observed

######  ADULTS  #####

ZZad <- myResultsSubset$sims.list$z
ZZad[!myResultsSubset$sims.list$age.cat %in% c(5) ]  <- 3 # Considers all individuals that are not 5 (adults) as dead

DensityCountriesRegions <- list()
n.years = 5
for(t in 1:n.years){
  DensityCountriesRegions[[t]] <- GetDensity_PD(
    sx = densityInputRegions$sx[,,t],# X COORDINATES
    sy =  densityInputRegions$sy[,,t],# Y COORDINATES
    z = ZZad[,,t],# Z 
    IDmx = densityInputRegions$habitat.id,
    aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
    regionID = densityInputRegions$regions.rgmx,
    returnPosteriorCells = F)
}#t

#mean density in each cell 
DensityCountriesRegions[[1]]$MeanCell

#plot()
ACDens<-list()
for(t in 1:n.years){
  ACDens[[t]] <- densityInputRegions$regions.r
  ACDens[[t]][] <- NA
  ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  plot(ACDens[[t]], main = paste(yearnames[t], "ADULTS"))
}

for(t in 1:n.years){
  summary_class_observed[,,t] <- DensityCountriesRegions[[t]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
}

all_class_observed$ad <- summary_class_observed


######  ADULT FEMALES  ###### 
# Sex: 0 Females; 1 Males
ZZadFEM <- myResultsSubset$sims.list$z
ZZadFEM[!myResultsSubset$sims.list$sex %in% c(0) ]  <- 3 # Considers all individuals that are not 0 (females) as dead
ZZadFEM[!myResultsSubset$sims.list$age.cat %in% c(5) ]  <- 3 # Considers all individuals that are not 5 (ADULTS) as dead

DensityCountriesRegions <- list()
n.years = 5
for(t in 1:n.years){
  DensityCountriesRegions[[t]] <- GetDensity_PD(
    sx = densityInputRegions$sx[,,t],# X COORDINATES
    sy =  densityInputRegions$sy[,,t],# Y COORDINATES
    z = ZZadFEM[,,t],# Z 
    IDmx = densityInputRegions$habitat.id,
    aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
    regionID = densityInputRegions$regions.rgmx,
    returnPosteriorCells = F)
}#t


#plot()
ACDens<-list()
for(t in 1:n.years){
  ACDens[[t]] <- densityInputRegions$regions.r
  ACDens[[t]][] <- NA
  ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  plot(ACDens[[t]], main = paste(yearnames[t], "FEMALES"))
}

for(t in 1:n.years){
  summary_class_observed[,,t] <- DensityCountriesRegions[[t]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
}

all_class_observed$adFemales <- summary_class_observed

######  ADULT MALES  ###### 
# Sex: 0 Females; 1 Males
ZZadMAL <- myResultsSubset$sims.list$z
ZZadMAL[!myResultsSubset$sims.list$sex %in% c(1) ]  <- 3 # Considers all individuals that are not 0 (females) as dead
ZZadMAL[!myResultsSubset$sims.list$age.cat %in% c(5) ]  <- 3 # Considers all individuals that are not 5 (ADULTS) as dead

DensityCountriesRegions <- list()
n.years = 5
for(t in 1:n.years){
  DensityCountriesRegions[[t]] <- GetDensity_PD(
    sx = densityInputRegions$sx[,,t],# X COORDINATES
    sy =  densityInputRegions$sy[,,t],# Y COORDINATES
    z = ZZadMAL[,,t],# Z 
    IDmx = densityInputRegions$habitat.id,
    aliveStates = 1,# WHICH Z STATE IS CONSIDERED ALIVE, E.G. IF MULTIPLE = C(1,2)
    regionID = densityInputRegions$regions.rgmx,
    returnPosteriorCells = F)
}#t


#plot()
ACDens<-list()
for(t in 1:n.years){
  ACDens[[t]] <- densityInputRegions$regions.r
  ACDens[[t]][] <- NA
  ACDens[[t]][!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  plot(ACDens[[t]], main = paste(yearnames[t], "MALES"))
}

for(t in 1:n.years){
  summary_class_observed[,,t] <- DensityCountriesRegions[[t]]$summary # 1 = habitat (state-space); 2 = buffer (sampling area)
}

all_class_observed$adMales <- summary_class_observed

## ---- Descriptive results ONLY OBSERVED: Explore sex ratio ----

######  ALL INDIVIDUALS ####

## Check if males and females sum the total number of individuals in sampling buffer
all_class_observed$males[2,1,5] + all_class_observed$females[2,1,5] == all_class_observed$all[2,1,5] # YES

## Sex ratio per year
sex.ratio.observed <- matrix(NA, nrow = n.years, ncol = 2)
rownames(sex.ratio.observed) <- yearnames
colnames(sex.ratio.observed) <- c("females", "males")

for(t in 1:n.years){
  sex.ratio.observed[t,1] <- all_class_observed$females[2,1,t]/all_class_observed$all[2,1,t]
  sex.ratio.observed[t,2] <- all_class_observed$males[2,1,t]/all_class_observed$all[2,1,t]
} # Sex structure of observed individuals remains stable. I think it has to do with the p of detection
# Males move more and more detected. More uncertainty in females, pushes up estimates?

## Total sex ratio: Very skewed to females! 63/36
apply(all_class_observed$males, c(1,2), sum)[2,1]/apply(all_class_observed$all, c(1,2), sum)[2,1]
apply(all_class_observed$females, c(1,2), sum)[2,1]/apply(all_class_observed$all, c(1,2), sum)[2,1]

## Total sex ratio in all ss: Also quite skewed to females: 68/32
apply(all_class_observed$males, c(1,2), sum)[3,1]/apply(all_class_observed$all, c(1,2), sum)[3,1]
apply(all_class_observed$females, c(1,2), sum)[3,1]/apply(all_class_observed$all, c(1,2), sum)[3,1]

#### Check if we leave outside a higher proportion of males than females (and that is why sex ratio is skewed)
apply(all_class_observed$males, c(1,2), sum)[1,1]/apply(all_class_observed$males, c(1,2), sum)[3,1] # Proportion of males outside
apply(all_class_observed$females, c(1,2), sum)[1,1]/apply(all_class_observed$females, c(1,2), sum)[3,1] # Proportion of males outside

prop.outside.observed <- sex.ratio
for(t in 1:n.years){
  prop.outside.observed[t,1] <- all_class_observed$females[1,1,t]/all_class_observed$females[3,1,t]
  prop.outside.observed[t,2] <- all_class_observed$males[1,1,t]/all_class_observed$males[3,1,t]
} 

# From the observed individuals, we leave 30% outside, while we leave only 7% of the females outside across years
# Males have a lower p, and that is why the model places them outside?

######  ADULTs  #####

## Check if adult females and males sum the total number of adults in sampling buffer
all_class_observed$adMales[2,1,5] + all_class_observed$adFemales[2,1,5] == all_class_observed$ad[2,1,5] # YES

## Sex ratio adults per year

sex.ratio.observed.adults <- matrix(NA, nrow = n.years, ncol = 2)
rownames(sex.ratio.observed.adults) <- yearnames
colnames(sex.ratio.observed.adults) <- c("females", "males")
for(t in 1:n.years){
  sex.ratio.observed.adults[t,1] <- all_class_observed$adFemales[2,1,t]/all_class_observed$ad[2,1,t]
  sex.ratio.observed.adults[t,2] <- all_class_observed$adMales[2,1,t]/all_class_observed$ad[2,1,t]
} # Sex structure is shifting along years

## Total sex ratio adults: Very skewed to females, even more than in all individuals! 75/25
apply(all_class_observed$adFemales, c(1,2), sum)[2,1]/apply(all_class_observed$ad, c(1,2), sum)[2,1]
apply(all_class_observed$adMales, c(1,2), sum)[2,1]/apply(all_class_observed$ad, c(1,2), sum)[2,1]

## Total sex ratio adults in all ss: Also quite skewed to females: 72/27

apply(all_class_observed$adMales, c(1,2), sum)[3,1]/apply(all_class_observed$ad, c(1,2), sum)[3,1]
apply(all_class_observed$adFemales, c(1,2), sum)[3,1]/apply(all_class_observed$ad, c(1,2), sum)[3,1]

## Check if we leave outside a higher proportion of males than females (and that is why sex ratio is skewed)

apply(all_class_observed$adMales, c(1,2), sum)[1,1]/apply(all_class_observed$adMales, c(1,2), sum)[3,1] # Proportion of males outside
apply(all_class_observed$adFemales, c(1,2), sum)[1,1]/apply(all_class_observed$adFemales, c(1,2), sum)[3,1] # Proportion of females outside

prop.outside.observed.adults <- sex.ratio
for(t in 1:n.years){
  prop.outside.observed.adults[t,1] <- all_class_observed$adFemales[1,1,t]/all_class_observed$adFemales[3,1,t]
  prop.outside.observed.adults[t,2] <- all_class_observed$adMales[1,1,t]/all_class_observed$adMales[3,1,t]
} # We leave a slight highest proportion of males outside, but not so many really (5% more)




## ---- Plot to see where AC of adult males fall ----

Tt <- 5
M.aug <- 300

# Put sxy in right format for plotting
sampmat2 <- do.call(rbind, nimOutputSXY) 
s.which <- grep('sxy', colnames(sampmat2)) # ASP: index columns all sxy (sampmat matrix)
sampmat2_sxy <- sampmat2[, s.which]
dim(sampmat2_sxy)
mean_sampmat <- colMeans(sampmat2_sxy) # I will plot the mean location over iterations

sxy <- array(NA, c(300, 2, 5))
for(t in 1:Tt){
  s.which.year <- grep(paste(t,"]", sep = ""), colnames(sampmat2_sxy)) # ASP: index columns all sxy (sampmat matrix)
  sxy[,,t] <- matrix(mean_sampmat[s.which.year] , M.aug, 2) 
}

dimnames(sxy)[[2]] <- c('x','y') 
sxy.uns <- scaleCoordsToHabitatGrid(coordsData = sxy,## this are your sxy
                                          coordsHabitatGridCenter = G,# this is your unscaled habitat (as you used when scaling the habitat/detector to the habitat. G?
                                          scaleToGrid = FALSE)$coordsDataScaled

######  ALL ADULT MALES ####

ZZadMAL <- myResultsSXYZ$sims.list$z
ZZadMAL[!myResultsSXYZ$sims.list$age.cat %in% c(5) ]  <- 0 # Considers all individuals that are not 5 (ADULTS) as dead
ZZadMAL[!myResultsSXYZ$sims.list$sex %in% c(1) ]  <- 0 # Considers all individuals that dont have a sex 1 (males) as dead

dim(ZZadMAL)
ZZadMAL_mean <- colMeans(ZZadMAL)
ZZadMAL_mean2 <- ifelse(ZZadMAL_mean > 0.5, 1, 0)

ZZadMAL_mean2[100:300]
setwd("D:/MargSalas/Oso/OPSCR_project/Results/Plots/model3.1")
pdf("AC_AdultMales_ALL.pdf",7,7)
par(mfrow = c(3,2))
for (t in 1:5){
  plot(Xbuf, main = yearnames[t])
  plot(Xbuf2, add = TRUE)
  p <- sxy.uns[which(ZZadMAL_mean2[,t] == 1),,t]
  sp <- SpatialPoints(p, proj4string=CRS(proj4string(Xbuf2)))
  points(sp, pch= 21, col = "red")
}
dev.off()

######  ONLY OBSERVED ADULT MALES ####

all_class_observed$adMales
ZZadMAL_obs <- myResultsSubset$sims.list$z
ZZadMAL_obs[!myResultsSubset$sims.list$sex %in% c(1) ]  <- 0 # Considers all individuals that are not 0 (females) as dead
ZZadMAL_obs[!myResultsSubset$sims.list$age.cat %in% c(5) ]  <- 0 # Considers all individuals that are not 5 (ADULTS) as dead

sxy.uns_obs <- sxy.uns[1:61,,]

dim(ZZadMAL)
ZZadMAL_obs_mean <- colMeans(ZZadMAL_obs)
ZZadMAL_obs_mean2 <- ifelse(ZZadMAL_obs_mean > 0.5, 1, 0)

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Plots/model3.1")
pdf("AC_AdultMales_obs.pdf",7,7)
par(mfrow = c(3,2))
for (t in 1:5){
  plot(Xbuf, main = yearnames[t])
  plot(Xbuf2, add = TRUE)
  p <- sxy.uns_obs[which(ZZadMAL_obs_mean2[,t] == 1),,t]
  sp <- SpatialPoints(p, proj4string=CRS(proj4string(Xbuf2)))
  points(sp, pch= 21, col = "red")
}
dev.off()

all_class$adMales[,,5]
all_class_observed$adMales[,,5]

######  ALL ADULT FEMALES ####

ZZadFEM <- myResultsSXYZ$sims.list$z
ZZadFEM[!myResultsSXYZ$sims.list$age.cat %in% c(5) ]  <- 0 # Considers all individuals that are not 5 (ADULTS) as dead
ZZadFEM[!myResultsSXYZ$sims.list$sex %in% c(0) ]  <- 0 # Considers all individuals that dont have a sex 1 (males) as dead

dim(ZZadFEM)
ZZadFEM_mean <- colMeans(ZZadFEM)
ZZadFEM_mean2 <- ifelse(ZZadFEM_mean > 0.5, 1, 0)

ZZadFEM_mean2[100:300]

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Plots/model3.1")
pdf("AC_AdultFemales_ALL.pdf",7,7)
par(mfrow = c(3,2))
for (t in 1:5){
  plot(Xbuf, main = yearnames[t])
  plot(Xbuf2, add = TRUE)
  p <- sxy.uns[which(ZZadFEM_mean2[,t] == 1),,t]
  sp <- SpatialPoints(p, proj4string=CRS(proj4string(Xbuf2)))
  points(sp, pch= 21, col = "red")
}
dev.off()

######  ONLY OBSERVED ADULT FEMALES ####

ZZadFEM_obs <- myResultsSubset$sims.list$z
ZZadFEM_obs[!myResultsSubset$sims.list$sex %in% c(0) ]  <- 0 # Considers all individuals that are not 0 (females) as dead
ZZadFEM_obs[!myResultsSubset$sims.list$age.cat %in% c(5) ]  <- 0 # Considers all individuals that are not 5 (ADULTS) as dead

sxy.uns_obs <- sxy.uns[1:61,,]

dim(ZZadFEM)
ZZadFEM_obs_mean <- colMeans(ZZadFEM_obs)
ZZadFEM_obs_mean2 <- ifelse(ZZadFEM_obs_mean > 0.5, 1, 0)

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Plots/model3.1")
pdf("AC_AdultFemales_obs.pdf",7,7)
par(mfrow = c(3,2))
for (t in 1:5){
  plot(Xbuf, main = yearnames[t])
  plot(Xbuf2, add = TRUE)
  p <- sxy.uns_obs[which(ZZadFEM_obs_mean2[,t] == 1),,t]
  sp <- SpatialPoints(p, proj4string=CRS(proj4string(Xbuf2)))
  points(sp, pch= 21, col = "red")
}
dev.off()

all_class$adFemales[,,5]
all_class_observed$adFemales[,,5]

## ---- Plot to see where AC fall: Only observed individuals known alive ----

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")
load("zObserved_yAgeDeaths.RData")

zdatAGE[is.na(zdatAGE)] <- 0
ageMatAug
sex

######  ALL INDIVIDUALS ####

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Plots/model3.1")
pdf("AC_observed_alive.pdf",7,7)
par(mfrow = c(3,2))
for (t in 1:5){
  plot(Xbuf, main = yearnames[t])
  plot(Xbuf2, add = TRUE)
  p <- sxy.uns[which(zdatAGE[,t] == 1),,t]
  sp <- SpatialPoints(p, proj4string=CRS(proj4string(Xbuf2)))
  points(sp, pch= 21, col = "red")
  
  over(sp,Xbuf2)
}
dev.off()

# Identify which ones are outside
o <- list()
w <- list()
for (t in 1:5){
  which.alive <- which(zdatAGE[,t] == 1)
  p <- sxy.uns[which.alive,,t]
  sp <- SpatialPoints(p, proj4string=CRS(proj4string(Xbuf2)))
  which.in <- which(is.na(over(sp,Xbuf2)))
  w$age <- ageMatAug[which.alive[which.in],t]
  w$sex <- sex[which.alive[which.in]]
  o[[t]] <- w
} # All are adult and subadult males
