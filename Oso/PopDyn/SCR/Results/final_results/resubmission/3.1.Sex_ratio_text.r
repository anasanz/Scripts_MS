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
Xbuf2 <- readOGR("Buffer_8500_traps_sxyObs_resub.shp") # This sampling buffer includes AC of observed individuals a bit outside the trapping array

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

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/RESUBMISSION/3-3.4_allparams")
load("myResults_RESUB_3-3.4_sxy.RData")

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


