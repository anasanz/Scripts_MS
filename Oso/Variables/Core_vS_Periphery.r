
## -------------------------------------------------
##                SPATIAL NUCLEUS vS PERIPHERY
## ------------------------------------------------- 

# POST-HOC ANALYSIS

# Identify the spatial nucleus (core) of the population and the periphery, to estimate the age structure in both areas
# The age structure will be estimated as the number of individuals of each age class in nucleus and periphery

# The core and periphery need to occupy the sampling area (not the state space), so I need to check how to define it
# 1º trial: Buffer of female average dispersal distance around the core polygons (origin) for the distance to core covariate

rm(list = ls())

library(rgdal)
library(dplyr)
library(rgeos)
library(MCMCvis)
library(rgdal)
library(nimbleSCR)


stderr <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

# Load Dispersal distances
setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
d <- read.csv("disp_distance.csv", header = TRUE, row.names = NULL, sep = ",")
mean_dist <- d %>% # Calculate dispersal distance by sex
  group_by(Sex) %>%
  summarise(
    count = n(),
    mean = mean(Distance, na.rm = TRUE),
    sd = sd(Distance, na.rm = TRUE),
    se = stderr(Distance, na.rm = TRUE)
  )

distfem <- mean_dist$mean[1]

# Load core polygon
setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/DistCore/corepol")
corepol <- readOGR("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/DistCore/corepol/corepol.shp", "corepol")
corepol@data$ID <- 1
corepol@data <- corepol@data[2]

# Load state space
setwd("D:/MargSalas/Oso/Datos/GIS/Countries")
Xbuf <- readOGR("Buffer_statespace.shp")
Xbuf2 <- readOGR("Buffer_8500_traps.shp") # Load sampling area (where we estimate abundance)
proj4string(Xbuf2) <- proj4string(corepol)

## ---- 1. Buffer of disp. distances arounf core area----
nuc <- gBuffer(corepol, width = distfem, byid = TRUE)

plot(Xbuf)
plot(Xbuf2, col = "grey", border = NA, add = TRUE)
plot(corepol, col = "red", add = TRUE)
plot(nuc, border = "red", add = TRUE)

# To compare core area vS periphery, I need to create the periphery polygon, extracting the core area.
# I will test the age structure defining core vS periphery in 2 different ways

## ---- 1.1. Core vS Trap buffer ----

# Nucleous: Polygon used to define core area in DistCore variable (corepol)
# Periphery: Sampling buffer (extracting the corepol)

Xbuf2_peri1 <- gDifference(Xbuf2, corepol)

plot(Xbuf2, col = "grey", border = NA)
plot(corepol, col = "red", add = TRUE)
plot(Xbuf2_peri1, col = "green", add = TRUE)

## ---- 1.2. Intersection Disp distance Nucleous and Core vS Trap buffer ----

# Nucleous: Polygon used to define core area in DistCore variable, buffered by female dispersal distance (nuc)
# Periphery: Sampling buffer (extracting the nuc)

Xbuf2_peri2 <- gDifference(Xbuf2, nuc)

plot(Xbuf2, col = "grey", border = NA)
plot(nuc, col = "red", add = TRUE)
plot(Xbuf2_peri2, col = "blue", add = TRUE)

## ---- 2. Estimation of age structure in nuc vS periphery, different definitions from - to + restrictive ----

# Load data

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams")
load("myResults_3-3.1_param.RData")
summary(nimOutput)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")
load("habcoord.RData")

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1") # Load results
load("myResults_3-3.1_sxy.RData")

dimnames(myResultsSXYZ$sims.list$sxy)[[3]] <- c('x','y') # Unscale the sxy coordinates
myResultsSXYZ$sims.list$sxy <- scaleCoordsToHabitatGrid(coordsData = myResultsSXYZ$sims.list$sxy,## this are your sxy
                                                        coordsHabitatGridCenter = G,# this is your unscaled habitat (as you used when scaling the habitat/detector to the habitat. G?
                                                        scaleToGrid = FALSE)$coordsDataScaled

# Function for estimating age structure in two polygons: nucleous and periphery

age_st <- function(nucleus = corepol, periphery = Xbuf2_peri1, sampling_buffer = Xbuf2){
  
  
  # 1. Total abundance in buffer area (sampling area)
  
  # Matrix to store abundance in the buffer each iteration and year
  NIn_trapBuf <- matrix(NA,nrow=dim(myResultsSXYZ$sims.list$z)[1],ncol=dim(myResultsSXYZ$sims.list$z)[3]) # nrow = iterations, ncol = year
  
  for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
    for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
      
      which.alive <- which(myResultsSXYZ$sims.list$z[ite,,t]==1) # Select only the individuals alive (z=1)
      
      which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
      
      sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(sampling_buffer))) # CONVERT SXY TO SPATIAL POINTS 
      
      which.In <- over(sp, sampling_buffer) # Check which ones are in the buffer
      
      NIn_trapBuf[ite,t] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
    }
  }
  
  # 2. Abundance of each age category in both polygons
  
  zones <- c(nucleus, periphery)
  
  ## CUBS
  ZZcubs <- myResultsSXYZ$sims.list$z
  ZZcubs[!myResultsSXYZ$sims.list$age.cat %in% c(1,2) ]  <- 3 # Considers all individuals that are not 1 and 2 as dead
  
  cub_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], 2)) # nrow = iterations, ncol = year, dim3 = nuc/periphery
  
  
  for(s in 1:2){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(ZZcubs[ite,,t]==1) # Select only the individuals alive (z=1)
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        cub_trapBuf[ite,t,s] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}
  
  
  ## SUBADULTS 
  
  ZZsub <- myResultsSXYZ$sims.list$z
  ZZsub[!myResultsSXYZ$sims.list$age.cat %in% c(3,4) ]  <- 3 # Considers all individuals that are not 3 and 4 (SUBADULTS) as dead
  
  subad_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], 2))
  
  for(s in 1:2){
  for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
    for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
      
      which.alive <- which(ZZsub[ite,,t]==1) # Select only the individuals alive (z=1)
      
      which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
      
      sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
      which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
      
      subad_trapBuf[ite,t,s] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
    }}}
  
  
  ## ADULTS
  
  ZZad <- myResultsSXYZ$sims.list$z
  ZZad[!myResultsSXYZ$sims.list$age.cat %in% c(5) ]  <- 3 # Considers all individuals that are not 3 and 4 (SUBADULTS) as dead
  
  ad_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], 2))
  
  for(s in 1:2){
  for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
    for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
      
      which.alive <- which(ZZad[ite,,t]==1) # Select only the individuals alive (z=1)
      
      which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
      
      sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
      which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
      
      ad_trapBuf[ite,t,s] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
    }}}
  
  ## ADULT FEMALES
  
  ZZadFEM <- myResultsSXYZ$sims.list$z
  ZZadFEM[!myResultsSXYZ$sims.list$age.cat %in% c(5) ]  <- 3 # Considers all individuals that are not 5 (ADULTS) as dead
  ZZadFEM[!myResultsSXYZ$sims.list$sex %in% c(0) ]  <- 3 # Considers all individuals that dont have a sex 0 (females) as dead
  
  adFEM_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], 2))
  
  for(s in 1:2){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(ZZad[ite,,t]==1) # Select only the individuals alive (z=1)
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        adFEM_trapBuf[ite,t,s] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}
  
  ## ADULT MALES
  
  ZZadMAL <- myResultsSXYZ$sims.list$z
  ZZadMAL[!myResultsSXYZ$sims.list$age.cat %in% c(5) ]  <- 3 # Considers all individuals that are not 5 (ADULTS) as dead
  ZZadMAL[!myResultsSXYZ$sims.list$sex %in% c(1) ]  <- 3 # Considers all individuals that dont have a sex 1 (males) as dead
  
  adMAL_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], 2))
  
  for(s in 1:2){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(ZZad[ite,,t]==1) # Select only the individuals alive (z=1)
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        adMAL_trapBuf[ite,t,s] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}
  
  
  # 3. Estimate abundance proportion of each age class 
  
  prop_years <- list()
  
  prop <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1],5,2)) # nrow = iterations, ncol = c(cub, subad, ad, adFEM, adMAL), z = Core/Periphery
  
  for(s in 1:2){
    for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      prop[ite, 1, s] <- cub_trapBuf[ite,t,s]/NIn_trapBuf[ite,t]
      prop[ite, 2, s] <- subad_trapBuf[ite,t,s]/NIn_trapBuf[ite,t]
      prop[ite, 3, s] <- ad_trapBuf[ite,t,s]/NIn_trapBuf[ite,t]
      prop[ite, 3, s] <- adFEM_trapBuf[ite,t,s]/NIn_trapBuf[ite,t]
      prop[ite, 3, s] <- adMAL_trapBuf[ite,t,s]/NIn_trapBuf[ite,t]
      
    }
    prop_years[[t]] <- prop
  }}
  
  return(results)
  
}


  

## ---- 2.2. Disp distance Nucleous vS whole periphery ----

