## ------------------------------------------------------------------------------------
##                      SPATIAL NUCLEUS vS PERIPHERY
##  Definition BIG CORE: Polygon female dispersal distance around original core
## ------------------------------------------------------------------------------------

# POST-HOC ANALYSIS

# Identify the spatial nucleus (core) of the population and the periphery, to estimate the age structure in both areas
# The age structure will be estimated as the number of individuals of each age class in nucleus and periphery


rm(list = ls())

## ---- Load stuff ----

library(rgdal)
library(dplyr)
library(rgeos)
library(MCMCvis)
library(rgdal)
library(nimbleSCR)
library(raster)

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

# Load distCore (ONLY FOR VISUALIZATION PURPOSES HERE)
setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
logDistcore <- raster("logDistcore_hrbear.tif")

# Load data
setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams")
load("myResults_3-3.1_param.RData")

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")
load("habcoord.RData")

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams") # Load results
load("myResults_3-3.1_sxy.RData")

dimnames(myResultsSXYZ$sims.list$sxy)[[3]] <- c('x','y') # Unscale the sxy coordinates
myResultsSXYZ$sims.list$sxy <- scaleCoordsToHabitatGrid(coordsData = myResultsSXYZ$sims.list$sxy,## this are your sxy
                                                        coordsHabitatGridCenter = G,# this is your unscaled habitat (as you used when scaling the habitat/detector to the habitat. G?
                                                        scaleToGrid = FALSE)$coordsDataScaled

## ---- 1. Define core and periphery (big core polygon) ----

# Buffer of disp. distances arounf core area
nuc <- gBuffer(corepol, width = distfem, byid = TRUE)

# Intersection Disp distance Nucleous and Core (BIG CORE) vS Trap buffer 
# Nucleous: Polygon used to define core area in DistCore variable, buffered by female dispersal distance and cropped so that it falls inside buffer (nuc2)
# Periphery: Sampling buffer (extracting the nuc): Xbuf2_peri2
nuc2 <- crop(nuc, Xbuf2)
Xbuf2_peri2 <- gDifference(Xbuf2, nuc2)

## ---- 2. Estimation of age structure in nuc vS periphery per sex and year ----

age_st_full_sex <- function(nucleus, periphery){
  
  zones <- c(nucleus, periphery)
  
  # 2. Abundance of each age category and sex in each zone
  
  ##FEMALES
  
  ## FEMALE 1
  
  ZZFEM1 <- myResultsSXYZ$sims.list$z
  ZZFEM1[!myResultsSXYZ$sims.list$age.cat %in% c(1) ]  <- 3 # Considers all individuals that are not 1 and 2 as dead
  ZZFEM1[!myResultsSXYZ$sims.list$sex %in% c(0) ]  <- 3 # Considers all individuals that dont have a sex 0 (females) as dead
  
  
  FEM1_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], length(zones))) # nrow = iterations, ncol = year, dim3 = Zone
  
  for(s in 1:length(zones)){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(ZZFEM1[ite,,t]==1) # Select only the individuals alive (z=1)
        
        if(length(which.alive) == 0) next
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        if(length(which.alive) == 1){
          which.aliveSXY <- data.frame(x = which.aliveSXY[1], y = which.aliveSXY[2])
        }
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        FEM1_trapBuf[ite,t,s] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}
  
  ## FEMALE 2
  
  ZZFEM2 <- myResultsSXYZ$sims.list$z
  ZZFEM2[!myResultsSXYZ$sims.list$age.cat %in% c(2) ]  <- 3 # Considers all individuals that are not 1 and 2 as dead
  ZZFEM2[!myResultsSXYZ$sims.list$sex %in% c(0) ]  <- 3 # Considers all individuals that dont have a sex 0 (females) as dead
  
  
  FEM2_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], length(zones))) # nrow = iterations, ncol = year, dim3 = Zone
  
  for(s in 1:length(zones)){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(ZZFEM2[ite,,t]==1) # Select only the individuals alive (z=1)
        
        if(length(which.alive) == 0) next
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        if(length(which.alive) == 1){
          which.aliveSXY <- data.frame(x = which.aliveSXY[1], y = which.aliveSXY[2])
        }
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        FEM2_trapBuf[ite,t,s] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}
  
  
  ## FEMALE 3
  
  ZZFEM3 <- myResultsSXYZ$sims.list$z
  ZZFEM3[!myResultsSXYZ$sims.list$age.cat %in% c(3) ]  <- 3 # Considers all individuals that are not 1 and 2 as dead
  ZZFEM3[!myResultsSXYZ$sims.list$sex %in% c(0) ]  <- 3 # Considers all individuals that dont have a sex 0 (females) as dead
  
  
  FEM3_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], length(zones))) # nrow = iterations, ncol = year, dim3 = Zone
  
  for(s in 1:length(zones)){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(ZZFEM3[ite,,t]==1) # Select only the individuals alive (z=1)
        
        if(length(which.alive) == 0) next
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        if(length(which.alive) == 1){
          which.aliveSXY <- data.frame(x = which.aliveSXY[1], y = which.aliveSXY[2])
        }
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        FEM3_trapBuf[ite,t,s] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}
  
  ## FEMALE 4
  
  ZZFEM4 <- myResultsSXYZ$sims.list$z
  ZZFEM4[!myResultsSXYZ$sims.list$age.cat %in% c(4) ]  <- 3 # Considers all individuals that are not 1 and 2 as dead
  ZZFEM4[!myResultsSXYZ$sims.list$sex %in% c(0) ]  <- 3 # Considers all individuals that dont have a sex 0 (females) as dead
  
  
  FEM4_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], length(zones))) # nrow = iterations, ncol = year, dim3 = Zone
  
  for(s in 1:length(zones)){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(ZZFEM4[ite,,t]==1) # Select only the individuals alive (z=1)
        
        if(length(which.alive) == 0) next
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        if(length(which.alive) == 1){
          which.aliveSXY <- data.frame(x = which.aliveSXY[1], y = which.aliveSXY[2])
        }
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        FEM4_trapBuf[ite,t,s] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}
  
  ## FEMALE 5
  
  ZZFEM5 <- myResultsSXYZ$sims.list$z
  ZZFEM5[!myResultsSXYZ$sims.list$age.cat %in% c(5) ]  <- 3 # Considers all individuals that are not 1 and 2 as dead
  ZZFEM5[!myResultsSXYZ$sims.list$sex %in% c(0) ]  <- 3 # Considers all individuals that dont have a sex 0 (females) as dead
  
  
  FEM5_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], length(zones))) # nrow = iterations, ncol = year, dim3 = Zone
  
  for(s in 1:length(zones)){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(ZZFEM5[ite,,t]==1) # Select only the individuals alive (z=1)
        
        if(length(which.alive) == 0) next
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        if(length(which.alive) == 1){
          which.aliveSXY <- data.frame(x = which.aliveSXY[1], y = which.aliveSXY[2])
        }
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        FEM5_trapBuf[ite,t,s] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}
  
  
  ##MALES
  
  ## MALE 1
  
  ZZMAL1 <- myResultsSXYZ$sims.list$z
  ZZMAL1[!myResultsSXYZ$sims.list$age.cat %in% c(1) ]  <- 3 # Considers all individuals that are not 1 and 2 as dead
  ZZMAL1[!myResultsSXYZ$sims.list$sex %in% c(1) ]  <- 3 # Considers all individuals that dont have a sex 0 (MALEs) as dead
  
  
  MAL1_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], length(zones))) # nrow = iterations, ncol = year, dim3 = Zone
  
  for(s in 1:length(zones)){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(ZZMAL1[ite,,t]==1) # Select only the individuals alive (z=1)
        
        if(length(which.alive) == 0) next
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        if(length(which.alive) == 1){
          which.aliveSXY <- data.frame(x = which.aliveSXY[1], y = which.aliveSXY[2])
        }
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        MAL1_trapBuf[ite,t,s] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}
  
  ## MALE 2
  
  ZZMAL2 <- myResultsSXYZ$sims.list$z
  ZZMAL2[!myResultsSXYZ$sims.list$age.cat %in% c(2) ]  <- 3 # Considers all individuals that are not 1 and 2 as dead
  ZZMAL2[!myResultsSXYZ$sims.list$sex %in% c(1) ]  <- 3 # Considers all individuals that dont have a sex 0 (MALEs) as dead
  
  
  MAL2_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], length(zones))) # nrow = iterations, ncol = year, dim3 = Zone
  
  for(s in 1:length(zones)){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(ZZMAL2[ite,,t]==1) # Select only the individuals alive (z=1)
        
        if(length(which.alive) == 0) next
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        if(length(which.alive) == 1){
          which.aliveSXY <- data.frame(x = which.aliveSXY[1], y = which.aliveSXY[2])
        }
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        MAL2_trapBuf[ite,t,s] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}
  
  
  ## MALE 3
  
  ZZMAL3 <- myResultsSXYZ$sims.list$z
  ZZMAL3[!myResultsSXYZ$sims.list$age.cat %in% c(3) ]  <- 3 # Considers all individuals that are not 1 and 2 as dead
  ZZMAL3[!myResultsSXYZ$sims.list$sex %in% c(1) ]  <- 3 # Considers all individuals that dont have a sex 0 (MALEs) as dead
  
  
  MAL3_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], length(zones))) # nrow = iterations, ncol = year, dim3 = Zone
  
  for(s in 1:length(zones)){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(ZZMAL3[ite,,t]==1) # Select only the individuals alive (z=1)
        
        if(length(which.alive) == 0) next
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        if(length(which.alive) == 1){
          which.aliveSXY <- data.frame(x = which.aliveSXY[1], y = which.aliveSXY[2])
        }
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        MAL3_trapBuf[ite,t,s] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}
  
  ## MALE 4
  
  ZZMAL4 <- myResultsSXYZ$sims.list$z
  ZZMAL4[!myResultsSXYZ$sims.list$age.cat %in% c(4) ]  <- 3 # Considers all individuals that are not 1 and 2 as dead
  ZZMAL4[!myResultsSXYZ$sims.list$sex %in% c(1) ]  <- 3 # Considers all individuals that dont have a sex 0 (MALEs) as dead
  
  
  MAL4_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], length(zones))) # nrow = iterations, ncol = year, dim3 = Zone
  
  for(s in 1:length(zones)){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(ZZMAL4[ite,,t]==1) # Select only the individuals alive (z=1)
        
        if(length(which.alive) == 0) next
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        if(length(which.alive) == 1){
          which.aliveSXY <- data.frame(x = which.aliveSXY[1], y = which.aliveSXY[2])
        }
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        MAL4_trapBuf[ite,t,s] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}
  
  ## MALE 5
  
  ZZMAL5 <- myResultsSXYZ$sims.list$z
  ZZMAL5[!myResultsSXYZ$sims.list$age.cat %in% c(5) ]  <- 3 # Considers all individuals that are not 1 and 2 as dead
  ZZMAL5[!myResultsSXYZ$sims.list$sex %in% c(1) ]  <- 3 # Considers all individuals that dont have a sex 0 (MALEs) as dead
  
  
  MAL5_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], length(zones))) # nrow = iterations, ncol = year, dim3 = Zone
  
  for(s in 1:length(zones)){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(ZZMAL5[ite,,t]==1) # Select only the individuals alive (z=1)
        
        if(length(which.alive) == 0) next
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        if(length(which.alive) == 1){
          which.aliveSXY <- data.frame(x = which.aliveSXY[1], y = which.aliveSXY[2])
        }
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        MAL5_trapBuf[ite,t,s] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}
  
  results <- list(FEM1 = FEM1_trapBuf, FEM2 = FEM2_trapBuf, FEM3 = FEM3_trapBuf, FEM4 = FEM4_trapBuf, FEM5 = FEM5_trapBuf,
                  MAL1 = MAL1_trapBuf, MAL2 = MAL2_trapBuf, MAL3 = MAL3_trapBuf, MAL4 = MAL4_trapBuf, MAL5 = MAL5_trapBuf) 
  return(results)
  
}
ageSt_bear_full <- age_st_full_sex(nucleus = nuc2, periphery = Xbuf2_peri2)

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams")
save(ageSt_bear_full, file = "NFullbuffer.RData")

## ---- 3. Plots ----

rm(list = ls())

library(nimbleSCR)
library(nimble)
library(rgdal)
library(R.utils)


#source("D:/MargSalas/Scripts_MS/Functions/plot.violins3.r")
source("D:/MargSalas/Scripts_MS/Functions/DoScale.r")
source("D:/MargSalas/Scripts_MS/Functions/PlotViolinsHoriz.r")

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams")
load("NFullbuffer.RData")

lim.mal <- max(unlist(lapply(ageSt_bear_full[6:10], max, na.rm = TRUE)))
lim.fem <- max(unlist(lapply(ageSt_bear_full[1:5], max, na.rm = TRUE)))

for(n in 1:length(ageSt_bear_full)){
  if(all(complete.cases(ageSt_bear_full[[n]])) == FALSE){ # delete rows with NA if there are
    delete <- unique(c(which(!complete.cases(ageSt_bear_full[[n]][,,1])), which(!complete.cases(ageSt_bear_full[[n]][,,2]))))
    ageSt_bear_full[[n]] <- ageSt_bear_full[[n]][-delete,,]
  }
}

ageSt_bear_full_F <- ageSt_bear_full[1:5]
ageSt_bear_full_M <- ageSt_bear_full[6:10]

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Plots/model3.1")
pdf("Core_Periphery_AgeFullSex.pdf", 5, 5)

par(mfrow = c(3,2)
    ,mar = c(1,1,1.3,1),
    oma = c(2,3,2,2))

#abline(h=c(0:15))
#colZn4 <- c("#9970ab", "#a6dba0")

for (t in 1:5){
  
  plot(-1000,xlim=c(-lim.mal,lim.fem),ylim=c(0,(length(ageSt_bear_full)/2)+1))
  
  abline(v=0)
  
for(i in 1:length(ageSt_bear_full_F)){
  
  
  plot.violinsHoriz(dat.list = list(ageSt_bear_full_F[[i]][,t,1]) ,x = i,at = i-0.2,
                    
                    add = T,horizontal = T,col="#9970ab",border.col = "#9970ab", alpha = 0.5)
  
  plot.violinsHoriz(dat.list =list(-ageSt_bear_full_M[[i]][,t,1]) ,x = i,at = i-0.2,
                    
                    add = T,horizontal = T,col="#9970ab",border.col = "#9970ab", alpha = 0.5)
  
  
  plot.violinsHoriz(dat.list =list(ageSt_bear_full_F[[i]][,t,2]) ,x = i,at = i+0.2,
                    
                    add = T,horizontal = T,col="#a6dba0",border.col = "#a6dba0", alpha = 0.5)
  
  
  plot.violinsHoriz(dat.list =list(-ageSt_bear_full_M[[i]][,t,2]) ,x = i,at = i+0.2,
                    
                    add = T,horizontal = T,col="#a6dba0",border.col = "#a6dba0", alpha = 0.5)
  
}
}


dev.off()







NAgeFIn <- rmultinom(100,50,c(0.3,0.40,0.10,0.1,0.1))

NAgeMIn <- rmultinom(100,50,c(0.4,0.30,0.20,0.1,0.0))

NAgeFOut <- rmultinom(100,20,c(0.3,0.40,0.10,0.1,0.1))

NAgeMOut <- rmultinom(100,20,c(0.4,0.30,0.20,0.1,0.0))




plot(-1000,xlim=c(-max(NAgeFIn),max(NAgeMIn)),ylim=c(0,nrow(NAgeFIn)+2))

abline(v=0)

abline(h=c(0:15))



for(i in 1:nrow(NAgeFIn)){
  
  
  
  plot.violinsHoriz(dat.list =list(NAgeFIn[i,]) ,x = i,at = i-0.2,
                    
                    add = T,horizontal = T,col="red",alpha = 0.5)
  
  
  
  plot.violinsHoriz(dat.list =list(-NAgeMIn[i,]) ,x = i,at = i-0.2,
                    
                    add = T,horizontal = T,col="red",alpha = 0.5)
  
  
  
  plot.violinsHoriz(dat.list =list(NAgeFOut[i,]) ,x = i,at = i+0.2,
                    
                    add = T,horizontal = T,col="blue",alpha = 0.5)
  
  
  
  plot.violinsHoriz(dat.list =list(-NAgeMOut[i,]) ,x = i,at = i+0.2,
                    
                    add = T,horizontal = T,col="blue",alpha = 0.5)
  
  
}

