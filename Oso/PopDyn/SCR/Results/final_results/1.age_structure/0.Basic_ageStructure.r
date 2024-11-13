
## -------------------------------------------------
##              Age structure basic info
## ------------------------------------------------- 

library(MCMCvis)
library(rgdal)
library(nimbleSCR)

source("D:/MargSalas/Scripts_MS/Functions/plot.violins3.r")
source("D:/MargSalas/Scripts_MS/Functions/DoScale.r")
source("D:/MargSalas/Scripts_MS/Functions/PlotViolinsHoriz.r")

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")
load("habcoord.RData")

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL")
load("myResults_3-3.1_sxy.RData")

# Unscale the sxy coordinates

dimnames(myResultsSXYZ$sims.list$sxy)[[3]] <- c('x','y')
myResultsSXYZ$sims.list$sxy <- scaleCoordsToHabitatGrid(coordsData = myResultsSXYZ$sims.list$sxy,## this are your sxy
                                                        coordsHabitatGridCenter = G,# this is your unscaled habitat (as you used when scaling the habitat/detector to the habitat. G?
                                                        scaleToGrid = FALSE)$coordsDataScaled


setwd("D:/MargSalas/Oso/Datos/GIS/Countries")
Xbuf <- readOGR("Buffer_statespace.shp")
Xbuf2 <- readOGR("Buffer_8500_traps_sxyObs.shp") # This sampling buffer includes AC of observed individuals a bit outside the trapping array

# Matrix to store abundance in the buffer each iteration and year
NIn_trapBuf <- matrix(NA,nrow=dim(myResultsSXYZ$sims.list$z)[1],ncol=dim(myResultsSXYZ$sims.list$z)[3]) # nrow = iterations, ncol = year


for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
  for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
    
    which.alive <- which(myResultsSXYZ$sims.list$z[ite,,t]==1) # Select only the individuals alive (z=1)
    
    which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
    
    sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(Xbuf2))) # CONVERT SXY TO SPATIAL POINTS 
    
    which.In <- over(sp, Xbuf2) # Check which ones are in the buffer
    
    NIn_trapBuf[ite,t] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
  }
}

## ---- 3. AGE STRUCTURE ----
# 1. Estimate abundance in buffer of the traps of each age category
# 2. Estimate proportion of total abundance in the sampling buffer (NIn_trapBuf)

## ---- 3.1. CUBS ----
ZZcubs <- myResultsSXYZ$sims.list$z
ZZcubs[!myResultsSXYZ$sims.list$age.cat %in% c(1,2) ]  <- 3 # Considers all individuals that are not 1 and 2 as dead

cub_trapBuf <- matrix(NA,nrow=dim(myResultsSXYZ$sims.list$z)[1],ncol=dim(myResultsSXYZ$sims.list$z)[3]) # nrow = iterations, ncol = year


for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
  for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
    
    which.alive <- which(ZZcubs[ite,,t]==1) # Select only the individuals alive (z=1)
    
    which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
    
    sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(Xbuf2))) # CONVERT SXY TO SPATIAL POINTS 
    which.In <- over(sp, Xbuf2) # Check which ones are in the buffer
    
    cub_trapBuf[ite,t] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
  }}


## ---- 3.2. SUBADULTS ----

ZZsub <- myResultsSXYZ$sims.list$z
ZZsub[!myResultsSXYZ$sims.list$age.cat %in% c(3,4) ]  <- 3 # Considers all individuals that are not 3 and 4 (SUBADULTS) as dead

subad_trapBuf <- matrix(NA,nrow=dim(myResultsSXYZ$sims.list$z)[1],ncol=dim(myResultsSXYZ$sims.list$z)[3]) # nrow = iterations, ncol = year

for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
  for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
    
    which.alive <- which(ZZsub[ite,,t]==1) # Select only the individuals alive (z=1)
    
    which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
    
    sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(Xbuf2))) # CONVERT SXY TO SPATIAL POINTS 
    which.In <- over(sp, Xbuf2) # Check which ones are in the buffer
    
    subad_trapBuf[ite,t] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
  }}


## ---- 3.3. ADULTS ----

ZZad <- myResultsSXYZ$sims.list$z
ZZad[!myResultsSXYZ$sims.list$age.cat %in% c(5) ]  <- 3 # Considers all individuals that are not 3 and 4 (SUBADULTS) as dead

ad_trapBuf <- matrix(NA,nrow=dim(myResultsSXYZ$sims.list$z)[1],ncol=dim(myResultsSXYZ$sims.list$z)[3]) # nrow = iterations, ncol = year

for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
  for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
    
    which.alive <- which(ZZad[ite,,t]==1) # Select only the individuals alive (z=1)
    
    which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
    
    sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(Xbuf2))) # CONVERT SXY TO SPATIAL POINTS 
    which.In <- over(sp, Xbuf2) # Check which ones are in the buffer
    
    ad_trapBuf[ite,t] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
  }}


## ---- 3.4. Estimate abundance proportion of each age class ----

prop_years <- list()

prop <- matrix(NA,nrow=dim(myResultsSXYZ$sims.list$z)[1],ncol=3) # nrow = iterations, ncol = age categories

ite = 1
t = 1

for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
  for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
    prop[ite, 1] <- cub_trapBuf[ite,t]/NIn_trapBuf[ite,t]
    prop[ite, 2] <- subad_trapBuf[ite,t]/NIn_trapBuf[ite,t]
    prop[ite, 3] <- ad_trapBuf[ite,t]/NIn_trapBuf[ite,t]
  }
  prop_years[[t]] <- prop
}

basic_ageSt <- list(cub_trapBuf, subad_trapBuf, ad_trapBuf, prop_years)

#setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams")
#save(basic_ageSt, file = "Nbuffer_BASIC_newSize.RData")

f <- lapply(prop_years, colMeans)
merged_matrix <- do.call(rbind, f)
colMeans(merged_matrix) # OVERALL AGE STRUCTURE RESULTS

## ---- 3.5. Plot ----

source("D:/PhD/MyScripts_PhD/Ch. 2-3/Ch. 3/Results/Functions/plot.violins3.r")
years <- c(2017, 2018, 2019, 2020, 2021)

setwd("D:/MargSalas/Oso/Results/Plots/model3.1")
pdf("age_structure.pdf", 9, 7)

par(mfrow = c(2,3))

for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
  
  plot(1, ylim = c(0.5, ncol(prop_years[[t]])+0.5), 
       xlim = c(0,1), 
       type ="n", yaxt="n", 
       #xaxt="n", 
       xlab = " ", ylab = "", main = years[t],
       cex.axis = 0.8)
  
  axis(2, c(1:ncol(prop_years[[t]])), labels = c("Cub","Subadult","Adult"), las = 2, cex.axis = 1)
  
  for(i in 1:ncol(prop_years[[t]])){
    plot.violins3(list(prop_years[[t]][ ,i]),
                  x = i,
                  at = i,
                  violin.width = 0.2,
                  plot.ci = 0.95,
                  col = c("darksalmon"),
                  add = T,
                  alpha = 0.3,
                  scale.width = FALSE,
                  border.col = "black",
                  horizontal = TRUE)}}

dev.off()

# Plot bonito

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Results_section/Plots")
pdf("SI_age_structure.pdf", 6, 10)

par(mfrow = c(1,1))

plot(1, ylim = c(0, 19 + 0.5), 
     xlim = c(0,1), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "",
     cex.axis = 0.8)

#axis(2, c(1:20), las = 2, cex.axis = 1)
axis(2, c(1:20), las = 2, cex.axis = 1, labels = c(2017, 2018, 2019, 2020, 2021),
     at = c(2,6,10,14,18))

cubs_at <- c(1,5,9,13,17)
subad_at <- c(2,6,10,14,18)
ad_at <- c(3,7,11,15,19)
ats <- list(cubs_at, subad_at, ad_at)

polygon(x = c(0,0,1,1), y = c(4,8,8,4), col = adjustcolor("grey", alpha.f = 0.5), border = NA)
polygon(x = c(0,0,1,1), y = c(12,16,16,12), col = adjustcolor("grey", alpha.f = 0.5), border = NA)

#color_cat1 <- c("lightgoldenrod1", "darkolivegreen4", "darkslategrey")
#color_cat2 <- c("lightgoldenrod1", "darkgoldenrod1", "sienna4") # no mucho
#color_cat3 <-c("aquamarine1", "aquamarine4","darkslategrey")
color_cat4 <- c("aquamarine1","lightseagreen", "darkslategrey")
#color_cat5 <- c("lightgoldenrod1", "darksalmon", "sienna4") # no mucho

for(i in 1:ncol(prop_years[[1]])){ # Look into cubs
  
  for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
    
    plot.violins3(list(prop_years[[t]][ ,i]),
                  x = i,
                  at = ats[[i]][t],
                  violin.width = 0.4,
                  plot.ci = 0.95,
                  col = color_cat4[i],
                  add = T,
                  alpha = 0.3,
                  scale.width = FALSE,
                  border.col = "black",
                  horizontal = TRUE)
  }}

legend("topright", legend = c("Juveniles", "Subadults", "Adults"), fill = rev(color_cat4), border = NA)

dev.off()


