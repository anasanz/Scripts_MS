
## -------------------------------------------------
##                 Plot predictions
## ------------------------------------------------- 

rm(list = ls())

library(nimbleSCR)
library(nimble)
library(rgdal)

source("D:/MargSalas/Scripts_MS/Functions/plot.violins3.r")
source("D:/MargSalas/Scripts_MS/Functions/DoScale.r")


# Load results from model (to get 5th year)

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL")
load("myResults_3-3.1_param.RData")
sampmat1 <- do.call(rbind, nimOutput)

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL")
load("myResults_3-3.1_sxy.RData")
sampmat2 <- do.call(rbind, nimOutputSXY)

sampmat <- cbind(sampmat1, sampmat2)

# Constants:
M.aug <- 300 # Augmented individuals estimation model
Tt <- 5 # Nyears estimation model (2017:2021)

M.new <- 700 # New augmentation limit future prediction
t.new <- 5 # Extra years future prediction

# Load results from predictions

# ** IMPORTANT
# - The files proj_pcr.core and proj_pcr.all come from the projections of 5000 random iterations (vector "itera")
#   with a PCR estimated for all individuals (per capita recruitment per males and females)
# - The files proj_pcr.core.fem and proj_pcr.all.fem come from the projections of 5000 random iterations (vector "itera")
#   with a PCR estimated for ONLY FEMALES (per capita recruitment per female)


setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL/Predictions/ALLiter")
#load("proj_pcr.core1-5000.RData")
#z.proj.core <- z.proj.core2
#sxy.proj.core <- sxy.proj.core2
#age.cat.proj.core <- age.cat.proj.core2
load("proj_pcr.core.RData")
### I use for the final predictions the abundance projections from the pcr estimated as per capita recruitment / FEMALE


#load("proj_pcr.all1-5000.RData")
#z.proj.all <- z.proj.all2
#sxy.proj.all <- sxy.proj.all2
#age.cat.proj.all <- age.cat.proj.all2
load("proj_pcr.all.RData")


# Load buffer core area and habitat grid (to subset in sampling area)

setwd("D:/MargSalas/Oso/Datos/GIS/Countries")
Xbuf2 <- readOGR("Buffer_8500_traps.shp")

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")
load("habcoord.RData")


## ---- 1. Join all results (past+future) in a single matrix with simlist format ----

# I add another dimension for the number of predictions
ndim <- 2 # 1 with pcr.core, 2 with pcr.all

sxy.allyears <- array(NA, c(dim(z.proj.core)[1], M.new, 2, Tt + t.new, ndim)) # Dataframes all years to store
z.allyears <- age.cat.allyears <- array(NA, c(dim(z.proj.core)[1], M.new, Tt + t.new, ndim))

# Different number of iterations. Predictions only from iter 1 to 5000, but it is in order.
# So the first 5000 values of the matrix of past years are equivalent to the 5000 iterations of the projection.

# Fill in data from past years

for(n in 1:ndim){ # Same for both dimensions of predictions (past years remain the same)
  for(t in 1:Tt){
    for(ite in 1:dim(z.proj.core)[1]){
      # 300 augmented individuals. Store in first five years
      sxy.allyears[ite,1:M.aug,,t,n] <- myResultsSXYZ$sims.list$sxy[ite,,,t] # Only 300 augmented individuals from past
      z.allyears[ite,1:M.aug,t,n] <- myResultsSXYZ$sims.list$z[ite,,t]
      age.cat.allyears[ite,1:M.aug,t,n] <- myResultsSXYZ$sims.list$age.cat[ite,,t]
    }}}

# Fill in data of future years
## PREDICTION 1: PCR ESTIMATED IN THE CORE (Dimension 1)
for(t in 1:t.new){ # 
  for(ite in 1:dim(z.proj.core)[1]){
    # 700 augmented individuals. We take year 5 from results estimation model. Store in 5 last years
    sxy.allyears[ite,,,(t+Tt),1] <- sxy.proj.core[ite,,,(t+1)] 
    z.allyears[ite,,(t+Tt),1] <- z.proj.core[ite,,(t+1)]
    age.cat.allyears[ite,,(t+Tt),1] <- age.cat.proj.core[ite,,(t+1)]
  }}

## PREDICTION 2: PCR ESTIMATED IN ALL STATE SPACE (Dimension 2)
for(t in 1:t.new){ # 
  for(ite in 1:dim(z.proj.all)[1]){
    # 700 augmented individuals. We take year 5 from results estimation model. Store in 5 last years
    sxy.allyears[ite,,,(t+Tt),2] <- sxy.proj.all[ite,,,(t+1)] 
    z.allyears[ite,,(t+Tt),2] <- z.proj.all[ite,,(t+1)]
    age.cat.allyears[ite,,(t+Tt),2] <- age.cat.proj.all[ite,,(t+1)]
  }}

## ---- 2. Age categories ----
# Do this for pcr estimated within core buffer area (Prediction 1)

# Unscale sxy

dimnames(sxy.allyears)[[3]] <- c('x','y')
sxy.allyears.uns <- scaleCoordsToHabitatGrid(coordsData = sxy.allyears,## this are your sxy
                                             coordsHabitatGridCenter = G,# this is your unscaled habitat (as you used when scaling the habitat/detector to the habitat. G?
                                             scaleToGrid = FALSE)$coordsDataScaled

## ---- 2.1. CUBS ----
ZZcubs <- z.allyears[,,,1]
ZZcubs[!age.cat.allyears[,,,1] %in% c(1,2) ]  <- 3 # Considers all individuals that are not 1 and 2 as dead

cub_trapBuf <- matrix(NA,nrow=dim(ZZcubs)[1],ncol=dim(ZZcubs)[3]) # nrow = iterations, ncol = year
cub_allss <- matrix(NA,nrow=dim(ZZcubs)[1],ncol=dim(ZZcubs)[3])

for(ite in 1:dim(ZZcubs)[1]){
  for(t in 1:dim(ZZcubs)[3]){
    
    which.alive <- which(ZZcubs[ite,,t]==1) # Select only the individuals alive (z=1)
    #if (length(which.alive) == 0) next # ite 348 year 10 says that there are no cubs alive!
    if(ite == 4141 & t == 10) next
    
    which.aliveSXY <- sxy.allyears.uns[ite,which.alive,,t,1] # Retrieve the activity center for those individuals
    
    sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(Xbuf2))) # CONVERT SXY TO SPATIAL POINTS 
    which.In <- over(sp, Xbuf2) # Check which ones are in the buffer
    
    cub_trapBuf[ite,t] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
    cub_allss[ite,t] <- length(which.alive) # The sum of the points alive is the abundance that year and iteration in the state space. Store
    
  }}


## ---- 2.2. SUBADULTS ----

ZZsub <- z.allyears[,,,1]
ZZsub[!age.cat.allyears[,,,1] %in% c(3,4) ]  <- 3 # Considers all individuals that are not 3 and 4 (SUBADULTS) as dead

subad_trapBuf <- matrix(NA,nrow=dim(ZZsub)[1],ncol=dim(ZZsub)[3]) # nrow = iterations, ncol = year
subad_allss <- matrix(NA,nrow=dim(ZZcubs)[1],ncol=dim(ZZcubs)[3])

for(ite in 1:dim(ZZsub)[1]){
  for(t in 1:dim(ZZsub)[3]){
    
    which.alive <- which(ZZsub[ite,,t]==1) # Select only the individuals alive (z=1)
    
    which.aliveSXY <- sxy.allyears.uns[ite,which.alive,,t,1] # Retrieve the activity center for those individuals
    
    sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(Xbuf2))) # CONVERT SXY TO SPATIAL POINTS 
    which.In <- over(sp, Xbuf2) # Check which ones are in the buffer
    
    subad_trapBuf[ite,t] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
    subad_allss[ite,t] <- length(which.alive)
    }}


## ---- 2.3. ADULTS ----

ZZad <- z.allyears[,,,1]
ZZad[!age.cat.allyears[,,,1] %in% c(5) ]  <- 3 # Considers all individuals that are not 3 and 4 (SUBADULTS) as dead

ad_trapBuf <- matrix(NA,nrow=dim(ZZad)[1],ncol=dim(ZZad)[3]) # nrow = iterations, ncol = year
ad_allss <- matrix(NA,nrow=dim(ZZad)[1],ncol=dim(ZZad)[3])

for(ite in 1:dim(ZZad)[1]){
  for(t in 1:dim(ZZad)[3]){
    
    which.alive <- which(ZZad[ite,,t]==1) # Select only the individuals alive (z=1)
    
    which.aliveSXY <- sxy.allyears.uns[ite,which.alive,,t,1] # Retrieve the activity center for those individuals
    
    sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(Xbuf2))) # CONVERT SXY TO SPATIAL POINTS 
    which.In <- over(sp, Xbuf2) # Check which ones are in the buffer
    
    ad_trapBuf[ite,t] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
    ad_allss[ite,t] <- length(which.alive)
    }}

# Check NA (there must be a couple)
cub_trapBuf[which(!complete.cases(cub_trapBuf)), ] <- cub_trapBuf[50, ] # Just to try
cub_allss[which(!complete.cases(cub_allss)), ] <- cub_allss[50, ] # Just to try

## ---- 3. Estimate abundance proportion of each age class ----

prop_years <- list()
prop <- array(NA,c(dim(ZZad)[1],3,2)) # nrow = iterations, ncol = age categories, ndim = 1: n subset in buffer/ 2: n all ss

for(t in 1:dim(ZZad)[3]){
  for(ite in 1:dim(ZZad)[1]){
    prop[ite, 1, 1] <- cub_trapBuf[ite,t]/sum(cub_trapBuf[ite,t], subad_trapBuf[ite,t], ad_trapBuf[ite,t])
    prop[ite, 2, 1] <- subad_trapBuf[ite,t]/sum(cub_trapBuf[ite,t], subad_trapBuf[ite,t], ad_trapBuf[ite,t])
    prop[ite, 3, 1] <- ad_trapBuf[ite,t]/sum(cub_trapBuf[ite,t], subad_trapBuf[ite,t], ad_trapBuf[ite,t])
    
    prop[ite, 1, 2] <- cub_allss[ite,t]/sum(cub_allss[ite,t], subad_allss[ite,t], ad_allss[ite,t])
    prop[ite, 2, 2] <- subad_allss[ite,t]/sum(cub_allss[ite,t], subad_allss[ite,t], ad_allss[ite,t])
    prop[ite, 3, 2] <- ad_allss[ite,t]/sum(cub_allss[ite,t], subad_allss[ite,t], ad_allss[ite,t])
  }
  prop_years[[t]] <- prop
}

df.core <- matrix(NA,nrow = 10, ncol = 3)
for (t in 1:10){
  df.core[t,] <- apply(prop_years[[t]][,,1], c(2), FUN = mean)
}

df.ss <- matrix(NA,nrow = 10, ncol = 3)
for (t in 1:10){
  df.ss[t,] <- apply(prop_years[[t]][,,2], c(2), FUN = mean)
}




## ---- 3.5. Plot ----

## ---- 3.5.1. Inside core area ----

source("D:/PhD/MyScripts_PhD/Ch. 2-3/Ch. 3/Results/Functions/plot.violins3.r")

years <- c(2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026)

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Plots/model3.1/Predictions")
#pdf("age_structure_pred_core_pcr.fem.pdf", 7, 7)
pdf("age_structure_pred_core_pcr.allind2.pdf", 7, 7)

par(mfrow = c(1,2),
    mar = c(3,3,3,2))

# Years estimation

plot(1, ylim = c(0, 19 + 0.5), 
     xlim = c(0,0.8), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Estimation",
     cex.axis = 0.8)

#axis(2, c(1:20), las = 2, cex.axis = 1)
axis(2, labels = years[1:5], lwd.ticks = 0,at = c(2,6,10,14,18), cex.axis = 1)

cubs_at <- c(1,5,9,13,17)
subad_at <- c(2,6,10,14,18)
ad_at <- c(3,7,11,15,19)
ats <- list(cubs_at, subad_at, ad_at)

#polygon(x = c(0,0,1,1), y = c(0,4,4,0), col = adjustcolor("grey", alpha.f = 0.5), border = NA)

color_cat1 <- c("lightgoldenrod1", "darkolivegreen4", "darkslategrey")
color_cat2 <- c("lightgoldenrod1", "darkgoldenrod1", "sienna4") # no mucho
color_cat3 <-c("aquamarine1", "aquamarine4","darkslategrey")
color_cat4 <- c("aquamarine1","lightseagreen", "darkslategrey")
color_cat5 <- c("lightgoldenrod1", "darksalmon", "sienna4") # no mucho

for(i in 1:ncol(prop_years[[1]])){ # Look into cubs
  
  for(t in 1:5){
    
    plot.violins3(list(prop_years[[t]][ ,i,1]),
                  x = i,
                  at = ats[[i]][t],
                  violin.width = 0.4,
                  plot.ci = 0.95,
                  col = color_cat1[i],
                  add = T,
                  alpha = 0.3,
                  scale.width = FALSE,
                  border.col = "black",
                  horizontal = TRUE)}}

# Years prediction

plot(1, ylim = c(0, 19 + 0.5), 
     xlim = c(0,0.8), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Prediction",
     cex.axis = 0.8)

#axis(2, c(1:20), las = 2, cex.axis = 1)
axis(2, labels = years[6:10], lwd.ticks = 0,at = c(2,6,10,14,18), cex.axis = 1)

for(i in 1:ncol(prop_years[[1]])){ # Look into cubs
  
  for(t in 1:5){
    
    plot.violins3(list(prop_years[[t+5]][ ,i,1]),
                  x = i,
                  at = ats[[i]][t],
                  violin.width = 0.4,
                  plot.ci = 0.95,
                  col = color_cat1[i],
                  add = T,
                  alpha = 0.3,
                  scale.width = FALSE,
                  border.col = "black",
                  horizontal = TRUE)}}

dev.off()

## ---- 3.5.1. All state space ----

source("D:/PhD/MyScripts_PhD/Ch. 2-3/Ch. 3/Results/Functions/plot.violins3.r")

years <- c(2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026)

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Plots/model3.1/Predictions")
#pdf("age_structure_pred_allss_pcr.fem.pdf", 7, 7)
pdf("age_structure_pred_allss_pcr.pcr.allind2.pdf", 7, 7)

par(mfrow = c(1,2),
    mar = c(3,3,3,2))

# Years estimation

plot(1, ylim = c(0, 19 + 0.5), 
     xlim = c(0,0.8), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Estimation",
     cex.axis = 0.8)

#axis(2, c(1:20), las = 2, cex.axis = 1)
axis(2, labels = years[1:5], lwd.ticks = 0,at = c(2,6,10,14,18), cex.axis = 1)

cubs_at <- c(1,5,9,13,17)
subad_at <- c(2,6,10,14,18)
ad_at <- c(3,7,11,15,19)
ats <- list(cubs_at, subad_at, ad_at)

#polygon(x = c(0,0,1,1), y = c(0,4,4,0), col = adjustcolor("grey", alpha.f = 0.5), border = NA)

color_cat1 <- c("lightgoldenrod1", "darkolivegreen4", "darkslategrey")
color_cat2 <- c("lightgoldenrod1", "darkgoldenrod1", "sienna4") # no mucho
color_cat3 <-c("aquamarine1", "aquamarine4","darkslategrey")
color_cat4 <- c("aquamarine1","lightseagreen", "darkslategrey")
color_cat5 <- c("lightgoldenrod1", "darksalmon", "sienna4") # no mucho

for(i in 1:ncol(prop_years[[1]])){ # Look into cubs
  
  for(t in 1:5){
    
    plot.violins3(list(prop_years[[t]][ ,i,2]),
                  x = i,
                  at = ats[[i]][t],
                  violin.width = 0.4,
                  plot.ci = 0.95,
                  col = color_cat1[i],
                  add = T,
                  alpha = 0.3,
                  scale.width = FALSE,
                  border.col = "black",
                  horizontal = TRUE)}}

# Years prediction

plot(1, ylim = c(0, 19 + 0.5), 
     xlim = c(0,0.8), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Prediction",
     cex.axis = 0.8)

#axis(2, c(1:20), las = 2, cex.axis = 1)
axis(2, labels = years[6:10], lwd.ticks = 0,at = c(2,6,10,14,18), cex.axis = 1)

for(i in 1:ncol(prop_years[[1]])){ # Look into cubs
  
  for(t in 1:5){
    
    plot.violins3(list(prop_years[[t+5]][ ,i,2]),
                  x = i,
                  at = ats[[i]][t],
                  violin.width = 0.4,
                  plot.ci = 0.95,
                  col = color_cat1[i],
                  add = T,
                  alpha = 0.3,
                  scale.width = FALSE,
                  border.col = "black",
                  horizontal = TRUE)}}

dev.off()


