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

#   - For sc0, sc1, sc2
M.new <- 700 # New augmentation limit future prediction
t.new <- 5 # Extra years future prediction

#   - For sc3
t.sc2 <- 3 # For sc3, the first three years (2022,2023,2024) are the projections from sc2
t.new2 <- 2 # Extra years sc3
M.aug2 <- 700 # Augmented individuals from first projection model
M.new2 <- 1000 # All augmented individuals projection model sc3 (700 previous + 300 new)

# Load results from predictions


setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL/Predictions/ALLiter_allscin")

load("proj_pcr.all.fem.sc0.RData")
load("proj_pcr.all.fem.sc1.RData")
load("proj_pcr.all.fem.sc2.RData")
load("proj_pcr.all.fem.sc3.RData")
load("proj_pcr.all.fem.sc4.1.RData")
load("proj_pcr.all.fem.sc4.2.RData")


# Load buffer core area and habitat grid (to subset in sampling area)

setwd("D:/MargSalas/Oso/Datos/GIS/Countries")
#Xbuf2 <- readOGR("Buffer_8500_traps.shp")
Xbuf2 <- readOGR("Buffer_8500_traps_sxyObs.shp") # This sampling buffer includes AC of observed individuals a bit outside the trapping array


setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")
load("habcoord.RData")

## ---- 1. Join all results (past+future) in a single matrix with simlist format ----

# I add another dimension for the number of predictions
ndim <- 6 # Sc0, Sc1, Sc2, Sc3, SC4.1, SC4.2

sxy.allyears <- array(NA, c(dim(z.proj.all.sc0)[1], M.new2, 2, Tt + t.new, ndim)) # Dataframes all years to store
z.allyears <- age.cat.allyears <- array(NA, c(dim(z.proj.all.sc0)[1], M.new2, Tt + t.new, ndim))

# Fill in data from past years
# index of iterations taken are in "itera" object
setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL/Predictions/ALLiter_allscin") 
load("itera.RData") # Itera taken from 20iter_sc3in

for(n in 1:ndim){ # Same for both dimensions of predictions (past years remain the same)
  for(t in 1:Tt){
    for(ite in 1:dim(z.proj.all.sc0)[1]){
      # 300 augmented individuals. Store in first five years
      sxy.allyears[ite,1:M.aug,,t,n] <- myResultsSXYZ$sims.list$sxy[itera[ite],,,t] # Only 300 augmented individuals from past
      z.allyears[ite,1:M.aug,t,n] <- myResultsSXYZ$sims.list$z[itera[ite],,t]
      age.cat.allyears[ite,1:M.aug,t,n] <- myResultsSXYZ$sims.list$age.cat[itera[ite],,t]
    }}}

# Fill in data of future years

## PREDICTION 1: SC0 (Dimension 1) -> Normal projection
for(t in 1:t.new){ # 
  for(ite in 1:dim(z.proj.all.sc0)[1]){
    # 700 augmented individuals. We take year 5 from results estimation model. Store in 5 last years
    sxy.allyears[ite,1:M.aug2,,(t+Tt),1] <- sxy.proj.all.sc0[ite,,,(t+1)] 
    z.allyears[ite,1:M.aug2,(t+Tt),1] <- z.proj.all.sc0[ite,,(t+1)]
    age.cat.allyears[ite,1:M.aug2,(t+Tt),1] <- age.cat.proj.all.sc0[ite,,(t+1)]
  }}

## PREDICTION 2: SC1 (Dimension 2) -> Death of 5% of the females in year 2021 INSIDE buffer
for(t in 1:t.new){ # 
  for(ite in 1:dim(z.proj.all.sc0)[1]){
    sxy.allyears[ite,1:M.aug2,,(t+Tt),2] <- sxy.proj.all.sc1[ite,,,(t+1)] 
    z.allyears[ite,1:M.aug2,(t+Tt),2] <- z.proj.all.sc1[ite,,(t+1)]
    age.cat.allyears[ite,1:M.aug2,(t+Tt),2] <- age.cat.proj.all.sc1[ite,,(t+1)]
  }}

## PREDICTION 3: SC2 (Dimension 3) <- Death of 10% of the females in year 2021 INSIDE buffer
for(t in 1:t.new){ # 
  for(ite in 1:dim(z.proj.all.sc0)[1]){
    sxy.allyears[ite,1:M.aug2,,(t+Tt),3] <- sxy.proj.all.sc2[ite,,,(t+1)] 
    z.allyears[ite,1:M.aug2,(t+Tt),3] <- z.proj.all.sc2[ite,,(t+1)]
    age.cat.allyears[ite,1:M.aug2,(t+Tt),3] <- age.cat.proj.all.sc2[ite,,(t+1)]
  }}

## PREDICTION 4: SC3 (Dimension 4) <- Death of 20% of the females in year 2021 INSIDE BUFFER
for(t in 1:t.new){ # 
  for(ite in 1:dim(z.proj.all.sc0)[1]){
    # The projection is stored from year 2022 (t+1 in sxy.proj.all.sc2) -> 2021 is stored from data model
    sxy.allyears[ite,1:M.aug2,,(t+Tt),4] <- sxy.proj.all.sc3[ite,,,(t+1)] 
    z.allyears[ite,1:M.aug2,(t+Tt),4] <- z.proj.all.sc3[ite,,(t+1)]
    age.cat.allyears[ite,1:M.aug2,(t+Tt),4] <- age.cat.proj.all.sc3[ite,,(t+1)]
  }}

## PREDICTION 5: SC4.1 (Dimension 5) <- Death of 10% of the females in year 2021 (sc2) + reintroduction of 5 females in 2024

# This prediction is composed of 2 scenarios:
# 1. Sc2: 10 % of females removed in year 2021. 

for(t in 1:t.sc2){ # Years 2022, 2023, 2024 from sc2 
  for(ite in 1:dim(z.proj.all.sc0)[1]){
    sxy.allyears[ite,1:M.aug2,,(t+Tt),5] <- sxy.proj.all.sc2[ite,,,(t+1)] 
    z.allyears[ite,1:M.aug2,(t+Tt),5] <- z.proj.all.sc2[ite,,(t+1)]
    age.cat.allyears[ite,1:M.aug2,(t+Tt),5] <- age.cat.proj.all.sc2[ite,,(t+1)]
  }}

# 2. Sc4.1: 5 females are re-introduced in year 2024

for(t in 1:t.new2){ # Years 2022, 2023, 2024 from sc2 
  for(ite in 1:dim(z.proj.all.sc0)[1]){
    sxy.allyears[ite,,,(t+Tt+t.sc2),5] <- sxy.proj.all.sc4.1[ite,,,(t+1)] 
    z.allyears[ite,,(t+Tt+t.sc2),5] <- z.proj.all.sc4.1[ite,,(t+1)]
    age.cat.allyears[ite,,(t+Tt+t.sc2),5] <- age.cat.proj.all.sc4.1[ite,,(t+1)]
  }}

## PREDICTION 6: SC4.2 (Dimension 6) <- Death of 20% of the females in year 2021 (sc3) + reintroduction of 5 females in 2024

# This prediction is composed of 2 scenarios:
# 1. Sc3: 20 % of females removed in year 2021. 

for(t in 1:t.sc2){ # Years 2022, 2023, 2024 from sc2 
  for(ite in 1:dim(z.proj.all.sc0)[1]){
    sxy.allyears[ite,1:M.aug2,,(t+Tt),6] <- sxy.proj.all.sc3[ite,,,(t+1)] 
    z.allyears[ite,1:M.aug2,(t+Tt),6] <- z.proj.all.sc3[ite,,(t+1)]
    age.cat.allyears[ite,1:M.aug2,(t+Tt),6] <- age.cat.proj.all.sc3[ite,,(t+1)]
  }}

# 2. Sc4.2: 5 females are re-introduced in year 2024

for(t in 1:t.new2){ # Years 2022, 2023, 2024 from sc2 
  for(ite in 1:dim(z.proj.all.sc0)[1]){
    sxy.allyears[ite,,,(t+Tt+t.sc2),6] <- sxy.proj.all.sc4.2[ite,,,(t+1)] 
    z.allyears[ite,,(t+Tt+t.sc2),6] <- z.proj.all.sc4.2[ite,,(t+1)]
    age.cat.allyears[ite,,(t+Tt+t.sc2),6] <- age.cat.proj.all.sc4.2[ite,,(t+1)]
  }}



## ---- 2. Subset abundance in buffer ----

# Unscale sxy

dimnames(sxy.allyears)[[3]] <- c('x','y')
sxy.allyears.uns <- scaleCoordsToHabitatGrid(coordsData = sxy.allyears,## this are your sxy
                                             coordsHabitatGridCenter = G,# this is your unscaled habitat (as you used when scaling the habitat/detector to the habitat. G?
                                             scaleToGrid = FALSE)$coordsDataScaled

# Array to store abundance in the buffer each iteration and year

#NIn <- matrix(NA,nrow = dim(z.allyears)[1], ncol=dim(z.allyears)[3]) # nrow = iterations, ncol = year
NIn <- array(NA, c(dim(z.allyears)[1], dim(z.allyears)[3], ndim))

for(n in 1:ndim){
  for(ite in 1:dim(z.allyears)[1]){
    for(t in 1:dim(z.allyears)[3]){
      
      which.alive <- which(z.allyears[ite,,t,n]==1) # Select only the individuals alive (z=1)
      
      which.aliveSXY <- sxy.allyears.uns[ite,which.alive,,t,n] # Retrieve the activity center for those individuals
      
      sp <- SpatialPoints(which.aliveSXY, proj4string=CRS(proj4string(Xbuf2))) # CONVERT SXY TO SPATIAL POINTS 
      
      which.In <- over(sp, Xbuf2) # Check which ones are in the buffer
      
      NIn[ite,t,n] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
      
    }
  }
}

#average number of individuals without the buffer for each year 
colMeans(NIn[,,1])
colMeans(NIn[,,2])
colMeans(NIn[,,3])
colMeans(NIn[,,4])
colMeans(NIn[,,5])
colMeans(NIn[,,6])



# Sum of individuals alive in total each year (without buffer)
colMeans(apply(z.allyears[,,,1],c(1,3),function(x) sum(x==1, na.rm = TRUE)))
colMeans(apply(z.allyears[,,,2],c(1,3),function(x) sum(x==1, na.rm = TRUE)))
colMeans(apply(z.allyears[,,,3],c(1,3),function(x) sum(x==1, na.rm = TRUE)))
colMeans(apply(z.allyears[,,,4],c(1,3),function(x) sum(x==1, na.rm = TRUE)))
colMeans(apply(z.allyears[,,,5],c(1,3),function(x) sum(x==1, na.rm = TRUE)))
colMeans(apply(z.allyears[,,,6],c(1,3),function(x) sum(x==1, na.rm = TRUE)))


## ---- Plot ----


setwd("D:/MargSalas/Oso/OPSCR_project/Results/Plots/model3.1/Predictions")
pdf("prediction_abundance_pcr.all.fem_DemScenALLIter_allscin.pdf", 7, 5)

# Adjust according to the number of scenarios

Nscenarios <- dim(NIn)[3]
offset <- seq(0,2,length.out = Nscenarios)
at.sc0 <- c(1,2,3,4,5,6,9,12,15,18)


plot(1, ylim = c(-10,max(NIn)+50), 
     xlim = c(0.5, at.sc0[10] + max(offset) + 0.5) , 
     type ="n", 
     #yaxt="n", 
     xaxt="n", 
     xlab = " ", ylab = "", main = "Abundance",
     cex.axis = 0.8)

axis(1, c(1:ncol(NIn)), labels = c(2017:2026), 
     at = c(1,2,3,4,5, at.sc0[6:10] + max(offset)/2 ),las = 2, cex.axis = 1)

# First five years (present)
for (i in 1:5){
  plot.violins3(list(NIn[ ,i,1]),
                x = i,
                at = i,
                violin.width = 0.3,
                plot.ci = 0.95,
                col = c("wheat4"),
                add = T,
                alpha = 0.8,
                scale.width = FALSE,
                border.col = "wheat4",
                horizontal = FALSE)}

# Scenarios
colSc <- c("yellow4", "yellow3","lightgoldenrod3","lightgoldenrod1", "khaki", "antiquewhite")
for(s in 1:Nscenarios){
  for (i in 6:10){
    plot.violins3(list(NIn[ ,i,s]),
                  x = i,
                  at = at.sc0[i]+offset[s],
                  violin.width = 0.2,
                  plot.ci = 0.95,
                  col = colSc[s],
                  add = T,
                  alpha = 0.8,
                  scale.width = FALSE,
                  border.col = colSc[s],
                  horizontal = FALSE)
  }
}

legend("topleft", inset=c(0,0), legend = c("Sc0(Normal)", "Sc1(-5%F in 2021)", "Sc2(-10%F in 2021)", "Sc3(-20%F in 2021)", "Sc4.1 (-10%F in 2021 + 5ind in 2024)", "Sc4.2 (-20%F in 2021 + 5ind in 2024)" ),
       fill = colSc, border = NA, cex = 0.8)

dev.off()


