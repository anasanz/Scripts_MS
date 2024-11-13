rm(list = ls())

library(tidyverse)
library(sf)
library(oSCR)
library(raster)

setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Functions")
#source("extract.rast.multip.r")


setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Oso/PopDyn/SCR/Effort/France/Syst_Opport")
load("edf2017_2019_fr.RData")
load("tdf2017_fr.RData")


# We found in a preliminary analysis of the whole data set that sigma was very high! 
# exp(9.68) = 16058.6

# We remove Nere and Goiat, which are two individuals moving a lot that don't represent
# the rest of the population.

edf <- edf[-which(edf$ind %in% c("Néré", "Goiat")), ]
edf <- edf[which(edf$session == 1), ] # Keep only detections of 2017 (year 1)

## ---- Prepare format the data for analysis using the data2oscr() ----

# Important!! Include "site" as a trap covariate, as it denotes the type of trap:
# Hair snag alone (hair), a camera trap combine with hair snag (both), or an opportunistic trap

rbs.data <- data2oscr(edf = edf,
                      #the EDF
                      sess.col = 1,
                      #session column
                      id.col = 2,
                      #individual column
                      occ.col = 3,
                      #occasion column
                      trap.col = 4, # BOARD
                      #trap column
                      sex.col = 5,
                      #sex column (optional)
                      tdf = list(tdf2017),
                      K = c(7), # Same number of occasions per season
                      #occasion vector
                      ntraps = c(900),
                      #no. traps vector
                      trapcov.names = c("effort","site"), #covariate names 
                      tdf.sep = "/")  #char used separator
#sex.nacode = "U")   #unknown sex code

# The one we need for analysis using oSCR is the scrFrame
rbs.sf <- rbs.data$scrFrame
rbs.sf$caphist
## ---- Create state space ----

# Parameters of Maelis: 25 km and 5 km2 resolution
rbs.ss <- make.ssDF(scrFrame = rbs.sf, buffer = 25000, res = 5000)

## ---- Extract raster covariate in density ----

#setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/Variables_hrscale")
#forest <- raster("forest_hrbear.tif")
#dem <- raster("dem_hrbear.tif")
#slope <- raster("slope_hrbear.tif")
#rough <- raster("rough_hrbear.tif")

#rasters <- list(forest, dem, slope, rough)

#rbs.ss_cov <- extract.rast.multip(rbs.ss, rasters, cov.name = c("forest", "dem", "slope", "rough"), func = mean)

#setwd("~/Scripts_MS/Oso/PopDyn/SCR/Effort/France")
setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Oso/PopDyn/SCR/Effort/France")
#save(rbs.ss_cov, file = "rbs.ss_cov_SystOpport.RData")
load("rbs.ss_cov_SystOpport.RData")

rbs.ss_cov <- rbs.ss_cov[1]

forest_mean <- mean(rbs.ss_cov[[1]]$forest)
forest_sd <- sd(rbs.ss_cov[[1]]$forest)
rbs.ss_cov[[1]]$forest_sc <- (rbs.ss_cov[[1]]$forest - forest_mean) / forest_sd

rbs.ss_cov <- list(rbs.ss_cov[[1]])

## ---- Local evaluation parameter ----

# 1. Backtramsform sigma value obtained by Maelis to do local evaluation -

# sigma <- sqrt(1433.3/3.14159265359)/5.99
# sigma_m <- sigma*1000

# We use 3*sigma for the trim parameter to make sure that it covers all used area
# locev <- sigma_m*3

#2. In our data, to start with: set local evaluation to 60000 so that is very big and doesn't constrain
# locev <- 60000

# 3. After some results; Local evaluation parameter to 3*sigma (3*6640 = ~20000)
locev <- 20000

## ---- Model fitting ----

# M1 # This ran faster! 1h15min
if(1 == 2){ 
  
  m1 <- oSCR.fit(model = list(D ~ 1,    #density
                              p0 ~ 1,   #detection 
                              sig ~ 1), #space use 
                 scrFrame = rbs.sf, ssDF = rbs.ss,
                 trimS = locev) 
  
  m2 <- oSCR.fit(model = list(D ~ 1,    #density
                              p0 ~ effort + site,   #detection 
                              sig ~ 1), #space use 
                 scrFrame = rbs.sf, ssDF = rbs.ss_cov,
                 trimS = locev)
  
  m3 <- oSCR.fit(model = list(D ~ 1,    #density
                              p0 ~ effort + site + sex,   #detection 
                              sig ~ 1), #space use 
                 scrFrame = rbs.sf, ssDF = rbs.ss_cov,
                 trimS = locev)
  
  m4 <- oSCR.fit(model = list(D ~ 1,    #density
                              p0 ~ effort + site,   #detection 
                              sig ~ sex), #space use 
                 scrFrame = rbs.sf, ssDF = rbs.ss_cov,
                 trimS = locev)
  
  m5 <- oSCR.fit(model = list(D ~ forest,    #density
                              p0 ~ effort + site,   #detection 
                              sig ~ sex), #space use 
                 scrFrame = rbs.sf, ssDF = rbs.ss_cov,
                 trimS = locev)
  
  setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Oso/PopDyn/SCR/Data/oSCR/Results/3.SystOpport2017only")
  save(m5, file = "m5.RData")
  
  m5.sc <- oSCR.fit(model = list(D ~ forest_sc,    #density
                              p0 ~ effort + site,   #detection 
                              sig ~ sex), #space use 
                 scrFrame = rbs.sf, ssDF = rbs.ss_cov,
                 trimS = locev)
  
  setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Oso/PopDyn/SCR/Data/oSCR/Results/3.SystOpport2017only")
  save(m5.sc, file = "m5.sc.RData")
  

  load("m5.RData")
  # Density intercept is log(per pixel) so total N = exp(intercept)*[# statespace points]
  Nhat <- sum(exp(m5.sc$outStats[8,"mle"] + m5.sc$outStats[9,"mle"]*rbs.ss_cov[[1]]$forest_sc))
  
  
  m6 <- oSCR.fit(model = list(D ~ dem,    #density
                              p0 ~ effort + site,   #detection 
                              sig ~ sex), #space use 
                 scrFrame = rbs.sf, ssDF = rbs.ss_cov,
                 trimS = locev)
  
  m7 <- oSCR.fit(model = list(D ~ slope,    #density
                              p0 ~ effort + site,   #detection 
                              sig ~ sex), #space use 
                 scrFrame = rbs.sf, ssDF = rbs.ss_cov,
                 trimS = locev)
  
  m8 <- oSCR.fit(model = list(D ~ rough,    #density
                              p0 ~ effort + site,   #detection 
                              sig ~ sex), #space use 
                 scrFrame = rbs.sf, ssDF = rbs.ss_cov,
                 trimS = locev)
  
  setwd("~/Scripts_MS/Oso/PopDyn/SCR/Data/oSCR/Results/2.SystOpport")
  save(m8, file = "m8.RData")
  
  m9 <- oSCR.fit(model = list(D ~ rough + forest,    #density
                              p0 ~ effort + site,   #detection 
                              sig ~ sex), #space use 
                 scrFrame = rbs.sf, ssDF = rbs.ss_cov,
                 trimS = locev)
  
  setwd("~/Scripts_MS/Oso/PopDyn/SCR/Data/oSCR/Results/2.SystOpport")
  save(m9, file = "m9.RData")
  
  m10 <- oSCR.fit(model = list(D ~ forest_sc,    #density
                              p0 ~ 1,   #detection 
                              sig ~ 1), #space use 
                 scrFrame = rbs.sf, ssDF = rbs.ss_cov,
                 trimS = locev)
  
  setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Oso/PopDyn/SCR/Data/oSCR/Results/3.SystOpport2017only")
  save(m10, file = "m10.RData")
  
  Nhat <- sum(exp(m10$outStats[3,"mle"] + m10$outStats[4,"mle"]*rbs.ss_cov[[1]]$forest_sc))
  sig <- exp(m10$outStats[2,"mle"]) # meters
  
}

## ---- Model results ----

path <- "D:/MargSalas/Scripts_MS/Scripts_MS/Oso/PopDyn/SCR/Data/oSCR/Results/3.SystOpport2017only"
files <- list.files(path=path)
setwd(path)
results <- sapply(files, function(x) mget(load(x)), simplify = TRUE) 
fl <- fitList.oSCR(results, drop.asu=T, rename=TRUE) #rename=T adds sesnible model names

# SAME RESULTS AS MAELIS, ROUGHNESS AND FOREST AS BEST COVARIATES

results[[1]]

