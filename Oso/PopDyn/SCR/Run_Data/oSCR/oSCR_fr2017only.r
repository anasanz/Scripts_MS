## -------------------------------------------------
##                          oSCR
##                  French 2017 - 2019
## ------------------------------------------------- 

rm(list = ls())

library(tidyverse)
library(sf)
library(oSCR)
library(raster)

setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Functions")
source("extract.rast.multip.r")


setwd("D:/MargSalas/Oso/Datos/Effort/France")
load("edf2017_2019_fr.RData")
load("tdf2017_fr.RData")
load("tdf2018_fr.RData")
load("tdf2019_fr.RData")


# We found in a preliminary analysis of the whole data set that sigma was very high! 
# exp(9.68) = 16058.6
rbs.sf$mmdm

# We remove Nere and Goiat, which are two individuals moving a lot that don't represent
# the rest of the population.
edf <- edf[-which(edf$ind %in% c("Néré", "Goiat")), ]
edf <- edf[which(edf$session == 1), ] # Keep only detections of 2017 (year 1)

## ---- Prepare format the data for analysis using the data2oscr() ----

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
                      ntraps = c(654),
                      #no. traps vector
                      trapcov.names = c("effort"), #covariate names 
                      tdf.sep = "/")  #char used separator
#sex.nacode = "U")   #unknown sex code

# The one we need for analysis using oSCR is the scrFrame
rbs.sf <- rbs.data$scrFrame

## ---- Create state space ----

# Parameters of Maelis: 25 km and 5 km2 resolution
rbs.ss <- make.ssDF(scrFrame = rbs.sf, buffer = 25000, res = 5000)

## ---- Extract raster covariate in density ----

setwd("D:/MargSalas/Oso/Datos/Effort/France")
#save(rbs.ss_cov, file = "rbs.ss_cov.RData")
load("rbs.ss_cov.RData")

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
  
  m1_no_Nere_Goiat <- oSCR.fit(model = list(D ~ 1,    #density
                                            p0 ~ 1,   #detection 
                                            sig ~ 1), #space use 
                               scrFrame = rbs.sf, ssDF = rbs.ss,
                               trimS = locev) 
  
  setwd("D:/MargSalas/Oso/SCR/oSCR")
  load("m1_noNG.RData")
  exp(8.801) # 6640.882 More reasonable although not perfect
  rbs.sf$mmdm
  
  
  m2 <- oSCR.fit(model = list(D ~ forest_sc,    #density
                                            p0 ~ 1,   #detection 
                                            sig ~ 1), #space use 
                               scrFrame = rbs.sf, ssDF = rbs.ss_cov,
                               trimS = locev)
  
  m2
  
  
  m3_no_NG <- oSCR.fit(model = list(D ~ session,    #density
                                    p0 ~ 1,   #detection 
                                    sig ~ 1), #space use 
                       scrFrame = rbs.sf, ssDF = rbs.ss_cov,
                       trimS = locev)
  
  m4_no_NG <- oSCR.fit(model = list(D ~ session,    #density
                                    p0 ~ effort + session,   #detection 
                                    sig ~ 1), #space use 
                       scrFrame = rbs.sf, ssDF = rbs.ss_cov,
                       trimS = locev)
  
  m5_no_NG <- oSCR.fit(model = list(D ~ session,    #density
                                    p0 ~ effort + session,   #detection 
                                    sig ~ sex), #space use 
                       scrFrame = rbs.sf, ssDF = rbs.ss_cov,
                       trimS = locev)
  
  m6_no_NG <- oSCR.fit(model = list(D ~ session,    #density
                                    p0 ~ effort + session + sex,   #detection 
                                    sig ~ sex), #space use 
                       scrFrame = rbs.sf, ssDF = rbs.ss_cov,
                       trimS = locev)
  
  m7_no_NG_fullmae <- oSCR.fit(model = list(D ~ session,    #density
                                            p0 ~ effort + session + sex,   #detection 
                                            sig ~ session + sex), #space use 
                               scrFrame = rbs.sf, ssDF = rbs.ss,
                               trimS = locev)
}

setwd("D:/MargSalas/Oso/SCR/oSCR/1.no_NG")
save(m7_no_NG_fullmae, file = "m7_no_NG_fullmae")

# Load models

path <- "D:/MargSalas/Oso/SCR/oSCR/1.no_NG"
files <- list.files(path=path)
setwd(path)
results <- sapply(files, function(x) mget(load(x)), simplify = TRUE) 
fl <- fitList.oSCR(results, drop.asu=T, rename=TRUE) #rename=T adds sesnible model names

# Results for now:
#     - No difference per session in density (?) Intercept significant (2017): In m3 (without) and m4 (with effort)

results[[7]]


