

## -------------------------------------------------
##     Presence/Absence of columbid species 
## ------------------------------------------------- 


rm(list=ls())

library(rjags)
library(jagsUI)
library(dplyr)
library(stringr)


## -------------------------------------------------
##                 STDEC
## ------------------------------------------------- 

setwd("D:/Otros/Tórtola/Data/Columbid")

tor <- read.csv("stdec_ds_ready_02_21.csv", sep = ",")
tor[,1] <- "STDEC"

yrs <- unique(tor$Year) 
nyrs <- length(yrs)

all.sites <- unique(tor$Site)
max.sites <- length(all.sites)

m <- matrix(NA, nrow = length(all.sites), ncol = nyrs)
rownames(m) <- all.sites
colnames(m) <- yrs

# Add counts > 0
tor_yes <- tor[which(tor$count == 1), ]
count <- aggregate(X ~ Year + Site, FUN = length, data = tor_yes)

for (i in 1:nrow(count)){
  m[which(rownames(m) %in% count$Site[i]), which(colnames(m) %in% count$Year[i])] <- count$X[i]
}
m2 <- m

# Add absences (0)
tor_no <- tor[which(tor$count == 0), ]
tor_no2 <- aggregate(count ~ site_year + Site + Year, FUN = sum, data = tor_no)

for (i in 1:nrow(tor_no2)){
  value.m <- m[which(rownames(m) %in% tor_no2$Site[i]), which(colnames(m) %in% tor_no2$Year[i])]
  if(is.na(value.m)) { # If there is a detection already in any section, is present
  m[which(rownames(m) %in% tor_no2$Site[i]), which(colnames(m) %in% tor_no2$Year[i])] <- tor_no2$count[i]
  }}

# Convert abundance into presence/absence
m2 <- ifelse(m > 0, 1,0)

# Select analyzed transects for STTUR
setwd("D:/Otros/Tórtola/Results/Study2/Model_results_1.1/STTUR")
tr_tor <- read.csv("results_analysis1.csv", sep = ",")

m2_subset <- m2[which(rownames(m2) %in% tr_tor$transect_ID), ]

# setwd("D:/Otros/Tórtola/Data/Columbid")
# write.csv(m2_subset, file = "presenceAbsence_stdec.csv")

## -------------------------------------------------
##                 COPAL
## ------------------------------------------------- 

setwd("D:/Otros/Tórtola/Data/Columbid")

tor <- read.csv("copal_ds_ready_02_21.csv", sep = ",")
tor[,1] <- "STDEC"

yrs <- unique(tor$Year) 
nyrs <- length(yrs)

all.sites <- unique(tor$Site)
max.sites <- length(all.sites)

m <- matrix(NA, nrow = length(all.sites), ncol = nyrs)
rownames(m) <- all.sites
colnames(m) <- yrs

# Add counts > 0
tor_yes <- tor[which(tor$count == 1), ]
count <- aggregate(X ~ Year + Site, FUN = length, data = tor_yes)

for (i in 1:nrow(count)){
  m[which(rownames(m) %in% count$Site[i]), which(colnames(m) %in% count$Year[i])] <- count$X[i]
}
m2 <- m

# Add absences (0)
tor_no <- tor[which(tor$count == 0), ]
tor_no2 <- aggregate(count ~ site_year + Site + Year, FUN = sum, data = tor_no)

for (i in 1:nrow(tor_no2)){
  value.m <- m[which(rownames(m) %in% tor_no2$Site[i]), which(colnames(m) %in% tor_no2$Year[i])]
  if(is.na(value.m)) { # If there is a detection already in any section, is present
    m[which(rownames(m) %in% tor_no2$Site[i]), which(colnames(m) %in% tor_no2$Year[i])] <- tor_no2$count[i]
  }}

# Convert abundance into presence/absence
m2 <- ifelse(m > 0, 1,0)

# Select analyzed transects for STTUR
setwd("D:/Otros/Tórtola/Results/Study2/Model_results_1.1/STTUR")
tr_tor <- read.csv("results_analysis1.csv", sep = ",")

m2_subset_copal <- m2[which(rownames(m2) %in% tr_tor$transect_ID), ]

#setwd("D:/Otros/Tórtola/Data/Columbid")
#write.csv(m2_subset_copal, file = "presenceAbsence_copal.csv")

# Transects not done should be the same
is.na(m2_subset_copal) == is.na(m2_subset)

(tr_tor[,c(5:24)] == 0) == is.na(m2_subset)
