rm(list=ls())



###################################################################################

### Check if deleting observations the uncertainty is similar. Talking with Cyril, 
# this doesn't make sense because we don't know in reality how individuals are grouped.
# The problem is there is heterogeneity in detection, and some individuals may have higher
# prob. detection if they are in a group. (although this is hard, I dont take that into account
# either in my papers
# it would work in farmdindis though (puting the group
# observations as independent). But it is too much work .

# So this script is incomplete

library(dplyr)
library(tidyr)

# Load analyzed transects (all species)
setwd("D:/Otros/Tórtola/Results/Study2/Model_results_1.1/STTUR")
tr <- read.csv("pInc_bEst_sd_transects_converged.csv")


# 1. Check transects with more detections (taking into account the analyzed years)

setwd("D:/Otros/Tórtola/Data")
m <- read.csv("m_counts_all_transects.csv", sep = ",")
m$tr2 <- substr(m$X, 1, nchar(m$X)-2)
colnames(m) <- c("X", 2002:2021, "tr2")

transect <- tr$transect_ID

yearly_det <- as.data.frame(matrix(NA, nrow = length(transect), ncol = 22))
rownames(yearly_det) <- transect
colnames(yearly_det) <- c(2002:2021,"nyears_trend", "detections")

for (xxx in 1:length(transect)){
  
  data_transect <- m[which(m$tr2 %in% transect[xxx]), ]
  
  # Year: Min and max year (if there are years with NA in between, those are considered within the time series)
  
  data_transect_nona <- data_transect[, colSums(is.na(data_transect)) != nrow(data_transect)] # Delete columns with NA
  
  min_year <- min(as.numeric(colnames(data_transect_nona)[-which(colnames(data_transect_nona) %in% c("X", "tr2"))]))
  max_year <- max(as.numeric(colnames(data_transect_nona)[-which(colnames(data_transect_nona) %in% c("X", "tr2"))]))
  
  calc <- data_transect_nona[,-which(colnames(data_transect_nona) %in% c("X", "tr2"))]
  calc_years <- apply(calc, 2, sum) 
  
  for (x in 1:length(calc_years)){
    yearly_det[which(rownames(yearly_det) %in% transect[xxx]),names(calc_years)[x]] <- calc_years[x]
  }
  
  yearly_det[which(rownames(yearly_det) %in% transect[xxx]),"nyears_trend"] <- max_year - min_year
  yearly_det[which(rownames(yearly_det) %in% transect[xxx]),"detections"] <- sum(calc)
  
}

yearly_det$transect <- rownames(yearly_det)
yearly_det <- arrange(yearly_det, desc(detections))

# Check if it is the transects that were hard to converge
setwd("D:/Otros/Tórtola/Results/Study2/Model_results_1.1/STTUR")
load("1.1TortoData_noconverge_2002_2021.RData")
load("1.1TortoData_noconverge2_9e5iter_2002_2021.RData")
load("1.1TortoData_noconverge3(added)_9e5iter_2002_2021.RData")

hard_conv <- c(no_converge, no_converge2, no_converge3)
yearly_det$transect[1:10] %in% hard_conv # NONE
tr_test <- yearly_det$transect[1:10]

# 2. Check Average cluster size of sttur in farmdindis

setwd("D:/Otros/Adrien/Send data") # From 2010-2020
clusdat <- openxlsx::read.xlsx("Data_DS_10_20.xlsx")
clusdat <- clusdat[which(clusdat$Especie == "STTUR"),]
mean(clusdat$Nombre)

# 3. Divide sample size by average sample size

mat <- matrix(NA, nrow = 1, ncol = ncol(m))

for (i in 1:length(tr_test)){
  
  reduced_sample <- yearly_det$detections[yearly_det$transect %in% tr_test[i]]/mean(clusdat$Nombre) # detections/clust size: sample size to end up
  delete <- round(yearly_det$detections[yearly_det$transect %in% tr_test[i]] - reduced_sample) 
  
  mat_trans1 <- m[which(m$tr2 %in% tr_test[i]),]
  mat_trans1 <- mat_trans1[, -which(colnames(mat_trans1) %in% c("X", "tr2"))]
  mat_trans1 <- as.matrix(mat_trans1)
  mat_trans2 <- mat_trans1
  
  mat_trans = mat_trans1
  
  samp <- function(N) {
    id <- which(mat_trans>1 & !is.na(mat_trans[]))
    id <- sample(id,size = N)
    mat_trans[id] <- mat_trans[id]-1
    return(mat_trans)
  }
  
}


# Delete observations
# Script checking model output, deleting observations and checking again