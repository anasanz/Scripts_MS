
## -------------------------------------------------
##        Explore individual capture histories
## ------------------------------------------------- 

rm(list = ls())

library(dplyr)
library(tidyr)

setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Functions")
source("capt_hist_bear.r")

# Load monitoring data
setwd("D:/MargSalas/Oso/Datos/Tablas_finales")
os <- read.csv("Data_os_96_20.csv", header = TRUE, row.names = NULL)
os <- os[,-1] 

os_id <- os[which(os$Confirmed_Individual != "Indetermined"), ] # Only identified
os_id$Method[which(os_id$Method == "sampling_station")] <- "Sampling_station"
os_id$Obs_type[which(os_id$Obs_type == "hair")] <- "Hair"

# Explore number of observations per method
os_nobs <- data.frame(os_id %>% group_by(Method, Obs_type) %>% summarise(n()))
sum_method <- t(spread(os_nobs, Obs_type, n..))
sum_method[is.na(sum_method)] <- 0

sum_method2 <- data.frame(matrix(data = as.numeric(sum_method[c(3:nrow(sum_method)), ]), ncol = ncol(sum_method), nrow = nrow(sum_method)-2))
colnames(sum_method2) <- sum_method[1, ]
rownames(sum_method2) <- rownames(sum_method)[-c(1,2)]
setwd("D:/MargSalas/Oso/Exp_analysis")
#write.csv(sum_method2, file = "sum_method.csv")

# Transects? are suposed to be only in France. Mostly data on hair
os_transect <- os_id[which(os_id$Method == "Transect" & os_id$Obs_type == "Hair"), ]
os_transect[which(os_transect$Country == "Spain"), ]
nrow(os_transect[which(os_transect$Country == "France"), ])

os_nobs_transect <- data.frame(os_id %>% group_by(Obs_type, Country, Year) %>% summarise(n()))
t(spread(os_nobs_transect, Obs_type, n..))

## -------------------------------------------------
##        GENERATE CAPTURE HISTORIES PER METHOD 
##  Relevant methods: 
##  Sampling station -> Photo/Video + Hair
##  Transect -> Only French, Hair es lo que más tiene
##  Radiotracking -> VHF + GPS
## ------------------------------------------------- 

CamSampStCH <- capt_hist_bear(data = os_id, 
                   method = "Sampling_station", 
                   obs_type = c("Photo", "Photo/Video", "Video"))

HairSampStCH <- capt_hist_bear(data = os_id, 
                              method = "Sampling_station", 
                              obs_type = c("Hair"))

HairTransCH <- capt_hist_bear(data = os_id, 
                         method = "Transect", 
                         obs_type = c("Hair"))

GPSCH <- capt_hist_bear(data = os_id, 
                         method = "Radiotracking", 
                         obs_type = c("GPS", "VHF"))

## ---- Number of individuals identified the last years ----

os_id_years <- os_id[which(os_id$Year %in% c(2015:2020)), ]

CamSampStCH <- capt_hist_bear(data = os_id_years, 
                              method = "Sampling_station", 
                              obs_type = c("Photo", "Photo/Video", "Video"))

HairSampStCH <- capt_hist_bear(data = os_id_years, 
                               method = "Sampling_station", 
                               obs_type = c("Hair"))

HairTransCH <- capt_hist_bear(data = os_id_years, 
                              method = "Transect", 
                              obs_type = c("Hair"))

