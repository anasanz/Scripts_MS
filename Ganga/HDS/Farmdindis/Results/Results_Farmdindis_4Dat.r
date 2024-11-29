rm(list=ls())

## ---- Results ----

setwd("D:/MargSalas/Ganga/Results/HDS/Farmdindis/Model_results")
load("2.2.Dat_HDS_trendmodel_lam[hq]_sigHR.RData")

out$summary

# Load hq areas

setwd("D:/MargSalas/Ganga/Data/FarmdindisDS")
hq_area <- read.csv(file = "HQ_area.csv", sep = ";")

area_transect <- 500*1000 # m2
average_clus2022 <- 2.5