## -------------------------------------------------
##           Run HDS model for PTS with DATA
##           Abundance estimation in 2022
##          Informative prior in detection
## ------------------------------------------------- 

rm(list=ls())

library(rjags)
library(jagsUI)

setwd("D:/MargSalas/Ganga/Data/FarmdindisDS")
dat <- read.csv(file = "Data_HDS_Farmdindis.csv")

dat <- dat[which(dat$Year %in% 2022 & dat$Count > -1), ] # Select data in 2022 (transects sampled)

## ---- 1. Survey characteristics ----

all.sites <- unique(dat$transectID)
nSites <- length(all.sites)

strip.width <- 500 				# strip half-width
dist.breaks <- c(0,25,50,100,200,500)
int.w <- diff(dist.breaks) # width of distance categories (v)
midpt <- diff(dist.breaks)/2+dist.breaks[-5]
nG <- length(dist.breaks)-1	

## ---- 2. Detection component: Sig = mu.sig (obs) + b*temp ---- 
## ---- 2.1. Data for informative prior ----

#setwd("D:/Otros/Ganga/Trend_HDS_model_ch2/6. 10_20_AF/10_20") 
#load("PTALC_HDS_HQ_10_20.RData") # RESULTS MODEL 2010-2020
#summary <- as.data.frame(as.matrix(out$summary))

# Intercept of sigma (Mean and SD of observer random effect)
mu.sig <- 3.672591 # summary$mean[rownames(summary) %in% "mu.sig"]
sig.sig <- 0.2941518 #summary$mean[rownames(summary) %in% "sig.sig"]

log.sig <- rnorm(1, mu.sig, sig.sig)

# Effect of temperature  (model beta coefficient)
mu.bTemp <- 0.08144117 # summary$mean[rownames(summary) %in% "bTemp.sig"]
sig.bTemp <- 0.02336605 # summary$sd[rownames(summary) %in% "bTemp.sig"]

## ---- 2.2. Temperature covariate ----

temp <- dat$Temp
temp_mean <- mean(temp)
temp_sd <- sd(temp)
temp_sc <- (temp - temp_mean) / temp_sd

# ---- 3. Abundance component: lam = alpha + b*HQ ----
## ---- 3.1. HQ covariate ----

setwd("D:/MargSalas/Ganga/Data/FarmdindisDS")


hq <- 
hq_mean <- mean(hq)
hq_sd <- sd(hq)
hq_sc <- (hq - hq_mean) / hq_sd


