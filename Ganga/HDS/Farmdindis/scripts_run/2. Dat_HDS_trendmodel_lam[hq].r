## -------------------------------------------------
##      Run HDS trend model for PTS with DATA
##                2010 - 2022
## Same model as 6.CompareMethods_GANGA_10_20_AF_HQ
##    but include habitat quality cov in lamda
## ------------------------------------------------- 

rm(list=ls())

library(rjags)
library(jagsUI)
library(dplyr)


setwd("D:/MargSalas/Ganga/Data/FarmdindisDS")
dat <- read.csv(file = "Data_HDS_Farmdindis.csv")

# Delete AF39, only 2 observations and doesn't appear in transect list
dat <- dat[-which(dat$transectID == "AF39"), ]

dat_det <- dat[which(dat$Count >0), ] # Only presences

## ---- 1. Survey characteristics ----

all.sites <- unique(dat$transectID)
nSites <- length(all.sites)

strip.width <- 500 				# strip half-width
dist.breaks <- c(0,25,50,100,200,500)
int.w <- diff(dist.breaks) # width of distance categories (v)
midpt <- diff(dist.breaks)/2+dist.breaks[-5]
nG <- length(dist.breaks)-1

yrs <- c(2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022) # I HAVE TO CONVERT THIS FROM 0-12 (but nyrs is still 13)
nyrs <- length(yrs)

## ---- 2. Detection data ----

m <- matrix(NA, nrow = length(all.sites), ncol = nyrs)
rownames(m) <- all.sites
colnames(m) <- yrs

# Add counts > 0
count <- aggregate(Species ~ Year + transectID, FUN = length, data = dat_det)

for (i in 1:nrow(count)){
  m[which(rownames(m) %in% count$transectID[i]), which(colnames(m) %in% count$Year[i])] <- count$Species[i]
}

# Add absences
dat_abs <- dat[which(dat$Count == 0), ] # Only presences

for (i in 1:nrow(dat_abs)){
  m[which(rownames(m) %in% dat_abs$transectID[i]), which(colnames(m) %in% dat_abs$Year[i])] <- dat_abs$Count[i]
}

# Check if m = NA fit with dat$Count == -1
dat_notmade <- dat[which(dat$Count == -1), ]
nrow(dat_notmade) == length(m[is.na(m)])

# Only to check: Count of individuals per year
count.year <- colSums(m,na.rm = TRUE)

# Count of individuals per year corrected by cluster size
average_clus <- mean(dat_det$Count) # TO INCLUDE IN THE MODEL
count.year_clus <- count.year*average_clus

## ---- 3. Detection and abundance covariates ----

# Year
yrs2 <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12) # To make it as a continuous variable, otherwise it doesnt work
year <- matrix(NA, nrow = length(all.sites), ncol = nyrs)
rownames(year) <- all.sites
colnames(year) <- yrs

for (i in 1:nyrs){
  year[ ,which(colnames(year) %in% yrs[i])] <- rep(yrs2[i], nSites)
}

# Observer 
### DATA DAVID ###

# Habitat quality

hq_dat2011 <- read.csv(file = "HQvariable2011.csv") 
hq_dat2021 <- read.csv(file = "HQvariable2021.csv") 
hq_dat <- left_join(hq_dat2011, hq_dat2021, "transectID")

hq <- matrix(NA, nrow = length(all.sites), ncol = nyrs)
rownames(hq) <- all.sites
colnames(hq) <- yrs

for (i in 1:nrow(hq_dat)) {
  hq[which(rownames(hq) %in% hq_dat$transectID[i]), 1:6] <- hq_dat$WeightedQuality2011[i]
  hq[which(rownames(hq) %in% hq_dat$transectID[i]), 7:13] <- hq_dat$WeightedQuality2021[i]
}

hq_mean <- mean(hq)
hq_sd <- sd(hq)
hq_sc <- (hq - hq_mean) / hq_sd

# ---- Specify data in JAGS format ----

# Distance class and ind
nind <- nrow(dat_det)
dclass <- dat_det$Banda

m  # Counts per year and site

# Number of sites sampled per year
nsites_year <- apply(m, 2, function(x) sum(complete.cases(x)))

# Co-variates
yrs <- 1:13 
year_number <- 0:12

# Matrix with observers
ob <- matrix(as.numeric(factor(obs)), nrow = max.sites, ncol = nyrs) # JAGS doesn't accept categorical variables
unique(factor(ob))
obs_id <- unique(factor(ob))[-3]
ob[which(is.na(ob))] <- sample(obs_id, length(which(is.na(ob))), replace = TRUE) # No NA in covariate

nobs <- length(unique(factor(ob)))


# Index for random effects
site <- c(1:nSites)
year <- c(1:nyrs)

sitesYears <- NULL
for (i in 1:nyrs){
  sitesYears <- c(sitesYears,c(1:length(all.sites)))}

# Fixed index to map dclass onto site and year 
# For the index, create a matrix m where NA are 0 (because I need the same length)

m_index <- m
m_index[which(is.na(m_index))] <- 0

site.dclass <- year.dclass <- NULL

for (t in 1:nyrs){ # sites has to be nested on years because dclass first indexes the sites on the same year
  for (j in 1:nSites){
    site.dclass <- c(site.dclass, rep(j, m_index[j,t]))
    year.dclass <- c(year.dclass, rep(t, m_index[j,t]))
  } }

