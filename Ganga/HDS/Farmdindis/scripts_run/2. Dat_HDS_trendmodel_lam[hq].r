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

# Delete AF39, only 1 observation and doesn't appear in transect list
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
median_clust <- median(dat_det$Count)

hist(dat_det$Count)
abline(v = average_clus, col = "red")
abline(v = median_clust, col = "blue") # Better by the median, the mean will be really pushed by big occasional groups

count.year_clus <- count.year*median_clust

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
 
# Format
obs <- matrix(NA, nrow = length(all.sites), ncol = nyrs)
rownames(obs) <- all.sites
colnames(obs) <- yrs

# Add observers for fields with counts > 0
for (i in 1:nrow(dat)){
  obs[which(rownames(obs) %in% dat$transectID[i]), which(colnames(obs) %in% dat$Year[i])] <- dat$Observer[i]
}
obs[which(obs == "")] <- NA

# Temperature
temp <- matrix(NA, nrow = length(all.sites), ncol = nyrs)
rownames(temp) <- all.sites
colnames(temp) <- yrs

# Add temper for fields with counts > 0
for (i in 1:nrow(dat)){
  temp[which(rownames(temp) %in% dat$transectID[i]), which(colnames(temp) %in% dat$Year[i])] <- dat$Temp[i]
}


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
ob <- matrix(as.numeric(factor(obs)), nrow = nSites, ncol = nyrs) # JAGS doesn't accept categorical variables
unique(factor(ob))
obs_id <- unique(factor(ob))[-3]
ob[which(is.na(ob))] <- sample(obs_id, length(which(is.na(ob))), replace = TRUE) # No NA in covariate

nobs <- length(unique(factor(ob)))

# Matrix with temperature (put random values where NA)
unique(factor(temp))
temp_id <- unique(factor(temp))[-9]
temp[which(is.na(temp))] <- sample(temp_id, length(which(is.na(temp))), replace = TRUE) # No NA in covariate
temp[which(temp == 0 | temp == 1)] <- 10 # It must be errors

temp_mean <- mean(temp)
temp_sd <- sd(temp)
temp_sc <- (temp - temp_mean) / temp_sd


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

# ---- Compile data for JAGS model ----

data1 <- list(nyears = nyrs, nsites = nSites, nG=nG, int.w=int.w, strip.width = strip.width, midpt = midpt, db = dist.breaks,
              year.dclass = year.dclass, site.dclass = site.dclass, y = m, nind=nind, dclass=dclass,
              hqCov = hq_sc, tempCov = temp_sc, ob = ob, nobs = nobs, year1 = year_number, site = site, year_index = yrs)

## ---- Inits ----

Nst <- m + 1
inits <- function(){list(mu.sig = runif(1, log(30), log(50)), sig.sig = runif(1), sig.sig.year = runif(1),
                         mu.lam.site = runif(1), sig.lam.site = 0.2, sig.lam.year = 0.3, 
                         bYear.lam = runif(1), bHQ = runif(1), 
                         N = Nst)} 
## ---- Params ----

params <- c( "mu.sig", "sig.sig", "log.sigma.year", "bTemp.sig",
             "mu.lam.site", "sig.lam.site", "sig.lam.year", "bYear.lam", "log.lambda.year", "bHQ",  # Save year effect
             "popindex", "sd", "rho", "w", "lam.tot",'Bp.Obs', 'Bp.N')


## ---- MCMC settings ----

nc <- 3 ; ni <- 700000 ; nb <- 100000 ; nt <- 5

## ---- Run model ----

setwd("D:/MargSalas/Scripts_MS/Ganga/HDS/Farmdindis/Model")
#setwd("~/Scripts_MS/Ganga/HDS/Farmdindis/Model")
source("2.HDS_trendmodel_lam[hq].r")

# With jagsUI 
out <- jags(data1, inits, params, "2.HDS_trendmodel_lam[hq].txt", n.chain = nc,
            n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

summary <- out$summary
print(out)

#setwd("~/Model_results")
#save(out, file = "2.Dat_HDS_trendmodel_lam[hq].RData") # 60000 iter, 4 thining

## ---- Results ----

setwd("D:/MargSalas/Ganga/Results/HDS/Model_results")
load("2.Dat_HDS_trendmodel_lam[hq].RData")
summary <- out$summary

