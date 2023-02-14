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
inits <- function(){list(mu.sig = runif(1, log(30), log(50)), sig.sig = runif(1), sig.sig.year = runif(1), b = runif(1),
                         mu.lam.site = runif(1), sig.lam.site = 0.2, sig.lam.year = 0.3, 
                         bYear.lam = runif(1), bHQ = runif(1), 
                         N = Nst)} 
## ---- Params ----

params <- c( "mu.sig", "sig.sig", "log.sigma.year", "bTemp.sig", "b",
             "mu.lam.site", "sig.lam.site", "sig.lam.year", "bYear.lam", "log.lambda.year", "bHQ",  # Save year effect
             "popindex", "sd", "rho", "w", "lam.tot",'Bp.Obs', 'Bp.N')


## ---- MCMC settings ----

nc <- 3 ; ni <- 700000 ; nb <- 100000 ; nt <- 5

## ---- Run model ----

setwd("D:/MargSalas/Scripts_MS/Ganga/HDS/Farmdindis/Model")
#setwd("~/Scripts_MS/Ganga/HDS/Farmdindis/Model")
source("2.2.HDS_trendmodel_lam[hq]_sigHR.r")

# With jagsUI 
out <- jags(data1, inits, params, "2.2.HDS_trendmodel_lam[hq]_sigHR.txt", n.chain = nc,
            n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

summary <- out$summary
print(out)

setwd("~/Model_results")
save(out, file = "2.2.Dat_HDS_trendmodel_lam[hq]_sigHR.RData") # 60000 iter, 4 thining

## ---- Results ----

setwd("D:/MargSalas/Ganga/Results/HDS/Model_results")
load("2.2.Dat_HDS_trendmodel_lam[hq]_sigHR.RData")

# Load hq areas

setwd("D:/MargSalas/Ganga/Data/FarmdindisDS")
hq_area <- read.csv(file = "HQ_area.csv")

area_transect <- 500*1000 # m2

mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

## ---- 1. Predictions from posterior distribution ----

## ---- 1.1. Ignoring w ----

plot(density(out$sims.list$popindex[,13]))
max(out$sims.list$popindex[,13])
mean(out$sims.list$popindex[,13])
min(out$sims.list$popindex[,13])
# It is very skewed....

# Parameters to predict abundance in each hq zone
mu.site <- out$sims.list$mu.lam.site
random.year.2022 <- out$sims.list$log.lambda.year[,13]
bYear.lam <- out$sims.list$bYear.lam
year1 <- 12
bHQ <- out$sims.list$bHQ

hqzones <- c("hq1", "hq2", "hq3")

ab <- data.frame(matrix(NA, nrow = length(hqzones), ncol = 6))
rownames(ab) <- hqzones
colnames(ab) <- c("mean_lambda", "density", "abundance", "total_abundance", "total_abundance_lci", "total_abundance_uci")

for (i in 1:length(hqzones)) {
  
  lambda <- exp(mu.site + random.year.2022 + bYear.lam*year1 + bHQ * i) # Expected abundance
  
  mean.pred  <- mean(lambda)
  dens_obs <- density(lambda)
  mode.pred <- dens_obs$x[dens_obs$y == max(dens_obs$y)]
  
  plot(density(lambda), main = hqzones[i])
  abline(v = mean.pred, col = "red")
  abline(v = mode.pred, col = "blue")
  
  
  # Mean
  
  #mean.pred  <- mean(lambda)
  #
  #dens <- mean.pred/area_transect
  #abundance <- dens*hq_area$x[i]
  #total_abundance <- abundance*average_clus
  
  # Mode
  
  dens <- mode.pred/area_transect
  abundance <- dens*hq_area$x[i]
  total_abundance <- abundance*average_clus
  
  # LCI
  
  lci <- quantile(lambda,probs = 0.025)
  
  dens_lci <- lci/area_transect
  abundance_lci <- dens_lci*hq_area$x[i]
  total_abundance_lci <- abundance_lci*average_clus
  
  # UCI
  
  uci <- quantile(lambda,probs = 0.975)
  
  dens_uci <- uci/area_transect
  abundance_uci <- dens_uci*hq_area$x[i]
  total_abundance_uci <- abundance_uci*average_clus
  
  #ab[i,1] <- mean.pred 
  ab[i,1] <- mode.pred 
  ab[i,2] <- dens
  ab[i,3] <- abundance
  ab[i,4] <- total_abundance
  ab[i,5] <- total_abundance_lci
  ab[i,6] <- total_abundance_uci
}

ta <- colSums(ab[,c(3:6)]) 

# Abundance can be estimated
# 1. Multiplying by average cluster size
average_clus # Result: 135.74363 [ 17.49625 - 293.34180]
# 2. Multiplying by median cluster size
median_clust # Result: 67.761811 [ 8.733949 - 146.433183]

## ---- 1.2. Including w ----

dim(out$sims.list$w) # Take last year

w.2022.sites <- out$sims.list$w[,,13] # Mean of w accross sites per iteration
w.2022 <- rowMeans(w.2022.sites)

ab_w <- data.frame(matrix(NA, nrow = length(hqzones), ncol = 6))
rownames(ab_w) <- hqzones
colnames(ab_w) <- c("mean_lambda", "density", "abundance", "total_abundance", "total_abundance_lci", "total_abundance_uci")

for (i in 1:length(hqzones)) {
  
  lambda <- exp(mu.site + random.year.2022 + bYear.lam*year1 + bHQ * i + w.2022) # Expected abundance
  
  # Mean
  
  mean.pred  <- mean(lambda)
  
  dens <- mean.pred/area_transect
  abundance <- dens*hq_area$x[i]
  total_abundance <- abundance*average_clus
  
  # LCI
  
  lci <- quantile(lambda,probs = 0.025)
  
  dens_lci <- lci/area_transect
  abundance_lci <- dens_lci*hq_area$x[i]
  total_abundance_lci <- abundance_lci*average_clus
  
  # UCI
  
  uci <- quantile(lambda,probs = 0.975)
  
  dens_uci <- uci/area_transect
  abundance_uci <- dens_uci*hq_area$x[i]
  total_abundance_uci <- abundance_uci*average_clus
  
  ab_w[i,1] <- mean.pred 
  ab_w[i,2] <- dens
  ab_w[i,3] <- abundance
  ab_w[i,4] <- total_abundance
  ab_w[i,5] <- total_abundance_lci
  ab_w[i,6] <- total_abundance_uci
}

ta_w <- colSums(ab_w[,c(3:6)]) 

# It doesn't change much in any case

## ---- 2. Predictions from summary ----

# Parameters to predict abundance in each hq zone
mu.site <- summary[which(rownames(summary) %in% "mu.lam.site"), 1]
random.year.2022 <- summary[which(rownames(summary) %in% "log.lambda.year[13]"), 1]
bYear.lam <- summary[which(rownames(summary) %in% "bYear.lam"), 1]
year1 <- 0:12
bHQ <- summary[which(rownames(summary) %in% "bHQ"), 1]


## ---- 2.1. Predict ignoring w ----

# Model: lambda[j,t] <- exp(log.lambda.site[site[j]] + log.lambda.year[year_index[t]] + bYear.lam*year1[t] + bHQ*hqCov[j,t] + w[j,t])

# Zone 1 (example on how to calculate it)
lambda1 <- exp(mu.site + random.year.2022 + bYear.lam*year1[13] + bHQ * 1)
dens1 <- lambda1/area_transect
abundance1 <- dens1*hq_area$x[1]
total_abundance1 <- abundance1*average_clus

lambda2 <- exp(mu.site + random.year.2022 + bYear.lam*year1[13] + bHQ * 2)
dens2 <- lambda2/area_transect
abundance2 <- dens2*hq_area$x[2]
total_abundance2 <- abundance2*average_clus

lambda3 <- exp(mu.site + random.year.2022 + bYear.lam*year1[13] + bHQ * 3)
dens3 <- lambda3/area_transect
abundance3 <- dens3*hq_area$x[3]
total_abundance3 <- abundance3*average_clus

ta <- total_abundance1 + total_abundance2 + total_abundance3


# Results in table to export

hqzones <- c("hq1", "hq2", "hq3")

ab <- data.frame(matrix(NA, nrow = length(hqzones), ncol = 4))
rownames(ab) <- hqzones
colnames(ab) <- c("lambda", "density", "abundance", "total_abundance")

for (i in 1:length(hqzones)) {
  ab[i,1] <- exp(mu.site + random.year.2022 + bYear.lam*year1[13] + bHQ * i) # Expected abundance
  ab[i,2] <- ab[i,1]/area_transect # Average transect density
  ab[i,3] <- ab[i,2]*hq_area$x[i] # 
  ab[i,4] <- ab[i,3]*average_clus
}

total_abundance <- sum(ab[,4])

## ---- 2.2. Predict taking into account w ----

# w in a year t depends on the overdispersion and ac that year, and the previous one 
# Explore how it changes each year
wvec <- summary[grep("w", rownames(summary)), 1]
w <- matrix(wvec, nrow = nSites, ncol = nyrs, byrow = FALSE)

setwd("D:/MargSalas/Ganga/Data/FarmdindisDS")
save(w, file = "w.RData")

library("RColorBrewer")

plot(colMeans(w) , type = "l")

plot(0, xlim = c(0,13), ylim = c(-1,1))
cl <- rainbow(nSites)
for (i in 1:nSites){
  m <- lm(w[i,] ~ c(1:13))
  points(w[i,] ~ c(1:13), col = cl[i])
  abline(m, col = cl[i])
}

# They don't follow a clear pattern...take the mean site w of the last year

# Zone 1 (example on how to calculate it)
lambda1 <- exp(mu.site + random.year.2022 + bYear.lam*year1[13] + bHQ * 1 + colMeans(w)[13])
dens1 <- lambda1/area_transect
abundance1 <- dens1*hq_area$x[1]
total_abundance1 <- abundance1*average_clus

lambda2 <- exp(mu.site + random.year.2022 + bYear.lam*year1[13] + bHQ * 2 + colMeans(w)[13])
dens2 <- lambda2/area_transect
abundance2 <- dens2*hq_area$x[2]
total_abundance2 <- abundance2*average_clus

lambda3 <- exp(mu.site + random.year.2022 + bYear.lam*year1[13] + bHQ * 3 + colMeans(w)[13])
dens3 <- lambda3/area_transect
abundance3 <- dens3*hq_area$x[3]
total_abundance3 <- abundance3*average_clus

ta_w <- total_abundance1 + total_abundance2 + total_abundance3

# It doesn't change much, difference of 1 individual