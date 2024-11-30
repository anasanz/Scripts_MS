## -------------------------------------------------
##      Run HDS trend model for PTS with DATA
##                2010 - 2022
## Same model as 6.CompareMethods_GANGA_10_20_AF_HQ
##    but include habitat quality cov in lamda
## ------------------------------------------------- 

## THIS IS THE GOOD FINAL MODEL USED (RESULTS BELOW)

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

## ---- Check cluster size ----
# Average and median cluster size

average_clus <- mean(dat_det$Count) # TO INCLUDE IN THE MODEL
median_clust <- median(dat_det$Count)

#Remove outliers
clus_no_out <- dat_det$Count[-which(dat_det$Count >= 10)]
hist(dat_det$Count)
average_clus_no_out <- mean(clus_no_out)

#hist(dat_det$Count)
#abline(v = average_clus, col = "red")
#abline(v = median_clust, col = "blue") # Better by the median, the mean will be really pushed by big occasional groups

# Just to check: Average cluster size in 2022
dat_det2022 <- dat_det[which(dat_det$Year %in% 2022), ]
average_clus2022 <- mean(dat_det2022$Count) # TO INCLUDE IN THE MODEL

# Finally, I will extrapolate to the median
medgroup <- aggregate(Count ~ Year, data = dat_det, median)
median(medgroup$Count)
?aggregate



# Count of individuals per year corrected by cluster size

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
setwd("D:/MargSalas/Ganga/Data/FarmdindisDS")

hq_dat2011 <- read.csv(file = "HQvariable2011_Farmdindis.csv") 
hq_dat2021 <- read.csv(file = "HQvariable2021_Farmdindis.csv") 
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
colSums(m, na.rm = TRUE)
sum(colSums(m, na.rm = TRUE))

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
             "popindex", "lam.tot",'Bp.Obs', 'Bp.N')


## ---- MCMC settings ----

nc <- 3 ; ni <- 700000 ; nb <- 100000 ; nt <- 5

## ---- Run model ----

setwd("D:/MargSalas/Scripts_MS/Ganga/HDS/Farmdindis/Model")
#setwd("~/Scripts_MS/Ganga/HDS/Farmdindis/Model")
source("2.2.HDS_trendmodel_lam[hq]_sigHR_noW.r")

# With jagsUI 
out <- jags(data1, inits, params, "2.2.HDS_trendmodel_lam[hq]_sigHR_noW.txt", n.chain = nc,
            n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

summary <- out$summary
print(out)

setwd("D:/MargSalas/Ganga/Results/HDS/Farmdindis/Model_results")
save(out, file = "2.2.Dat_HDS_trendmodel_lam[hq]_sigHR_noW.RData") # 60000 iter, 4 thining

## ---- Results ----

setwd("D:/MargSalas/Ganga/Results/HDS/Farmdindis/Model_results")
load("2.2.Dat_HDS_trendmodel_lam[hq]_sigHR.RData")

# Load hq areas

setwd("D:/MargSalas/Ganga/Data/FarmdindisDS")
hq_area <- read.csv(file = "HQ_area.csv", sep = ";")

area_transect <- 500*1000 # m2


## ---- 1. Predictions from posterior distribution ----

## ---- 1.1. Ignoring w ----

#plot(density(out$sims.list$popindex[,13]))
#abline(v = mean(out$sims.list$popindex[,13]), col = "blue")
# It is very skewed....

# Parameters to predict abundance in each hq zone
mu.site <- out$sims.list$mu.lam.site # Mean = summary -1.532456
random.year.2022 <- out$sims.list$log.lambda.year[,13] # Mean = summary -0.2761982
bYear.lam <- out$sims.list$bYear.lam # Mean = summary -0.02565323
year1 <- 12
bHQ <- out$sims.list$bHQ # Mean = summary

#par(mfrow = c(2,2))
#plot(density(out$sims.list$mu.lam.site))
#abline(v = mean(out$sims.list$mu.lam.site), col = "blue")
#
#plot(density(out$sims.list$log.lambda.year[,13])) ### THIS ONE IS THE ONE THAT SEEMS SKEWED?
#abline(v = mean(out$sims.list$log.lambda.year[,13]), col = "blue")
#
#plot(density(out$sims.list$bYear.lam))
#abline(v = mean(out$sims.list$bYear.lam), col = "blue")
#
#plot(density(out$sims.list$bHQ))
#abline(v = mean(out$sims.list$bHQ), col = "blue")


hqzones <- c("hq1", "hq2", "hq3")

ab <- data.frame(matrix(NA, nrow = length(mu.site), ncol = 4))
colnames(ab) <- c(hqzones, "total")

de <- data.frame(matrix(NA, nrow = length(mu.site), ncol = 4))
colnames(de) <- c(hqzones, "total")

for (i in 1:length(hqzones)) {
  
  lambda <- exp(mu.site + random.year.2022 + bYear.lam*year1 + bHQ * i) # Expected abundance
  
  dens <- lambda/area_transect
  abundance <- dens*hq_area$x[i]
  total_abundance <- abundance*average_clus2022
  ab[,i] <- total_abundance
  
  densHA <- lambda/(area_transect/10000)
  
  de[,i] <- densHA*average_clus2022
  
}

ab[,4] <- rowSums(ab[,c(1:3)])
de[,4] <- rowSums(de[,c(1:3)])

# All observations

dens_obs1 <- density(ab[,4]) 
mode_ab <- dens_obs1$x[dens_obs1$y == max(dens_obs1$y)]
mean_ab <- mean(ab[,4])

# 95% CI (excludes the 2.5% of obs with lower and higher values)

lci <- quantile(ab[,4],probs = 0.025) 
uci <- quantile(ab[,4],probs = 0.975)

# Observations of value lower than the upper ci 474.58 (97.5%)

dens_obs2 <- density(ab[,4][which(ab[,4]<uci)]) 

# 90% CI (excludes the 5% of obs with lower and higher values)

lci2 <- quantile(ab[,4],probs = 0.05)
uci2 <- quantile(ab[,4],probs = 0.95)

# 85% CI (excludes the 7.5% of obs with lower and higher values)

lci3 <- quantile(ab[,4],probs = 0.075)
uci3 <- quantile(ab[,4],probs = 0.925)


## ---- 1.1.1. Plot ----

setwd("D:/MargSalas/Ganga/Results/HDS/Plots")
pdf("Abundance_estimate.pdf", 7, 9)

par(mfrow = c(2,1),
    mar = c(3.2,3,2,1))
plot(dens_obs1, main = "Full posterior distribution", col = "darkcyan", lwd = 1.2, xlab = " ", ylab = " ", bty = "n", axes = FALSE)
polygon(c(dens_obs1$x, 0), c(dens_obs1$y, 0), col="darkcyan", border = "darkcyan") # ?? I still don't know if its right
polygon(c(dens_obs1$x[which(dens_obs1$x > lci3 & dens_obs1$x < uci3)], uci3, lci3), 
        c(dens_obs1$y[which(dens_obs1$x > lci3 & dens_obs1$x < uci3)], 0, 0), 
        col=adjustcolor("yellow", alpha = 0.5),
        border = adjustcolor("yellow", alpha = 0.5)) 

axis(1, pos = 0, tck = -0.05, cex.axis = 0.9)
axis(2, pos = -0.5, tck = -0.05, cex.axis = 0.9)

mtext("Abundance (N)", side = 1, line = 1.7)
mtext("Probability", side = 2, line = 1.7)


plot(dens_obs2, main = "Posterior distribution (0 - 97.5% CI)", col = "darkcyan",lwd = 1.5, xlab = "", ylab = "", bty = "n", axes = FALSE)
polygon(c(dens_obs2$x, 0), c(dens_obs2$y, 0), col="darkcyan", border = "darkcyan",lwd = 1.5) # ?? I still don't know if its right
polygon(c(dens_obs2$x[which(dens_obs2$x > lci3 & dens_obs2$x < uci3)], uci3, lci3), 
        c(dens_obs2$y[which(dens_obs2$x > lci3 & dens_obs2$x < uci3)], 0, 0), 
        col=adjustcolor("yellow", alpha = 0.5),
        border = "yellow", lwd = 1.5) 

axis(1, pos = 0, tck = -0.05, cex.axis = 0.9)
axis(2, pos = -0.5, tck = -0.05, cex.axis = 0.9)

mtext("Abundance (N)", side = 1, line = 1.7)
mtext("Probability", side = 2, line = 1.7)

abline(v = mean_ab, col = "darkslategrey", lwd = 2)
abline(v = mode_ab, col = "darkslategrey", lwd = 2)

segments(x0=lci3,y0=0,x1=lci3,y1=dens_obs2$y[41],col="yellow", lwd = 1.5)
segments(x0=uci3,y0=0,x1=uci3,y1=dens_obs2$y[270],col="yellow", lwd = 1.5)

text(x = 82, y = 0.0007, labels = "Mean:\n85 ind", adj = 0, pos = 4, col = "darkslategrey", cex = 1, font = 2)
text(x = 34, y = 0.0007, labels = "Mode:\n37 ind", adj = 0, pos = 4, col = "darkslategrey", cex = 1, font = 2)


dev.off()


## ---- 1.1.1. Table ----

hqzones <- c("hq1", "hq2", "hq3")

results <- data.frame(matrix(NA, nrow = length(hqzones)+1, ncol = 6))
rownames(results) <- c(hqzones,"total")
colnames(results) <- c("abundance_mean","abundance_lci", "abundance_uci", "density_mean","density_lci", "density_uci")

for (i in 1:(length(hqzones)+1)) {
  results[i,1] <- round(mean(ab[,i]),2)
  results[i,2] <- round(quantile(ab[,i],probs = 0.075),2)
  results[i,3] <- round(quantile(ab[,i],probs = 0.925),2)
  results[i,4] <- mean(de[,i])
  results[i,5] <- quantile(de[,i],probs = 0.075)
  results[i,6] <- quantile(de[,i],probs = 0.925)
}

results[4,1] <- paste("Mean = ", round(mean_ab, 2), ", Mode = ", round(mode_ab, 2), sep = "")

setwd("D:/MargSalas/Ganga/Results/HDS")
write.csv(results, file = "resultsHDS_2022.csv")

## --------- HDI interval ----

library(HDInterval)

hdi(ab[,4], credMass = 0.85)

## ---- 1.2. Including w ----

summary <- out$summary 
# w in a year t depends on the overdispersion and ac that year, and the previous one 
# Explore how it changes each year
wvec <- summary[grep("w", rownames(summary)), 1]
w <- matrix(wvec, nrow = nSites, ncol = nyrs, byrow = FALSE)

#setwd("D:/MargSalas/Ganga/Data/FarmdindisDS")
#save(w, file = "w.RData")

library("RColorBrewer")

plot(colMeans(w) , type = "l")

plot(0, xlim = c(0,13), ylim = c(-1,1))
cl <- rainbow(nSites)
for (i in 1:nSites){
  m <- lm(w[i,] ~ c(1:13))
  points(w[i,] ~ c(1:13), col = cl[i])
  abline(m, col = cl[i])
}

dim(out$sims.list$w) # Take last year

w.2022.sites <- out$sims.list$w[,,13] # Mean of w accross sites per iteration
w.2022 <- rowMeans(w.2022.sites)

hqzones <- c("hq1", "hq2", "hq3")

ab <- data.frame(matrix(NA, nrow = length(mu.site), ncol = 4))
colnames(ab) <- c(hqzones, "total")

for (i in 1:length(hqzones)) {
  
  lambda <- exp(mu.site + random.year.2022 + bYear.lam*year1 + bHQ * i + w.2022) # Expected abundance
  
  dens <- lambda/area_transect
  abundance <- dens*hq_area$x[i]
  total_abundance <- abundance*average_clus2022
  ab[,i] <- total_abundance
}

ab[,4] <- rowSums(ab[,c(1:3)])

# All observations

dens_obs1 <- density(ab[,4]) 
mode_ab <- dens_obs1$x[dens_obs1$y == max(dens_obs1$y)]
mean_ab <- mean(ab[,4])

# 95% CI (excludes the 2.5% of obs with lower and higher values)

lci <- quantile(ab[,4],probs = 0.025) 
uci <- quantile(ab[,4],probs = 0.975)

#### CONCLUSION: IT DOESN'T CHANGE MUCH, SO I KEEP THE MODEL PREDICTION WITHOUT W


## ---- 2. Predictions from posterior distribution ALL YEARS ----

# Parameters to predict abundance in each hq zone
mu.site <- out$sims.list$mu.lam.site # Mean = summary -1.532456
random.year <- out$sims.list$log.lambda.year # Mean = summary -0.2761982
bYear.lam <- out$sims.list$bYear.lam # Mean = summary -0.02565323
year1 <- 0:12
bHQ <- out$sims.list$bHQ # Mean = summary

hqzones <- c("hq1", "hq2", "hq3")

ab_years <- list()
lam_years <- list()

ab <- data.frame(matrix(NA, nrow = length(mu.site), ncol = 4))
colnames(ab) <- c(hqzones, "total")

lam <- data.frame(matrix(NA, nrow = length(mu.site), ncol = 4))
colnames(lam) <- c(hqzones, "total")

for (t in 1:length(year1)){
  for (i in 1:length(hqzones)) {
    
    lambda <- exp(mu.site + random.year[,t] + bYear.lam*year1[t] + bHQ * i) # Expected abundance
    
    dens <- lambda/area_transect
    abundance <- dens*hq_area$x[i]
    total_abundance <- abundance*average_clus2022
    ab[,i] <- total_abundance
    ab_years[[t]] <- ab
    
    lam[,i] <- lambda
    lam_years[[t]] <- lam
    
  }}

total_ab_years <- list()
total_lam_years <- list()


for (t in 1:length(year1)){
  ab <- ab_years[[t]]
  ab[,4] <- rowSums(ab[,c(1:3)])
  total_ab_years[[t]] <- ab
  
  lam <- lam_years[[t]] 
  lam[,4] <- rowSums(lam[,c(1:3)])
  total_lam_years[[t]] <- lam
}

### Results total abundance

results_allyears <- data.frame(matrix(NA, nrow = length(year1), ncol = 4))

yrs <- c(2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022)
rownames(results_allyears) <- yrs
colnames(results_allyears) <- c("Mean", "Mode", "Low 85% BCI", "Upper 85% BCI")

for (t in 1:length(year1)){
  
  dens_obs1 <- density(total_ab_years[[t]][,4]) 
  mode_ab <- dens_obs1$x[dens_obs1$y == max(dens_obs1$y)]
  mean_ab <- mean(total_ab_years[[t]][,4])
  
  # 85% CI (excludes the 7.5% of obs with lower and higher values)
  lci3 <- quantile(total_ab_years[[t]][,4],probs = 0.075)
  uci3 <- quantile(total_ab_years[[t]][,4],probs = 0.925)
  
  results_allyears[t,1] <- mean_ab
  results_allyears[t,2] <- mode_ab
  results_allyears[t,3] <- lci3
  results_allyears[t,4] <- uci3
}

### Results expected abundance per transect

results_allyears2 <- data.frame(matrix(NA, nrow = length(year1), ncol = 4))

yrs <- c(2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022)
rownames(results_allyears2) <- yrs
colnames(results_allyears2) <- c("Mean", "Mode", "Low 85% BCI", "Upper 85% BCI")

for (t in 1:length(year1)){
  
  dens_obs1 <- density(total_lam_years[[t]][,4]) 
  mode_ab <- dens_obs1$x[dens_obs1$y == max(dens_obs1$y)]
  mean_ab <- mean(total_lam_years[[t]][,4])
  
  # 85% CI (excludes the 7.5% of obs with lower and higher values)
  lci3 <- quantile(total_lam_years[[t]][,4],probs = 0.075)
  uci3 <- quantile(total_lam_years[[t]][,4],probs = 0.925)
  
  results_allyears2[t,1] <- mean_ab
  results_allyears2[t,2] <- mode_ab
  results_allyears2[t,3] <- lci3
  results_allyears2[t,4] <- uci3
}

setwd("D:/MargSalas/Ganga/Results/HDS")
write.csv(results_allyears[6:13,], file = "resultsHDS_totalab_2015-2022.csv")
write.csv(results_allyears, file = "resultsHDS_totalab_allyears.csv")
write.csv(results_allyears2, file = "resultsHDS_expectedlam_allyears.csv")

## Plot

setwd("D:/MargSalas/Ganga/Results/HDS/Plots")
pdf("Full_posterior_allyears.pdf", 7, 9)

par(mfrow = c(3,3),
    mar = c(3.2,3,2,1),
    oma = c(3,3,3,1))

for (t in 6:13){
  
  dens_obs1 <- density(total_ab_years[[t]][,4]) 
  mode_ab <- dens_obs1$x[dens_obs1$y == max(dens_obs1$y)]
  mean_ab <- mean(total_ab_years[[t]][,4])
  
  # 85% CI (excludes the 7.5% of obs with lower and higher values)
  lci3 <- quantile(total_ab_years[[t]][,4],probs = 0.075)
  uci3 <- quantile(total_ab_years[[t]][,4],probs = 0.925)
  
  plot(dens_obs1, main = yrs[t], col = "darkcyan", lwd = 1.2, xlab = " ", ylab = " ", bty = "n", axes = FALSE)
  polygon(c(dens_obs1$x, 0), c(dens_obs1$y, 0), col="darkcyan", border = "darkcyan") # ?? I still don't know if its righ
  polygon(c(dens_obs1$x[which(dens_obs1$x > lci3 & dens_obs1$x < uci3)], uci3, lci3), 
          c(dens_obs1$y[which(dens_obs1$x > lci3 & dens_obs1$x < uci3)], 0, 0), 
          col=adjustcolor("yellow", alpha = 0.5),
          border = adjustcolor("yellow", alpha = 0.5)) 
  axis(1, pos = 0, tck = -0.05, cex.axis = 0.9)
  axis(2, pos = -0.5, tck = -0.05, cex.axis = 0.9)
}

mtext("Full posterior distribution", side = 3, line = 1, outer = TRUE)
mtext("Abundance (N)", side = 1, line = 1, outer = TRUE)
mtext("Probability", side = 2, line = 1, outer = TRUE)

dev.off()

## ---- Save parameter estimates ----

sum <- out$summary
sum <- sum[which(rownames(sum) %in% c("mu.sig", "sig.sig", "bTemp.sig", "b", "mu.lam.site", "bHQ", "sd", "rho", "Bp.Obs", "Bp.N")), c(1,2,3,7)]
sum <- round(sum,3)

setwd("D:/MargSalas/Ganga/Results/HDS")
write.csv(sum, file = "model_estimates_2.2.HR.csv")

