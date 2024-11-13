## -------------------------------------------------
##           Run HDS model for PTS with Data from Specific transects
##           Abundance estimation in 2022
##          Informative prior in detection
##          Informative prior in habitat quality
## ------------------------------------------------- 

# WORK IN PROGRESS
# I did this very fast and the abundance estimates are very high, so I would need to check:
# - That habitat quality variable is well calculated
# - That this model of including prior information for abundance and detection is OKAY
# - That data processing and model are well set up


rm(list=ls())

library(rjags)
library(jagsUI)


setwd("D:/MargSalas/Ganga/Data/SpecificPTS")
dat <- read.csv(file = "Data_HDS_Specific.csv")

dat_det <- dat[which(dat$Count >0), ] # Only presences

## ---- 1. Survey characteristics ----

all.sites <- unique(dat$transectID)
nSites <- length(all.sites)

strip.width <- 400 				# strip half-width
dist.breaks <- c(0,25,50,100,200,400)
int.w <- diff(dist.breaks) # width of distance categories (v)
midpt <- diff(dist.breaks)/2+dist.breaks[-5]
nG <- length(dist.breaks)-1	

## ---- 2. Detection data ----

m <- matrix(0, nrow = length(all.sites), ncol = 1)
rownames(m) <- all.sites

# Add counts > 0
count <- aggregate( Species ~ transectID, FUN = length, data = dat_det) # Only one detection per transect so its the same as detections, but I leave this line to bear in mind that you have to sum the number of rows per transect

for (i in 1:nrow(count)){
  m[which(rownames(m) %in% count$transectID[i]), 1] <- count$Species[i]
}


# Only to check: Count of individuals
count.ind <- colSums(m)
# Count of individuals corrected by cluster size
average_clus <- mean(dat_det$Count) # TO INCLUDE IN THE MODEL
count.ind_clus <- count.ind*average_clus

average_clus2 <- mean(dat_det$Count[-3]) # Remove outlier from 19 individuals


## ---- 3. Detection component: Sig = mu.sig (obs) + b*temp ---- 

## ---- 3.1. Data for informative prior ----

# FOR NOW I TAKE THE INFORMATIVE PRIOR FROM FARMDINDIS MODEL

# setwd("D:/Otros/Ganga/Trend_HDS_model_ch2/6. 10_20_AF/10_20") 
# load("PTALC_HDS_HQ_10_20.RData") # RESULTS MODEL 2010-2020
# summary <- as.data.frame(as.matrix(out$summary))

# Intercept of sigma (Mean and SD of observer random effect)
mu.sig <- 3.672591 # summary$mean[rownames(summary) %in% "mu.sig"]
sig.sig <- 0.2941518 #summary$mean[rownames(summary) %in% "sig.sig"]

#########################################################################
### RS: If you want to fix sigma to the mean from the trend model, just use
### mu.sig directly
### If you want to set an informative prior, provide mu.sig and sig.sig
### as data to be used as parameters of the prior for log.sig inside the 
### model

#log.sig <- rnorm(1, mu.sig, sig.sig)

# ---- 3. Abundance component: lam = alpha + b*HQ ----
## ---- 3.1. HQ covariate ----

hq_dat <- read.csv(file = "HQvariable2021_Specific.csv") # TO REMOVE TRANSECTS OF HQ = 0

hq <- matrix(NA, nrow = length(all.sites), ncol = 1)
rownames(hq) <- all.sites

for (i in 1:nrow(hq_dat)){
  hq[which(rownames(hq) %in% dat$transectID[i]), ] <- hq_dat$WeightedQuality[i]
}

hq <- as.vector(hq)
hq_mean <- mean(hq)
hq_sd <- sd(hq)
hq_sc <- (hq - hq_mean) / hq_sd

## ---- 3.2. HQ prior ----

# We have very scarce observations that will not allow estimating the habitat quality
# so we build an informative prior on bHQ from the farmdindis model with years 2010-2022

# setwd("D:/MargSalas/Ganga/Results/HDS/Farmdindis/Model_results")
# load("2.2.Dat_HDS_trendmodel_lam[hq]_sigHR.RData")
# summary <- as.data.frame(as.matrix(out$summary))

mu.bHQ <- 0.125796 # summary$mean[rownames(summary) %in% "bHQ"]
sig.bHQ <- 0.2112865 # summary$sd[rownames(summary) %in% "bHQ"]

# ---- 4. Specify data in JAGS format ----

# Distance class and ind
nind <- nrow(dat_det)
dclass <- dat_det$Bin

m <- as.vector(m)
# Counts per site

# Fixed index to map dclass onto site and year 

# site.dclass <- NULL
# 
# for (j in 1:nSites){
#     site.dclass <- c(site.dclass, rep(j, m[j]))}

# ---- 5. Compile data for JAGS model ----

data1 <- list(nSites = nSites, mu.sig = mu.sig, sig.sig = sig.sig, mu.bHQ = mu.bHQ, sig.bHQ = sig.bHQ, hq = hq_sc,
              nG=nG, int.w=int.w, strip.width = strip.width, midpt = midpt, db = dist.breaks, 
              y = m, nind=nind, dclass=dclass
              #, site.dclass = site.dclass
              )

## ---- Inits ----

Nst <- m + 1
inits <- function(){list(alpha = runif(1), bHQ = runif(1),
                         N = Nst
)}

## ---- Params ----

params <- c("Ntotal",
            "alpha", "bHQ", "sigma")

## ---- MCMC settings ----

nc <- 3 ; ni <- 1500 ; nb <- 200 ; nt <- 2

## ---- Run model ----

setwd("D:/MargSalas/Scripts_MS/Ganga/HDS/SpecificPTS/Model")
source("1.HDS_sig[inf]_lam[hq[inf]].r")

# With jagsUI 
out <- jags(data1, inits, params, "1.HDS_sig[inf]_lam[hq[inf]].txt", n.chain = nc,
            n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

sum <- out$summary

## ---- Results ----

# Load hq areas

setwd("D:/MargSalas/Ganga/Data/FarmdindisDS")
hq_area <- read.csv(file = "HQ_area.csv", sep = ";")

area_transect <- 500*800 # m2


## ---- 1. Predictions from posterior distribution ----

# Parameters to predict abundance in each hq zone
alpha <- out$sims.list$alpha # Mean = summary -1.532456
bHQ <- out$sims.list$bHQ # Mean = summary

hqzones <- c("hq1", "hq2", "hq3")

ab <- data.frame(matrix(NA, nrow = length(alpha), ncol = 4))
colnames(ab) <- c(hqzones, "total")

de <- data.frame(matrix(NA, nrow = length(alpha), ncol = 4))
colnames(de) <- c(hqzones, "total")

for (i in 1:length(hqzones)) {
  
  lambda <- exp(alpha + bHQ * i) # Expected abundance
  
  dens <- lambda/area_transect
  abundance <- dens*hq_area$x[i]
  total_abundance <- abundance*average_clus2
  ab[,i] <- total_abundance
  
  densHA <- lambda/(area_transect/10000)
  
  de[,i] <- densHA*average_clus2
  
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

#setwd("D:/MargSalas/Ganga/Results/HDS/Plots")
#pdf("Abundance_estimate.pdf", 7, 9)

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


#dev.off()


















