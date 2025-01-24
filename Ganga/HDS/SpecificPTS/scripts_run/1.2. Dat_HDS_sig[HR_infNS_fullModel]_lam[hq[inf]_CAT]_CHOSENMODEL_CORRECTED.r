## -------------------------------------------------
##           Run HDS model for PTS with Data from Specific transects
##           Abundance estimation in 2022
##     Corrected version! VEbro included in intercept
##         Included HQ as categorical covariate
## ------------------------------------------------- 


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

## ---- 3. Detection component ---- 

## ---- 3.1. Variables (same as NS detection) ----

# HOUR OF DAY
# Mean and sd from NS model to center our variable
h_mean <- 0.46054
h_sd <- 0.2002978


h <- dat$time01 ## -> NORMAL VARIABLE (not included in the model)
h_sc <- (h - h_mean) / h_sd # 2. Centered on the variable used in NS (squared mean and sd) -> It should be the same of course
hSquare_sc <- h_sc^2 # 3. Squared

summary(h_sc) # Range compared with pepe's range: -1.155147, 2.145494

# Save variable hour to get the range for simulation:
# setwd("D:/MargSalas/Ganga/Data/SpecificPTS")
# save(h,file = "hour.RData")

# JULIAN DAY
# Mean and sd from NS model to center our variable
jd_mean <- 140.2123
jd_sd <- 23.35039

jd <- dat$jd ## -> NORMAL VARIABLE
jd_sc <- (jd - jd_mean) / jd_sd # 2. Centered on the variable used in NS
jdSquare_sc <- jd_sc^2
summary(jd_sc) # Value range from NS model is (-3.606462, 4.787403)

# Save jd to get the range for simulation:
# setwd("D:/MargSalas/Ganga/Data/SpecificPTS")
# save(jd, file = "jd.RData")

## ---- 3.2. Data for informative prior ----

# Prior from national survey model: Hazard rate
mu.b <- 0.997
sig.b <- 0.0438 # This is not the sd, is the se. But according to RS it doesn't matter in ML

# Prior info for coefficients of sigma

mu.alpha <- 3.87567
sig.alpha <- 0.0875

mu.hSquare <- -0.00781
sig.hSquare <- 0.0154

mu.jd <- -0.05391
sig.jd <- 0.0240

mu.jdSquare <- 0.01176
sig.jdSquare <- 0.0100

mu.Vebro <- 0.26297
sig.Vebro <- 0.0690

# CORRECTED: intercept includes the VEbro effect now
mu.alpha.new <- mu.alpha+mu.Vebro ##new alpha.sig
sig.alpha.new <- sqrt(sig.alpha^2 + sig.Vebro^2)
alpha.sig <- mu.alpha.new


# ---- 3. Abundance component: lam = alpha + b*HQ ----
## ---- 3.1. HQ covariate ----
setwd("D:/MargSalas/Ganga/Data/SpecificPTS")
hq_dat <- read.csv(file = "HQvariable2021_Specific.csv") # TO REMOVE TRANSECTS OF HQ = 0

# New rounding to convert to categorical
hq_dat$WeightedQuality2021_CAT2 <- 0
for (i in 1:nrow(hq_dat)){
  if(hq_dat$WeightedQuality2021[i] > 0 & hq_dat$WeightedQuality2021[i] <= 1){
    hq_dat$WeightedQuality2021_CAT2[i] <- 1
  } else if (hq_dat$WeightedQuality2021[i] > 1 & hq_dat$WeightedQuality2021[i] <= 2){
    hq_dat$WeightedQuality2021_CAT2[i] <- 2
  } else if (hq_dat$WeightedQuality2021[i] > 2){
    hq_dat$WeightedQuality2021_CAT2[i] <- 3
  }
}

hq <- matrix(NA, nrow = length(all.sites), ncol = 1)
rownames(hq) <- all.sites

for (i in 1:nrow(hq_dat)){
  hq[which(rownames(hq) %in% dat$transectID[i]), ] <- hq_dat$WeightedQuality2021_CAT2[i]
}


# Create dummy variable

hq.dummy <- matrix(1, nrow = length(all.sites), ncol = 3) # 3th dimension includes 3 arrays: one per level (intercept doesnt count)

# hq.dummy[,,1] =0  hq.dummy[,,2] =0  hq.dummy[,,3] =0  ==> (hq 0, cat 1): Intercept, no need to add
# hq.dummy[,,1] =1  hq.dummy[,,2] =0 hq.dummy[,,3] =0   ==> (hq 1, cat 2): Multiply b.hq1*array #1 
# hq.dummy[,,1] =0  hq.dummy[,,2] =1 hq.dummy[,,3] =0   ==> (hq 2, cat 3): Multiply b.hq2*array #2
# hq.dummy[,,1] =0  hq.dummy[,,2] =0 hq.dummy[,,3] =1   ==> (hq 3, cat 4): Multiply b.hq3*array #3

  tmp <- tmp1 <- tmp2 <- tmp3 <- hq
  
  # Dummy variable hq 1 (only 1 appear as 1)
  tmp1[tmp1[] %in% c(2,3)] <- 0
  tmp1[tmp1[] %in% c(1)] <- 1
  hq.dummy[,1] <- tmp1
  
  # Dummy variable hq2 (only de 2 appear as 1)
  tmp2[tmp2[] %in% c(1,3)] <- 0
  tmp2[tmp2[] %in% c(2)] <- 1
  hq.dummy[,2] <- tmp2
  
  # Dummy variable hq3 (only de 3 appear as 1)
  tmp3[tmp3[] %in% c(1,2)] <- 0
  tmp3[tmp3[] %in% c(3)] <- 1
  hq.dummy[,3] <- tmp3


# Name it as covariate
hqCov1 <- hq.dummy[,1]
hqCov2 <- hq.dummy[,2]
hqCov3 <- hq.dummy[,3]


 #Save hq for simulation:
 setwd("D:/MargSalas/Ganga/Data/SpecificPTS")
 save(hq,file = "hqCAT.RData")

## ---- 3.2. HQ prior ----

# We have very scarce observations that will not allow estimating the habitat quality
# so we build an informative prior on bHQ from the farmdindis model with years 2010-2022

#setwd("D:/MargSalas/Ganga/Results/HDS/Farmdindis/Model_results")
#load("3.Dat_HDS_trendmodel_lam[hq_weight_CAT2]_sigHR.RData")
#summary <- as.data.frame(as.matrix(out$summary))

mu.bHQ1 <- 0.9268979 # summary$mean[rownames(summary) %in% "bHQ1"]
sig.bHQ1 <- 0.5948143 # summary$sd[rownames(summary) %in% "bHQ1"]

mu.bHQ2 <- 1.083205 # summary$mean[rownames(summary) %in% "bHQ2"]
sig.bHQ2 <- 0.5536007 # summary$sd[rownames(summary) %in% "bHQ2"]

mu.bHQ3 <- 0.635831 # summary$mean[rownames(summary) %in% "bHQ3"]
sig.bHQ3 <- 0.5884766 # summary$sd[rownames(summary) %in% "bHQ3"]

# ---- 4. Specify data in JAGS format ----

# Distance class and ind
nind <- nrow(dat_det)
dclass <- dat_det$Bin

m <- as.vector(m)
# Counts per site

# Fixed index to map dclass onto site and year 

site.dclass <- NULL

for (j in 1:nSites){
  site.dclass <- c(site.dclass, rep(j, m[j]))}

# ---- 5. Compile data for JAGS model ----
####RS: changed input to mu.alpha.new, sig.alpha.new; removed Vebro

data1 <- list(nSites = nSites, 
              mu.bHQ1 = mu.bHQ1, sig.bHQ1 = sig.bHQ1, # Informative prior for habitat quality
              mu.bHQ2 = mu.bHQ2, sig.bHQ2 = sig.bHQ2, # Informative prior for habitat quality
              mu.bHQ3 = mu.bHQ3, sig.bHQ3 = sig.bHQ3, # Informative prior for habitat quality
              mu.alpha = mu.alpha.new, sig.alpha = sig.alpha.new, # Informative priors for sigma
              mu.hSquare = mu.hSquare, sig.hSquare = sig.hSquare, 
              mu.jd = mu.jd, sig.jd = sig.jd, 
              mu.jdSquare = mu.jdSquare, sig.jdSquare = sig.jdSquare, 
              mu.b = mu.b, sig.b = sig.b, # Informative prior for shape parameter
              nG=nG, int.w=int.w, strip.width = strip.width, midpt = midpt, db = dist.breaks, 
              y = m, nind=nind, dclass=dclass, site.dclass = site.dclass,
              hqCov1 = hqCov1, # Variables
              hqCov2 = hqCov2, # Variables
              hqCov3 = hqCov3, # Variables
              hSquare = hSquare_sc,
              jd = jd_sc,
              jdSquare = jdSquare_sc
)
## ---- Inits ----

Nst <- m + 1

inits <- function(){list(alpha = runif(1),
                         N = Nst,
                         alpha.sig = mu.alpha.new,
                         b.hSquare = mu.hSquare,
                         b.jd = mu.jd,
                         b.jdSquare = mu.jdSquare,
                         log.b = mu.b,
                         bHQ1 = mu.bHQ1,
                         bHQ2 = mu.bHQ2,
                         bHQ3 = mu.bHQ3)}

## ---- Params ----

params <- c("Ntotal",
            "alpha", "bHQ1", "bHQ2", "bHQ3", 
            #"sigma", 
            "log.b",
            "alpha.sig", "b.hSquare", "b.jd", "b.jdSquare")

## ---- MCMC settings ----

nc <- 3 ; ni <- 150000 ; nb <- 28000 ; nt <- 2

## ---- Run model ----

setwd("D:/MargSalas/Scripts_MS/Ganga/HDS/SpecificPTS/Model")
source("1.2.HDS_sig[HR_inf_fullModel]_lam[hq[inf]_CAT]_corrected.r")

# With jagsUI 
out <- jags(data1, inits, params, "1.2.HDS_sig[HR_inf_fullModel]_lam[hq[inf]_CAT]_corrected.txt", n.chain = nc,
            n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

setwd("D:/MargSalas/Ganga/Results/HDS/SpecificPTS/Model")
save(out,file = "out_1.2_5_corrected.RData")

sum <- out$summary

library(MCMCvis)

MCMCtrace(out, 
          params = params, 
          pdf = FALSE)

