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

dat <- dat[which(dat$Year %in% 2022 & dat$Count > -1), ] # Select data in 2022 (transects sampled, presences + absences)
dat_det <- dat[which(dat$Count >0), ] # Only presences

## ---- 1. Survey characteristics ----

all.sites <- unique(dat$transectID)
nSites <- length(all.sites)

strip.width <- 500 				# strip half-width
dist.breaks <- c(0,25,50,100,200,500)
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


## ---- 3. Detection component: Sig = mu.sig (obs) + b*temp ---- 

## ---- 3.1. Data for informative prior ----

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

# Effect of temperature  (model beta coefficient)
mu.bTemp <- 0.08144117 # summary$mean[rownames(summary) %in% "bTemp.sig"]
sig.bTemp <- 0.02336605 # summary$sd[rownames(summary) %in% "bTemp.sig"]

## ---- 3.2. Temperature covariate ----

temp <- matrix(NA, nrow = length(all.sites), ncol = 1)
rownames(temp) <- all.sites

for (i in 1:nrow(dat)){
  temp[which(rownames(temp) %in% dat$transectID[i]), ] <- dat$Temp[i]
}

temp <- as.vector(temp) # I don't scale it because the coefficient was obtained from the unscaled variable
#temp_mean <- mean(temp)
#temp_sd <- sd(temp)
#temp_sc <- (temp - temp_mean) / temp_sd

# ---- 3. Abundance component: lam = alpha + b*HQ ----
## ---- 3.1. HQ covariate ----

hq_dat <- read.csv(file = "HQvariable.csv") # TO REMOVE TRANSECTS OF HQ = 0


hq <- matrix(NA, nrow = length(all.sites), ncol = 1)
rownames(hq) <- all.sites

for (i in 1:nrow(hq_dat)){
  hq[which(rownames(hq) %in% dat$transectID[i]), ] <- hq_dat$WeightedQuality[i]
}

hq <- as.vector(hq)
hq_mean <- mean(hq)
hq_sd <- sd(hq)
hq_sc <- (hq - hq_mean) / hq_sd

# ---- 4. Specify data in JAGS format ----

# Distance class and ind
nind <- nrow(dat_det)
dclass <- dat_det$Banda

m <- as.vector(m)
# Counts per site

# Fixed index to map dclass onto site and year 

site.dclass <- NULL

for (j in 1:nSites){
    site.dclass <- c(site.dclass, rep(j, m[j]))}

# ---- 5. Compile data for JAGS model ----

data1 <- list(nSites = nSites, mu.sig = mu.sig, sig.sig = sig.sig, mu.bTemp = mu.bTemp, sig.bTemp = sig.bTemp, temp = temp, 
              nG=nG, int.w=int.w, strip.width = strip.width, midpt = midpt, db = dist.breaks, 
              y = m, nind=nind, dclass=dclass, site.dclass = site.dclass, hq = hq_sc)

## ---- Inits ----

Nst <- m + 1
inits <- function(){list(alpha = runif(1), bHQ = runif(1),
                         N = Nst
)}

## ---- Params ----

params <- c("Ntotal",
            "alpha", "bHQ")

## ---- MCMC settings ----

nc <- 3 ; ni <- 1500 ; nb <- 200 ; nt <- 2

## ---- Run model ----

setwd("D:/MargSalas/Scripts_MS/Ganga/HDS/Farmdindis/Model")
source("1.HDS_sig[inf]_lam[hq].r")

# With jagsUI 
out <- jags(data1, inits, params, "1.HDS_sig[inf]_lam[hq].txt", n.chain = nc,
            n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

sum <- out$summary

# Save for Rahel
#setwd("D:/MargSalas/Ganga/Data/FarmdindisDS")
#save(dat, hq_dat, file = "datagrouse.RData")

## ---- RESULTS ----

# Predict average density in each habitat quality
hq_cat <- data.frame(hq = c(1,2,3))

# Scale with mean and SD parameters of habitat quality variable
hq_cat_sc <- (hq_cat - hq_mean) / hq_sd

#### ESTIMATE ABUNDANCE PER AREA:

# 1. Predict average abundance with each habitat category
lam.pred <- exp(sum[rownames(sum) %in% "alpha", 1] +  sum[rownames(sum) %in% "bHQ", 1] * hq_cat_sc)

# 2. Average density per habitat category
area_transect <- 500*1000 # m2
dens.pred <- lam.pred/area_transect

# 3. Abundance per habitat category (except hq = 0)
Ahq1 <- 17552400 # area in m2
Ahq2 <- 15332400
Ahq3 <- 31892800

(n1 <- dens.pred[1,1] * Ahq1)
(n2 <-dens.pred[2,1] * Ahq2)
(n3 <- dens.pred[3,1] * Ahq3)

(n <- sum(n1, n2, n3))


#####################################################

# Here I try to predict transect-level predicted lambda from model, but it is very weird

hq_df <- as.data.frame(cbind(hq, hq_sc))
hq_df <- arrange(hq_df, hq_sc)


outall <- do.call(rbind,out$samples)

pred <- list()
for(i in 1:dim(outall)[1]){ 
  pred[[i]] <- exp(outall[i,"alpha"] + 
                     outall[i,"bHQ"]*hq_df$hq_sc) } # Pred contains the list of the prediction of sigma for each iteration #(one prediction line per iteration)

predall <- do.call(rbind,pred) # All predictions/iterations together in one data frame (where columns are the prediction per each predictor (hq) values)
# predall is the abundance per site for all iterations

lci <- uci <- mean.pred <- 0 

for(i in 1:length(hq_sc)){
  lci[i]  <- quantile(predall[,i],probs = 0.025) 
  uci[i]  <- quantile(predall[,i],probs = 0.975)
  mean.pred[i]  <- mean(predall[,i])
}
# Mean abundance per site

par(mfrow = c(1,1))
plot(-15, xlim=c(-2,2), ylim=c(min(lci), max(uci)), main = "Prediction hq", xlab = "hq", ylab = "Abundance") # area_SG_HA unscaled variable

points(mean.pred ~ hq_sc, type="p")

