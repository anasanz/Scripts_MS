## -------------------------------------------------
##           Run HDS model for PTS with Data from Specific transects
##           Abundance estimation in 2022
##     Corrected version! VEbro included in intercept
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

# Save hq for simulation:
# setwd("D:/MargSalas/Ganga/Data/SpecificPTS")
# save(hq,file = "hq.RData")

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

site.dclass <- NULL

for (j in 1:nSites){
  site.dclass <- c(site.dclass, rep(j, m[j]))}

# ---- 5. Compile data for JAGS model ----
####RS: changed input to mu.alpha.new, sig.alpha.new; removed Vebro

data1 <- list(nSites = nSites, 
              mu.bHQ = mu.bHQ, sig.bHQ = sig.bHQ, # Informative prior for habitat quality
              mu.alpha = mu.alpha.new, sig.alpha = sig.alpha.new, # Informative priors for sigma
              mu.hSquare = mu.hSquare, sig.hSquare = sig.hSquare, 
              mu.jd = mu.jd, sig.jd = sig.jd, 
              mu.jdSquare = mu.jdSquare, sig.jdSquare = sig.jdSquare, 
              mu.b = mu.b, sig.b = sig.b, # Informative prior for shape parameter
              nG=nG, int.w=int.w, strip.width = strip.width, midpt = midpt, db = dist.breaks, 
              y = m, nind=nind, dclass=dclass, site.dclass = site.dclass,
              hq = hq_sc, # Variables
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
                         bHQ = mu.bHQ)}

## ---- Params ----

params <- c("Ntotal",
            "alpha", "bHQ", 
            #"sigma", 
            "log.b",
            "alpha.sig", "b.hSquare", "b.jd", "b.jdSquare")

## ---- MCMC settings ----

nc <- 3 ; ni <- 150000 ; nb <- 28000 ; nt <- 2

## ---- Run model ----

setwd("D:/MargSalas/Scripts_MS/Ganga/HDS/SpecificPTS/Model")
source("1.2.HDS_sig[HR_inf_fullModel]_lam[hq[inf]]_corrected.r")

# With jagsUI 
out <- jags(data1, inits, params, "1.2.HDS_sig[HR_inf_fullModel]_lam[hq[inf]]_corrected.txt", n.chain = nc,
            n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

setwd("D:/MargSalas/Ganga/Results/HDS/SpecificPTS/Model")
save(out,file = "out_1.2_5_corrected.RData")

sum <- out$summary

library(MCMCvis)

MCMCtrace(out, 
          params = params, 
          pdf = FALSE)

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

#text(x = 82, y = 0.0007, labels = "Mean:\n1797 ind", adj = 0, pos = 4, col = "darkslategrey", cex = 1, font = 2)
#text(x = 34, y = 0.0007, labels = "Mode:\n238 ind", adj = 0, pos = 4, col = "darkslategrey", cex = 1, font = 2)


#dev.off()

mode_ab 
mean_ab 


