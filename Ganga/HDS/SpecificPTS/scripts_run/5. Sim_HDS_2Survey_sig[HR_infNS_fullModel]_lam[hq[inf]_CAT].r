## -------------------------------------------------
##          Data simulation HDS model for PTS
##               Specific HDS model 
##           Simulate repeated surveys
##           HQ as categorical covariate on lambda
## ------------------------------------------------- 

rm(list=ls())

library(rjags)
library(jagsUI)
library(plyr)

# Define HR detection function

g <- function(x, sig, b) 1 - exp(-(x/sig)^-b)

### Covariates ###
nSites <- 22			# number of line transect surveys
set.seed(2) 
hq <- rbinom(nSites,3,0.5)

#### SIMULATE DATASETS ####

surv <- c(10) # Number of repeated surveys

niter <- 100

# To store
# For each iteration:
dlist <- list() # To store the data
rlist <- list() # To store the results

Dlist <- list()
Rlist <- list() # To store the results

ntot <- as.data.frame(matrix(NA, nrow = niter, ncol = length(surv)))
colnames(ntot) <- as.character(surv)

cv <- as.data.frame(matrix(NA, nrow = niter, ncol = length(surv)))
colnames(cv) <- as.character(surv)

acu <- as.data.frame(matrix(NA, nrow = niter, ncol = length(surv)))
colnames(acu) <- as.character(surv)

for (sss in 1:length(surv)){
  
  for (iter in 1:niter){
    
    ## ---- 1. Survey characteristics ----
    nSurveys <- surv[[sss]]
    nSites <- 22			# number of line transect surveys
    site <- 1:nSites # Variable site
    
    strip.width <- 400 				# SAME THAN IN NATIONAL CENSUS
    dist.breaks <- c(0,100,200,300, 400)
    int.w <- diff(dist.breaks) # width of distance categories (v)
    midpt <- diff(dist.breaks)/2+dist.breaks[-5]
    nG <- length(dist.breaks)-1	
    
    ## ---- 2. Detection ----
    
    ### Parameters ###
    
    # Parameter values determined by informative priors
    # Same model than the national PTS monitoring program (informative priors in all)
    
    ## RS: In all of the below, I would use just the mean as input
    ##     parameters: so use mu.alpha as alpha.sig etc
    ##     no need to generate a new random number using the uncertainty
    ##     as we are saying, if these values are correct and we re-sample the population
    ##     when do we get good estimates
    
    mu.hSquare <- -0.00781
    sig.hSquare <- 0.0154
    b.hSquare <- mu.hSquare #rnorm(1, mu.hSquare, sig.hSquare)
    
    mu.jd <- -0.05391
    sig.jd <- 0.0240
    b.jd <- mu.jd #rnorm(1, mu.jd, sig.jd)
    
    mu.jdSquare <- 0.01176
    sig.jdSquare <- 0.0100
    b.jdSquare <- mu.jdSquare #rnorm(1, mu.jdSquare, sig.jdSquare)
    
    mu.alpha <- 3.87567
    sig.alpha <- 0.0875
    
    mu.Vebro <- 0.26297
    sig.Vebro <- 0.0690
    
    mu.alpha.new <- mu.alpha+mu.Vebro ##new alpha.sig
    sig.alpha.new <- sqrt(sig.alpha^2 + sig.Vebro^2)
    alpha.sig <- mu.alpha.new
    
    mu.b <- 0.997 # Beta
    sig.b <- 0.0438 
    beta <- mu.b # rnorm(1, mu.b, sig.b)
    
    ### Covariates ###
    
    # HOUR OF DAY
    
    # Mean and sd from NS model to center our variable
    h_mean <- 0.46054
    h_sd <- 0.2002978
    
    setwd("D:/MargSalas/Ganga/Data/SpecificPTS")
    load("hour.RData") # From the real data
    
    h <- matrix(sample(h, nSites * nSurveys, replace = TRUE), nrow = nSites) # 1. Simulate variable (using values from the real data, more realistic)
    h_sc <- (h - h_mean) / h_sd # 2. Centered on the variable used in NS (squared mean and sd) -> It should be the same of course
    hSquare_sc <- h_sc^2 # 3. Squared
    
    # JULIAN DAY
    
    # Mean and sd from NS model to center our variable
    jd_mean <- 140.2123
    jd_sd <- 23.35039
    
    setwd("D:/MargSalas/Ganga/Data/SpecificPTS")
    load("jd.RData") # 1. Simulate variable from real data
    
    jdmat <- matrix(NA, ncol = nSurveys, nrow = nSites)
    for (i in 2:ncol(jdmat)){
      jdmat[,1] <- jd
      jdmat[,i] <- jdmat[,i-1] + 7} # To mimic weekly surveys (which is not completely realistic, but the best I can think off now)
    
    jd_sc <- (jdmat - jd_mean) / jd_sd # 2. Centered on the variable used in NS
    jdSquare_sc <- jd_sc^2 #3. Squared
    
    
    ### Sigma ###
    
    ## RS: removed Vebro effect
    
    sigma <- array(NA, dim = c(nSites, nSurveys))
    for (s in 1:nSurveys) {
      sigma[, s] <- exp(alpha.sig + b.hSquare * hSquare_sc[, s] +
                          b.jd * jd_sc[, s] + b.jdSquare * jdSquare_sc[, s])}
    
    ## ---- Abundance ----
    
    ### Parameters ###
    # informative prior on bHQ from the farmdindis model with years 2010-2022
    mu.bHQ1 <- 0.9268979 
    sig.bHQ1 <- 0.5948143 
    
    mu.bHQ2 <- 1.083205 
    sig.bHQ2 <- 0.5536007 
    
    mu.bHQ3 <- 0.635831 
    sig.bHQ3 <- 0.5884766
    
    bHQ1 <- mu.bHQ1 #rnorm(1, mu.bHQ, sig.bHQ)
    bHQ2 <- mu.bHQ2 #rnorm(1, mu.bHQ, sig.bHQ)
    bHQ3 <- mu.bHQ3 #rnorm(1, mu.bHQ, sig.bHQ)
    
    alpha <- 0.6 # Intercept
    
    # Create dummy variable
    
    hq.dummy <- matrix(1, nrow = nSites, ncol = 3) # 3th dimension includes 3 arrays: one per level (intercept doesnt count)
    
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
    
    lambda <- exp(alpha + bHQ1*hqCov1 + bHQ2*hqCov2 + bHQ3*hqCov3)
    
    ### RS: First data check - are total abundances reasonable?
    N <- rpois(nSites,lambda) # Abundance per site
    N.tot <- sum(N) # Total number of individuals
    
    # ---- Simulate continuous distance data ----
    
    y <- array(0, dim = c(nSites, nG, nSurveys))
    
    for(j in 1:nSites) {
      
      if(N[j] == 0) next
      
      for(s in 1:nSurveys) {
        
        # Distance from observer to the individual
        d <- runif(N[j], 0, strip.width) 		# Uniform distribution of animals
        
        # Simulates one distance for each individual in the site (N[j])
        p <- g(x=d, sig=sigma[j,s], b = beta)   		# Detection probability. Sigma is site-time specific
        
        seen <- rbinom(N[j], 1, p)
        
        if(all(seen == 0))
          next
        
        d1 <- d[seen==1] 				# The distance data for seen individuals
        counts <- table(cut(d1, dist.breaks, include.lowest=TRUE))
        y[j,,s] <- counts
      } }
    
    y.sum.sites <- apply(y, c(1, 3), sum)
    
    # ---- Convert data to JAGS format ----
    
    nind <- sum(y.sum.sites)
    
    # Vector with counts
    yLong <- y.sum.sites
    
    # Distance category
    site.dclass <- dclass <- survey.dclass <- NULL
    
    for(j in 1:nSites){ ## HERE see with the model if I map first into sites and then into surveys or viceversa
      for (s in 1:nSurveys) {
        if (yLong[j,s] == 0) # Refers for the ditance classes to the list with years and bins
          next 
        site.dclass <- c(site.dclass, rep(j, yLong[j,s]))
        survey.dclass <- c(survey.dclass, rep(s, yLong[j, s]))
        for (k in 1:nG){
          dclass <- c(dclass, rep(k, y[j, k, s]))	# Distance category index
        } } }
    
    
    # ---- Compile data for JAGS model ----
    
    ####RS: changed input to mu.alpha.new, sig.alpha.new; removed Vebro
    
    data1 <- list(nSites = nSites, nSurveys = nSurveys,
                  mu.bHQ1 = mu.bHQ1, sig.bHQ1 = sig.bHQ1, # Informative prior for habitat quality
                  mu.bHQ2 = mu.bHQ2, sig.bHQ2 = sig.bHQ2, # Informative prior for habitat quality
                  mu.bHQ3 = mu.bHQ3, sig.bHQ3 = sig.bHQ3, # Informative prior for habitat quality
                  mu.alpha = mu.alpha.new, sig.alpha = sig.alpha.new, # Informative priors for sigma
                  mu.hSquare = mu.hSquare, sig.hSquare = sig.hSquare, 
                  mu.jd = mu.jd, sig.jd = sig.jd, 
                  mu.jdSquare = mu.jdSquare, sig.jdSquare = sig.jdSquare, 
                  mu.b = mu.b, sig.b = sig.b, # Informative prior for shape parameter
                  nG=nG, int.w=int.w, strip.width = strip.width, midpt = midpt, db = dist.breaks, 
                  y = yLong, nind=nind, dclass=dclass, site.dclass = site.dclass, survey.dclass = survey.dclass,
                  hqCov1 = hqCov1, # Variables
                  hqCov2 = hqCov2, # Variables
                  hqCov3 = hqCov3, # Variables
                  hSquare = hSquare_sc,
                  jd = jd_sc,
                  jdSquare = jdSquare_sc)
    
    ## ---- Inits ----
    
    Nst <- apply(yLong,1,max) + 1
    
    ### Inits for all parameters, especially in a simulation
    ## you want to start the chains off at a reasonable location
    ### Given that you use informative priors, you can reasonably start params
    ### at the prior mean
    
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
                "log.b",
                "b",
                "alpha.sig", "b.hSquare", "b.jd", "b.jdSquare")
    
    ## ---- Save the data ----
    
    ntot[iter,sss] <- N.tot
    dlist[[iter]] <- data1
    
    #### RUN MODEL ####
    
    ## ---- MCMC settings ----
    nc <- 3 ; ni <- 15000 ; nb <- 5000 ; nt <- 1
    
    ## ---- Run model ----
    setwd("D:/MargSalas/Scripts_MS/Ganga/HDS/SpecificPTS/Model")
    source("5.HDS_2Survey_sig[HR_inf_fullModel]_lam[hq[inf]_CAT].r")
    
    # With jagsUI 
    out <- jags(data1, inits, params, "5.HDS_2Survey_sig[HR_inf_fullModel]_lam[hq[inf]_CAT].txt", n.chain = nc,
                n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

    rlist[[iter]] <- out
    cv[iter,sss] <- (out$sd$N/out$mean$N)*100
    acu[iter,sss] <- out$mean$N - N.tot
  }
  Dlist[[sss]] <- dlist # To store the data
  Rlist[[sss]] <- rlist
}

# Save
setwd("D:/MargSalas/Ganga/Results/HDS/SpecificPTS/Model/Simulation_nSurveys")

save(Dlist, file = "data10_cat.RData")
save(Rlist, file = "results10_cat.RData")
save(cv, file = "cv10_cat.RData")
save(acu, file = "acu10_cat.RData")
save(ntot, file = "ntot10_cat.RData")


## ---- Results ----

# Load
setwd("D:/MargSalas/Ganga/Results/HDS/SpecificPTS/Model/Simulation_nSurveys")
load("data3-5_cat.RData")
load("results3-5_cat.RData")
load("ntot3-5_cat.RData")
load("cv3-5_cat.RData")
load("acu3-5_cat.RData")

Dlist1 <- Dlist
Rlist1 <- Rlist
cv1 <- cv
acu1 <- acu
ntot1 <- ntot

setwd("D:/MargSalas/Ganga/Results/HDS/SpecificPTS/Model/Simulation_nSurveys")
load("data8_cat.RData")
load("results8_cat.RData")
load("ntot8_cat.RData")
load("cv8_cat.RData")
load("acu8_cat.RData")

Dlist2 <- Dlist
Rlist2 <- Rlist
cv2 <- cv
acu2 <- acu
ntot2 <- ntot

setwd("D:/MargSalas/Ganga/Results/HDS/SpecificPTS/Model/Simulation_nSurveys")
load("data10_cat.RData")
load("results10_cat.RData")
load("ntot10_cat.RData")
load("cv10_cat.RData")
load("acu10_cat.RData")

Dlist3 <- Dlist
Rlist3 <- Rlist
cv3 <- cv
acu3 <- acu
ntot3 <- ntot

Dlist <- c(Dlist1, Dlist2, Dlist3)
Rlist <- c(Rlist1, Rlist2, Rlist3)
cv <- cbind(cv1, cv2, cv3)
acu <- cbind(acu1, acu2, acu3)
ntot <- cbind(ntot1, ntot2, ntot3)

# Save for ZENODO
setwd("D:/MargSalas/Ganga/Results/HDS/SpecificPTS/Model/Simulation_nSurveys")
save(Dlist, file = "data_cat.RData")
save(Rlist, file = "results_cat.RData")
save(cv, file = "cv_cat.RData")
save(acu, file = "acu_cat.RData")
save(ntot, file = "ntot_cat.RData")

# Summary statistics

means <- colMeans(cv)
cilo <- apply(cv, 2, quantile, probs = 0.025, na.rm = TRUE)
cihi <- apply(cv, 2, quantile, probs = 0.975, na.rm = TRUE)
plot_cv <- data.frame(group = as.numeric(1:4), occasions = as.factor(colnames(cv)), means, cilo, cihi)

means <- colMeans(acu)
cilo <- apply(acu, 2, quantile, probs = 0.025, na.rm = TRUE)
cihi <- apply(acu, 2, quantile, probs = 0.975, na.rm = TRUE)
plot_acu <- data.frame(group = as.numeric(1:4), occasions = as.factor(colnames(acu)), means, cilo, cihi)


# Import data from GMR model (to add point of mean abundance and CV)

setwd("D:/MargSalas/Scripts_MS/Ganga/CR/model_results")
load("out_sim_M0_psex.RData")

mean_ab_original <- out_sim_M0_psex$mean$N
mean_ab_original_5Up <- mean_ab_original + (mean_ab_original*0.05)
mean_ab_original_5Low <- mean_ab_original - (mean_ab_original*0.05)

cv_original <- (out_sim_M0_psex$sd$N/out_sim_M0_psex$mean$N)*100

# Plot

par(mfrow = c(2,1),
    mar = c(3, 4.1, 1, 2.1),
    par(oma = c(2, 1, 2, 3) + 0.1))

#CV
plot(1, xlim = c(1,4), ylim = c(10, 170), axes = FALSE, xlab = " ", ylab = "CV")
axis(side = 1, at = c(1:4), labels = colnames(cv))
axis(side = 2, at = c(10, 40, 70, 100, 130, 160))
points(plot_cv$means, pch = 19)
arrows(plot_cv$group, plot_cv$cilo, plot_cv$group, plot_cv$cihi, 
       length=0.1, angle=90, code=3)
lines(plot_cv$group, plot_cv$means, col = "#35978f")
abline(h = cv_original, col = "grey")


#ACURACY
plot(500, xlim = c(1,4), ylim = c(-50, +400), axes = FALSE, xlab = " ", ylab = "Accuracy")
axis(side = 1, at = c(1:4), labels = colnames(cv))
axis(side = 2)
points(plot_acu$means, pch = 19)
arrows(plot_acu$group, plot_acu$cilo, plot_acu$group, plot_acu$cihi, 
       length=0.1, angle=90, code=3)
lines(plot_acu$group, plot_acu$means, col = "#35978f")

mtext("Number of occasions", side = 1, line = 0, outer = TRUE)

######## The mean of accuracy is not 0, which means it is very biased, possibilities: 

## The distribution of abundance is quite skewed. 
## So estimate the posterior mode as an estimate of bias

acu_mode <- as.data.frame(matrix(NA, nrow = niter, ncol = 4))
colnames(acu_mode) <- colnames(acu)
niter = 100
for (s in 1:4){
  for (ite in 1:niter){
    dens_obs1 <- density(Rlist[[s]][[ite]]$sims.list$Ntotal)
    mode_ab <- dens_obs1$x[dens_obs1$y == max(dens_obs1$y)]
    acu_mode_ab <- mode_ab - ntot[ite,s]
    acu_mode[ite,s] <- acu_mode_ab
  }
}

means <- colMeans(acu_mode)
cilo <- apply(acu_mode, 2, quantile, probs = 0.025, na.rm = TRUE)
cihi <- apply(acu_mode, 2, quantile, probs = 0.975, na.rm = TRUE)
plot_acu_mode <- data.frame(group = as.numeric(1:4), occasions = as.factor(colnames(acu_mode)), means, cilo, cihi)

# Plot
plot(500, xlim = c(1,4), ylim = c(-40, +70), axes = FALSE, xlab = " ", ylab = "Accuracy")
axis(side = 1, at = c(1:4), labels = colnames(acu_mode))
axis(side = 2, at = c(-40, -20, 0, 20, 40, 60))
points(plot_acu_mode$means, pch = 19)
arrows(plot_acu_mode$group, plot_acu_mode$cilo, plot_acu_mode$group, plot_acu_mode$cihi, 
       length=0.1, angle=90, code=3)
lines(plot_acu_mode$group, plot_acu_mode$means, col = "#35978f")

mtext("Number of occasions", side = 1, line = 0, outer = TRUE)


## ---- Plot final ----

#setwd("D:/MargSalas/Ganga/Results/Plots/Increase_occasions")
#pdf(file = "Increase_oc_cat.pdf", 7,7)

par(mfrow = c(2,1),
    mar = c(3, 4.1, 1, 2.1),
    par(oma = c(2, 1, 2, 3) + 0.1))

#CV
plot(1, xlim = c(1,4), ylim = c(10, 170), axes = FALSE, xlab = " ", ylab = "CV")
axis(side = 1, at = c(1:4), labels = plot_cv$occasions)
axis(side = 2, at = c(10, 40, 70, 100, 130, 160))
points(plot_cv$means, pch = 19)
arrows(plot_cv$group, plot_cv$cilo, plot_cv$group, plot_cv$cihi, 
       length=0.1, angle=90, code=3)
lines(plot_cv$group, plot_cv$means, col = "#35978f")
abline(h = cv_original, col = "grey", lwd = 2)

#ACCURACY (MODE)
plot(500, xlim = c(1,4), ylim = c(-40, +70), axes = FALSE, xlab = " ", ylab = "Accuracy")
axis(side = 1, at = c(1:4), labels = plot_acu_mode$occasions)
axis(side = 2, at = c(-40, -20, 0, 20, 40, 60))
points(plot_acu_mode$means, pch = 19)
arrows(plot_acu_mode$group, plot_acu_mode$cilo, plot_acu_mode$group, plot_acu_mode$cihi, 
       length=0.1, angle=90, code=3)
lines(plot_acu_mode$group, plot_acu_mode$means, col = "#35978f")

mtext("Number of surveys", side = 1, line = 0, outer = TRUE)

#dev.off()

## ---- Plot final BOXPLOT ----

setwd("D:/MargSalas/Ganga/Results/Plots/Increase_occasions")
pdf(file = "Increase_oc_cat_boxplot.pdf", 7,7)

par(mfrow = c(2,1),
    mar = c(3, 4.1, 1, 2.1),
    par(oma = c(2, 1, 2, 3) + 0.1))

#CV
boxplot(cv, pch = 19, cex = 0.5, medlwd = 1.3, col = adjustcolor("#253494", alpha.f = 0.5), axes = FALSE, xlab = " ", ylab = "CV")
axis(side = 1, at = c(1:4), labels = plot_cv$occasions)
axis(side = 2, at = c(10, 40, 70, 100, 130, 160))
points(plot_cv$means, pch = 18, col = "white")
abline(h = cv_original, col = "red", lwd = 1.3)

#ACCURACY (MODE)
boxplot(acu_mode, pch = 19, cex = 0.5,  medlwd = 1.3, col = adjustcolor("#253494", alpha.f = 0.5), axes = FALSE, xlab = " ", ylab = "Bias")
axis(side = 1, at = c(1:4), labels = plot_acu_mode$occasions)
axis(side = 2, at = c(-40, 0, 40, 80, 120, 160, 200))
points(plot_acu_mode$means, pch = 18, col = "white")

mtext("Number of surveys", side = 1, line = 0, outer = TRUE)

dev.off()

