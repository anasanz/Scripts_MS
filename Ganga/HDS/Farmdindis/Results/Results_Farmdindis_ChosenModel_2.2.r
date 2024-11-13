
rm(list=ls())

## ---- Results ----

setwd("D:/MargSalas/Ganga/Results/HDS/Farmdindis/Model_results")
load("2.2.Dat_HDS_trendmodel_lam[hq]_sigHR.RData")

out$summary

# Load hq areas

setwd("D:/MargSalas/Ganga/Data/FarmdindisDS")
hq_area <- read.csv(file = "HQ_area.csv", sep = ";")

area_transect <- 500*1000 # m2
average_clus2022 <- 2.5

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

setwd("D:/MargSalas/Ganga/Results/HDS/Farmdindis/Model_results")
save(ab, de, file = "2.2.Procesed.RData")

# All observations

dens_obs1 <- density(ab[,4]) 
mode_ab <- dens_obs1$x[dens_obs1$y == max(dens_obs1$y)]
mean_ab <- mean(ab[,4])

# 95% CI (excludes the 2.5% of obs with lower and higher values)

lci <- quantile(ab[,4],probs = 0.025) 
uci <- quantile(ab[,4],probs = 0.975)

# Observations of value lower than the upper ci (97.5%)

dens_obs2 <- density(ab[,4][which(ab[,4]<uci)]) 

# 90% CI (excludes the 5% of obs with lower and higher values)

lci2 <- quantile(ab[,4],probs = 0.05)
uci2 <- quantile(ab[,4],probs = 0.95)

# 85% CI (excludes the 7.5% of obs with lower and higher values)

lci3 <- quantile(ab[,4],probs = 0.075)
uci3 <- quantile(ab[,4],probs = 0.925)

# Estimate Coefficient of Variation (CV) = posterior sd / posterior mean
(sd(ab[,4])/mean(ab[,4]))*100


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

# High density interval: All points within this interval have a higher 
#  probability density than points outside the interval

library(HDInterval)
# 95%
(c(lci,uci))
hdi(ab[,4], credMass = 0.95)

# 90%
(c(lci2,uci2))
hdi(ab[,4], credMass = 0.90)

# 85%
(c(lci3,uci3))
hdi(ab[,4], credMass = 0.85)

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
sum <- sum[which(rownames(sum) %in% c("mu.sig", "sig.sig", "bTemp.sig", "b", "mu.lam.site", "sig.lam.site", "sig.lam.year", "bYear.lam", "bHQ", "sd", "rho", "Bp.Obs", "Bp.N")), c(1,2,3,7)]
sum <- round(sum,3)

setwd("D:/MargSalas/Ganga/Results/HDS")
write.csv(sum, file = "model_estimates_2.2.HR.csv")

