
#setwd("D:/MargSalas/Ganga/Results/HDS/Model_results")
#load("2.2.Dat_HDS_trendmodel_lam[hq]_sigHR.RData")
#
## Load hq areas
#
#setwd("D:/MargSalas/Ganga/Data/FarmdindisDS")
#hq_area <- read.csv(file = "HQ_area.csv")
#
#area_transect <- 500*1000 # m2


## ---- 1. Predictions from posterior distribution ----

## ---- 1.1. Ignoring w ----

#plot(density(out$sims.list$popindex[,13]))
#abline(v = mean(out$sims.list$popindex[,13]), col = "blue")
# It is very skewed....

## Parameters to predict abundance in each hq zone
#mu.site <- out$sims.list$mu.lam.site # Mean = summary -1.532456
#random.year.2022 <- out$sims.list$log.lambda.year[,13] # Mean = summary -0.2761982
#bYear.lam <- out$sims.list$bYear.lam # Mean = summary -0.02565323
#year1 <- 12
#bHQ <- out$sims.list$bHQ # Mean = summary
#average_clus2022 <- 2.5

# Save for Rahel
#save(mu.site,random.year.2022,bYear.lam,year1,bHQ,average_clus2022,hq_area, file = "Dat_Rahel.RData")

setwd("D:/MargSalas/Ganga/Data/FarmdindisDS")
load("Dat_Rahel.RData")
area_transect <- 500*1000 # m2


par(mfrow = c(2,2))
plot(density(mu.site))
abline(v = mean(mu.site), col = "blue")

plot(density(random.year.2022)) ### THIS ONE IS THE ONE THAT SEEMS SKEWED?
abline(v = mean(random.year.2022), col = "blue")

plot(density(bYear.lam))
abline(v = mean(bYear.lam), col = "blue")

plot(density(bHQ))
abline(v = mean(bHQ), col = "blue")


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

## --------- HDI interval ----

library(HDInterval)

hdi(ab[,4], credMass = 0.85)
