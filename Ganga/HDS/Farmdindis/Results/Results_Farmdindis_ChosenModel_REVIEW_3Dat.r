rm(list=ls())

## ---- Results ----

setwd("D:/MargSalas/Ganga/Results/HDS/Farmdindis/Model_results")
load("3.Dat_HDS_trendmodel_lam[hq_weight_CAT2]_sigHR.RData")

out$summary

# Load hq areas

setwd("D:/MargSalas/Ganga/Data/FarmdindisDS")
hq_area <- read.csv(file = "HQ_area.csv", sep = ";")

area_transect <- 500*1000 # m2

## -------------------------------------------------
##                 AVERAGE CLUST 2022
## ------------------------------------------------- 

average_clus2022 <- 2.5

## ---- 1. Predictions from posterior distribution ----

## ---- 1.1. Ignoring w ----

# Parameters to predict abundance in each hq zone
mu.site <- out$sims.list$mu.lam.site # Mean = summary -1.532456
random.year.2022 <- out$sims.list$log.lambda.year[,13] # Mean = summary -0.2761982
bYear.lam <- out$sims.list$bYear.lam # Mean = summary -0.02565323
year1 <- 12
bHQ1 <- out$sims.list$bHQ1 # Mean = summary
bHQ2 <- out$sims.list$bHQ2 # Mean = summary
bHQ3 <- out$sims.list$bHQ3 # Mean = summary


# Dummy covariate habitat quality to predict
hqzones <- c("hq1", "hq2", "hq3")
hq_predict <- data.frame(hq1 = c(1,0,0), hq2 = c(0,1,0), hq3 = c(0,0,1))

ab <- data.frame(matrix(NA, nrow = length(mu.site), ncol = 4))
colnames(ab) <- c(hqzones, "total")

de <- data.frame(matrix(NA, nrow = length(mu.site), ncol = 4))
colnames(de) <- c(hqzones, "total")

for (i in 1:length(hqzones)) {
  
  lambda <- exp(mu.site + random.year.2022 + bYear.lam*year1 + 
                  bHQ1 * hq_predict[1,i] + bHQ2 * hq_predict[2,i] + bHQ3 * hq_predict[3,i]) # Expected abundance
  
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
save(ab, de, file = "3Dat.Procesed_review.RData")

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

setwd("D:/MargSalas/Ganga/Results/Plots")
pdf("Abundance_estimate_revision.pdf", 7, 9)

par(mfrow = c(2,1),
    mar = c(3.2,3,2,1))
plot(dens_obs1, main = "Full posterior distribution", col = "darkcyan", lwd = 1.2, xlab = " ", ylab = " ", bty = "n", axes = FALSE)
polygon(c(dens_obs1$x, 0), c(dens_obs1$y, 0), col="darkcyan", border = "darkcyan") # ?? I still don't know if its right
polygon(c(dens_obs1$x[which(dens_obs1$x > lci & dens_obs1$x < uci)], uci, lci), 
        c(dens_obs1$y[which(dens_obs1$x > lci & dens_obs1$x < uci)], 0, 0), 
        col=adjustcolor("yellow", alpha = 0.5),
        border = adjustcolor("yellow", alpha = 0.5)) 

axis(1, pos = 0, tck = -0.05, cex.axis = 0.9)
axis(2, pos = -0.5, tck = -0.05, cex.axis = 0.9)

mtext("Abundance (N)", side = 1, line = 1.7)
mtext("Probability", side = 2, line = 1.7)


plot(dens_obs2, main = "Posterior distribution (0 - 97.5% CI)", col = "darkcyan",lwd = 1.5, xlab = "", ylab = "", bty = "n", axes = FALSE)
polygon(c(dens_obs2$x, 0), c(dens_obs2$y, 0), col="darkcyan", border = "darkcyan",lwd = 1.5) # ?? I still don't know if its right
polygon(c(dens_obs2$x[which(dens_obs2$x > lci & dens_obs2$x < uci)], uci, lci), 
        c(dens_obs2$y[which(dens_obs2$x > lci & dens_obs2$x < uci)], 0, 0), 
        col=adjustcolor("yellow", alpha = 0.5),
        border = "yellow", lwd = 1.5) 

axis(1, pos = 0, tck = -0.05, cex.axis = 0.9)
axis(2, pos = -0.5, tck = -0.05, cex.axis = 0.9)

mtext("Abundance (N)", side = 1, line = 1.7)
mtext("Probability", side = 2, line = 1.7)

abline(v = mean_ab, col = "darkslategrey", lwd = 2)
abline(v = mode_ab, col = "darkslategrey", lwd = 2)

segments(x0=lci,y0=0,x1=lci,y1=dens_obs2$y[41],col="yellow", lwd = 1.2)
segments(x0=uci,y0=0,x1=uci,y1=dens_obs2$y[270],col="yellow", lwd = 1.2)

text(x = mean_ab - 1.5, y = 0.0009, labels = "Mean:\n50 ind", adj = 0, pos = 4, col = "darkslategrey", cex = 1, font = 2)
text(x = mode_ab - 1.5, y = 0.0009, labels = "Mode:\n38 ind", adj = 0, pos = 4, col = "darkslategrey", cex = 1, font = 2)


dev.off()


## ---- 1.1.1. Table ----

hqzones <- c("hq1", "hq2", "hq3")

results <- data.frame(matrix(NA, nrow = length(hqzones)+1, ncol = 6))
rownames(results) <- c(hqzones,"total")
colnames(results) <- c("abundance_mean","abundance_lci", "abundance_uci", "density_mean","density_lci", "density_uci")

for (i in 1:(length(hqzones)+1)) {
  results[i,1] <- round(mean(ab[,i]),2)
  results[i,2] <- round(quantile(ab[,i],probs = 0.025),2)
  results[i,3] <- round(quantile(ab[,i],probs = 0.975),2)
  results[i,4] <- mean(de[,i])
  results[i,5] <- quantile(de[,i],probs = 0.025)
  results[i,6] <- quantile(de[,i],probs = 0.975)
}

results[4,1] <- paste("Mean = ", round(mean_ab, 2), ", Mode = ", round(mode_ab, 2), sep = "")

setwd("D:/MargSalas/Ganga/Results/HDS/Revision")
write.csv(results, file = "resultsHDS_2022_revision.csv")

## ---- Save parameter estimates ----

sum <- out$summary
sum <- sum[which(rownames(sum) %in% c("mu.sig", "sig.sig", "bTemp.sig", "b", "mu.lam.site", "sig.lam.site", "sig.lam.year", "bYear.lam",
                                      "bHQ1", "bHQ2", "bHQ3", "Bp.Obs", "Bp.N")), c(1,2,3,7)]
sum <- round(sum,3)

setwd("D:/MargSalas/Ganga/Results/HDS/Revision")
write.csv(sum, file = "model_estimates_3Dat.HR.csv")



#### FInalmente, he decidido no extrapolar abundance todos los años y usar el average cluster size, por lo que elimino resto del script
# (mirar script no revision para la extrapolacion de todos los años)