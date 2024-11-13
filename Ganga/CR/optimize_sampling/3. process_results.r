

## -------------------------------------------------
##                 Process results 
##          Analysis to improve cost-efficiency
## ------------------------------------------------- 

rm(list=ls())

library(jagsUI)
library(rjags)

# Import data from original model (to add point of mean abundance and CV)

setwd("D:/MargSalas/Scripts_MS/Ganga/CR/model_results")
load("out_sim_M0_psex.RData")

mean_ab_original <- out_sim_M0_psex$mean$N
mean_ab_original_5Up <- mean_ab_original + (mean_ab_original*0.05)
mean_ab_original_5Low <- mean_ab_original - (mean_ab_original*0.05)

cv_original <- (out_sim_M0_psex$sd$N/out_sim_M0_psex$mean$N)*100

## ---- 1.1. Subsampling a % ----

setwd("D:/MargSalas/Scripts_MS/Ganga/CR/optimize_sampling/results")
load("results_per2.RData")

ab_sub <- as.data.frame(matrix(NA, nrow = 100, ncol = 6))
colnames(ab_sub) <- c("50%", "75%", "80%", "90%", "95%", "100%")
max(ab_sub)

cv_sub <- as.data.frame(matrix(NA, nrow = 100, ncol = 6))
colnames(cv_sub) <- c("50%", "75%", "80%", "90%", "95%", "100%")

for(p in 1:ncol(ab_sub)){
  for (i in 1:nrow(ab_sub)){
    ab_sub[i,p] <- results_per[[p]][[i]]$mean$N
    cv_sub[i,p] <- (results_per[[p]][[i]]$sd$N/results_per[[p]][[i]]$mean$N)*100
  }
}

# Summary statistics

means <- colMeans(ab_sub)
cilo <- apply(ab_sub, 2, quantile, probs = 0.025)
cihi <- apply(ab_sub, 2, quantile, probs = 0.975)
plot_sum_ab_sub <- data.frame(group = as.numeric(1:6), subsample = as.factor(colnames(ab_sub)), means, cilo, cihi)

means <- colMeans(cv_sub)
cilo <- apply(cv_sub, 2, quantile, probs = 0.025)
cihi <- apply(cv_sub, 2, quantile, probs = 0.975)
plot_cv_sub <- data.frame(group = as.numeric(1:6), subsample = as.factor(colnames(cv_sub)), means, cilo, cihi)


# Plot

setwd("D:/MargSalas/Ganga/Results/Plots/Optimize")
pdf(file = "1. subsample.pdf", 7,7)

par(mfrow = c(2,1),
    mar = c(3, 4.1, 1, 2.1),
    par(oma = c(2, 1, 2, 3) + 0.1))
#Abundance
plot(1, xlim = c(1,6), ylim = c(50, 150), axes = FALSE, xlab = " ", ylab = "N")
axis(side = 1, at = c(1:6), labels = colnames(ab_sub))
axis(side = 2)
points(plot_sum_ab_sub$means, pch = 19)
arrows(plot_sum_ab_sub$group, plot_sum_ab_sub$cilo, plot_sum_ab_sub$group, plot_sum_ab_sub$cihi, 
       length=0.1, angle=90, code=3)
abline(h = mean_ab_original, col = "grey")
abline(h = mean_ab_original_5Up, col = "grey", lty = 2)
abline(h = mean_ab_original_5Low, col = "grey", lty = 2)


#CV
plot(1, xlim = c(1,7), ylim = c(10, 30), axes = FALSE, xlab = " ", ylab = "CV")
axis(side = 1, at = c(1:6), labels = colnames(ab_sub))
axis(side = 2)
points(plot_cv_sub$means, pch = 19)
arrows(plot_cv_sub$group, plot_cv_sub$cilo, plot_cv_sub$group, plot_cv_sub$cihi, 
       length=0.1, angle=90, code=3)
lines(plot_cv_sub$group, plot_cv_sub$means, col = "#35978f")


mtext("Subsample", side = 1, line = 0, outer = TRUE)

dev.off()

## ---- 1.2. Removing occasions ----

setwd("D:/MargSalas/Scripts_MS/Ganga/CR/optimize_sampling/results")
load("results_oc_remove.RData")

ab_oc <- as.data.frame(matrix(NA, nrow = 10, ncol = 3))
colnames(ab_oc) <- c("1", "2", "3")

cv_oc <- as.data.frame(matrix(NA, nrow = 10, ncol = 3))
colnames(cv_oc) <- c("1", "2", "3")

for (i in 1:5){ # For occasion 1
  ab_oc[i,1] <- results_oc_remove[[1]][[i]]$mean$N
  cv_oc[i,1] <- (results_oc_remove[[1]][[i]]$sd$N/results_oc_remove[[1]][[i]]$mean$N)*100
}

for(p in 2:3){
  for (i in 1:nrow(ab_oc)){
    ab_oc[i,p] <- results_oc_remove[[p]][[i]]$mean$N
    cv_oc[i,p] <- (results_oc_remove[[p]][[i]]$sd$N/results_oc_remove[[p]][[i]]$mean$N)*100
  }
}



# Summary statistics (adding the original value)

means <- c(mean_ab_original, apply(ab_oc, 2, mean, na.rm = TRUE))
cilo <- c(NA, apply(ab_oc, 2, quantile, probs = 0.025, na.rm = TRUE))
cihi <- c(NA, apply(ab_oc, 2, quantile, probs = 0.975, na.rm = TRUE))
plot_sum_ab_oc <- data.frame(group = as.numeric(1:4), means, cilo, cihi)

means <- c(cv_original, apply(cv_oc, 2, mean, na.rm = TRUE))
cilo <- c(NA, apply(cv_oc, 2, quantile, probs = 0.025, na.rm = TRUE))
cihi <- c(NA, apply(cv_oc, 2, quantile, probs = 0.975, na.rm = TRUE))
plot_cv_oc <- data.frame(group = as.numeric(1:4), means, cilo, cihi)


# Plot

setwd("D:/MargSalas/Ganga/Results/Plots/Optimize")
pdf(file = "2. rem_occasions.pdf", 7,7)

par(mfrow = c(2,1),
    mar = c(3, 4.1, 1, 2.1),
    par(oma = c(2, 1, 2, 3) + 0.1))
#Abundance
plot(1, xlim = c(1,5), ylim = c(50, 150), axes = FALSE, xlab = " ", ylab = "N")
axis(side = 1, at = c(1:4), labels = c("0", colnames(ab_oc)))
axis(side = 2)
points(plot_sum_ab_oc$means, pch = 19)
arrows(plot_sum_ab_oc$group, plot_sum_ab_oc$cilo, plot_sum_ab_oc$group, plot_sum_ab_oc$cihi, 
       length=0.1, angle=90, code=3)
abline(h = mean_ab_original, col = "grey")
abline(h = mean_ab_original_5Up, col = "grey", lty = 2)
abline(h = mean_ab_original_5Low, col = "grey", lty = 2)


#CV
plot(1, xlim = c(1,5), ylim = c(10, 35), axes = FALSE, xlab = " ", ylab = "CV")
axis(side = 1, at = c(1:4), labels = c("0", colnames(ab_oc)))
axis(side = 2)
points(plot_cv_oc$means, pch = 19)
arrows(plot_cv_oc$group, plot_cv_oc$cilo, plot_cv_oc$group, plot_cv_oc$cihi, 
       length=0.1, angle=90, code=3)
lines(plot_cv_oc$group, plot_cv_oc$means, col = "#35978f")

mtext("Occasions removed", side = 1, line = 0, outer = TRUE)

dev.off()


## ---- FINAL PLOT (Results section) ----

setwd("D:/MargSalas/Ganga/Results/Plots/Optimize")
pdf(file = "Optimize_plot.pdf", 7,6)

par(mfrow = c(2,2),
    mar = c(3, 2.5, 1, 1),
    par(oma = c(2, 2, 2, 3) + 0.1))

#Abundance sub
plot(1, xlim = c(0.8,6), ylim = c(60, 145), axes = FALSE, xlab = " ", ylab = "")
axis(side = 1, at = c(1:6), labels = colnames(ab_sub))
axis(side = 2)
points(plot_sum_ab_sub$means, pch = 19)
arrows(plot_sum_ab_sub$group, plot_sum_ab_sub$cilo, plot_sum_ab_sub$group, plot_sum_ab_sub$cihi, 
       length=0.1, angle=90, code=3)
abline(h = mean_ab_original, col = "grey")
abline(h = mean_ab_original_5Up, col = "grey", lty = 2)
abline(h = mean_ab_original_5Low, col = "grey", lty = 2)
mtext("Abundance (N)", side = 2, line = 3)


# Abundance oc
plot(1, xlim = c(0.8,4), ylim = c(60, 145), axes = FALSE, xlab = " ", ylab = "")
axis(side = 1, at = c(1:4), labels = c("0", colnames(ab_oc)))
axis(side = 2)
points(plot_sum_ab_oc$means, pch = 19)
arrows(plot_sum_ab_oc$group, plot_sum_ab_oc$cilo, plot_sum_ab_oc$group, plot_sum_ab_oc$cihi, 
       length=0.1, angle=90, code=3)
abline(h = mean_ab_original, col = "grey")
abline(h = mean_ab_original_5Up, col = "grey", lty = 2)
abline(h = mean_ab_original_5Low, col = "grey", lty = 2)

# CV sub
plot(1, xlim = c(0.8,6), ylim = c(10, 35), axes = FALSE, xlab = " ", ylab = "")
axis(side = 1, at = c(1:6), labels = colnames(ab_sub))
axis(side = 2)
points(plot_cv_sub$means, pch = 19)
arrows(plot_cv_sub$group, plot_cv_sub$cilo, plot_cv_sub$group, plot_cv_sub$cihi, 
       length=0.1, angle=90, code=3)
lines(plot_cv_sub$group, plot_cv_sub$means, col = "#35978f")
mtext("Removed samples (%)", side = 1, line = 3)
mtext("Coefficient of Variation (CV)", side = 2, line = 3)


#CV oc
plot(1, xlim = c(0.8,4), ylim = c(10, 35), axes = FALSE, xlab = " ", ylab = " ")
axis(side = 1, at = c(1:4), labels = c("0", colnames(ab_oc)))
axis(side = 2)
points(plot_cv_oc$means, pch = 19)
arrows(plot_cv_oc$group, plot_cv_oc$cilo, plot_cv_oc$group, plot_cv_oc$cihi, 
       length=0.1, angle=90, code=3)
lines(plot_cv_oc$group, plot_cv_oc$means, col = "#35978f")
mtext("Removed occasions (nº)", side = 1, line = 3)

dev.off()
