
## -------------------------------------------------
##              Plot abundance estimates   
## ------------------------------------------------- 

rm(list=ls())

source("D:/MargSalas/Scripts_MS/Functions/plot.violins3.r")
source("D:/MargSalas/Scripts_MS/Functions/DoScale.r")
library(tidyverse)

## ---- CR ----

setwd("D:/MargSalas/Scripts_MS/Ganga/CR/model_results")
load("out_sim_M0_psex.RData")

out_sim_M0_psex$summary # For recapture probability for females (p[1]) and males (p[2])

# Estimate number of individuals of each sex

zzall <- out_sim_M0_psex$sims.list$z

zzfem <- out_sim_M0_psex$sims.list$z
zzfem[!out_sim_M0_psex$sims.list$sex %in% c(1)] <- 3 # Only females (set males as dead)

zzmal <- out_sim_M0_psex$sims.list$z
zzmal[!out_sim_M0_psex$sims.list$sex %in% c(0)] <- 3 # Only Males (set females as dead)

N <- matrix(NA, nrow = dim(out_sim_M0_psex$sims.list$z)[1], ncol = 3)

for(ite in 1:dim(out_sim_M0_psex$sims.list$z)[1]){
  
  which.alive.all <- which(zzall[ite,] == 1) # Select only the individuals alive (z=1)
  N[ite,1] <- length(which.alive.all)
  
  which.alive.fem <- which(zzfem[ite,] == 1) # Select only the individuals alive (z=1)
  N[ite,2] <- length(which.alive.fem) 
  
  which.alive.mal <- which(zzmal[ite,] == 1) # Select only the individuals alive (z=1)
  N[ite,3] <- length(which.alive.mal) 
}

colMeans(N)
lciN <- quantile(N[,1],probs = 0.025) 
uciN <- quantile(N[,1],probs = 0.975)

lciNf <- quantile(N[,2],probs = 0.025) 
uciNf <- quantile(N[,2],probs = 0.975)

lciNm <- quantile(N[,3],probs = 0.025) 
uciNm <- quantile(N[,3],probs = 0.975)

# Sex ratio
sr <- N[,3]/N[,2]
mean(sr)
(lcSR <- quantile(sr,probs = 0.025))
(ucSR <- quantile(sr,probs = 0.975))


## ---- Farmdindis ----

setwd("D:/MargSalas/Ganga/Results/HDS/Farmdindis/Model_results")
load("2.2.Procesed.RData")

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


## ---- Plot ----

setwd("D:/MargSalas/Ganga/Results/Plots")
pdf("Abundance_estimates2.pdf", 5, 7)
par(mfrow = c(1,1),
    mar = c(5,6,4,2) + 0.1)

plot(1, ylim = c(0.5, 5 +0.5), 
     xlim = c(0,299), 
     type ="n", yaxt="n", bty = "n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = " ",
     cex.axis = 0.8)


plot.violins3(list(ab[,4]),
              x = 1,
              at = 1,
              violin.width = 0.2,
              plot.ci = 0.85,
              col = adjustcolor("#253494", alpha.f = 0.4),
              add = T,
              alpha = 0.8,
              scale.width = FALSE,
              border.col = "#253494",
              horizontal = TRUE,
              median = FALSE)


points(x = mode_ab, y = 1,
       pch = 17,
       cex = 0.7,
       col = "white")


offs <- seq(0,1,length.out = 3)

for (i in 1:3){
  plot.violins3(list(N[,i]),
                x = i,
                at = 2 + offs[i] ,
                violin.width = 0.2,
                plot.ci = 0.95,
                col = adjustcolor("#35978f", alpha.f = 0.4),
                add = T,
                alpha = 0.8,
                scale.width = FALSE,
                border.col = "#35978f",
                horizontal = TRUE,
                median = FALSE)}

axis(2, c(1:length(params)), 
     labels = c("Total", "Total", "Females", 'Males'),
     at = c(1, 2 + offs[1], 2 + offs[2], 2 + offs[3]),
     las = 2, cex.axis = 1, tick = FALSE)

#ab2 <- c(round(mean(N[,1]),0), round(mean(N[,2]),0), round(mean(N[,3]),0))
#text(x = ab2 - 8, y = c(1,2,3), labels = ab2, cex = 0.8)
dev.off()


