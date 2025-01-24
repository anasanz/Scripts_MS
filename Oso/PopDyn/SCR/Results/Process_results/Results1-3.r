
## -------------------------------------------------
##                   RESULTS 1-3
##          Check results different models 
## - From single year models with density covariate (scripts 1)
##- To open models with age structure (scripts 3)
## ------------------------------------------------- 

rm(list = ls())

library(MCMCvis)
library(coda)
source("D:/PhD/MyScripts_PhD/Ch. 2-3/Ch. 3/Results/Functions/plot.violins3.r")
source("D:/PhD/MyScripts_PhD/Ch. 2-3/Ch. 3/Results/Functions/DoScale.r")
source("D:/PhD/MyScripts_PhD/Ch. 2-3/Ch. 3/Results/Functions/ProcessCodaOutput.R")

## -------------------------------------------------
##           SCRdenscov_singleyear      
## ------------------------------------------------- 

## ---- Forest ----

# Load and plot at the same time

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year")
pdf("Topographic.pdf")

par(mfrow = c(2,2))
plot(1, ylim = c(0.5, 3+0.5), 
     xlim = c(-1.5,2), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Forest",
     cex.axis = 0.8)
axis(2, c(1:3), labels = c("2017", "2018", "2019"), las = 2, cex.axis = 1)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/Forest")

## 2017
load("sampSCRdenscov2017.RData")
summ17_forest <- MCMCsummary(samp)
#trace17 <- MCMCtrace(samp17, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 1
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2018
load("sampSCRdenscov2018.RData")
summ18_forest <- MCMCsummary(samp)
#trace18 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 2
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2019
load("sampSCRdenscov2019.RData")
summ19_forest <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 3
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## ---- Elevation ----

# Load and plot at the same time
plot(1, ylim = c(0.5, 3+0.5), 
     xlim = c(-2, 2.5), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Elevation",
     cex.axis = 0.8)
axis(2, c(1:3), labels = c("2017", "2018", "2019"), las = 2, cex.axis = 1)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/Elevation")

## 2017
load("sampSCRdenscov2017.RData")
summ17_elevation <- MCMCsummary(samp)
#trace17 <- MCMCtrace(samp17, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 1
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2018
load("sampSCRdenscov2018.RData")
summ18_elevation <- MCMCsummary(samp)
#trace18 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 2
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2019
load("sampSCRdenscov2019.RData")
summ19_elevation <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 3
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)


## ---- Roughness ----

# Load and plot at the same time
plot(1, ylim = c(0.5, 3+0.5), 
     xlim = c(-1.5, 3.5), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Roughness",
     cex.axis = 0.8)
axis(2, c(1:3), labels = c("2017", "2018", "2019"), las = 2, cex.axis = 1)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/Roughness")

## 2017
load("sampSCRdenscov2017.RData")
summ17_roughness <- MCMCsummary(samp)
#trace17 <- MCMCtrace(samp17, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 1
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2018
load("sampSCRdenscov2018.RData")
summ18_roughness <- MCMCsummary(samp)
#trace18 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 2
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2019
load("sampSCRdenscov2019.RData")
summ19_roughness <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 3
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## ---- Slope ----

# Load and plot at the same time
plot(1, ylim = c(0.5, 3+0.5), 
     xlim = c(-1.5, 4.5), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Slope",
     cex.axis = 0.8)
axis(2, c(1:3), labels = c("2017", "2018", "2019"), las = 2, cex.axis = 1)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/Slope")

## 2017
load("sampSCRdenscov2017.RData")
summ17_slope <- MCMCsummary(samp)
#trace17 <- MCMCtrace(samp17, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 1
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2018
load("sampSCRdenscov2018.RData")
summ18_slope <- MCMCsummary(samp)
#trace18 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 2
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2019
load("sampSCRdenscov2019.RData")
summ19_slope <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 3
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

dev.off()


## ---- Distcore ----

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year")
pdf("Distcore.pdf")

# Load and plot at the same time
plot(1, ylim = c(0.5, 3+0.5), 
     xlim = c(-1.5, 2), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Distance to core (log)",
     cex.axis = 0.8)
axis(2, c(1:3), labels = c("2017", "2018", "2019"), las = 2, cex.axis = 1)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/Distcore")

## 2017
load("sampSCRdenscov2017.RData")
summ17_distcore <- MCMCsummary(samp)
#trace17 <- MCMCtrace(samp17, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 1
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2018
load("sampSCRdenscov2018.RData")
summ18_distcore <- MCMCsummary(samp)
#trace18 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 2
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2019
load("sampSCRdenscov2019.RData")
summ19_distcore <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 3
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

dev.off()

## ---- obsDens ----

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year")
pdf("obsDens.pdf",9,6)

par(mfrow = c(1,2))

# Load and plot at the same time
plot(1, ylim = c(0.5, 3+0.5), 
     xlim = c(-0.5, 1), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Observation density (1996-2021)",
     cex.axis = 0.8)
axis(2, c(1:3), labels = c("2017", "2018", "2019"), las = 2, cex.axis = 1)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/obsDens")

## 2017
load("sampSCRdenscov2017.RData")
summ17_obsdens <- MCMCsummary(samp)
#trace17 <- MCMCtrace(samp17, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 1
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2018
load("sampSCRdenscov2018.RData")
summ18_obsdens <- MCMCsummary(samp)
#trace18 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 2
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2019
load("sampSCRdenscov2019.RData")
summ19_obsdens <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 3
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## ---- obsDens_preST ----

# Load and plot at the same time
plot(1, ylim = c(0.5, 3+0.5), 
     xlim = c(-0.5, 1), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Observation density (2010-2016)",
     cex.axis = 0.8)
axis(2, c(1:3), labels = c("2017", "2018", "2019"), las = 2, cex.axis = 1)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/obsDensPre")

## 2017
load("sampSCRdenscov2017.RData")
summ17_obsdens_preST <- MCMCsummary(samp)
#trace17 <- MCMCtrace(samp17, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 1
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2018
load("sampSCRdenscov2018.RData")
summ18_obsdens_preST <- MCMCsummary(samp)
#trace18 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 2
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2019
load("sampSCRdenscov2019.RData")
summ19_obsdens_preST <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 3
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

dev.off()

## ---- roads1 ----

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year")
pdf("Roads.pdf")

par(mfrow = c(2,2))
# Load and plot at the same time
plot(1, ylim = c(0.5, 3+0.5), 
     xlim = c(-7, 1), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "roads1",
     cex.axis = 0.8)
axis(2, c(1:3), labels = c("2017", "2018", "2019"), las = 2, cex.axis = 1)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/roads1")

## 2017
load("sampSCRdenscov2017.RData")
summ17_roads1 <- MCMCsummary(samp)
#trace17 <- MCMCtrace(samp17, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 1
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2018
load("sampSCRdenscov2018.RData")
summ18_roads1 <- MCMCsummary(samp)
#trace18 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 2
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2019
load("sampSCRdenscov2019.RData")
summ19_roads1 <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 3
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## ---- roads4 ----
# Load and plot at the same time
plot(1, ylim = c(0.5, 3+0.5), 
     xlim = c(-5, 1), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "roads4",
     cex.axis = 0.8)
axis(2, c(1:3), labels = c("2017", "2018", "2019"), las = 2, cex.axis = 1)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/roads4")

## 2017
load("sampSCRdenscov2017.RData")
summ17_roads4 <- MCMCsummary(samp)
#trace17 <- MCMCtrace(samp17, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 1
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2018
load("sampSCRdenscov2018.RData")
summ18_roads4 <- MCMCsummary(samp)
#trace18 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 2
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2019
load("sampSCRdenscov2019.RData")
summ19_roads4 <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 3
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## ---- roads5 ----
# Load and plot at the same time
plot(1, ylim = c(0.5, 3+0.5), 
     xlim = c(-2, 2.5), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "roads5",
     cex.axis = 0.8)
axis(2, c(1:3), labels = c("2017", "2018", "2019"), las = 2, cex.axis = 1)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/roads5")

## 2017
load("sampSCRdenscov2017.RData")
summ17_roads5 <- MCMCsummary(samp)
#trace17 <- MCMCtrace(samp17, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 1
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2018
load("sampSCRdenscov2018.RData")
summ18_roads5 <- MCMCsummary(samp)
#trace18 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 2
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2019
load("sampSCRdenscov2019.RData")
summ19_roads5 <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 3
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)


dev.off()

## -------------------------------------------------
##           SCRdenscov_singleyear
##             FinalData17-21
## ------------------------------------------------- 

# Data frame to store values of p0 accross models and years

covs <- c("forest", "elevation", "roughness", "slope", "logDistcore", "obsDens200m", "obsDens200m_preST",
          "roads1", "roads4", "roads5", "roads6")
yrs <- c("2017", "2018", "2019", "2020", "2021")

p0 <- data.frame(matrix(NA, nrow = length(yrs), ncol = length(covs)))
rownames(p0) <- yrs
colnames(p0) <- covs


## ---- Forest ----

# Load and plot at the same time

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_17_19")
pdf("Topographic.pdf")

par(mfrow = c(2,2))
plot(1, ylim = c(0.5, 5+0.5), 
     xlim = c(-1.5,2), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Forest",
     cex.axis = 0.8)
axis(2, c(1:5), labels = c("2017", "2018", "2019", "2020", "2021"), las = 2, cex.axis = 1)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_17_19/forest")

## 2017
load("sampSCRdenscov2017.RData")
summ17_forest <- MCMCsummary(samp)
#trace17 <- MCMCtrace(samp17, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 1
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2018
load("sampSCRdenscov2018.RData")
summ18_forest <- MCMCsummary(samp)
#trace18 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 2
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2019
load("sampSCRdenscov2019.RData")
summ19_forest <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 3
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2020
load("sampSCRdenscov2020.RData")
summ20_forest <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 4
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2021
load("sampSCRdenscov2021.RData")
summ21_forest <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 5
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

files <- list(summ17_forest, summ18_forest, summ19_forest, summ20_forest, summ21_forest)

for (i in 1:length(yrs)){
  p0[i,which(colnames(p0) %in% "forest")] <- paste(round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 1],3), 
                                                   "(",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 3],3),
                                                   "-",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 5],3),
                                                   ")",sep = "")
}


## ---- Elevation ----

# Load and plot at the same time
plot(1, ylim = c(0.5, 5+0.5), 
     xlim = c(-2, 2.5), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Elevation",
     cex.axis = 0.8)
axis(2, c(1:5), labels = c("2017", "2018", "2019", "2020", "2021"), las = 2, cex.axis = 1)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_17_19/dem")

## 2017
load("sampSCRdenscov2017.RData")
summ17_elevation <- MCMCsummary(samp)
#trace17 <- MCMCtrace(samp17, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 1
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2018
load("sampSCRdenscov2018.RData")
summ18_elevation <- MCMCsummary(samp)
#trace18 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 2
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2019
load("sampSCRdenscov2019.RData")
summ19_elevation <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 3
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2020
load("sampSCRdenscov2020.RData")
summ20_elevation <- MCMCsummary(samp)
#trace20 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 4
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2021
load("sampSCRdenscov2021.RData")
summ21_elevation <- MCMCsummary(samp)
#trace21 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 5
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

files <- list(summ17_elevation, summ18_elevation, summ19_elevation, summ20_elevation, summ21_elevation)

for (i in 1:length(yrs)){
  p0[i,which(colnames(p0) %in% "elevation")] <- paste(round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 1],3), 
                                                   "(",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 3],3),
                                                   "-",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 5],3),
                                                   ")",sep = "")
}


## ---- Roughness ----

# Load and plot at the same time
plot(1, ylim = c(0.5, 5+0.5), 
     xlim = c(-1.5, 3.5), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Roughness",
     cex.axis = 0.8)
axis(2, c(1:5), labels = c("2017", "2018", "2019", "2020", "2021"), las = 2, cex.axis = 1)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_17_19/rough")

## 2017
load("sampSCRdenscov2017.RData")
summ17_roughness <- MCMCsummary(samp)
#trace17 <- MCMCtrace(samp17, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 1
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2018
load("sampSCRdenscov2018.RData")
summ18_roughness <- MCMCsummary(samp)
#trace18 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 2
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2019
load("sampSCRdenscov2019.RData")
summ19_roughness <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 3
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2020
load("sampSCRdenscov2020.RData")
summ20_roughness <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 4
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2021
load("sampSCRdenscov2021.RData")
summ21_roughness <- MCMCsummary(samp)
#trace21 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 5
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

files <- list(summ17_roughness, summ18_roughness, summ19_roughness, summ20_roughness, summ21_roughness)

for (i in 1:length(yrs)){
  p0[i,which(colnames(p0) %in% "roughness")] <- paste(round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 1],3), 
                                                      "(",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 3],3),
                                                      "-",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 5],3),
                                                      ")",sep = "")
}

## ---- Slope ----

# Load and plot at the same time
plot(1, ylim = c(0.5, 5+0.5), 
     xlim = c(-1.5, 4.5), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Slope",
     cex.axis = 0.8)
axis(2, c(1:5), labels = c("2017", "2018", "2019", "2020", "2021"), las = 2, cex.axis = 1)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_17_19/slope")

## 2017
load("sampSCRdenscov2017.RData")
summ17_slope <- MCMCsummary(samp)
#trace17 <- MCMCtrace(samp17, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 1
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2018
load("sampSCRdenscov2018.RData")
summ18_slope <- MCMCsummary(samp)
#trace18 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 2
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2019
load("sampSCRdenscov2019.RData")
summ19_slope <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 3
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2020
load("sampSCRdenscov2020.RData")
summ20_slope <- MCMCsummary(samp)
#trace20 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 4
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2021
load("sampSCRdenscov2021.RData")
summ21_slope <- MCMCsummary(samp)
#trace21 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 5
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

dev.off()

files <- list(summ17_slope, summ18_slope, summ19_slope, summ20_slope, summ21_slope)

for (i in 1:length(yrs)){
  p0[i,which(colnames(p0) %in% "slope")] <- paste(round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 1],3), 
                                                      "(",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 3],3),
                                                      "-",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 5],3),
                                                      ")",sep = "")
}

## ---- Distcore ----

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_17_19")
pdf("Distcore.pdf")

# Load and plot at the same time
plot(1, ylim = c(0.5, 5+0.5), 
     xlim = c(-1.5, 2), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Distance to core (log)",
     cex.axis = 0.8)
axis(2, c(1:5), labels = c("2017", "2018", "2019", "2020", "2021"), las = 2, cex.axis = 1)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_17_19/logDistcore")

## 2017
load("sampSCRdenscov2017.RData")
summ17_distcore <- MCMCsummary(samp)
#trace17 <- MCMCtrace(samp17, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 1
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2018
load("sampSCRdenscov2018.RData")
summ18_distcore <- MCMCsummary(samp)
#trace18 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 2
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2019
load("sampSCRdenscov2019.RData")
summ19_distcore <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 3
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2020
load("sampSCRdenscov2020.RData")
summ20_distcore <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 4
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2021
load("sampSCRdenscov2021.RData")
summ21_distcore <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 5
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

dev.off()

files <- list(summ17_distcore, summ18_distcore, summ19_distcore, summ20_distcore, summ21_distcore)

for (i in 1:length(yrs)){
  p0[i,which(colnames(p0) %in% "logDistcore")] <- paste(round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 1],3), 
                                                  "(",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 3],3),
                                                  "-",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 5],3),
                                                  ")",sep = "")
}

## ---- obsDens ----

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_17_19")
pdf("obsDens.pdf",9,6)

par(mfrow = c(1,2))

# Load and plot at the same time
plot(1, ylim = c(0.5, 5+0.5), 
     xlim = c(-0.5, 1), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Observation density (1996-2021)",
     cex.axis = 0.8)
axis(2, c(1:5), labels = c("2017", "2018", "2019", "2020", "2021"), las = 2, cex.axis = 1)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_17_19/obsDens200m")

## 2017
load("sampSCRdenscov2017.RData")
summ17_obsdens <- MCMCsummary(samp)
#trace17 <- MCMCtrace(samp17, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 1
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2018
load("sampSCRdenscov2018.RData")
summ18_obsdens <- MCMCsummary(samp)
#trace18 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 2
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2019
load("sampSCRdenscov2019.RData")
summ19_obsdens <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 3
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2020
load("sampSCRdenscov2020.RData")
summ20_obsdens <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 4
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2021
load("sampSCRdenscov2021.RData")
summ21_obsdens <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 5
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

files <- list(summ17_obsdens, summ18_obsdens, summ19_obsdens, summ20_obsdens, summ21_obsdens)

for (i in 1:length(yrs)){
  p0[i,which(colnames(p0) %in% "obsDens200m")] <- paste(round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 1],3), 
                                                        "(",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 3],3),
                                                        "-",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 5],3),
                                                        ")",sep = "")
}

## ---- obsDens_preST ----

# Load and plot at the same time
plot(1, ylim = c(0.5, 5+0.5), 
     xlim = c(-0.5, 1), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Observation density (2010-2016)",
     cex.axis = 0.8)
axis(2, c(1:5), labels = c("2017", "2018", "2019", "2020", "2021"), las = 2, cex.axis = 1)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_17_19/obsDens200m_preST")

## 2017
load("sampSCRdenscov2017.RData")
summ17_obsdens_preST <- MCMCsummary(samp)
#trace17 <- MCMCtrace(samp17, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 1
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2018
load("sampSCRdenscov2018.RData")
summ18_obsdens_preST <- MCMCsummary(samp)
#trace18 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 2
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2019
load("sampSCRdenscov2019.RData")
summ19_obsdens_preST <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 3
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2020
load("sampSCRdenscov2020.RData")
summ20_obsdens_preST <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 4
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2021
load("sampSCRdenscov2021.RData")
summ21_obsdens_preST <- MCMCsummary(samp)
#trace21 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 5
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

dev.off()

files <- list(summ17_preST, summ18_preST, summ19_preST, summ20_preST, summ21_preST)

for (i in 1:length(yrs)){
  p0[i,which(colnames(p0) %in% "obsDens200m_preST")] <- paste(round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 1],3), 
                                                        "(",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 3],3),
                                                        "-",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 5],3),
                                                        ")",sep = "")
}


## ---- roads1 ----

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_17_19")
pdf("Roads.pdf")

par(mfrow = c(2,2))
# Load and plot at the same time
plot(1, ylim = c(0.5, 5+0.5), 
     xlim = c(-7, 1), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "roads1",
     cex.axis = 0.8)
axis(2, c(1:5), labels = c("2017", "2018", "2019", "2020", "2021"), las = 2, cex.axis = 1)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_17_19/roads1")

## 2017
load("sampSCRdenscov2017.RData")
summ17_roads1 <- MCMCsummary(samp)
#trace17 <- MCMCtrace(samp17, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 1
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2018
load("sampSCRdenscov2018.RData")
summ18_roads1 <- MCMCsummary(samp)
#trace18 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 2
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2019
load("sampSCRdenscov2019.RData")
summ19_roads1 <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 3
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2020
load("sampSCRdenscov2020.RData")
summ20_roads1 <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 4
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2021
load("sampSCRdenscov2021.RData")
summ21_roads1 <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 5
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)


files <- list(summ17_roads1, summ18_roads1, summ19_roads1, summ20_roads1, summ21_roads1)

for (i in 1:length(yrs)){
  p0[i,which(colnames(p0) %in% "roads1")] <- paste(round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 1],3), 
                                                              "(",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 3],3),
                                                              "-",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 5],3),
                                                              ")",sep = "")
}


## ---- roads6 ----
# Load and plot at the same time
plot(1, ylim = c(0.5, 5+0.5), 
     xlim = c(-2, 2.5), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "roads6",
     cex.axis = 0.8)
axis(2, c(1:5), labels = c("2017", "2018", "2019", "2020", "2021"), las = 2, cex.axis = 1)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_17_19/roads6")

## 2017
load("sampSCRdenscov2017.RData")
summ17_roads6 <- MCMCsummary(samp)
#trace17 <- MCMCtrace(samp17, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 1
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2018
load("sampSCRdenscov2018.RData")
summ18_roads6 <- MCMCsummary(samp)
#trace18 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 2
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2019
load("sampSCRdenscov2019.RData")
summ19_roads6 <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 3
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2020
load("sampSCRdenscov2020.RData")
summ20_roads6 <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 4
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2021
load("sampSCRdenscov2021.RData")
summ21_roads6 <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 5
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

files <- list(summ17_roads6, summ18_roads6, summ19_roads6, summ20_roads6, summ21_roads6)

for (i in 1:length(yrs)){
  p0[i,which(colnames(p0) %in% "roads6")] <- paste(round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 1],3), 
                                                   "(",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 3],3),
                                                   "-",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 5],3),
                                                   ")",sep = "")
}


## ---- roads4 ----
# Load and plot at the same time
plot(1, ylim = c(0.5, 5+0.5), 
     xlim = c(-5, 1), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "roads4",
     cex.axis = 0.8)
axis(2, c(1:5), labels = c("2017", "2018", "2019", "2020", "2021"), las = 2, cex.axis = 1)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_17_19/roads4")

## 2017
load("sampSCRdenscov2017.RData")
summ17_roads4 <- MCMCsummary(samp)
#trace17 <- MCMCtrace(samp17, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 1
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2018
load("sampSCRdenscov2018.RData")
summ18_roads4 <- MCMCsummary(samp)
#trace18 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 2
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2019
load("sampSCRdenscov2019.RData")
summ19_roads4 <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 3
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2020
load("sampSCRdenscov2020.RData")
summ20_roads4 <- MCMCsummary(samp)
#trace20 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 4
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2021
load("sampSCRdenscov2021.RData")
summ21_roads4 <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 5
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

files <- list(summ17_roads4, summ18_roads4, summ19_roads4, summ20_roads4, summ21_roads4)

for (i in 1:length(yrs)){
  p0[i,which(colnames(p0) %in% "roads4")] <- paste(round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 1],3), 
                                                   "(",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 3],3),
                                                   "-",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 5],3),
                                                   ")",sep = "")
}


## ---- roads5 ----
# Load and plot at the same time
plot(1, ylim = c(0.5, 5+0.5), 
     xlim = c(-2, 2.5), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "roads5",
     cex.axis = 0.8)
axis(2, c(1:5), labels = c("2017", "2018", "2019", "2020", "2021"), las = 2, cex.axis = 1)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_17_19/roads5")

## 2017
load("sampSCRdenscov2017.RData")
summ17_roads5 <- MCMCsummary(samp)
#trace17 <- MCMCtrace(samp17, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 1
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2018
load("sampSCRdenscov2018.RData")
summ18_roads5 <- MCMCsummary(samp)
#trace18 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 2
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2019
load("sampSCRdenscov2019.RData")
summ19_roads5 <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 3
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2020
load("sampSCRdenscov2020.RData")
summ20_roads5 <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 4
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2021
load("sampSCRdenscov2021.RData")
summ21_roads5 <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 5
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

dev.off()

files <- list(summ17_roads5, summ18_roads5, summ19_roads5, summ20_roads5, summ21_roads5)

for (i in 1:length(yrs)){
  p0[i,which(colnames(p0) %in% "roads5")] <- paste(round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 1],3), 
                                                   "(",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 3],3),
                                                   "-",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 5],3),
                                                   ")",sep = "")
}

setwd("D:/MargSalas/Oso/SCR/Exp_analysis")
#openxlsx::write.xlsx(p0, 'p0_singleyear_allcovs.xlsx')



## -------------------------------------------------
##           SCRdenscov_singleyear
##             FinalData17-21
##           p0 = effort + trap + b / Sig[sex]
## ------------------------------------------------- 

# Data frame to store values of p0 accross models and years

covs <- c("forest", "slope", "logDistcore", "roads1")

yrs <- c("2017", "2018", "2019", "2020", "2021")

p0 <- data.frame(matrix(NA, nrow = length(yrs), ncol = length(covs)))
rownames(p0) <- yrs
colnames(p0) <- covs

## ---- Forest ----

# Load and plot at the same time

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_covsp0")
pdf("Topographic.pdf")

par(mfrow = c(2,2))
plot(1, ylim = c(0.5, 5+0.5), 
     xlim = c(-1.5,2), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Forest",
     cex.axis = 0.8)
axis(2, c(1:5), labels = c("2017", "2018", "2019", "2020", "2021"), las = 2, cex.axis = 1)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_covsp0/Results_server_cyril/forest")

## 2017
load("myResults_2017.RData")
samp <- nimOutput
summ17_forest <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

MCMCtrace(samp,
          ind = TRUE,
          pdf = TRUE)

i = 1
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2018
load("myResults_2018.RData")
samp <- nimOutput
summ18_forest <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

MCMCtrace(samp,
          ind = TRUE,
          pdf = TRUE)

i = 2
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2019
load("myResults_2019.RData")
samp <- nimOutput
summ19_forest <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

MCMCtrace(samp,
          ind = TRUE,
          pdf = TRUE)

i = 3
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2020
load("myResults_2020.RData")
samp <- nimOutput
summ20_forest <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

MCMCtrace(samp,
          ind = TRUE,
          pdf = TRUE)


i = 4
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2021
load("myResults_2021.RData")
samp <- nimOutput
summ21_forest <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

MCMCtrace(samp,
          ind = TRUE,
          pdf = TRUE)


i = 5
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

files <- list(summ17_forest, summ18_forest, summ19_forest, summ20_forest, summ21_forest)

for (i in 1:length(yrs)){
  p0[i,which(colnames(p0) %in% "forest")] <- paste(round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 1],3), 
                                                   "(",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 3],3),
                                                   "-",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 5],3),
                                                   ")",sep = "")
}

## ---- Slope ----

# Load and plot at the same time

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_covsp0")
pdf("Topographic.pdf")

par(mfrow = c(2,2))
plot(1, ylim = c(0.5, 5+0.5), 
     xlim = c(-1.5,2), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Forest",
     cex.axis = 0.8)
axis(2, c(1:5), labels = c("2017", "2018", "2019", "2020", "2021"), las = 2, cex.axis = 1)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_covsp0/Results_server_cyril/slope")

## 2017
load("myResults_2017.RData")
samp <- nimOutput
summ17_slope <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

MCMCtrace(samp,
          ind = TRUE,
          pdf = TRUE)

i = 1
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2018
load("myResults_2018.RData")
samp <- nimOutput
summ18_slope <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

MCMCtrace(samp,
          ind = TRUE,
          pdf = TRUE)

i = 2
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2019
load("myResults_2019.RData")
samp <- nimOutput
summ19_slope <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

MCMCtrace(samp,
          ind = TRUE,
          pdf = TRUE)

i = 3
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2020
load("myResults_2020.RData")
samp <- nimOutput
summ20_slope <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

MCMCtrace(samp,
          ind = TRUE,
          pdf = TRUE)


i = 4
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2021
load("myResults_2021.RData")
samp <- nimOutput
summ21_slope <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

MCMCtrace(samp,
          ind = TRUE,
          pdf = TRUE)


i = 5
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

files <- list(summ17_forest, summ18_forest, summ19_forest, summ20_forest, summ21_forest)

for (i in 1:length(yrs)){
  p0[i,which(colnames(p0) %in% "forest")] <- paste(round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 1],3), 
                                                   "(",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 3],3),
                                                   "-",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 5],3),
                                                   ")",sep = "")
}

## ---- logDistcore ----

# Load and plot at the same time

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_covsp0")
pdf("Topographic.pdf")

par(mfrow = c(2,2))
plot(1, ylim = c(0.5, 5+0.5), 
     xlim = c(-1.5,2), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Forest",
     cex.axis = 0.8)
axis(2, c(1:5), labels = c("2017", "2018", "2019", "2020", "2021"), las = 2, cex.axis = 1)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_covsp0/Results_server_cyril/logDistcore")

## 2017
load("myResults_2017.RData")
samp <- nimOutput
summ17_logDistcore <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

MCMCtrace(samp,
          ind = TRUE,
          pdf = TRUE)

i = 1
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2018
load("myResults_2018.RData")
samp <- nimOutput
summ18_logDistcore <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

MCMCtrace(samp,
          ind = TRUE,
          pdf = TRUE)

i = 2
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2019
load("myResults_2019.RData")
samp <- nimOutput
summ19_logDistcore <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

MCMCtrace(samp,
          ind = TRUE,
          pdf = TRUE)

i = 3
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2020
load("myResults_2020.RData")
samp <- nimOutput
summ20_logDistcore <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

MCMCtrace(samp,
          ind = TRUE,
          pdf = TRUE)


i = 4
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2021
load("myResults_2021.RData")
samp <- nimOutput
summ21_logDistcore <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

MCMCtrace(samp,
          ind = TRUE,
          pdf = TRUE)


i = 5
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

files <- list(summ17_forest, summ18_forest, summ19_forest, summ20_forest, summ21_forest)

for (i in 1:length(yrs)){
  p0[i,which(colnames(p0) %in% "forest")] <- paste(round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 1],3), 
                                                   "(",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 3],3),
                                                   "-",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 5],3),
                                                   ")",sep = "")
}

## ---- roads1 ----

# Load and plot at the same time

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_covsp0")
pdf("Topographic.pdf")

par(mfrow = c(2,2))
plot(1, ylim = c(0.5, 5+0.5), 
     xlim = c(-1.5,2), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Forest",
     cex.axis = 0.8)
axis(2, c(1:5), labels = c("2017", "2018", "2019", "2020", "2021"), las = 2, cex.axis = 1)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_covsp0/Results_server_cyril/roads1")

## 2017
load("myResults_2017.RData")
samp <- nimOutput
summ17_roads1 <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

MCMCtrace(samp,
          ind = TRUE,
          pdf = TRUE)

i = 1
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2018
load("myResults_2018.RData")
samp <- nimOutput
summ18_roads1 <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

MCMCtrace(samp,
          ind = TRUE,
          pdf = TRUE)

i = 2
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2019
load("myResults_2019.RData")
samp <- nimOutput
summ19_roads1 <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

MCMCtrace(samp,
          ind = TRUE,
          pdf = TRUE)

i = 3
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2020
load("myResults_2020.RData")
samp <- nimOutput
summ20_roads1 <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

MCMCtrace(samp,
          ind = TRUE,
          pdf = TRUE)


i = 4
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2021
load("myResults_2021.RData")
samp <- nimOutput
summ21_roads1 <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

MCMCtrace(samp,
          ind = TRUE,
          pdf = TRUE)


i = 5
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

files <- list(summ17_forest, summ18_forest, summ19_forest, summ20_forest, summ21_forest)

for (i in 1:length(yrs)){
  p0[i,which(colnames(p0) %in% "forest")] <- paste(round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 1],3), 
                                                   "(",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 3],3),
                                                   "-",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 5],3),
                                                   ")",sep = "")
}


## -------------------------------------------------
##           SCRdenscov_singleyear
##             FinalData17-21
##           p0 = effort + trap + b / NO sig[sex]
##          I tried without sig sex because they didnt converge
## ------------------------------------------------- 

# Data frame to store values of p0 accross models and years

covs <- c("forest", "slope", "logDistcore", "roads1")

yrs <- c("2017", "2018", "2019", "2020", "2021")

p0 <- data.frame(matrix(NA, nrow = length(yrs), ncol = length(covs)))
rownames(p0) <- yrs
colnames(p0) <- covs

## ---- Forest ----

# Load and plot at the same time

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_covsp0")
pdf("Topographic.pdf")

par(mfrow = c(2,2))
plot(1, ylim = c(0.5, 5+0.5), 
     xlim = c(-1.5,2), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Forest",
     cex.axis = 0.8)
axis(2, c(1:5), labels = c("2017", "2018", "2019", "2020", "2021"), las = 2, cex.axis = 1)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_covsp0/Results_server_cyril_covsp0/forest")

## 2017
load("myResults_2017.RData")
samp <- nimOutput
summ17_forest <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

MCMCtrace(samp,
          ind = TRUE,
          pdf = TRUE)

i = 1
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2018
load("myResults_2018.RData")
samp <- nimOutput
summ18_forest <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

MCMCtrace(samp,
          ind = TRUE,
          pdf = TRUE)

i = 2
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2019
load("myResults_2019.RData")
samp <- nimOutput
summ19_forest <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

MCMCtrace(samp,
          ind = TRUE,
          pdf = TRUE)

i = 3
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2020
load("myResults_2020.RData")
samp <- nimOutput
summ20_forest <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

MCMCtrace(samp,
          ind = TRUE,
          pdf = TRUE)


i = 4
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2021
load("myResults_2021.RData")
samp <- nimOutput
summ21_forest <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

MCMCtrace(samp,
          ind = TRUE,
          pdf = TRUE)


i = 5
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

files <- list(summ17_forest, summ18_forest, summ19_forest, summ20_forest, summ21_forest)

for (i in 1:length(yrs)){
  p0[i,which(colnames(p0) %in% "forest")] <- paste(round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 1],3), 
                                                   "(",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 3],3),
                                                   "-",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 5],3),
                                                   ")",sep = "")
}

## ---- Slope ----

# Load and plot at the same time

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_covsp0")
pdf("Topographic.pdf")

par(mfrow = c(2,2))
plot(1, ylim = c(0.5, 5+0.5), 
     xlim = c(-1.5,2), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Forest",
     cex.axis = 0.8)
axis(2, c(1:5), labels = c("2017", "2018", "2019", "2020", "2021"), las = 2, cex.axis = 1)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_covsp0/Results_server_cyril_covsp0/slope")

## 2017
load("myResults_2017.RData")
samp <- nimOutput
summ17_slope <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

MCMCtrace(samp,
          ind = TRUE,
          pdf = TRUE)

i = 1
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2018
load("myResults_2018.RData")
samp <- nimOutput
summ18_slope <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

MCMCtrace(samp,
          ind = TRUE,
          pdf = TRUE)

i = 2
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2019
load("myResults_2019.RData")
samp <- nimOutput
summ19_slope <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

MCMCtrace(samp,
          ind = TRUE,
          pdf = TRUE)

i = 3
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2020
load("myResults_2020.RData")
samp <- nimOutput
summ20_slope <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

MCMCtrace(samp,
          ind = TRUE,
          pdf = TRUE)


i = 4
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2021
load("myResults_2021.RData")
samp <- nimOutput
summ21_slope <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

MCMCtrace(samp,
          ind = TRUE,
          pdf = TRUE)


i = 5
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

files <- list(summ17_forest, summ18_forest, summ19_forest, summ20_forest, summ21_forest)

for (i in 1:length(yrs)){
  p0[i,which(colnames(p0) %in% "forest")] <- paste(round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 1],3), 
                                                   "(",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 3],3),
                                                   "-",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 5],3),
                                                   ")",sep = "")
}

## ---- logDistcore ----

# Load and plot at the same time

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_covsp0/Results_server_cyril_covsp0")
pdf("Distcore.pdf")

# Load and plot at the same time
plot(1, ylim = c(0.5, 5+0.5), 
     xlim = c(-1.5, 2), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Distance to core (log)",
     cex.axis = 0.8)
axis(2, c(1:5), labels = c("2017", "2018", "2019", "2020", "2021"), las = 2, cex.axis = 1)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_covsp0/Results_server_cyril_covsp0/logDistcore")

## 2017
load("myResults_2017.RData")
samp <- nimOutput
summ17_logDistcore <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

#MCMCtrace(samp,
#          ind = TRUE,
#          pdf = FALSE)

i = 1
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2018
load("myResults_2018.RData")
samp <- nimOutput
summ18_logDistcore <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

#MCMCtrace(samp,
#          ind = TRUE,
#          pdf = FALSE)

i = 2
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2019
load("myResults_2019.RData")
samp <- nimOutput
summ19_logDistcore <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

#MCMCtrace(samp,
#          ind = TRUE,
#          pdf = FALSE)

i = 3
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2020
load("myResults_2020.RData")
samp <- nimOutput
summ20_logDistcore <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

#MCMCtrace(samp,
#          ind = TRUE,
#          pdf = FALSE)


i = 4
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2021
load("myResults_2021.RData")
samp <- nimOutput
summ21_logDistcore <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

#MCMCtrace(samp,
#          ind = TRUE,
#          pdf = FALSE)


i = 5
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

dev.off()

files <- list(summ17_forest, summ18_forest, summ19_forest, summ20_forest, summ21_forest)

for (i in 1:length(yrs)){
  p0[i,which(colnames(p0) %in% "forest")] <- paste(round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 1],3), 
                                                   "(",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 3],3),
                                                   "-",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 5],3),
                                                   ")",sep = "")}


## ---- roads1 ----

# Load and plot at the same time

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_covsp0")
pdf("Topographic.pdf")

par(mfrow = c(2,2))
plot(1, ylim = c(0.5, 5+0.5), 
     xlim = c(-1.5,2), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Forest",
     cex.axis = 0.8)
axis(2, c(1:5), labels = c("2017", "2018", "2019", "2020", "2021"), las = 2, cex.axis = 1)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_covsp0/Results_server_cyril_covsp0/roads1")

## 2017
load("myResults_2017.RData")
samp <- nimOutput
summ17_roads1 <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

MCMCtrace(samp,
          ind = TRUE,
          pdf = TRUE)

i = 1
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2018
load("myResults_2018.RData")
samp <- nimOutput
summ18_roads1 <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

MCMCtrace(samp,
          ind = TRUE,
          pdf = TRUE)

i = 2
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2019
load("myResults_2019.RData")
samp <- nimOutput
summ19_roads1 <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

MCMCtrace(samp,
          ind = TRUE,
          pdf = TRUE)

i = 3
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2020
load("myResults_2020.RData")
samp <- nimOutput
summ20_roads1 <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

MCMCtrace(samp,
          ind = TRUE,
          pdf = TRUE)


i = 4
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2021
load("myResults_2021.RData")
samp <- nimOutput
summ21_roads1 <- MCMCsummary(samp)

out <- ProcessCodaOutput(samp)

MCMCtrace(samp,
          ind = TRUE,
          pdf = TRUE)


i = 5
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

files <- list(summ17_forest, summ18_forest, summ19_forest, summ20_forest, summ21_forest)

for (i in 1:length(yrs)){
  p0[i,which(colnames(p0) %in% "forest")] <- paste(round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 1],3), 
                                                   "(",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 3],3),
                                                   "-",round(files[[i]][which(rownames(files[[i]]) %in% "p0"), 5],3),
                                                   ")",sep = "")
}


## -------------------------------------------------
##                  openSCRdenscov      
## ------------------------------------------------- 

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/2.openSCRdenscov")
load("sampOpenSCR.RData")
sampOpen <- samp2
summOpen1 <- MCMCsummary(sampOpen)
MCMCtrace(sampOpen)

write.csv(summOpen, file = "summOpen.csv")

summOpen

## -------------------------------------------------
##        openSCRdenscov + different trap arrays    
## ------------------------------------------------- 

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/2.openSCRdenscov")
load("sampOpenSCR_diftraps.RData")

##remove NAs
inn<-colnames(samp[[1]])
remm<-pmatch(c("R[1]", "pc.gam[1]"), inn)
samp2<-lapply(samp, function(x)x[,-remm])

summOpen3 <- MCMCsummary(samp2)
MCMCtrace(samp2)

## -------------------------------------------------
##        openSCRdenscov + Age structure + different trap arrays    
## ------------------------------------------------- 

library(MCMCvis)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/3.openSCRdenscov_Age")
load("sampOpenSCR_diftraps_age.RData")

summOpen1 <- MCMCsummary(samp)
#MCMCtrace(samp)

# Same model but monitoring B and Abundance each year and age category
setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/3.openSCRdenscov_Age")
load("sampOpenSCR_diftraps_age2.RData")

summOpen2 <- MCMCsummary(samp)
#MCMCtrace(samp)

## -------------------------------------------------
##        openSCRdenscov + Age structure + different trap arrays 
##                      2017 - 2021
## ------------------------------------------------- 

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/3.openSCRdenscov_Age/2021")
load("sampOpenSCR_diftraps_age_2021.RData")

out.list<- list()
out.list[[1]] <- as.mcmc(chain_output[[1]])
out.list[[2]] <- as.mcmc(chain_output[[2]])
out.list[[3]] <- as.mcmc(chain_output[[3]])

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)
out$sims.list$N

# Convergence

library(MCMCvis)

MCMCtrace(out$sims.list,
          ind = TRUE,
          pdf = TRUE, 
          open_pdf = FALSE, 
          filename = 'OpenSCR_diftraps_age_2021_Traceplot', 
          wd = 'D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/3.openSCRdenscov_Age/2021')

out$mean

summOpen1 <- MCMCsummary(out)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/3.openSCRdenscov_Age/2021")
pdf("openSCRdenscov_Age.pdf", 7, 7)


par(mfrow = c(2,2),
    oma = c(2,4,2,1),
    mar = c(3,3,2,3))

# Detection probability

params_p <- c("p.cub", "p.sub", 'p.ad')

plot(1, ylim = c(0.5, length(params_p)+0.5), 
     xlim = c(0,0.1), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Detection probability",
     cex.axis = 0.8)

axis(2, c(1:length(params_p)), labels = c("Cub(1+2)", "Subad(3+4)", "Adult(>4)"), las = 2, cex.axis = 1)

for (i in 1:length(params_p)){
  plot.violins3(list(out$sims.list[names(out$sims.list) %in% params_p[i]][[1]]),
                x = i,
                at = i,
                violin.width = 0.2,
                plot.ci = 0.95,
                col = c("purple"),
                add = T,
                alpha = 0.3,
                scale.width = FALSE,
                border.col = "black",
                horizontal = TRUE)}

# Survival

params_phi <- c("phi.cub", "phi.sub", 'phi.ad')

plot(1, ylim = c(0.5, length(params_phi)+0.5), 
     xlim = c(0.3,1), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Survival probability",
     cex.axis = 0.8)

axis(2, c(1:length(params_phi)), labels = c("Cub(1+2)", "Subad(3+4)", "Adult(>4)"), las = 2, cex.axis = 1)

for (i in 1:length(params_phi)){
  plot.violins3(list(out$sims.list[names(out$sims.list) %in% params_phi[i]][[1]]),
                x = i,
                at = i,
                violin.width = 0.2,
                plot.ci = 0.95,
                col = c("purple"),
                add = T,
                alpha = 0.3,
                scale.width = FALSE,
                border.col = "black",
                horizontal = TRUE)}

# Age structure

plot(1, ylim = c(0.5, ncol(out$sims.list$piAGE)+0.5), 
     xlim = c(0,1), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Age classes probabilities",
     cex.axis = 0.8)

axis(2, c(1:ncol(out$sims.list$piAGE)), labels = c("Cub(1)","Cub(2)","Subad(3)","Subad(4)","Adult(>4)"), las = 2, cex.axis = 1)

ncol(out$sims.list$piAGE)

for (i in 1:ncol(out$sims.list$piAGE)){
  plot.violins3(list(out$sims.list$piAGE[ ,i]),
                x = i,
                at = i,
                violin.width = 0.2,
                plot.ci = 0.95,
                col = c("purple"),
                add = T,
                alpha = 0.3,
                scale.width = FALSE,
                border.col = "black",
                horizontal = TRUE)}

# Abundance

plot(1, ylim = c(0.5, ncol(out$sims.list$N)+0.5), 
     xlim = c(30,max(out$sims.list$N)), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Abundance",
     cex.axis = 0.8)

axis(2, c(1:ncol(out$sims.list$N)), labels = c(2017:2021), las = 2, cex.axis = 1)


for (i in 1:ncol(out$sims.list$N)){
  plot.violins3(list(out$sims.list$N[ ,i]),
                x = i,
                at = i,
                violin.width = 0.4,
                plot.ci = 0.95,
                col = c("purple"),
                add = T,
                alpha = 0.3,
                scale.width = FALSE,
                border.col = "black",
                horizontal = TRUE)}
dev.off()

# Abundance

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/3.openSCRdenscov_Age/2021")
pdf("openSCRdenscov_Age_Abundance.pdf", 5, 5)

plot(1, ylim = c(30,max(out$sims.list$N)) , 
     xlim = c(0.5, ncol(out$sims.list$N)+0.5), 
     type ="n", yaxt="n", xaxt="n",
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Abundance",
     cex.axis = 0.8)

axis(2, at = seq(30,max(out$sims.list$N), by = 20), labels = seq(30,max(out$sims.list$N), by = 20),las = 2, cex.axis = 1)
axis(1, at = c(1:5), labels = 2017:2021,las = 2, cex.axis = 1)


for (i in 1:ncol(out$sims.list$N)){
  plot.violins3(list(out$sims.list$N[ ,i]),
                x = i,
                at = i,
                violin.width = 0.4,
                plot.ci = 0.95,
                col = c("purple"),
                add = T,
                alpha = 0.3,
                scale.width = FALSE,
                border.col = "black",
                horizontal = FALSE)}
dev.off()


## -------------------------------------------------
##        DistCore ALL (single models + open model) 
##                      2017 - 2021
## ------------------------------------------------- 


setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/3.openSCRdenscov_Age/2021")
pdf("Distcore.pdf")

# Load and plot at the same time
plot(1, ylim = c(0.5, 6 +0.5), 
     xlim = c(-1.5, 2), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Distance to core (log)",
     cex.axis = 0.8)
axis(2, c(1:6), labels = c("2017", "2018", "2019", "2020", "2021", "OPSCR"), las = 2, cex.axis = 1)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/1.SCRdenscov_year/FinalData_17_19/logDistcore")

## 2017
load("sampSCRdenscov2017.RData")
summ17_distcore <- MCMCsummary(samp)
#trace17 <- MCMCtrace(samp17, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 1
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2018
load("sampSCRdenscov2018.RData")
summ18_distcore <- MCMCsummary(samp)
#trace18 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 2
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2019
load("sampSCRdenscov2019.RData")
summ19_distcore <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 3
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2020
load("sampSCRdenscov2020.RData")
summ20_distcore <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 4
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

## 2021
load("sampSCRdenscov2021.RData")
summ21_distcore <- MCMCsummary(samp)
#trace19 <- MCMCtrace(samp, pdf = FALSE)

out.list<- list()
out.list[[1]] <- as.mcmc(samp$chain1)
out.list[[2]] <- as.mcmc(samp$chain2)
out.list[[3]] <- as.mcmc(samp$chain3)

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 5
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)

# OpenSCR

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/3.openSCRdenscov_Age/2021")
load("sampOpenSCR_diftraps_age_2021.RData")

out.list<- list()
out.list[[1]] <- as.mcmc(chain_output[[1]])
out.list[[2]] <- as.mcmc(chain_output[[2]])
out.list[[3]] <- as.mcmc(chain_output[[3]])

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 6
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("green"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "darkgreen",
              horizontal = TRUE)


dev.off()


## -------------------------------------------------
##                  DistCore open model
##                      2017 - 2021
## ------------------------------------------------- 

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/3.openSCRdenscov_Age/2021")
pdf("Distcore_open.pdf")

par(mfrow = c(2,2),
    oma = c(2,4,2,1),
    mar = c(3,3,2,3))

# Load and plot at the same time
plot(1, ylim = c(0.5, 1 +0.5), 
     xlim = c(-1, 1), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = "Beta coefficient", ylab = "", main = "Distance to core (log)",
     cex.axis = 0.8)
axis(2, c(1:1), labels = c(" "), las = 2, cex.axis = 1)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/3.openSCRdenscov_Age/2021")
load("sampOpenSCR_diftraps_age_2021.RData")

out.list<- list()
out.list[[1]] <- as.mcmc(chain_output[[1]])
out.list[[2]] <- as.mcmc(chain_output[[2]])
out.list[[3]] <- as.mcmc(chain_output[[3]])

out.list <- as.mcmc.list(out.list)

out <- ProcessCodaOutput(out.list)

i = 1
plot.violins3(list(out$sims.list[names(out$sims.list) %in% "beta.dens"][[1]]),
              x = i,
              at = i,
              violin.width = 0.2,
              plot.ci = 0.95,
              col = c("purple"),
              add = T,
              alpha = 0.3,
              scale.width = FALSE,
              border.col = "black",
              horizontal = TRUE)

dev.off()


## -------------------------------------------------
##        openSCRdenscov + Age structure + different trap arrays 
##          2017 - 2021 FINAL DATA (Cross check Elena)
## ------------------------------------------------- 

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021")
load("sampOpenSCR_diftraps_age_final1721.RData")

out.list<- list()
out.list[[1]] <- as.mcmc(chain_output[[1]])
out.list[[2]] <- as.mcmc(chain_output[[2]])
out.list[[3]] <- as.mcmc(chain_output[[3]])

out.list <- as.mcmc.list(out.list)
out <- ProcessCodaOutput(out.list)

summOpen1 <- MCMCsummary(out.list)


# Convergence

library(MCMCvis)

MCMCtrace(out$sims.list,
          ind = TRUE,
          pdf = TRUE, 
          open_pdf = FALSE, 
          filename = 'OpenSCR_diftraps_age_2021_Traceplot', 
          wd = 'D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/3.openSCRdenscov_Age/2021')

out$mean


setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/3.openSCRdenscov_Age/2021")
pdf("openSCRdenscov_Age_final1721.pdf", 7, 7)


par(mfrow = c(2,2),
    oma = c(2,4,2,1),
    mar = c(3,3,2,3))

# Detection probability

params_p <- c("p.cub", "p.sub", 'p.ad')

plot(1, ylim = c(0.5, length(params_p)+0.5), 
     xlim = c(0,0.1), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Detection probability",
     cex.axis = 0.8)

axis(2, c(1:length(params_p)), labels = c("Cub(1+2)", "Subad(3+4)", "Adult(>4)"), las = 2, cex.axis = 1)

for (i in 1:length(params_p)){
  plot.violins3(list(out$sims.list[names(out$sims.list) %in% params_p[i]][[1]]),
                x = i,
                at = i,
                violin.width = 0.2,
                plot.ci = 0.95,
                col = c("purple"),
                add = T,
                alpha = 0.3,
                scale.width = FALSE,
                border.col = "black",
                horizontal = TRUE)}

# Survival

params_phi <- c("phi.cub", "phi.sub", 'phi.ad')

plot(1, ylim = c(0.5, length(params_phi)+0.5), 
     xlim = c(0.3,1), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Survival probability",
     cex.axis = 0.8)

axis(2, c(1:length(params_phi)), labels = c("Cub(1+2)", "Subad(3+4)", "Adult(>4)"), las = 2, cex.axis = 1)

for (i in 1:length(params_phi)){
  plot.violins3(list(out$sims.list[names(out$sims.list) %in% params_phi[i]][[1]]),
                x = i,
                at = i,
                violin.width = 0.2,
                plot.ci = 0.95,
                col = c("purple"),
                add = T,
                alpha = 0.3,
                scale.width = FALSE,
                border.col = "black",
                horizontal = TRUE)}

# Age structure

plot(1, ylim = c(0.5, ncol(out$sims.list$piAGE)+0.5), 
     xlim = c(0,1), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Age classes probabilities",
     cex.axis = 0.8)

axis(2, c(1:ncol(out$sims.list$piAGE)), labels = c("Cub(1)","Cub(2)","Subad(3)","Subad(4)","Adult(>4)"), las = 2, cex.axis = 1)

ncol(out$sims.list$piAGE)

for (i in 1:ncol(out$sims.list$piAGE)){
  plot.violins3(list(out$sims.list$piAGE[ ,i]),
                x = i,
                at = i,
                violin.width = 0.2,
                plot.ci = 0.95,
                col = c("purple"),
                add = T,
                alpha = 0.3,
                scale.width = FALSE,
                border.col = "black",
                horizontal = TRUE)}

# Abundance

plot(1, ylim = c(0.5, ncol(out$sims.list$N)+0.5), 
     xlim = c(30,max(out$sims.list$N)), 
     type ="n", yaxt="n", 
     #xaxt="n", 
     xlab = " ", ylab = "", main = "Abundance",
     cex.axis = 0.8)

axis(2, c(1:ncol(out$sims.list$N)), labels = c(2017:2021), las = 2, cex.axis = 1)


for (i in 1:ncol(out$sims.list$N)){
  plot.violins3(list(out$sims.list$N[ ,i]),
                x = i,
                at = i,
                violin.width = 0.4,
                plot.ci = 0.95,
                col = c("purple"),
                add = T,
                alpha = 0.3,
                scale.width = FALSE,
                border.col = "black",
                horizontal = TRUE)}
dev.off()

## -------------------------------------------------
##        openSCRdenscov + Age structure + different trap arrays 
##                EFFORT Covariates
##               2017 - 2021 FINAL DATA 
##                 Time: 3.31 days
## ------------------------------------------------- 

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Run_Data/Nimble/Results/3.openSCRdenscov_Age/2021")
load("Results_Model3-2.RData")

out.list<- list()
out.list[[1]] <- as.mcmc(chain_output[[1]])
out.list[[2]] <- as.mcmc(chain_output[[2]])
out.list[[3]] <- as.mcmc(chain_output[[3]])

out.list <- as.mcmc.list(out.list)

# There are NA and the function ProcessCodaOutput doesnt work.
# I need to substitute them as deleting the columns with NA don't work

na <- function(x){which(!complete.cases(x))}
lapply(out.list, na)
lapply(out.list, function(x){print(x[1,])}) # Its pc-gam[1] and R[1], set to 0

for (i in 1:3){
  out.list[[i]][1,'R[1]'] <- 0
  out.list[[i]][1,'pc.gam[1]'] <- 0
}

out2 <- ProcessCodaOutput(out.list)

MCMCtrace(out.list,   
          ind = TRUE,
          pdf = TRUE)

summ_EffortTrapBh_sigSex_age <- MCMCsummary(out.list)

