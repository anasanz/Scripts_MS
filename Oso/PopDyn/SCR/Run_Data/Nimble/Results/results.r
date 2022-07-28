
## -------------------------------------------------
##          Check results different models
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


## -------------------------------------------------
##                  openSCRdenscov      
## ------------------------------------------------- 

setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Oso/PopDyn/SCR/Data/Nimble/Results/openSCRdenscov")
load("sampOpenSCR.RData")
sampOpen <- samp2
summOpen <- MCMCsummary(sampOpen)
MCMCtrace(sampOpen)

write.csv(summOpen, file = "summOpen.csv")

summOpen

