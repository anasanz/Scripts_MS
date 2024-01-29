
## -------------------------------------------------
##        PRELIMINARY RESULTS STTUR + COPAL
## ------------------------------------------------- 

## Update on status of study:
# - We ran species trend models in 56 transects
#   -> For STTUR, 46 models converged
#   -> For COPAL, 53 models converged

rm(list=ls())
 
library(dplyr)

# Load analyzed transects (all species)
setwd("D:/Otros/Tórtola/Results/Study2/Model_results_1.1/STTUR")
load("1.1TortoData_transects_0221.RData") 

## ---- 1. STTUR ----

## ---- 1.1. Identify the files that I will keep for the results (files that converged) ----

load("1.1TortoData_noconverge_2002_2021.RData") # Did not converge with 5e5 iter
converge <- transect[-which(transect %in% no_converge)] # Transects that converged with 5e5 iter
converge1_files <-paste("1.1TortoData_0221_", converge,".RData", sep = "") 

load("1.1TortoData_noconverge2_9e5iter_2002_2021.RData") # Did not converge with 9e5 iter
no_converge2 <- c(no_converge2,"169") # Añadir 169 que no consigo correr, error
converge2 <- transect[-which(transect %in% c(converge, no_converge2))] # Transects that converged with 9e5 iter (removing the ones added already (converge))
converge2_files <-paste("1.1TortoData_0221_", converge2, "_9e5iter", ".RData", sep = "")

load("1.1TortoData_transects_0221_ADD.RData") # ADDED transects (less restrictive)
load("1.1TortoData_noconverge3(added)_9e5iter_2002_2021.RData") # Did not converge with 9e5 iter
converge3 <- transect_add[-which(transect_add %in% c(no_converge3))]
converge3_files <-paste("1.1TortoData_0221_", converge3, "_9e5iter", ".RData", sep = "")

converge_all <- c(converge, converge2, converge3)
files <- c(converge1_files, converge2_files, converge3_files)

converged_files <- data.frame(transect_ID = converge_all, file = files)

## Load the ones that converged properly from rerun, to substitute them (I will keep these as priority)
load("1.1TortoData_transects_0221_RERUN.RData") # Rerun transects
load("1.1TortoData_noconverge4_15e5iter_2002_2021.RData") # Did not converge even rerunning
no_converge4 <- c(no_converge4,"169") # Añadir 169 que no consigo correr, error
converge4 <- rerun[-which(rerun %in% no_converge4)] # Transects that converged with 15e5 iter
converge4_files <-paste("1.1TortoData_0221_", converge4, "_15e5iter", ".RData", sep = "")

converged_files$file[which(converged_files$transect_ID %in% converge4)] <- converge4_files # Substitute, keep the ones that converged better with more iter

converged_files2 <- converged_files # To add info of trend and estimate together

## ---- Probability of increase in trend ----

converged_files$p_increasing <- NA

for (i in 1:nrow(converged_files)){
  load(converged_files$file[i])
  df.outall <- as.data.frame(out$sims.list)
  total.samples <- nrow(df.outall)
  increasing <- df.outall$bYear.lam[which(df.outall$bYear.lam > 0)]
  prob_increasing <- length(increasing)/total.samples
  converged_files$p_increasing[i] <- prob_increasing
}

converged_files$p_increasing <- round(converged_files$p_increasing,4)
converged_files <- arrange(converged_files,p_increasing)

#write.csv(converged_files, file = "pInc_transects_converged.csv")

## ---- Join info of trend and estimate together for data exploration ----
 
converged_files2$p_increasing <- NA
converged_files2$b_estimate <- NA
converged_files2$sd <- NA

for (i in 1:nrow(converged_files2)){
  load(converged_files2$file[i])
  df.outall <- as.data.frame(out$sims.list)
  total.samples <- nrow(df.outall)
  # Prob increasing
  increasing <- df.outall$bYear.lam[which(df.outall$bYear.lam > 0)]
  prob_increasing <- length(increasing)/total.samples
  converged_files2$p_increasing[i] <- prob_increasing
  # Estimate
  converged_files2$b_estimate[i] <- mean(df.outall$bYear.lam)
  converged_files2$sd[i] <- sd(df.outall$bYear.lam)
  }

converged_files2$p_increasing <- round(converged_files2$p_increasing,4)
converged_files2$b_estimate <- round(converged_files2$b_estimate,4)
converged_files2$sd <- round(converged_files2$sd,4)

converged_files2 <- arrange(converged_files2,p_increasing)

write.csv(converged_files2, file = "pInc_bEst_sd_transects_converged.csv")


## ---- COPAL ----

# List of file names to load (depending on how many iter)

setwd("C:/Users/anasa/OneDrive/deepthought/Results/Otros/Tortola/Study2/Model1.1/Copal_2002_2021")
load("1.1CopalData_noconverge_2002_2021.RData") # Did not converge with 5e5 iter
converge <- transect[-which(transect %in% no_converge)] # Transects that converged with 5e5 iter 
converge1_files <-paste("1.1CopalData_0221_", converge,".RData", sep = "") 

setwd("C:/Users/anasa/OneDrive/deepthought/Results/Otros/Tortola/Study2/Model1.1/Copal_2002_2021")
load("1.1CopalData_noconverge2_9e5iter_2002_2021.RData")
no_converge2 <- c(no_converge2,"169") # Añadir 169 que no consigo correr, error
converge2 <- transect[-which(transect %in% c(converge, no_converge2))] # Transects that converged with 9e5 iter (removing the ones added already (converge))
converge2_files <-paste("1.1CopalData_0221_", converge2, "_9e5iter", ".RData", sep = "")

converge_all <- c(converge, converge2)
files <- c(converge1_files, converge2_files)

# From this I would need to see other ways to reach convergence, increasing iter wouln't work

transectCOPAL <- transect[-which(transect %in% no_converge2)] 

## VALUES PROBABILITY OF INCREASE IN TREND

trendCOPAL <- data.frame(Site = converge_all, p_increasing_COPAL = NA, name_file = NA)

setwd("D:/Deepthought/Results/Otros/Tortola/Study2/Model1.1/2002_2021/COPAL")

for (i in 1:length(files)){
  load(files[i])
  df.outall <- as.data.frame(out$sims.list)
  total.samples <- nrow(df.outall)
  increasing <- df.outall$bYear.lam[which(df.outall$bYear.lam > 0)]
  prob_increasing <- length(increasing)/total.samples
  trendCOPAL[i,2] <- prob_increasing
  trendCOPAL[i,3] <- files[i]
}

trendCOPAL$p_increasing_COPAL <- round(trendCOPAL$p_increasing_COPAL,3)
trendCOPAL$label_COPAL <- paste("PI =", trendCOPAL$p_increasing_COPAL)
trendCOPAL <- arrange(trendCOPAL,p_increasing_COPAL)

setwd("D:/Otros/Tórtola/Results/Study2/Model_results")
#write.csv(trendCOPAL,"COPAL_trend1.1.csv")

## ---- Join ----

transectSTTUR[which(transectSTTUR %in% transectCOPAL == FALSE)] 
  # Transectos de STTUR que no están en COPAL. Super prioritario que converjan!!
  # De momento, los quito de la lista para unir las dos tablas

trends <- left_join(trendCOPAL, trendSTTUR, by = "Site")
trends <- trends[complete.cases(trends), ]
setwd("D:/Otros/Tórtola/Results/Study2/Model_results")
#write.csv(trends,"trends_COPAL_STTUR_1.1.csv")

## ---- Preliminary exploration ----

cor(trends$p_increasing_COPAL, trends$p_increasing_STTUR) 
m <- lm(p_increasing_STTUR ~ p_increasing_COPAL, data = trends)
summary(m)
# There is no correlation, no effects in linear model

plot(trends$p_increasing_STTUR ~ trends$p_increasing_COPAL, ylab = "Prob. increasing STTUR", xlab = "Prob. increasing COPAL", pch = 19)
abline(m, col = "red")

