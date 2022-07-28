
## -------------------------------------------------
##        PRELIMINARY RESULTS STTUR + COPAL
## ------------------------------------------------- 

## Update on status of study:
# - We ran species trend models in 56 transects
#   -> For STTUR, 46 models converged
#   -> For COPAL, X models converged

rm(list=ls())
 
library(dplyr)

# Load analyzed transects (all species)
setwd("C:/Users/anasa/OneDrive/deepthought/Results/Otros/Tortola/Study2/Model1.1/2002_2021")
load("1.1TortoData_transects_0221.RData") 

## ---- STTUR ----

# List of file names to load (depending on how many iter)

setwd("D:/Deepthought/Results/Otros/Tortola/Study2/Model1.1/2002_2021/STTUR")
load("1.1TortoData_noconverge_2002_2021.RData") # Did not converge with 5e5 iter
converge <- transect[-which(transect %in% no_converge)] # Transects that converged with 5e5 iter 
converge1_files <-paste("1.1TortoData_0221_", converge,".RData", sep = "") 

setwd("D:/Deepthought/Results/Otros/Tortola/Study2/Model1.1/2002_2021/STTUR")
load("1.1TortoData_noconverge2_9e5iter_2002_2021.RData") # Did not converge with 9e5 iter
no_converge2 <- c(no_converge2,"169") # Añadir 169 que no consigo correr, error
converge2 <- transect[-which(transect %in% c(converge, no_converge2))] # Transects that converged with 9e5 iter (removing the ones added already (converge))
converge2_files <-paste("1.1TortoData_0221_", converge2, "9e5iter", ".RData", sep = "")

converge_all <- c(converge, converge2)
files <- c(converge1_files, converge2_files)

# From this, except the transecta 14 and 22 ([1,2]), which might work by increasing iterations,
# I would have to see other ways to reach convergence

transectSTTUR <- transect[-which(transect %in% no_converge2)] 

## VALUES PROBABILITY OF INCREASE IN TREND

trendSTTUR <- data.frame(Site = converge_all, p_increasing_STTUR = NA, name_file = NA)

setwd("D:/Deepthought/Results/Otros/Tortola/Study2/Model1.1/2002_2021/STTUR")

for (i in 1:length(files)){
  load(files[i])
  df.outall <- as.data.frame(out$sims.list)
  total.samples <- nrow(df.outall)
  increasing <- df.outall$bYear.lam[which(df.outall$bYear.lam > 0)]
  prob_increasing <- length(increasing)/total.samples
  trendSTTUR[i,2] <- prob_increasing
  trendSTTUR[i,3] <- files[i]
}

trendSTTUR$p_increasing_STTUR <- round(trendSTTUR$p_increasing_STTUR,3)
trendSTTUR$label_STTUR <- paste("PI =", trendSTTUR$p_increasing_STTUR)
trendSTTUR <- arrange(trendSTTUR,p_increasing_STTUR)

setwd("D:/Otros/Tórtola/Results/Study2/Model_results")
write.csv(trendSTTUR,"STTUR_trend1.1.csv")

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
write.csv(trendCOPAL,"COPAL_trend1.1.csv")

## ---- Join ----

transectSTTUR[which(transectSTTUR %in% transectCOPAL == FALSE)] 
  # Transectos de STTUR que no están en COPAL. Super prioritario que converjan!!
  # De momento, los quito de la lista para unir las dos tablas

trends <- left_join(trendCOPAL, trendSTTUR, by = "Site")
trends <- trends[complete.cases(trends), ]

## ---- Preliminary exploration ----

cor(trends$p_increasing_COPAL, trends$p_increasing_STTUR) 
m <- lm(p_increasing_STTUR ~ p_increasing_COPAL, data = trends)
summary(m)
# There is no correlation, no effects in linear model

plot(trends$p_increasing_STTUR ~ trends$p_increasing_COPAL, ylab = "Prob. increasing STTUR", xlab = "Prob. increasing COPAL", pch = 19)
abline(m, col = "red")

