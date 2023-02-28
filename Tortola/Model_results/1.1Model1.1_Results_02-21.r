
rm(list=ls())
library(dplyr)

source("D:/PhD/MyScripts_PhD/Ch. 2-3/Ch. 3/Results/Functions/plot.violins3.r")
source("D:/PhD/MyScripts_PhD/Ch. 2-3/Ch. 3/Results/Functions/DoScale.r")

setwd("D:/Otros/Tórtola/Results/Study2/Model_results_1.1/STTUR")
load("1.1TortoData_transects_0221.RData") # Load analyzed transects

## ---- Check convergence of structural parameters ----

# With with 5e5 iter

no_converge <- NULL
no_converge_1.1Rhat <- NULL

setwd("D:/Otros/Tórtola/Results/Study2/Model_results_1.1/STTUR")
for (i in 1:length(transect)){
  load(paste("1.1TortoData_0221_", transect[i],".RData", sep = ""))
  
  sum <-as.data.frame(out$summary[rownames(out$summary) 
                                  %in% c("mu.lam.year", "sig.lam.year", "bYear.lam", "mu.sig.year", "sig.sig.year", "brough.sig", "sd", "rho",'Bp.Obs', 'Bp.N'), ]) # Only structural parameters
  if (sum(ifelse(sum$Rhat > 1.30, 1, 0)) > 0 )
    no_converge <- c(no_converge, transect[i])
  
  if (sum(ifelse(sum$Rhat > 1.15, 1, 0)) > 0 )
    no_converge_1.1Rhat <- c(no_converge_1.1Rhat, transect[i])
  next 
}

rerun1 <- no_converge_1.1Rhat[which(no_converge_1.1Rhat %in% no_converge == FALSE)] # This I need to re-run with more iters to really reach convergence

#setwd("D:/Otros/Tórtola/Results/Study2/Model_results_1.1/STTUR")
#save(no_converge, file = "1.1TortoData_noconverge_2002_2021.RData") # Save to run at 500000

# --> Check how bad is convergence of 5e5 iter

setwd("D:/Otros/Tórtola/Results/Study2/Model_results_1.1/STTUR")
load("1.1TortoData_noconverge_2002_2021.RData")
list_no_converge <- list()

for (i in 1:length(no_converge)){
  load(paste("1.1TortoData_0221_", no_converge[i],".RData", sep = ""))
  
  list_no_converge[[i]] <-as.data.frame(out$summary[rownames(out$summary) 
                                                    %in% c("mu.lam.year", "sig.lam.year", "bYear.lam", "mu.sig.year", "sig.sig.year", "brough.sig", "sd", "rho",'Bp.Obs', 'Bp.N'), ]) # Only structural parameters
}

list_no_converge

# With 9e5 iter
setwd("D:/Otros/Tórtola/Results/Study2/Model_results_1.1/STTUR")
load("1.1TortoData_noconverge_2002_2021.RData")

no_converge <- no_converge[-10] # Transecto 10 (num 169) no esta! No se ha corrido bien

no_converge2 <- NULL
no_converge2_1.1Rhat <- NULL


for (i in 1:length(no_converge)){
  load(paste("1.1TortoData_0221_", no_converge[i], "9e5iter", ".RData", sep = ""))
  
  sum <-as.data.frame(out$summary[rownames(out$summary) 
                                  %in% c("mu.lam.year", "sig.lam.year", "bYear.lam", "mu.sig.year", "sig.sig.year", "brough.sig", "sd", "rho",'Bp.Obs', 'Bp.N'), ]) # Only structural parameters
  if (sum(ifelse(sum$Rhat > 1.30, 1, 0)) > 0 )
    no_converge2 <- c(no_converge2, no_converge[i])
  
  if (sum(ifelse(sum$Rhat > 1.15, 1, 0)) > 0 )
    no_converge2_1.1Rhat <- c(no_converge2_1.1Rhat, no_converge[i])
  next 
}

rerun2 <- no_converge2_1.1Rhat[which(no_converge2_1.1Rhat %in% no_converge2 == FALSE)] # This I need to re-run with more iters to really reach convergence
rerun2 <- c(rerun2, "169")

#setwd("D:/Otros/Tórtola/Results/Study2/Model_results_1.1/STTUR")
#save(no_converge2, file = "1.1TortoData_noconverge2_9e5iter_2002_2021.RData") # Save to run at 500000

# --> Check how bad is convergence of 9e5 iter

setwd("D:/Otros/Tórtola/Results/Study2/Model_results_1.1/STTUR")
load("1.1TortoData_noconverge2_9e5iter_2002_2021.RData")

list_no_converge2 <- list()

for (i in 1:length(no_converge2)){
  load(paste("1.1TortoData_0221_", no_converge2[i], "9e5iter", ".RData", sep = ""))
  
  list_no_converge2[[i]] <-as.data.frame(out$summary[rownames(out$summary) 
                                                    %in% c("mu.lam.year", "sig.lam.year", "bYear.lam", "mu.sig.year", "sig.sig.year", "brough.sig", "sd", "rho",'Bp.Obs', 'Bp.N'), ]) # Only structural parameters
}


## ---- Probability of increase in trend ----

# Create vector with name of files to be loaded in results (5e5 iter)

#setwd("S:/Results/Otros/Tortola/Study2/Model1.1")
#load("1.1TortoData_noconverge.RData")

converge <- transect[-which(transect %in% no_converge)] 
setwd("D:/Deepthought/Results/Otros/Tortola/Study2/Model1.1/2002_2021")
files <-paste("1.1TortoData_0221_", converge,".RData", sep = "") # remove the ones that did not converge

trend <- data.frame(Site = converge, p_increasing = NA, name_file = NA)

for (i in 1:length(files)){
  load(files[i])
  df.outall <- as.data.frame(out$sims.list)
  total.samples <- nrow(df.outall)
  increasing <- df.outall$bYear.lam[which(df.outall$bYear.lam > 0)]
  prob_increasing <- length(increasing)/total.samples
  trend[i,2] <- prob_increasing
  trend[i,3] <- files[i]
}

trend$p_increasing <- round(trend$p_increasing,3)
trend$label <- paste("PI =", trend$p_increasing)
trend <- arrange(trend,p_increasing)

setwd("D:/Otros/Tórtola/Results/Study2/Model_results")
write.csv(trend,"trend1.1_0221.csv")

## ---- Plot coefficients and probabilities ----

files_order <- trend$name_file

setwd("D:/Otros/Tórtola/Results/Study2/Plots")

pdf("trend_transects_model1.1_0221.pdf", width = 9, height = 7)

par(mfrow = c(1,2))

# Plot trends
plot(10, ylim = c(1, length(files_order)), 
     xlim = c(-5,6), 
     type ="n", yaxt="n", xlab = " ", ylab = "", main = "Population trend")
axis(2, c(1:length(files_order)), labels = trend$Site, las = 2, cex.axis = 0.7)
mtext("Beta coefficient", side = 1, line = 2.5)

setwd("D:/Deepthought/Results/Otros/Tortola/Study2/Model1.1/2002_2021")
for(i in 1:length(files_order)){
  load(files_order[i])
  
  plot.violins3(list(out$sims.list$bYear.lam),
                x = i,
                at = i,
                violin.width = 0.3,
                plot.ci = 0.95,
                col = "grey52",
                add = T,
                alpha = 0.3,
                scale.width = FALSE,
                border.col = "grey52",
                horizontal = TRUE) 
}
segments(0, -1, 0, 37.4, col = "red")
text(5, sort(seq(1:length(files_order))), labels = trend$label, cex = 0.7)

# Plot roughness
plot(10, ylim = c(1, length(files_order)), 
     xlim = c(-7,7), 
     type ="n", yaxt="n", xlab = " ", ylab = "", main = "Roughness(detection)")

axis(2, c(1:length(files_order)), labels = trend$Site, las = 2, cex.axis = 0.7)
mtext("Beta coefficient", side = 1, line = 2.5)

setwd("D:/Deepthought/Results/Otros/Tortola/Study2/Model1.1/2002_2021")
for(i in 1:length(files_order)){
  load(files_order[i])
  
  plot.violins3(list(out$sims.list$brough.sig),
                x = i,
                at = i,
                violin.width = 0.3,
                plot.ci = 0.95,
                col = "grey52",
                add = T,
                alpha = 0.3,
                scale.width = FALSE,
                border.col = "grey52",
                horizontal = TRUE) 
  
  # For indicating significant effects
  over0<- as.numeric(data.table::between(0, out$q2.5$brough.sig, out$q97.5$brough.sig))
  sig <- ifelse(over0 == 1, " ", "*")
  text(6, i, labels = sig, cex = 2)
  
}
segments(0, -1, 0, 37.4, col = "red")



dev.off()

# Roughness tiene efecto en la detectabilidad de 11/36 transectos: 8 efecto negativo, 3 efecto positivo?

## ---- Added transects, less restrictive: cHECK CONVERGENCE AND ADD ----

# With 9e5 iter

setwd("D:/Otros/Tórtola/Results/Study2/Model_results_1.1/STTUR")
load("1.1TortoData_transects_0221_ADD.RData") 

no_converge3 <- NULL
no_converge3_1.1Rhat <- NULL


setwd("D:/Otros/Tórtola/Results/Study2/Model_results_1.1/Models_added_transects")
for (i in 1:length(transect_add)){
  load(paste("1.1TortoData_0221_", transect_add[i],"_9e5iter.RData", sep = ""))
  
  sum <-as.data.frame(out$summary[rownames(out$summary) 
                                  %in% c("mu.lam.year", "sig.lam.year", "bYear.lam", "mu.sig.year", "sig.sig.year", "brough.sig", "sd", "rho",'Bp.Obs', 'Bp.N'), ]) # Only structural parameters
  if (sum(ifelse(sum$Rhat > 1.30, 1, 0)) > 0 )
    no_converge3 <- c(no_converge3, transect_add[i])
  
  if (sum(ifelse(sum$Rhat > 1.15, 1, 0)) > 0 )
    no_converge3_1.1Rhat <- c(no_converge3_1.1Rhat, transect_add[i])
  next 
}

rerun3 <- no_converge3_1.1Rhat[which(no_converge3_1.1Rhat %in% no_converge3 == FALSE)] # This I need to re-run with more iters to really reach convergence

setwd("D:/Otros/Tórtola/Results/Study2/Model_results_1.1/STTUR")
#save(no_converge3, file = "1.1TortoData_noconverge3(added)_9e5iter_2002_2021.RData") # Save to run at 500000

## Compile all transects to re-run with more iterations
rerun <- c(rerun1,rerun2,rerun3)
rerun <- c(rerun, "14", "22") # Add this 2 as Rhat doesnt arrive to 2, just to see if there is luck

setwd("D:/Otros/Tórtola/Results/Study2/Model_results_1.1/STTUR")
#save(rerun, file = "1.1TortoData_transects_0221_RERUN.RData") 

## ---- RE-RUN TRANSECTS WITH MORE ITERATIONS ----

setwd("D:/Otros/Tórtola/Results/Study2/Model_results_1.1/STTUR")
load("1.1TortoData_transects_0221_RERUN.RData")

no_converge4 <- NULL
no_converge4_1.1Rhat <- NULL

rerun <- rerun[-which(rerun == "169")]
setwd("D:/Otros/Tórtola/Results/Study2/Model_results_1.1/Models_rerun")
for (i in 1:length(rerun)){
  load(paste("1.1TortoData_0221_", rerun[i],"_15e5iter.RData", sep = ""))
  
  sum <-as.data.frame(out$summary[rownames(out$summary) 
                                  %in% c("mu.lam.year", "sig.lam.year", "bYear.lam", "mu.sig.year", "sig.sig.year", "brough.sig", "sd", "rho",'Bp.Obs', 'Bp.N'), ]) # Only structural parameters
  if (sum(ifelse(sum$Rhat > 1.30, 1, 0)) > 0 )
    no_converge4 <- c(no_converge4, rerun[i])
  
  if (sum(ifelse(sum$Rhat > 1.15, 1, 0)) > 0 )
    no_converge4_1.1Rhat <- c(no_converge4_1.1Rhat, rerun[i])
  next 
}

rerun4 <- no_converge4_1.1Rhat[which(no_converge4_1.1Rhat %in% no_converge4 == FALSE)]

setwd("D:/Otros/Tórtola/Results/Study2/Model_results_1.1/STTUR")
save(no_converge4, file = "1.1TortoData_noconverge4_15e5iter_2002_2021.RData") # Save to run at 500000



