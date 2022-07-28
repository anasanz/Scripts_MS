
rm(list=ls())
library(dplyr)

source("D:/PhD/MyScripts_PhD/Ch. 2-3/Ch. 3/Results/Functions/plot.violins3.r")
source("D:/PhD/MyScripts_PhD/Ch. 2-3/Ch. 3/Results/Functions/DoScale.r")

setwd("C:/Users/anasa/OneDrive/deepthought/Results/Otros/Tortola/Study2/Model1.1/2002_2021")
load("1.1TortoData_transects_0221.RData") # Load analyzed transects

## ---- Check convergence of structural parameters ----

# With with 5e5 iter

no_converge <- NULL

setwd("D:/Deepthought/Results/Otros/Tortola/Study2/Model1.1/2002_2021/COPAL")
for (i in 1:length(transect)){
  load(paste("1.1CopalData_0221_", transect[i],".RData", sep = ""))
  
  sum <-as.data.frame(out$summary[rownames(out$summary) 
                                  %in% c("mu.lam.year", "sig.lam.year", "bYear.lam", "mu.sig.year", "sig.sig.year", "brough.sig", "sd", "rho",'Bp.Obs', 'Bp.N'), ]) # Only structural parameters
  
  sum$Rhat[is.na(sum$Rhat)] <- 1 # Just for one observation, the rest converge
  if (sum(ifelse(sum$Rhat > 1.30, 1, 0)) > 0 )
    no_converge <- c(no_converge, transect[i])
  next 
}

setwd("C:/Users/anasa/OneDrive/deepthought/Results/Otros/Tortola/Study2/Model1.1/Copal_2002_2021")
save(no_converge, file = "1.1CopalData_noconverge_2002_2021.RData") # Save to run at 500000

# --> Check how bad is convergence of 5e5 iter

setwd("C:/Users/anasa/OneDrive/deepthought/Results/Otros/Tortola/Study2/Model1.1/Copal_2002_2021")
load("1.1CopalData_noconverge_2002_2021.RData")
list_no_converge <- list()

for (i in 1:length(no_converge)){
  load(paste("1.1CopalData_0221_", no_converge[i],".RData", sep = ""))
  
  list_no_converge[[i]] <-as.data.frame(out$summary[rownames(out$summary) 
                                                    %in% c("mu.lam.year", "sig.lam.year", "bYear.lam", "mu.sig.year", "sig.sig.year", "brough.sig", "sd", "rho",'Bp.Obs', 'Bp.N'), ]) # Only structural parameters
}

list_no_converge

# With 9e5 iter
setwd("C:/Users/anasa/OneDrive/deepthought/Results/Otros/Tortola/Study2/Model1.1/Copal_2002_2021")
load("1.1CopalData_noconverge_2002_2021.RData")

no_converge <- no_converge[-7] # Transecto 7 (num 169) no he podido (error por resolver)

no_converge2 <- NULL

setwd("D:/Deepthought/Results/Otros/Tortola/Study2/Model1.1/2002_2021/COPAL")

for (i in 1:length(no_converge)){
  load(paste("1.1CopalData_0221_", no_converge[i], "_9e5iter", ".RData", sep = ""))
  
  sum <-as.data.frame(out$summary[rownames(out$summary) 
                                  %in% c("mu.lam.year", "sig.lam.year", "bYear.lam", "mu.sig.year", "sig.sig.year", "brough.sig", "sd", "rho",'Bp.Obs', 'Bp.N'), ]) # Only structural parameters
  sum$Rhat[is.na(sum$Rhat)] <- 1 # Just for one observation, the rest converge
  if (sum(ifelse(sum$Rhat > 1.30, 1, 0)) > 0 )
    no_converge2 <- c(no_converge2, no_converge[i])
  next 
}

setwd("C:/Users/anasa/OneDrive/deepthought/Results/Otros/Tortola/Study2/Model1.1/Copal_2002_2021")
save(no_converge2, file = "1.1CopalData_noconverge2_9e5iter_2002_2021.RData") # Save to run at 500000

# --> Check how bad is convergence of 9e5 iter

setwd("C:/Users/anasa/OneDrive/deepthought/Results/Otros/Tortola/Study2/Model1.1/Copal_2002_2021")
load("1.1CopalData_noconverge2_9e5iter_2002_2021.RData")

list_no_converge2 <- list()

for (i in 1:length(no_converge2)){
  load(paste("1.1CopalData_0221_", no_converge2[i], "_9e5iter", ".RData", sep = ""))
  
  list_no_converge2[[i]] <-as.data.frame(out$summary[rownames(out$summary) 
                                                    %in% c("mu.lam.year", "sig.lam.year", "bYear.lam", "mu.sig.year", "sig.sig.year", "brough.sig", "sd", "rho",'Bp.Obs', 'Bp.N'), ]) # Only structural parameters
}

list_no_converge2

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