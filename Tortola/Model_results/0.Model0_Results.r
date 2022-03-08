
rm(list=ls())

library(dplyr)

source("D:/MyScripts/Ch. 2-3/Ch. 3/Results/Functions/plot.violins3.r")
source("D:/MyScripts/Ch. 2-3/Ch. 3/Results/Functions/DoScale.r")

setwd("S:/Results/Otros/Tortola/Study2/Model0")

load("0TortoData_transects.RData") # Load analyzed transects

## ---- Check convergence of structural parameters ----

# With with 2e5 iter

no_converge <- NULL

for (i in 1:length(transect)){
 load(paste("0TortoData", "_", transect[i], ".RData"))
  
  sum <-as.data.frame(out$summary[rownames(out$summary) 
                    %in% c("mu.lam.year", "sig.lam.year", "bYear.lam", "mu.sig.year", "sig.sig.year", "sd", "rho",'Bp.Obs', 'Bp.N'), ]) # Only structural parameters
  if (sum(ifelse(sum$Rhat > 1.15, 1, 0)) > 0 )
    no_converge <- c(no_converge, transect[i])
    next 
}

save(no_converge, file = "0TortoData_noconverge.RData") # Save to run at 500000

# With 5e5 iter

load("0TortoData_noconverge.RData")

no_converge2 <- NULL

for (i in 1:length(no_converge)){
  load(paste("0TortoData", "_", no_converge[i], "5e5iter", ".RData"))
  
  sum <-as.data.frame(out$summary[rownames(out$summary) 
                                  %in% c("mu.lam.year", "sig.lam.year", "bYear.lam", "mu.sig.year", "sig.sig.year", "sd", "rho",'Bp.Obs', 'Bp.N'), ]) # Only structural parameters
  if (sum(ifelse(sum$Rhat > 1.15, 1, 0)) > 0 )
    no_converge2 <- c(no_converge2, no_converge[i])
  next 
}

## ---- Probability of increase in trend ----

# Create vector with name of files to be loaded in results (both 2e5 and 5e5 iter)

converge <- transect[-which(transect %in% no_converge)] # Converged with 2e5 iter: 
converge1_files <-paste("0TortoData", "_", converge, ".RData") # remove the ones that did not converge

converge2 <- transect[-which(transect %in% c(converge, no_converge2))] # Converged with 5e5 iter: 
converge2_files <-paste("0TortoData", "_", converge2, "5e5iter", ".RData")  #remove the ones that did converged before and did not converge now

converge_all <- c(converge, converge2)
files <- c(converge1_files, converge2_files)

trend <- data.frame(Site = converge_all, p_increasing = NA, name_file = NA)

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
write.csv(trend,"trend0.csv")

## ---- Plot coefficients and probabilities ----

files_order <- trend$name_file

setwd("D:/Otros/Tórtola/Results/Study2/Plots")

pdf("trend_transects_model0.pdf", width = 5, height = 9)
plot(10, ylim = c(1, length(files_order)), 
     xlim = c(-3,5), 
     type ="n", yaxt="n", xlab = " ", ylab = "", main = "Population trend")
axis(2, c(1:length(files_order)), labels = trend$Site, las = 2, cex.axis = 0.9)
mtext("Beta coefficient", side = 1, line = 2.5)

setwd("S:/Results/Otros/Tortola/Study2/Model0")
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
segments(0, -1, 0, 38.4, col = "red")
text(4, sort(seq(1:38)), labels = trend$label, cex = 0.7)

dev.off()

