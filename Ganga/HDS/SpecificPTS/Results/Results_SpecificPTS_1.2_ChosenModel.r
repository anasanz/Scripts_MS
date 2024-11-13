
## -------------------------------------------------
##  Load and process results from full model Specific PTS
##   I ran it 5 times to show how different are results
## ------------------------------------------------- 

rm(list=ls())

## ---- Results ----

## LOAD all model runs 
# Set the path to the folder containing the RData files
folder_path <- "D:/MargSalas/Ganga/Results/HDS/SpecificPTS/Model"

# Get a list of all RData files in the folder
file_list <- list.files(path = folder_path, pattern = "\\.RData$", full.names = TRUE)

# Create a list to store the data from each file
data_list <- list()

# Loop through the files, read the data, and assign to different object names
for (i in seq_along(file_list)) {
  # Extract the file name without the path and extension
  file_name <- basename(file_list[i])
  object_name <- gsub("\\.RData$", "", file_name)
  
  # Load the data from the file and assign it to the corresponding object name
  load(file_list[i])
  data_list[[object_name]] <- out
}

# Store all results in lists to just visualize

summary <- list()
mean <- list()
mode <- list()

for(xxx in 1:length(data_list)){

out <- data_list[[xxx]]

summary[[xxx]] <- out$summary

#library(MCMCvis)

#MCMCtrace(out, 
#          params = params, 
#          pdf = FALSE)

# Load hq areas

setwd("D:/MargSalas/Ganga/Data/FarmdindisDS")
hq_area <- read.csv(file = "HQ_area.csv", sep = ";")

area_transect <- 500*800 # m2
average_clus2 <- 4.33 # From model script

## ---- 1. Predictions from posterior distribution ----

# Parameters to predict abundance in each hq zone
alpha <- out$sims.list$alpha # Mean = summary -1.532456
bHQ <- out$sims.list$bHQ # Mean = summary

hqzones <- c("hq1", "hq2", "hq3")

ab <- data.frame(matrix(NA, nrow = length(alpha), ncol = 4))
colnames(ab) <- c(hqzones, "total")

de <- data.frame(matrix(NA, nrow = length(alpha), ncol = 4))
colnames(de) <- c(hqzones, "total")

for (i in 1:length(hqzones)) {
  
  lambda <- exp(alpha + bHQ * i) # Expected abundance
  
  dens <- lambda/area_transect
  abundance <- dens*hq_area$x[i]
  total_abundance <- abundance*average_clus2
  ab[,i] <- total_abundance
  
  densHA <- lambda/(area_transect/10000)
  
  de[,i] <- densHA*average_clus2
  
}

ab[,4] <- rowSums(ab[,c(1:3)])
de[,4] <- rowSums(de[,c(1:3)])

# All observations

dens_obs1 <- density(ab[,4]) 
mode_ab <- dens_obs1$x[dens_obs1$y == max(dens_obs1$y)]
mean_ab <- mean(ab[,4])

mode[[xxx]] <- mode_ab 
mean[[xxx]] <- mean_ab 

# 95% CI (excludes the 2.5% of obs with lower and higher values)

lci <- quantile(ab[,4],probs = 0.025) 
uci <- quantile(ab[,4],probs = 0.975)

# Observations of value lower than the upper ci 474.58 (97.5%)

dens_obs2 <- density(ab[,4][which(ab[,4]<uci)]) 

# 90% CI (excludes the 5% of obs with lower and higher values)

lci2 <- quantile(ab[,4],probs = 0.05)
uci2 <- quantile(ab[,4],probs = 0.95)

# 85% CI (excludes the 7.5% of obs with lower and higher values)

lci3 <- quantile(ab[,4],probs = 0.075)
uci3 <- quantile(ab[,4],probs = 0.925)

## ---- 1.1.1. Plot ----

setwd("D:/MargSalas/Ganga/Results/HDS/SpecificPTS/Plots/Abundance")
#pdf(paste("Abundance_estimate_",names(data_list)[[xxx]], ".pdf", sep = ""), 7, 9)

par(mfrow = c(2,1),
    mar = c(3.2,3,2,1))
plot(dens_obs1, main = "Full posterior distribution", col = "darkcyan", lwd = 1.2, xlab = " ", ylab = " ", bty = "n", axes = FALSE)
polygon(c(dens_obs1$x, 0), c(dens_obs1$y, 0), col="darkcyan", border = "darkcyan") # ?? I still don't know if its right
polygon(c(dens_obs1$x[which(dens_obs1$x > lci3 & dens_obs1$x < uci3)], uci3, lci3), 
        c(dens_obs1$y[which(dens_obs1$x > lci3 & dens_obs1$x < uci3)], 0, 0), 
        col=adjustcolor("yellow", alpha = 0.5),
        border = adjustcolor("yellow", alpha = 0.5)) 

axis(1, pos = 0, tck = -0.05, cex.axis = 0.9)
axis(2, pos = -0.5, tck = -0.05, cex.axis = 0.9)

mtext("Abundance (N)", side = 1, line = 1.7)
mtext("Probability", side = 2, line = 1.7)


plot(dens_obs2, main = "Posterior distribution (0 - 97.5% CI)", col = "darkcyan",lwd = 1.5, xlab = "", ylab = "", bty = "n", axes = FALSE)
polygon(c(dens_obs2$x, 0), c(dens_obs2$y, 0), col="darkcyan", border = "darkcyan",lwd = 1.5) # ?? I still don't know if its right
polygon(c(dens_obs2$x[which(dens_obs2$x > lci3 & dens_obs2$x < uci3)], uci3, lci3), 
        c(dens_obs2$y[which(dens_obs2$x > lci3 & dens_obs2$x < uci3)], 0, 0), 
        col=adjustcolor("yellow", alpha = 0.5),
        border = "yellow", lwd = 1.5) 

axis(1, pos = 0, tck = -0.05, cex.axis = 0.9)
axis(2, pos = -0.5, tck = -0.05, cex.axis = 0.9)

mtext("Abundance (N)", side = 1, line = 1.7)
mtext("Probability", side = 2, line = 1.7)

abline(v = mean_ab, col = "darkslategrey", lwd = 2)
abline(v = mode_ab, col = "darkslategrey", lwd = 2)

segments(x0=lci3,y0=0,x1=lci3,y1=dens_obs2$y[41],col="yellow", lwd = 1.5)
segments(x0=uci3,y0=0,x1=uci3,y1=dens_obs2$y[270],col="yellow", lwd = 1.5)

#text(x = 82, y = 0.0007, labels = "Mean:\n1797 ind", adj = 0, pos = 4, col = "darkslategrey", cex = 1, font = 2)
#text(x = 34, y = 0.0007, labels = "Mode:\n238 ind", adj = 0, pos = 4, col = "darkslategrey", cex = 1, font = 2)


#dev.off()

}

for(xxx in 1:length(data_list)){
  print(summary[[xxx]])
  print(mean[[xxx]])
  print(mode[[xxx]])
}
