
############################################################################# 
#####                 2. Analysis with distance                        #####
#############################################################################

rm(list=ls())
library(Distance)

## ---- 1.1. FARMDINDIS ----

setwd("D:/MargSalas/Ganga/Data")
dat <- read.csv("farmdindis_distance.csv", sep = ",")

# Plot to see the distribution of the data
hist(dat$midpoint, breaks = c(0,25,50,100,200,500),
     main = "Detection curve", col = "grey", freq = FALSE, xlab = " ")

dat <- dat[ ,-c(1,5,8)]

# Units: To know conversion factor
convert_units("metre", "kilometre", "square kilometre")

# Run models
m1 <- ds(data = dat, key = "hr", convert_units = 0.001)
m2 <- ds(data = dat, key = "hn", convert_units = 0.001)

sum1 <- summary(m1) 
sum2 <- summary(m2) 

# Assess model fit (BOTH FIT)
gof_ds(m1) # p > 0.05 -> Good fit (accept the null hypothesis that observed = expected)
gof_ds(m2) # p > 0.05 -> Good fit (accept the null hypothesis that observed = expected)


# Write results 
#PLOTS
setwd("D:/MargSalas/Ganga/Results/HDS")
pdf("DF_Farmdindis.pdf",7,5)
par(mfrow = c(1,2))
plot(m1, nc = 5)
mtext("Hazard rate", side = 3, line = 1)
plot(m2, nc = 5, main = )
mtext("Half normal", side = 3, line = 1)
mtext("Detection function FARMDINDIS", side = 3,line = -1.5, outer = TRUE)
dev.off()

# DENSITY
D1 <- sum1$dht$individuals$D
D2 <- sum2$dht$individuals$D
D1$Model <- "HR"
D2$Model <- "HN"
D <- rbind(D1,D2)
D$Valor <- "Densidad"

# ABUNDANCE
A1 <- sum1$dht$individuals$N
A2 <- sum2$dht$individuals$N
A1$Model <- "HR"
A2$Model <- "HN"
A <- rbind(A1,A2)
A$Valor<- "Abundancia"

results <- rbind(D,A)


## ---- 1.2. SOCC ----

setwd("D:/MargSalas/Ganga/Data")
dat <- read.csv("socc_distance.csv", sep = ",")

# Plot to see the distribution of the data
hist(dat$midpoint, breaks = c(0,25,50,100,200,500),
     main = "Detection curve", col = "grey", freq = FALSE, xlab = " ")

dat <- dat[ ,-c(1,4,5,8)]

# Units: To know conversion factor; convert_units(distance_units, effort_units, area_units)
convert_units("metre", "kilometre", "square kilometre") 

# Run models
m1 <- ds(data = dat, key = "hr", convert_units = 0.001)
m2 <- ds(data = dat, key = "hn", convert_units = 0.001)

sum1 <- summary(m1) 
sum2 <- summary(m2) 

# Assess model fit (BOTH FIT)
gof_ds(m1) # p > 0.05 -> Good fit (accept the null hypothesis that observed = expected)
gof_ds(m2) # p > 0.05 -> Good fit (accept the null hypothesis that observed = expected)


# Write results 
#PLOTS
setwd("D:/MargSalas/Ganga/Results/HDS")
pdf("DF_Socc.pdf",7,5)
par(mfrow = c(1,2))
plot(m1, nc = 5)
mtext("Hazard rate", side = 3, line = 1)
plot(m2, nc = 5, main = )
mtext("Half normal", side = 3, line = 1)
mtext("Detection function FARMDINDIS", side = 3,line = -1.5, outer = TRUE)
dev.off()


# DENSITY
D1 <- sum1$dht$individuals$D
D2 <- sum2$dht$individuals$D
D1$Model <- "HR"
D2$Model <- "HN"
D <- rbind(D1,D2)
D$Valor <- "Densidad"

# ABUNDANCE
A1 <- sum1$dht$individuals$N
A2 <- sum2$dht$individuals$N
A1$Model <- "HR"
A2$Model <- "HN"
A <- rbind(A1,A2)
A$Valor<- "Abundancia"

results <- rbind(D,A)
