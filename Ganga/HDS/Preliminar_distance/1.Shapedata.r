############################################################################# 
#####  1. Get data into shape for DS analysis (for package DISTANCE)    #####
#############################################################################

rm(list=ls())

library(dplyr)

#UNITS
# Para que se mas facil obtener densidades en km2, utilizar el procedimiento de Distance:
# Detection distance in m, Effort in km y Area en km2
# Utilizar parametro convert.units = 0.001 en modelo (para que converta m a km)

## ---- 1.1. FARMDINDIS ----

setwd("D:/MargSalas/Ganga/Data")
dat <- read.csv('Transectes_Farmdindis_Ganga2022.csv', sep = ";", header = TRUE)
dat <- dat[,-c(3:9)]

# Change column names to make it fit with distance
names(dat)
colnames(dat)[which(colnames(dat) == "ï..Codi_seca")] <- "Region.Label"
colnames(dat)[which(colnames(dat) == "Num_transecte")] <- "Sample.Label"
colnames(dat)[which(colnames(dat) == "Num")] <- "size"
colnames(dat)[which(colnames(dat) == "Banda")] <- "Bin"

datds <- dat
datds <- datds[-which(is.na(datds$Bin)), ]
# Create columns bins
datds$distbegin <- NA
datds$distend <- NA
datds$midpoint <- NA # Only to plot in histogram

unique(datds$Bin)

for (i in 1:nrow(datds)){
  if (datds$Bin[i] == 1) {datds$distbegin[i] = 0 ; datds$distend[i] = 25; datds$midpoint[i] = 12.5}
  else if (datds$Bin[i] == 2) {datds$distbegin[i] = 25 ; datds$distend[i] = 50; datds$midpoint[i] = 37.5}
  else if (datds$Bin[i] == 3) {datds$distbegin[i] = 50 ; datds$distend[i] = 100; datds$midpoint[i] = 75}
  else if (datds$Bin[i] == 4) {datds$distbegin[i] = 100 ; datds$distend[i] = 200; datds$midpoint[i] = 150}
  else if (datds$Bin[i] == 5) {datds$distbegin[i] = 200 ; datds$distend[i] = 500; datds$midpoint[i] = 350}
}

# Add column area in km2
# Area AF = 7000 ha
# Area GR = 1 ha en modelo HDS?????

datds$Area <- NA
for (i in 1:nrow(datds)){
  if (datds$Region.Label[i] == "AF") {datds$Area[i] = 7000/100}
  else if (datds$Region.Label[i] == "GR") {datds$Area[i] = 1/100} }

# Add effort column (500 m length) in km
datds$Effort <- 500/1000

# Plot to see the distribution of the data
hist(datds$midpoint, breaks = c(0,25,50,100,200,500),
     main = "Detection curve", col = "grey", freq = FALSE, xlab = " ")

setwd("D:/MargSalas/Ganga/Data")
write.csv(datds, "farmdindis_distance.csv")

## ---- 1.2. SOCC ----

setwd("D:/MargSalas/Ganga/Data")
dat2 <- read.csv('Transectes_Específics_Ganga.csv', sep = ";", header = TRUE)
dat2 <- dat2[,-c(2:13,15,17:20)]

names(dat2)
colnames(dat2)[which(colnames(dat2) == "ï..ID.Transecte")] <- "Sample.Label"
colnames(dat2)[which(colnames(dat2) == "N")] <- "size"
colnames(dat2)[which(colnames(dat2) == "DistÃ.ncia")] <- "Dist"
dat2$Bin <- NA

datds2 <- dat2
datds2 <- datds2[-which(is.na(datds2$Dist)), ]

# Assign to distance bins
#1: 0-25
#2: 25-50
#3: 50-100

for (i in 1:nrow(datds2)){
  if (datds2$Dist[i] > 0 & datds2$Dist[i] <= 25) {datds2$Bin[i] = 1 }
  else if (datds2$Dist[i] > 25 & datds2$Dist[i] <= 50) {datds2$Bin[i] = 2 }
  else if (datds2$Dist[i] > 50 & datds2$Dist[i] <= 100) {datds2$Bin[i] = 3 }
}

# Create columns bins
datds2$distbegin <- NA
datds2$distend <- NA
datds2$midpoint <- NA # Only to plot in histogram


for (i in 1:nrow(datds2)){
  if (datds2$Bin[i] == 1) {datds2$distbegin[i] = 0 ; datds2$distend[i] = 25; datds2$midpoint[i] = 12.5}
  else if (datds2$Bin[i] == 2) {datds2$distbegin[i] = 25 ; datds2$distend[i] = 50; datds2$midpoint[i] = 37.5}
  else if (datds2$Bin[i] == 3) {datds2$distbegin[i] = 50 ; datds2$distend[i] = 100; datds2$midpoint[i] = 75}
}

#datds2$Area <- 7000/100 # Area in km2; But it is a lot!! What if we put as with fox? 3810 ha, or previous ganga anlysis? 1500
datds2$Area <- 1500/100 # As in previous analysis

# Add effort column (500 m length) in km
datds2$Effort <- 500/1000

# Add region column (Only Alfes)
datds2$Region.Label <- "AF"

# Plot to see the distribution of the data
hist(datds2$midpoint, breaks = c(0,25,50,100,200,500),
     main = "Detection curve", col = "grey", freq = FALSE, xlab = " ")

setwd("D:/MargSalas/Ganga/Data")
write.csv(datds2, "socc_distance.csv")
