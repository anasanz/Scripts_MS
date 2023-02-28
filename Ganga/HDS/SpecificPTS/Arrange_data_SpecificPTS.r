## -------------------------------------------------
##      Arrange data DS Specific PTS 2022
## -------------------------------------------------

rm(list=ls())

setwd("D:/MargSalas/Ganga/Data/SpecificPTS")
dat <- read.csv('Transectes_Específics_Ganga.csv', sep = ";", header = TRUE)

colnames(dat)[which(colnames(dat) == "ï..ID.Transecte")] <- "transectID"
colnames(dat)[which(colnames(dat) == "Hora.Inici")] <- "Start_time"
colnames(dat)[which(colnames(dat) == "Hora.Final")] <- "End_time"
colnames(dat)[which(colnames(dat) == "Observador")] <- "Observer"
colnames(dat)[which(colnames(dat) == "Vent")] <- "Wind"
colnames(dat)[which(colnames(dat) == "N")] <- "Count"
colnames(dat)[which(colnames(dat) == "DistÃ.ncia")] <- "Distance"

dat$Species <- "PTALC"
dat$Year <- 2022
dat$Effort <- 500

dat <- dat[,c(1,22,2,3,4,7,21,14,16,5,23)]

dat_det <- dat[which(dat$Count>0),]
## ---- Explore detection curve ----

hist(dat$Distance[which(!is.na(dat$Distance))], breaks = c(0,24,49,99,200,500),
     main = "Detection function", col = "grey", freq = FALSE, xlab = "Distance") 

