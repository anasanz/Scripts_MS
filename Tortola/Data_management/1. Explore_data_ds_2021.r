##################################################################################
############      EXPLORE DATA TORTOLA SOCC AMPLIADO      ########################
##################################################################################



library(rgdal)
library(sp)
library(raster)
library(dplyr)
library(tidyr)
library(splitstackshape)

# Base de datos socc ampliado (con distance sampling)

rm(list=ls())
setwd("D:/Otros/Tórtola/Data")
ds <- read.csv("Dades_SOCC_columbids_2021_v2.csv", sep = ";")
colnames(ds)[1] <- "Itinerari"
id_itinerari <- unique(ds$Itinerari)

# Select Tórtola
ds <- ds[which(ds$Especie == "STTUR"), ] # Tórtola

 n <- ds %>% group_by(AnySOCC) %>% summarise(n()) 
sum(n[1:18,2]) # There is 100 observations more than in the data base of 2019, but we move on
sum(n[1:20,2]) # There is 100 observations more than in the data base of 2019, but we move on

## ---- Exploración datos DS ----

# Problem: rows are not independent
# Check average group size detected during farmdindis
setwd("D:/PhD/Third chapter/Data")
farm <- read.csv("DataDS_ch3_15_19_READY_FIXED.csv") 
farm <- farm[which(farm$Species == "STTUR"), ]
observations <- xtabs(~ Count, farm) 
hist(farm$Count, breaks = 12, at = c(0,1,2,3,4,5,6,7,8,9,10,11,12))
mean(farm$Count)

farm2 <- farm[which(farm$Count != 12), ]
observations2 <- xtabs(~ Count, farm2) 
mean(farm2$Count)

# Proporción escuchadas/vistas
nrow(farm[which(farm$Obs_type == "V"), ])/nrow(farm)

hist(farm2$Banda)

# El tamaño medio de grupo detectado en farmdindis no es muy alto.
# No podremos estimar abundancias, y violar la asumcion de datos independencia
# Puede hacer que los estimates sean "overly precise". 
# Pero podemos dar valores de "Probabilidad de declive"

# Put in ds format (1 row per observation)
datds <- gather(ds, Banda, count, banda1:banda3)
datds <- datds[which(datds$count != 0), ]
datds <- datds[ ,-c(5,6,8)]

dat <- expandRows(datds, 'count')


# Explorar periodo: Elegir periodo 2 porque ya están todas

freq <- dat %>% count(AnySOCC, Periode)
dat <- dat[which(dat$Periode == 2), -3]

# Cambiar columnas

colnames(dat) <- c("Site", "Year", "Observer", "Section", "Bin")

# Crear nueva columna: Itinerari_seccion
dat$site_sec <- paste(dat$Site, dat$Section, sep = "_")
dat$site_year <- paste(dat$Site, dat$Year, sep = "_")
length(unique(dat$site_year)) # 881 transect-year

## ---- OBSERVER ----

obs <- unique(dat$Observer)

# Cuantos transectos ha hecho cada observador

data_obs <- dat %>%
  group_by(Observer) %>%
  summarise(n_transects = n_distinct(site_year))

freq_obs <- arrange(data_obs, n_transects)
freq_obs$Observer <- factor(freq_obs$Observer, levels = freq_obs$Observer)

barplot(freq_obs$n_transects ~ freq_obs$Observer, las = 2)

# Select observers that did less than 5 census
# I don't do it unless I'm gonna add observer as a co-variate

# obs_less5 <- freq_obs[which(freq_obs$n_transects < 5), ]
# obs_less5 <- unique(obs_less5$Observer)
# dat_less5 <- dat[which(dat$Observer %in% obs_less5), ] # Perdemos 354 observaciones, pero para meter observador como variable es necesario

# Remove
# dat <- dat[-which(dat$Observer %in% obs_less5), ] # I DON'T REMOVE HERE, CHECK WHEN THE ABSENCES ARE ADDED
# length(unique(dat$site_year)) # 633 transect-year (se pierden 134 transect-year)

# ---- DETECTION CURVE ---- #

dat$distance <- NA # Medium point of each bin except in bin 4

for (i in 1:nrow(dat)){
  if (dat$Bin[i] == "banda1") {dat$distance[i] = 12.5} # 0-25
  else if (dat$Bin[i] == "banda2") {dat$distance[i] = 62.5} # 25-100
  else  {dat$distance[i] = 300} # 100 - 500
}

hist(dat$distance,breaks = c(0,25,99,300), main = "Detection curve", col = "grey", freq = FALSE) 

# Detection curve per observer

obs <- unique(dat$Observer) 

par(mfrow = c())
for (i in 1:length(obs)){
  hist(dat$distance[which(dat$Observer %in% obs[i])],breaks = c(0,25,99,1000),
       main = paste(obs[i], "- Distances"), col = "grey", freq = FALSE)
}

# Some curves are not really good, which could be worse per transect:

sum <- dat %>% group_by(Site, Observer) %>% summarise()
sum_obs_transect <- sum %>% group_by(Site) %>% summarise(n())
par(mfrow = c(1,2))
plot(sum_obs_transect$`n()` ~ sum_obs_transect$Site, pch = 19)
hist(sum_obs_transect$`n()`) # Most have one observer, some 2
(length(sum_obs_transect$`n()`[which(sum_obs_transect$`n()` == 2)])/nrow(sum_obs_transect))*100 # 16%

# Not worthy to include observer, or maybe only for the 35 transects it could be categorical
# Look at the detection curves for those with more than 1

trans_obs <- sum_obs_transect$Site[which(sum_obs_transect$`n()` > 1)]

setwd("D:/Otros/Tórtola/Results/Study2/Plots")

pdf("Transects_observers_2002_2021.pdf", width = 9, height = 7)

for (i in 1:length(unique(trans_obs))) {
  dat_trans <- dat[dat$Site == trans_obs[i], ]
  nobs <- unique(dat_trans$Observer)
  if (length(nobs == 2)) {
    par(mfrow = c(1,2))
      for (o in 1:length(nobs)){
        hist(dat_trans$distance[which(dat_trans$Observer %in% nobs[o])],breaks = c(0,25,99,1000),
             main = paste("Observer",nobs[o]), col = "grey", freq = FALSE, xlab = "Distance")
      mtext(paste("Transect",trans_obs[i]), side = 3, line = -1.3, outer = TRUE, cex = 1.5)}
  } else if (length(nobs == 3)) {
    par(mfrow = c(2,2))
    for (o in 1:length(nobs)){
      hist(dat_trans$distance[which(dat_trans$Observer %in% nobs[o])],breaks = c(0,25,99,1000),
           main = paste("Observer",nobs[o]), col = "grey", freq = FALSE, xlab = "Distance")
      mtext(paste("Transect",trans_obs[i]), side = 3, line = -1.3, outer = TRUE, cex = 1.5)}
    plot(1)
  } else if(length(nobs == 4)) {
    par(mfrow = c(2,2))
    for (o in 1:length(nobs)){
      hist(dat_trans$distance[which(dat_trans$Observer %in% nobs[o])],breaks = c(0,25,99,1000),
           main = paste("Observer",nobs[o]), col = "grey", freq = FALSE, xlab = "Distance")
      mtext(paste("Transect",trans_obs[i]), side = 3, line = -1.3, outer = TRUE, cex = 1.5)}
  }
}
  dev.off()

# Detection curve per site - FOR STUDY 2 (ANALYSIS INDEPENDENT PER SITE)
# This is for the previous data 2002-2019

setwd("C:/Users/anasa/OneDrive/deepthought/Results/Otros/Tortola/Study2/Model0")
load("0TortoData_transects.RData") # Load analyzed transects

site <- transect 

setwd("D:/Otros/Tórtola/Results/Study2/Plots")

pdf("Explore_df_transects_2021.pdf", width = 9, height = 7)

par(mfrow = c(3,3))
for (i in 1:length(site)){
  hist(dat$distance[which(dat$Site %in% site[i])],breaks = c(0,25,99,500),
       main = paste(site[i], "- Distances"), col = "grey", freq = FALSE)
}

dev.off()

# The detection curve per each transect are acceptable to run transect-specific models

# ---- YEARS ---- #

data_years <- dat %>%
  group_by(Year) %>%
  summarise(n_transects = n_distinct(site_year))


# Save
setwd("D:/Otros/Tórtola/Data")
write.csv(dat, "tortola_ds_02_21.csv")
