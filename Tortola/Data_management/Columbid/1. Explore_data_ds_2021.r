##################################################################################
###      EXPLORE DATA SOCC AMPLIADO OF PREDICTOR COLUMBIDAE     #################
##################################################################################



library(rgdal)
library(sp)
library(raster)
library(dplyr)
library(tidyr)
library(splitstackshape)

## -------------------------------------------------
##   COPAL (Columba palumbus/Paloma Torcaz/Common Wood pigeon)
## ------------------------------------------------- 

# Load base de datos socc ampliado (con distance sampling)

rm(list=ls())

setwd("D:/Otros/Tórtola/Data")
ds <- read.csv("Dades_SOCC_columbids_2021_v2.csv", sep = ";")
colnames(ds)[1] <- "Itinerari"
id_itinerari <- unique(ds$Itinerari)
unique(ds$Especie)

# Select specie
ds <- ds[which(ds$Especie == "COPAL"), ] 

n <- ds %>% group_by(AnySOCC) %>% summarise(n()) 
sum(n[1:20,2]) 

## ---- Exploración datos DS ----

# Problem: rows are not independent
# Check average group size detected during farmdindis
setwd("D:/PhD/Third chapter/Data")
farm <- read.csv("DataDS_ch3_15_19_READY_FIXED.csv") 
farm <- farm[which(farm$Species == "COPAL"), ]
observations <- xtabs(~ Count, farm) 
hist(farm$Count, breaks = 12, at = c(0,1,2,3,4,5,6,7,8,9,10,11,12))
mean(farm$Count) # Higher than in tortola

farm2 <- farm[which(farm$Count != 12), ]
observations2 <- xtabs(~ Count, farm2) 
mean(farm2$Count)
# Proporción escuchadas/vistas
nrow(farm[which(farm$Obs_type == "V"), ])/nrow(farm)


# Put in ds format (1 row per observation)
datds <- gather(ds, Banda, count, banda1:banda3)
datds <- datds[which(datds$count != 0), ]
datds <- datds[ ,-c(5,6,8)]

dat <- expandRows(datds, 'count')


# Explorar periodo: Elegir periodo 2 porque ya están todas

freq <- dat %>% count(AnySOCC, Periode)
sum(freq[which(freq$Periode == 2),3]) # 2 has more observations, good to take it
sum(freq[which(freq$Periode == 1),3])

dat <- dat[which(dat$Periode == 2), -3]

# Cambiar columnas

colnames(dat) <- c("Site", "Year", "Observer", "Section", "Bin")

# Crear nueva columna: Itinerari_seccion
dat$site_sec <- paste(dat$Site, dat$Section, sep = "_")
dat$site_year <- paste(dat$Site, dat$Year, sep = "_")
length(unique(dat$site_year)) # 881 transect-year


# ---- DETECTION CURVE ---- #

dat$distance <- NA # Medium point of each bin except in bin 4

for (i in 1:nrow(dat)){
  if (dat$Bin[i] == "banda1") {dat$distance[i] = 12.5} # 0-25
  else if (dat$Bin[i] == "banda2") {dat$distance[i] = 62.5} # 25-100
  else  {dat$distance[i] = 300} # 100 - 500
}

hist(dat$distance,breaks = c(0,25,99,300), main = "Detection curve", col = "grey", freq = FALSE) 

# Detection curve per site - FOR STUDY 2 (ANALYSIS INDEPENDENT PER SITE)
#### DO THIS LATER ON??? IT'T TOO MANY AND I DON'T KNOW IF i WILL ANALYZE THEM ALL
unique(dat$Site)

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
setwd("D:/Otros/Tórtola/Data/Columbid")
#write.csv(dat, "copal_ds_02_21.csv")

## -------------------------------------------------
##   STDEC (Streptopelia decaocto/Tórtola turca/Eurasian collared dove)
## ------------------------------------------------- 

# Load data again

rm(list=ls())

setwd("D:/Otros/Tórtola/Data")
ds <- read.csv("Dades_SOCC_columbids_2021_v2.csv", sep = ";")
colnames(ds)[1] <- "Itinerari"
id_itinerari <- unique(ds$Itinerari)
unique(ds$Especie)

# Select specie
ds <- ds[which(ds$Especie == "STDEC"), ] 

n <- ds %>% group_by(AnySOCC) %>% summarise(n()) 
sum(n[1:20,2]) 

## ---- Exploración datos DS ----

# In farmdindis I don't have the information, so I can't check average group size
# But I guess it is less than COPAL

# Put in ds format (1 row per observation)
datds <- gather(ds, Banda, count, banda1:banda3)
datds <- datds[which(datds$count != 0), ]
datds <- datds[ ,-c(5,6,8)]

dat <- expandRows(datds, 'count')


# Explorar periodo: Elegir periodo 2 porque ya están todas

freq <- dat %>% count(AnySOCC, Periode)

sum(freq[which(freq$Periode == 2),3]) # 2 has more observations, good to take it
sum(freq[which(freq$Periode == 1),3])

dat <- dat[which(dat$Periode == 2), -3]

# Cambiar columnas

colnames(dat) <- c("Site", "Year", "Observer", "Section", "Bin")

# Crear nueva columna: Itinerari_seccion
dat$site_sec <- paste(dat$Site, dat$Section, sep = "_")
dat$site_year <- paste(dat$Site, dat$Year, sep = "_")
length(unique(dat$site_year)) # 881 transect-year


# ---- DETECTION CURVE ---- #

dat$distance <- NA # Medium point of each bin except in bin 4

for (i in 1:nrow(dat)){
  if (dat$Bin[i] == "banda1") {dat$distance[i] = 12.5} # 0-25
  else if (dat$Bin[i] == "banda2") {dat$distance[i] = 62.5} # 25-100
  else  {dat$distance[i] = 300} # 100 - 500
}

hist(dat$distance,breaks = c(0,25,99,300), main = "Detection curve", col = "grey", freq = FALSE) 

# Detection curve per site - FOR STUDY 2 (ANALYSIS INDEPENDENT PER SITE)
#### DO THIS LATER ON??? IT'T TOO MANY AND I DON'T KNOW IF i WILL ANALYZE THEM ALL
unique(dat$Site)

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
setwd("D:/Otros/Tórtola/Data/Columbid")
write.csv(dat, "stdec_ds_02_21.csv")

## -------------------------------------------------
##   COOEN (Columba Oenas/Paloma zurita/Stock dove)
## ------------------------------------------------- 

# Load data again

rm(list=ls())

setwd("D:/Otros/Tórtola/Data")
ds <- read.csv("Dades_SOCC_columbids_2021_v2.csv", sep = ";")
colnames(ds)[1] <- "Itinerari"
id_itinerari <- unique(ds$Itinerari)
unique(ds$Especie)

# Select specie
ds <- ds[which(ds$Especie == "COOEN"), ] 

n <- ds %>% group_by(AnySOCC) %>% summarise(n()) 
sum(n[1:19,2]) # Casi no hay!! Seguro que es porque es muy común y no la cuentan

## ---- Exploración datos DS ----

# In farmdindis I don't have the information, so I can't check average group size
# But I guess it is less than COPAL

# Put in ds format (1 row per observation)
datds <- gather(ds, Banda, count, banda1:banda3)
datds <- datds[which(datds$count != 0), ]
datds <- datds[ ,-c(5,6,8)]

dat <- expandRows(datds, 'count')


# Explorar periodo: Elegir periodo 2 porque ya están todas

freq <- dat %>% count(AnySOCC, Periode)
sum(freq[which(freq$Periode == 2),3])
sum(freq[which(freq$Periode == 1),3]) # 1 has 60 more observations, but better to take 2

dat <- dat[which(dat$Periode == 2), -3]

# Cambiar columnas

colnames(dat) <- c("Site", "Year", "Observer", "Section", "Bin")

# Crear nueva columna: Itinerari_seccion
dat$site_sec <- paste(dat$Site, dat$Section, sep = "_")
dat$site_year <- paste(dat$Site, dat$Year, sep = "_")
length(unique(dat$site_year)) # 881 transect-year


# ---- DETECTION CURVE ---- #

dat$distance <- NA # Medium point of each bin except in bin 4

for (i in 1:nrow(dat)){
  if (dat$Bin[i] == "banda1") {dat$distance[i] = 12.5} # 0-25
  else if (dat$Bin[i] == "banda2") {dat$distance[i] = 62.5} # 25-100
  else  {dat$distance[i] = 300} # 100 - 500
}

hist(dat$distance,breaks = c(0,25,99,300), main = "Detection curve", col = "grey", freq = FALSE) 

# Detection curve per site - FOR STUDY 2 (ANALYSIS INDEPENDENT PER SITE)
#### DO THIS LATER ON??? IT'T TOO MANY AND I DON'T KNOW IF i WILL ANALYZE THEM ALL
unique(dat$Site)

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
setwd("D:/Otros/Tórtola/Data/Columbid")
write.csv(dat, "cooen_ds_02_21.csv")


