
# FOR ANALYSIS FOR GANGA WITH BIN 5 FOR 2010 - 2020

rm(list=ls())

library(dplyr)
library(stringr)

setwd("D:/Otros/Ganga/Trend_HDS_model_ch2")
dat <- read.csv("Data_DS_10_20.csv", sep = ";")

dat$Especie <- as.character(dat$Especie)
dat <- dat[ ,-4] 
colnames(dat)[which(colnames(dat) == "ï..Transecte_detall_Id_transecte_detall")] <- "Id_transecte_detall" # To make it equal to 2018 and 2019



# ---- Column names ----
names(dat)
colnames(dat)[which(colnames(dat) == "Id_transecte_detall")] <- "Sample.Label"
colnames(dat)[which(colnames(dat) == "Codi_seca")] <- "Region.Label"
colnames(dat)[which(colnames(dat) == "Any")] <- "Year"
colnames(dat)[which(colnames(dat) == "Especie")] <- "Species"
colnames(dat)[which(colnames(dat) == "Nombre")] <- "Count"
colnames(dat)[which(colnames(dat) == "Sexe")] <- "Sex"
colnames(dat)[which(colnames(dat) == "Us")] <- "Crop_type"
colnames(dat)[which(colnames(dat) == "Tipus_observacio")] <- "Obs_type"
colnames(dat)[which(colnames(dat) == "Hora_inici")] <- "Start_time"
colnames(dat)[which(colnames(dat) == "Observador")] <- "Observer"
colnames(dat)[which(colnames(dat) == "Vent")] <- "Wind"
colnames(dat)[which(colnames(dat) == "Nuvolositat")] <- "Clouds"
colnames(dat)[which(colnames(dat) == "Temperatura")] <- "Temp"
dat$Effort <- 500

# ----- Create variable transectID, than matches with the code of the GIS layers (i.e., two digits: 09) ----

#1. Add a 0 before the transect number
for (i in 1:nrow(dat)){ 
  dat$Num_transecte[i] <- paste(0,dat$Num_transecte[i], sep = "")
}

#2. Keep only the last 2 digits (or 3 in the case of the transects that contain 100)

for (i in 1:nrow(dat)) { 
  cent <- substr(dat$Num_transecte[i], 4,4)
  cent <- as.numeric(cent) # NA if it doesnt have 4 digits
  if(is.na(cent)) { # if is NA (has 3 digits)
    dat$Num_transecte[i] <- str_sub(dat$Num_transecte[i], start = -2) # Keep the last 2
  } else { dat$Num_transecte[i] <- str_sub(dat$Num_transecte[i], start = -3)} # Otherwise, keep the last 3
}


# Create variable by pasting it
for (i in 1:nrow(dat)){ 
  dat$transectID[i] <- paste(dat$Region.Label[i],dat$Num_transecte[i], sep = "")
}

# ---- Distance ----
setwd("D:/PhD/Second chapter/Data")
band <- read.csv("Banda.csv", sep = ";")
colnames(band)[1] <- "Banda"

dat <- left_join(dat, band, by = "Banda") # Joind bands (distbegin/distend)

dat <- dat[which(!is.na(dat$Banda)), ]

#Delete banda 6 (flying)
dat <- dat[-which(dat$Banda == 6), ]

dat$distance <- NA # Medium point of each bin except in bin 4

for (i in 1:nrow(dat)){
  if (dat$Banda[i] == 1) {dat$distance[i] = 12.5}
  else if (dat$Banda[i] == 2) {dat$distance[i] = 37.5}
  else if (dat$Banda[i] == 3) {dat$distance[i] = 75}
  else if (dat$Banda[i] == 4) {dat$distance[i] = 150}
  else  {dat$distance[i] = 350} 
}

# Flying? Banda 6


# ---- Repeated observations ----
# There are few transects that have 2 census in the same year-season.
# Because I have joined the bin5 it doesn't work with the indexes from 2010-2018 without bin5. 
# But it works from when I didn't delete any observations

for (i in 1:nrow(dat)){ 
  dat$T_Y[i] <- paste(dat$transectID[i],dat$Year[i], sep = "_")
}

trans <- dat[!duplicated(dat$Sample.Label), which(colnames(dat) %in% c("Sample.Label", "T_Y"))]
trans_rep <- trans[which(duplicated(trans$T_Y)), ]

# ----  Select only transects in ALFES and GRANJA (where are the animals) ----
dat <- dat[which(dat$Region.Label %in% c("AF")), ] # Better to analyze only ALFES


###################################################################################
# I DID THIS TO ANALYZE THE DATA OF GR BUT I FINALLY DIDNT DO IT IN 2019 SO I REMOVED IT (See script shape data 2019)

#dat <- dat[which(dat$Region.Label %in% c("AF", "GR")), ]

# For ALFES: Remove transects out of ZEPA (DONT NOW (IM GONNA CALCULATE DENSITY INSIDE AS A NEW CATEGORY))
#dat <- dat[-which(dat$transectID %in% c("AF31", "AF35", "AF36", "AF38", "AF39", "AF41",
#                                        "AF42", "AF43")), ]

# For GRANJA: Select the transects where they are always
#dat <- dat[which(dat$transectID %in% c("GR08", "GR10", "GR11", "GR12") | dat$Region.Label == "AF"), ]

###################################################################################

# AND THIS WAS HERE FROM WHEN I WAS ANALYZING THE WHOLE COMMUNITY:

# These are the ones without bin 5
#rem <- c(122, 181, 165, 110, 197, 125, 178, 160, 170, 192, 253, 
#1357, 804, 809, 243, 1203, 711)

# These are the ones from before, that work now as well:
# DUPLICATES: The sample label changed when removing observations related bin4. The previous
# label (data not modified is listed in blue)
# In some its because census were repeated in january, april and may. Take the ones of late April/May (In AL):
# Remove sample.label: 122, 198, 178, 114, 216, 131, 194, 172, 184, 210, 281
# In others, 2 of the same season
# Remove sample.label: 1505, 1737, 1744. Take the ones I could modify
# In others, different weather conditions. Take good coditions
# Remove sample.label: 268, 1350, 1594
# Comparison with the new data (modified is in excelfile DataDS_comparedup)
# Remove the one repeated from 2019: 

############################################################################

# SINCE I CHOOSE ALREADY THE ONES OF ALFES THIS DOESNT MAKE SENSE ANYMORE

#rem <- c(122, 198, 178, 114, 216, 131, 194, 172, 184, 210, 281,
#         1505, 1737, 1744, 268, 1350, 1594)
#dat <- dat[-which(dat$Sample.Label %in% rem), ] # No repeated observations
#
#trans2 <- dat[!duplicated(dat$Sample.Label), which(colnames(dat) %in% c("Sample.Label", "T_Y"))]
#trans_rep <- trans2[which(duplicated(trans2$T_Y)), ]


# ---- Remove transects that are irrigated (and therefore have very different conditions) ----
irri <- read.csv("TransecteAnyReg.csv", sep = ";")
colnames(irri)[2] <- "Num_transecte"
colnames(irri)[1] <- "Region.Label"

# CREATE TRANSECT ID VARIABLE
# Add a 0 before the transect number
for (i in 1:nrow(irri)){ 
  irri$Num_transecte[i] <- paste(0,irri$Num_transecte[i], sep = "")}
# Keep only the last 2 digits 
library(stringr)
for (i in 1:nrow(irri)){ 
  irri$Num_transecte[i] <- str_sub(irri$Num_transecte[i], start = -2)}
# Create variable by pasting it
for (i in 1:nrow(irri)){ 
  irri$transectID[i] <- paste(irri$Region.Label[i],irri$Num_transecte[i], sep = "")}

# Remove the ones irrigated all years
irri_all <- irri$transectID[which(irri$Regadio == 1)]
dat[which(dat$transectID %in% irri_all), ] # There is none! So I remove this because it has 0 obs
#dat <- dat[-which(dat$transectID %in% irri_all), ]

# Remove the ones irrigated the year it changed (report remove = 1 for the ones to remove)

irri_change <- irri[which(!is.na(irri$X1er.año.cambio)), ]
irri_change_ID <- irri_change$transectID
irri_change_year <- irri_change$X1er.año.cambio
dat$remove <- NA

for (i in 1:nrow(dat)){
  if (sum(dat$transectID[i] == irri_change_ID)>0) { # For the transects that changed from irrigation
    tmp_change <- irri_change_year[which(irri_change_ID == dat$transectID[i])] # Year of change
    
    if(dat$Year[i] >= tmp_change){
      dat[i,which(colnames(dat) %in% "remove")] <- 1 # Data from that year gets a 1 (to be removed)
    }}
}

# Remove the ones (1) and column remove

dat <- dat[-which(dat$remove == 1), ]
dat <- dat[ ,-which(colnames(dat) %in% "remove")]



# FIX OBSERVATION CO-VARIATES
# Temperature: mistakes typing
dat$Temp[which(dat$Temp == 0)] <- 10
dat$Temp[which(dat$Temp == 100)] <- 10

# Na (cojo el valor de el transecto anterior realizado o algo fiable)
unique(dat$Temp)
dat[which(is.na(dat$Temp)), ] 
dat$Temp[which(dat$T_Y == "AF09_2018")] <- 20

##################################################################
# JOIN WITH HABITAT QUALITY CATEGORY
# Spatial join

library(rgdal)
library(sf)

setwd("D:/Otros/Ganga/Trend_HDS_model_ch2")
tr <- readOGR("D:/PhD/Third chapter/GIS", "Trans_2018_EPSG23031") 
pa <- readOGR("D:/Otros/Ganga/Trend_HDS_model_ch2", "PA") 
trAF <- tr[grep("AF", tr$Codi), ]

trAF_sf <- st_as_sf(trAF)
pa_sf <- st_as_sf(pa)

j <- st_join(trAF_sf, pa_sf)

# Check duplicates to keep the one with the majority in the zone X

j <- j[,c(2,3,4,25)]

# Export and modify excel:
#write.csv(j, "transect_habquality_ganga.csv")
########################################################################

# Load manually reclassified
setwd("D:/Otros/Ganga/Trend_HDS_model_ch2")
hq <- read.csv("transect_habquality_ganga_reclas.csv", sep = ";")
colnames(hq)[3] <- "transectID"

dat2 <- left_join(dat,hq)


#write.csv(dat2,"DataDS_ready_para_Ganga_hq_10_20.csv") 

##################################################

# Calculate areas habitat quality AF
pa <- readOGR("D:/PhD/Otros/Ganga/Trend_HDS_model_ch2", "HQ_AF")

library(rgeos)

area <- as.vector(gArea(pa,byid = TRUE)) #Area
pa@data <- cbind(pa@data,area)
df <- as.data.frame(pa@data)
zona <- unique(df$ZONA)

area_zona <- aggregate(df$a, by = list(df$ZONA), FUN = sum)
area_zona_HA <- cbind(area_zona, area_zona$x/10000) # En HA

