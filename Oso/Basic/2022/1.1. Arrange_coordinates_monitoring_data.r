## -------------------------------------------------
##     Create DB with unified coordinate systems 
##           Geographic coordinates (WGS84)
## ------------------------------------------------- 

rm(list = ls())

library(sp)
library(rgdal)
library(stringr)
library(dplyr)

setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
os <- read.csv("Seguiment_Ossos_Pirineus_1996_2021_Pre-coordinates.csv", header = TRUE, row.names = NULL, sep = ";")
Encoding(os$Remarks) <- "UTF-8" # Change encoding for weird characters (Â imported from excel)
Encoding(os$X) <- "UTF-8" 
os$Y <- as.character(os$Y)
Encoding(os$Y) <- "UTF-8"
Encoding(os$Municipality) <- "UTF-8"
Encoding(os$Site) <- "UTF-8"
Encoding(os$Database) <- "UTF-8"
Encoding(os$Region) <- "UTF-8"
colnames(os)[1] <- "Confirmed_Individual"

os$X <- str_trim(os$X, side = c("both")) # Remove white spaces
os$Y <- str_trim(os$Y, side = c("both"))

os$Probable_Individual <- str_trim(os$Probable_Individual, side = c("both")) # Remove white spaces
os$Probable_Individual[which(os$Probable_Individual == "")] <- "Indetermined" # Need that all white spaces are "Indetermined"

# Remove duplicates 

os <- os[-which(duplicated(os)), ]


## ---- First: Join sex and age manually ----

# Load info table

setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
info <- readxl::read_xlsx("D:/MargSalas/Oso/Datos/Tablas_finales/2022/Info_individuals_2021.xlsx", sheet = 1)
info <- info[,c(4,5,8)]
colnames(info)[1] <- "Probable_Individual"

# To check that all names match between info and os
# os$Probable_Individual[which(os$Probable_Individual != "" & os$Probable_Individual != "Indetermined" & os$Probable_Individual %in% info$Probable_Individual == FALSE)]
# which(os$Confirmed_Individual != "" & os$Confirmed_Individual != "Indetermined" & os$Confirmed_Individual %in% info$Probable_Individual == FALSE)
# Except BRUNA in Probable_individuals (we don't know who is), all match

os <- left_join(os,info,by = "Probable_Individual")

# SEX
for (i in 1:nrow(os)){    # I could do with ifelse function but it doesn't work I don't know why
  if (is.na(os$Sex.y[i])){
    os$Sex.y2[i] <- os$Sex.x[i]
  } else {
    os$Sex.y2[i] <- os$Sex.y[i]
  }
}
os <- os[,-c(3,30)]
os <- os[ ,c(1:2,30,3:29)]
colnames(os)[3] <- "Sex"

# AGE
os$Year <- as.numeric(os$Year)
os$Year_birth <- as.numeric(os$Year_birth)
os <- os %>% mutate(Age2 = Year - Year_birth)

for (i in 1:nrow(os)){    # I could do with ifelse function but it doesn't work I don't know why
  if (is.na(os$Age2[i])){
    os$Age2[i] <- os$Age[i]
  } else {
    os$Age2[i] <- os$Age2[i]
  }
}

os <- os[,-c(4,30)]
os <- os[,c(1:3,29,4:28)]
colnames(os)[4] <- "Age"

## ---- SECOND: Prepare to arrange coordinates ----

os$X <- as.numeric(os$X)
os$Y <- as.numeric(os$Y)

os$X[os$Coordinate_system == "Lat_Long"]
os$Y[os$Coordinate_system == "Lat_Long"]

os_na <- os[which(is.na(os$X)), ] # Extract NA to join later
os_na$x_long <- NA
os_na$y_lat <- NA

os <- os[-which(is.na(os$X)), ] # Remove NA to convert coordinates

CRS_geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # Transformed coordinate system (geographic WGS 84)

## ---- Lat_Long ----
# Toponims from french data (obtained in Geoportail) - EPSG 3857 Proyección Google Mercator

os_latlong <- os[which(os$Coordinate_system == "Lat_Long"), ]
coordinates(os_latlong) <- os_latlong[,c("X","Y")] # Spatial object

# os_latlong@proj4string <- CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs")
# When we assign this projection the transformation does not work
# The difference with WGS 84 geographic is very small, and all the Lat_Long coordinates have
# a precision 2-3 anyway 
# So put the same coordinates and consider them as WGS 84 geographic

os_latlong@proj4string <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

os_latlong_TRANSF <- os_latlong
os_latlong_TRANSF$x_long <- os_latlong$X
os_latlong_TRANSF$y_lat <- os_latlong$Y

## ---- Lambert II ----
# French data - EPSG:27572 NTF (Paris) / Lambert zone II

os_lambert <- os[which(os$Coordinate_system == "LambertII"), ]
coordinates(os_lambert) <- os_lambert[,c("X","Y")] # Spatial object
os_lambert@proj4string <- CRS("+proj=lcc +lat_1=46.8 +lat_0=46.8 +lon_0=0 +k_0=0.99987742 +x_0=600000 +y_0=2200000 +a=6378249.2 +b=6356515 +towgs84=-168,-60,320,0,0,0,0 +pm=paris +units=m +no_defs")

os_lambert_TRANSF <- spTransform(os_lambert, CRS_geo)
os_lambert_TRANSF$x_long <- coordinates(os_lambert_TRANSF)[,1]
os_lambert_TRANSF$y_lat <- coordinates(os_lambert_TRANSF)[,2]

## ---- ETRS 89_31N ----

os_etrs31N <- os[which(os$Coordinate_system == "ETRS89_31N"), ]
coordinates(os_etrs31N) <- os_etrs31N[,c("X","Y")] # Spatial object
os_etrs31N@proj4string <- CRS("+proj=utm +zone=31 +ellps=GRS80 +units=m +no_defs") # Define coordinate system ETRS 89/31N

os_etrs31N_TRANSF <- spTransform(os_etrs31N, CRS_geo)
os_etrs31N_TRANSF$x_long <- coordinates(os_etrs31N_TRANSF)[,1]
os_etrs31N_TRANSF$y_lat <- coordinates(os_etrs31N_TRANSF)[,2]

## ---- ETRS 89_30N ----

os_etrs30N <- os[which(os$Coordinate_system == "ETRS89_30N"), ]
coordinates(os_etrs30N) <- os_etrs30N[,c("X","Y")] # Spatial object
os_etrs30N@proj4string <- CRS("+proj=utm +zone=30 +ellps=GRS80 +units=m +no_defs") # Define coordinate system ETRS 89/30N

os_etrs30N_TRANSF <- spTransform(os_etrs30N, CRS_geo)
os_etrs30N_TRANSF$x_long <- coordinates(os_etrs30N_TRANSF)[,1]
os_etrs30N_TRANSF$y_lat <- coordinates(os_etrs30N_TRANSF)[,2]

## ---- ED50  ----
# We don't know to which zone each of them belong, so we assign it spacially 
# Looking at arcgis it looks like: 
# ED50 X > 500000 = ED50 30N
#      X < 500000 = ED50 31N

os_ed50 <- os[which(os$Coordinate_system == "ED50"), ]
coordinates(os_ed50) <- os_ed50[,c("X","Y")] # Spatial object

## ED50_30N

os_ed5030N <- os_ed50[which(os_ed50$X > 500000), ]
os_ed5030N$Coordinate_system <- "ED50_30N"

os_ed5030N@proj4string <- CRS("+proj=utm +zone=30 +ellps=intl +units=m +no_defs") # Define coordinate system ED50

os_ed5030N_TRANSF <- spTransform(os_ed5030N, CRS_geo)
os_ed5030N_TRANSF$x_long <- coordinates(os_ed5030N_TRANSF)[,1]
os_ed5030N_TRANSF$y_lat <- coordinates(os_ed5030N_TRANSF)[,2]

## ED50_31N

os_ed5031N <- os_ed50[which(os_ed50$X < 500000), ]
os_ed5031N$Coordinate_system <- "ED50_31N"

os_ed5031N@proj4string <- CRS("+proj=utm +zone=31 +ellps=intl +units=m +no_defs") # Define coordinate system ED50

os_ed5031N_TRANSF <- spTransform(os_ed5031N, CRS_geo)
os_ed5031N_TRANSF$x_long <- coordinates(os_ed5031N_TRANSF)[,1]
os_ed5031N_TRANSF$y_lat <- coordinates(os_ed5031N_TRANSF)[,2]

## ---- WGS84-UTM  ----
# We don't know to which zone each of them belong, so we assign it spacially 
# Looking at arcgis it looks like all WGS84 are 31N

os_wgs84 <- os[which(os$Coordinate_system == "WGS 84"), ]
os_wgs84$Coordinate_system <- "WGS84_UTM_31N"
coordinates(os_wgs84) <- os_wgs84[,c("X","Y")] # Spatial object
os_wgs84@proj4string <- CRS("+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") # Define coordinate system wgs84

os_wgs84_TRANSF <- spTransform(os_wgs84, CRS_geo)
os_wgs84_TRANSF$x_long <- coordinates(os_wgs84_TRANSF)[,1]
os_wgs84_TRANSF$y_lat <- coordinates(os_wgs84_TRANSF)[,2]


## ---- Join and arrange all ----

os_spatial <- rbind(os_latlong_TRANSF, os_lambert_TRANSF, os_etrs31N_TRANSF, os_etrs30N_TRANSF, 
                 os_ed5030N_TRANSF, os_ed5031N_TRANSF, os_wgs84_TRANSF)

os_data <- rbind(os_spatial@data, os_na)

# Arrange by date
os_data$Date_register_formated <- os_data$Date_register # Create column with date format to sort data frame by date
os_data$Date_register_formated <- as.Date(os_data$Date_register_formated, format = "%d/%m/%Y")
os_data <- os_data[,c(1:8,32,9:31)]
os_data <- arrange(os_data, Year, Date_register_formated) # Arrange by year and date, to keep in order the registers without a specific day
os_data <- os_data[,-9] # Remove column with date format and keep the original one (with non-accurate dates)

# Sort out columns
os_data <- os_data[,c(1:12,30,31,13:29)]


# Assign number to the observations

os_data$ID_obs <- seq(1:nrow(os_data))
os_data <- os_data[,c(32,1:31)]

setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
write.csv(os_data, "Seguiment_Ossos_Pirineus_1996_2021.csv")

