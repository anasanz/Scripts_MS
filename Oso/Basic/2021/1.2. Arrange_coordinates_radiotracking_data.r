## -------------------------------------------------
##     Create DB with unified coordinate systems 
##           Geographic coordinates (WGS84)
## ------------------------------------------------- 

rm(list = ls())

library(sp)
library(rgdal)
library(stringr)
library(dplyr)

###+ NOTE: ADD  COLUMNS ON WHETHER THE INDIVIDUAL IS SUPPOSSED TO BE WITH CUBS BEFORE LOADING THE FILE:
### Radiotracking_ossos_1996_2020_info.csv
### Also added here the observations of Hvala radiotracking 2009 and 2010 from monitoring table

setwd("D:/MargSalas/Oso/Datos/GPS")
os <- read.csv("Radiotracking_ossos_1996_2020_Pre-coordinates_Info_cubs.csv", header = TRUE, row.names = NULL, sep = ";")

os$X <- as.numeric(os$X)
os$Y <- as.numeric(os$Y)

CRS_geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # Transformed coordinate system (geographic WGS 84)

## ---- Lat_Long ----

os_latlong <- os[which(os$Coordinate_system == "Lat_Long"), ]
coordinates(os_latlong) <- os_latlong[,c("X","Y")] # Spatial object

os_latlong@proj4string <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

os_latlong_TRANSF <- os_latlong
os_latlong_TRANSF$x_long <- os_latlong$X
os_latlong_TRANSF$y_lat <- os_latlong$Y

## ---- Lambert II ----
# French data - EPSG:27572 NTF (Paris) / Lambert zone II

os_lambert <- os[which(os$Coordinate_system == "Lambert II"), ]
coordinates(os_lambert) <- os_lambert[,c("X","Y")] # Spatial object
os_lambert@proj4string <- CRS("+proj=lcc +lat_1=46.8 +lat_0=46.8 +lon_0=0 +k_0=0.99987742 +x_0=600000 +y_0=2200000 +a=6378249.2 +b=6356515 +towgs84=-168,-60,320,0,0,0,0 +pm=paris +units=m +no_defs")

os_lambert_TRANSF <- spTransform(os_lambert, CRS_geo)
os_lambert_TRANSF$x_long <- coordinates(os_lambert_TRANSF)[,1]
os_lambert_TRANSF$y_lat <- coordinates(os_lambert_TRANSF)[,2]


## ---- ED50  ----

## ED50_31N

os_ed5031N <- os[which(os$Coordinate_system == "ED50_31N"), ]
coordinates(os_ed5031N) <- os_ed5031N[,c("X","Y")] # Spatial object

os_ed5031N@proj4string <- CRS("+proj=utm +zone=31 +ellps=intl +units=m +no_defs") # Define coordinate system ED50

os_ed5031N_TRANSF <- spTransform(os_ed5031N, CRS_geo)
os_ed5031N_TRANSF$x_long <- coordinates(os_ed5031N_TRANSF)[,1]
os_ed5031N_TRANSF$y_lat <- coordinates(os_ed5031N_TRANSF)[,2]

## ---- WGS84-UTM  ----

os_wgs84 <- os[which(os$Coordinate_system == "WGS84_UTM_31N"), ]
coordinates(os_wgs84) <- os_wgs84[,c("X","Y")] # Spatial object
os_wgs84@proj4string <- CRS("+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") # Define coordinate system wgs84

os_wgs84_TRANSF <- spTransform(os_wgs84, CRS_geo)
os_wgs84_TRANSF$x_long <- coordinates(os_wgs84_TRANSF)[,1]
os_wgs84_TRANSF$y_lat <- coordinates(os_wgs84_TRANSF)[,2]


## ---- Join and arrange all ----

os_spatial <- rbind(os_latlong_TRANSF, os_lambert_TRANSF, 
                 os_ed5031N_TRANSF, os_wgs84_TRANSF)

os_data <- os_spatial@data

# Arrange by date
os_data$Date_register_formated <- os_data$Date_GMT # Create column with date format to sort data frame by date
os_data$Date_register_formated <- as.Date(os_data$Date_register_formated, format = "%d/%m/%Y")
os_data <- arrange(os_data, Bear_name, Date_register_formated) # Arrange by year and date, to keep in order the registers without a specific day
os_data <- os_data[,-c(1,2,24)] # Remove column with date format and keep the original one (with non-accurate dates)

# Sort out columns
os_data <- os_data[,c(1:10,20,21,11:19)]

# Assign number to the observations

os_data$ID_obs <- seq(1:nrow(os_data))
os_data <- os_data[,c(22,1:21)]

# Save. In this one we will add the info on whether the individuals had cubs
setwd("D:/MargSalas/Oso/Datos/GPS")
write.csv(os_data, "Radiotracking_ossos_1996_2020_1.csv") 

