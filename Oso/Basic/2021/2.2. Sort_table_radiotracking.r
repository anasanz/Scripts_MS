
## -------------------------------------------------
##                 Sort out table
## ------------------------------------------------- 

rm(list = ls())

library(sp)
library(rgdal)
library(lubridate)
library(dplyr)
 

setwd("D:/MargSalas/Oso/Datos/GPS")
os <- read.csv("Radiotracking_ossos_1996_2020_1.csv", header = TRUE, row.names = NULL)
os <- os[,-c(1)]

## ---- Add year column ----

os$Date_register_formated <- os$Date_GMT
os$Date_register_formated <- as.Date(os$Date_register_formated, format = "%d/%m/%Y")
os$Year <- year(ymd(os$Date_register_formated)) 

os <- os[,c(1:3,24,4:22)]


## ---- Sex, Age and Age class ----

# Load info 

setwd("D:/MargSalas/Oso/Datos")
info <- read.csv("Info_individuals.csv", header = TRUE, row.names = NULL, sep = ";")
info <- info[,c(4,5,8)]
colnames(info)[1] <- "Bear_name"

os <- left_join(os,info,by = "Bear_name")

# AGE
os$Year <- as.numeric(os$Year)
os$Year_birth <- as.numeric(os$Year_birth)
os <- os %>% mutate(Age = Year - Year_birth)


os <- os[,c(1:2,24,26,4:6,3,7:23)]


os$Age2 <- os$Age # To convert to numeric to do classification (will produce NA)
os$Age2 <- as.numeric(os$Age2)

for (i in 1:nrow(os)) { 
  if(is.na(os$Age2[i])) next
  if(os$Age2[i] > 4) {
    os$Age_class[i] <- "Adult" } else if (os$Age2[i] <= 4 & os$Age2[i] > 1.5) {
    os$Age_class[i] <- "Subadult" } else if (os$Age2[i] == 1.5 | os$Age2[i] == 1) {
    os$Age_class[i] <- "Cub2" }
}

os <- os[ ,c(1:4,27,5:25)]

## ---- Add region and country ----

coordinates(os) <- os[,c("x_long","y_lat")] # Spatial object
os@proj4string <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# Load layer with regions

map <- readOGR(dsn = "D:/MargSalas/Oso/Datos/GIS/Countries", layer = "clip_pyros2")
d <- over(os,map) # Overlay to see where each point fall

os@data <- cbind(os@data, d)
colnames(os@data)[which(colnames(os@data) %in% c("NAME_0", "NAME_1"))] <- c("Country", "Region")

os@data <- os@data[,c(1:17,27:28,18:26)]

# Export GIS layer with coordinates and good attribute table (last version)
#writeOGR(os, "D:/MargSalas/Oso/Datos/GPS/Radiotracking_layer", "Radiotracking_coordinates", driver = "ESRI Shapefile")

os_data <- os@data
os_data <- arrange(os_data, ID_obs) 


setwd("D:/MargSalas/Oso/Datos/GPS")
write.csv(os_data, "Radiotracking_ossos_1996_2020_taula_final.csv")


