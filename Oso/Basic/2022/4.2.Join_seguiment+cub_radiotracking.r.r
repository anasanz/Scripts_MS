

## -------------------------------------------------
##                    Explore data
## ------------------------------------------------- 


rm(list = ls())

library(sp)
library(rgdal)
library(dplyr)
library(tidyr)
library(viridis) 
library(RColorBrewer)

# Load monitoring data
setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
os <- readxl::read_xlsx("D:/MargSalas/Oso/Datos/Tablas_finales/2022/Seguiment_Ossos_Pirineus_1996_2021_taula_final_cubLocations.xlsx", sheet = 1)

# To check 
#checks <- os[sample(nrow(os), 60), ]
#write.csv(checks, "check.csv")

os <- os[-which(is.na(os$X)), ] # Remove NA
os <- os[-which(os$Coordinates.precision == 3 | os$Coordinates.precision == 2), ] # Remove coordinates precision 2 and 3

# Load radiotracking data
setwd("D:/MargSalas/Oso/Datos/GPS")
os_gps <- read.csv("Radiotracking_ossos_1996_2020_taula_final.csv", header = TRUE, row.names = NULL)
os_gps <- os_gps[,-1]
os_gps$Bear_name[which(os_gps$Bear_name == "Melba")] <- "Mellba"

## ---- JOIN SPATIAL DATA FROM BOTH DATASETS ----

# Keep column names common to both
colnames(os)
names_os <- c("Confirmed_Individual", "Sex", "Age","Age_class","Year","Date",
              "X","Y","x_long","y_lat","Coordinate_system","Country","Region","With_cubs_estimated","N_cubs_estimated",
              "ID_cub1","ID_cub2","ID_cub3", "Method", "Obs_type", "Remarks")
unique(os$Method)

os <- os[,colnames(os) %in% names_os]

colnames(os_gps)
names_os_gps <- c("Bear_name","Sex","Age","Age_class","Year","Date_GMT","Tracking_system","Coordinate_system",
                  "X","Y", "x_long","y_lat","Country","Region","With_cubs_estimated", "N_cubs_estimated", "ID_cub1",
                  "ID_cub2", "ID_cub3", "Remarks")
os_gps <- os_gps[,colnames(os_gps) %in% names_os_gps]
colnames(os_gps)[c(1,6,7)] <- c("Confirmed_Individual", "Date", "Obs_type") # Change column names as in monitoring db
os_gps$Method <- "Radiotracking"
os_gps <- os_gps[ ,c(1:6,9:12,8,13:19,21,7,20)] # Change column order as in monitoring db

# Join
os_all <- rbind(os,os_gps)

# Arrange by date in case we plot it by colors
os_all$Date_register_formated <- os_all$Date 
os_all$Date_register_formated <- as.Date(os_all$Date_register_formated, format = "%d/%m/%Y")

os_all <- arrange(os_all, Date_register_formated, Year) 

# Set coordinates
coordinates(os_all) <- os_all[,c("x_long","y_lat")] # Spatial object
os_all@proj4string <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

os <- os_all # Call it os to make it fit with code script 3.1

# Save dataset with ALL observations (monitoring+radiotracking)
setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
write.csv(os@data, file = "Data_os_96_21_cubLocations.csv")

