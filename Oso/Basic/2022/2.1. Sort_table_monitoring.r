
## -------------------------------------------------
##                 Sort out table
## ------------------------------------------------- 

rm(list = ls())

library(sp)
library(rgdal)
library(lubridate)
library(dplyr)
library(adehabitatHR)

setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
os <- read.csv("Seguiment_Ossos_Pirineus_1996_2021_2.csv", header = TRUE, row.names = NULL)
os <- os[,-c(1)] 

## ---- Differenciate cubs in critic period ----

# Date column to differenciate month

os$Date_register_formated <- os$Date
os$Date_register_formated <- as.Date(os$Date_register_formated, format = "%d/%m/%Y")
os$month <- month(ymd(os$Date_register_formated)) 

# For female with cubs, re-classify the ones that have a cub of less than 6 month (critic period)

os$With_cubs_estimated_new <- os$With_cubs_estimated

for (i in 1:nrow(os)) {
  if(is.na(os$With_cubs_estimated[i])) next
  if(os$With_cubs_estimated[i] == 1) { # First year cubs
    if(is.na(os$month[i])) next
    if(os$month[i] < 7) { # Before July
      os$With_cubs_estimated_new[i] <- "<6month" # Place <6 month if it has cubs (estimated number) that are below 6 months
    }}}

os <- os[,c(1:22,37,23:36)] # Check variable is well calculated
os <- os[ ,-c(22)]# Delete extra-columns and sort table
colnames(os)[colnames(os) %in% c("With_cubs_estimated_new")] <- "With_cubs_estimated"


## ---- Age class ----

os$Age2 <- os$Age # To convert to numeric to do classification (will produce NA)
os$Age2 <- as.numeric(os$Age2)

for (i in 1:nrow(os)) { 
  if(is.na(os$Age2[i]) | is.na(os$month[i])) next
  if(os$Age2[i] > 4) {
    os$Age_class[i] <- "Adult" } else if (os$Age2[i] <= 4 & os$Age2[i] > 1.5) {
    os$Age_class[i] <- "Subadult" } else if (os$Age2[i] == 1.5 | os$Age2[i] == 1) {
    os$Age_class[i] <- "Cub2" } else if (os$Age2[i] == 0 & os$month[i] > 7) {
    os$Age_class[i] <- "Cub1" } else if (os$Age2[i] == 0 & os$month[i] < 7) {
    os$Age_class[i] <- "Cub0" }
}

# Determine age class when month is missing

os_manage <- os[which(is.na(os$month) & !is.na(os$Age2)), ] # I will add age classes on this one
os_removed <- os[-which(is.na(os$month) & !is.na(os$Age2)), ] # I remove them and I will add them later when fixed

for (i in 1:nrow(os_manage)) { 
  if(os_manage$Age2[i] > 4) {
    os_manage$Age_class[i] <- "Adult" } else if (os_manage$Age2[i] <= 4 & os_manage$Age2[i] > 1.5) {
      os_manage$Age_class[i] <- "Subadult" } else if (os_manage$Age2[i] == 1.5 | os_manage$Age2[i] == 1) {
        os_manage$Age_class[i] <- "Cub2" } else if (os_manage$Age2[i] == 0) {os_manage$Age_class[i] <- "Cub1" }
  }

os <- rbind(os_removed, os_manage)

os <- os[,-c(35:37)] # Delete extra-columns

## ---- Add region and country ----

os_na <- os[which(is.na(os$X)), ] # Extract NA to join later
os <- os[-which(is.na(os$X)), ] # Remove NA to convert coordinates

coordinates(os) <- os[,c("x_long","y_lat")] # Spatial object
os@proj4string <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# Load layer with regions
os@data <- os@data[ ,-which(colnames(os@data) %in% c("Country", "Region"))]

map <- readOGR(dsn = "D:/MargSalas/Oso/Datos/GIS/Countries", layer = "clip_pyros2")
Encoding(map$NAME_1) <- "UTF-8"
d <- over(os,map) # Overlay to see where each point fall
os@data <- cbind(os@data, d)
colnames(os@data)[which(colnames(os@data) %in% c("NAME_0", "NAME_1"))] <- c("Country", "Region")

os@data <- os@data[,c(1:17,33:34,18:32)]

# Export GIS layer with coordinates and good attribute table (last version)
#writeOGR(os, "D:/MargSalas/Oso/Datos/GIS/2022/Seguiment_GIS_layer", "Seguiment_Ossos_Pirineus_1996_2021_coordinates_final", driver = "ESRI Shapefile")

os_data <- rbind(os@data, os_na)
os_data <- arrange(os_data, ID_obs) 


setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
#openxlsx::write.xlsx(os_data, 'Seguiment_Ossos_Pirineus_1996_2021_taula_final.xlsx')


## ---- Create GIS layers Females with cubs of < 6 month and 1 year ----
## ONLY MONITORING DATA!!

# Column of critical observations (only cubs or females with cubs of less than 6 month)

os$Female_cubs_critic <- ifelse(os$Age_class == "Cub0" | 
                                       os$With_cubs_estimated == "<6month",
                                     1,0)
os$Female_cubs_critic[is.na(os$Female_cubs_critic)] <- 0


# Column of observations of 1st year (cubs or females with cubs of 1 year, INCLUDING THE ONES OF < 6 MONTH)

os$Female_cubs_year1 <- ifelse(os$Age_class == "Cub0" | 
                                      os$Age_class == "Cub1" |
                                      os$With_cubs_estimated == "<6month" |
                                      os$With_cubs_estimated == "1"
                                      ,
                                     1,0)
os$Female_cubs_year1[is.na(os$Female_cubs_year1)] <- 0  


# Save GIS layers

os_critic <- os[which(os$Female_cubs_critic == 1), ]
os_year1 <- os[which(os$Female_cubs_year1 == 1), ]

writeOGR(os_critic, "D:/MargSalas/Oso/Datos/GIS/2022/Seguiment_GIS_layers_OCA", "Seguiment_Ossos_Pirineus_1996_2021_OCA_critic", driver = "ESRI Shapefile")
writeOGR(os_year1, "D:/MargSalas/Oso/Datos/GIS/2022/Seguiment_GIS_layers_OCA", "Seguiment_Ossos_Pirineus_1996_2021_OCA_1any", driver = "ESRI Shapefile")


# Plot maps

# Load basemap
map1 <- readOGR(dsn = "D:/MargSalas/Oso/Datos/GIS/Countries", layer = "clip_pyros2")

setwd("D:/MargSalas/Oso/Datos/Plots")

pdf("Seguiment1.pdf", 6.7,5.9)

par(mfrow = c(1,1),
    oma = c(2,2,4,1),
    mar = c(0,0.5,0,0))

plot(map1, col = "lightgrey", border = "grey")
mtext("Dades de seguiment de l'ós bru als Pirineus (1996-2020)", side = 3, line = -1, cex = 1.5)

points(os, col = "black", pch = 19, cex = 0.8)
points(os_year1, col = "orange", pch = 19, cex = 0.8)
points(os_critic, col = "red", pch = 19, cex = 0.8)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',
       legend = c("Totes les observacions", "Observacions sensibles (femelles amb cries de menys d'un any)", "Observacions crítiques (femelles amb cries de 6 mesos o menys)"), 
       col = c("black", "Orange", "red"), fill = c("black", "Orange", "red"), xpd = TRUE, 
       cex = 1, seg.len = 1, bty = 'n')

dev.off()
