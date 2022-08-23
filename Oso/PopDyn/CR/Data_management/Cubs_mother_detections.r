


setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
os <- read.csv("Data_os_96_21_cubLocations.csv", header = TRUE, row.names = NULL)
os <- os[,-1] 

os_id <- os[which(os$Confirmed_Individual != "Indetermined"), ] # Only identified
os_id <- os_id[-which(os_id$Confirmed_Individual == "Camille_AspeOuest"), ] # Remove because I don't know age

