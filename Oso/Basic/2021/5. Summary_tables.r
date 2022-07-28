
library(dplyr)
library(tidyr)

rm(list = ls())

# Load monitoring data
setwd("D:/Oso/Datos/Tablas_finales")
os <- read.csv("Seguiment_Ossos_Pirineus_1996_2020_taula_final.csv", header = TRUE, row.names = NULL)
os <- os[,-1] 

unique(os$Coordinate_system)

## ---- Summary table monitoring ----

# Create column for identified/unidentified

os$Identified <- 0
os$Identified <- ifelse(os$Confirmed_Individual != "Indetermined", os$Identified <- 1, os$Identified <- 0)

# Pool observations 1996-1999, 2000-2010, 2010-2020
os$pool_year <- "1996-1999"
for (i in 1:nrow(os)){
  if(os$Year[i] > 1999 & os$Year[i] < 2011){ os$pool_year[i] <- "2000-2010" } else if (os$Year[i] > 2010 ) {
    os$pool_year[i] <- "2010-2020" }
}


d <- as.data.frame(os %>% group_by(pool_year, Identified, Obs_type) %>% summarise(n()))
colnames(d)[4] <- "Number"
d <- spread(d, Obs_type, Number)

d[is.na(d)] <- 0

d <- t(d)

setwd("D:/Oso/Datos/Tablas_finales")
write.csv(d, "summary_monitoring.csv")

# Radiotracking

setwd("D:/Oso/Datos/GPS")
os_gps <- read.csv("Radiotracking_ossos_1996_2020_taula_final.csv", header = TRUE, row.names = NULL)
os_gps <- os_gps[,-1]

e <- as.data.frame(os_gps %>% group_by(Bear_name, Year) %>% summarise(n()))
colnames(e)[3] <- "Number"
e <- spread(e, Bear_name, Number)

e[is.na(e)] <- 0

setwd("D:/Oso/Datos/Tablas_finales")
write.csv(e, "summary_radiotracking.csv")
