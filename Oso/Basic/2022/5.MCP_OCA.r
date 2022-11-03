
## -------------------------------------------------
##          Analysis for Anna planella
##      MCP area females with cubs critic period
## ------------------------------------------------- 

rm(list = ls())

library(tidyverse)
library(adehabitatHR)
library(rgeos)

# Load database of monitoring + radiotracking
# - Without independent cubs observations (doesn't make sense, only interested in mothers)
# - Without natal_established locations (not interested in natal territories)

setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
os <- read.csv("Data_os_96_21.csv")

coordinates(os) <- os[,c("x_long","y_lat")] # Spatial object
os@proj4string <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# Column of critical observations (females with cubs of less than 6 month)

os$Female_cubs_critic <- ifelse(os$With_cubs_estimated == "<6month",1,0)
os$Female_cubs_critic[is.na(os$Female_cubs_critic)] <- 0


# Column of observations of 1st year (cubs or females with cubs of 1 year, INCLUDING THE ONES OF < 6 MONTH)

os$Female_cubs_year1 <- ifelse(os$With_cubs_estimated == "<6month" | os$With_cubs_estimated == "1",
                               1,0)
os$Female_cubs_year1[is.na(os$Female_cubs_year1)] <- 0  

os_critic <- os[which(os$Female_cubs_critic == 1), ]
os_year1 <- os[which(os$Female_cubs_year1 == 1), ]

#setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
#openxlsx::write.xlsx(os_critic@data, 'OCA_critic_allmethods.xlsx')
#openxlsx::write.xlsx(os_year1@data, 'OCA_1year_allmethods.xlsx')

# Calculate MCP females with cubs < 6 month

os_critic <- os_critic[which(os_critic$Age_class %in% c("Adult", "Subadult") & os_critic$Confirmed_Individual != "Indetermined" & os_critic$Sex %in% c("F")), c(2,6) ]
os_critic <- spTransform(os_critic, CRS('+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'))

# resumen num. de observaciones
sum <- os_critic@data %>% group_by(Confirmed_Individual, Year) %>%
  summarise(n())

# Hay muy pocas localizaciones por año para estimar un MCP anual por bicho, no tiene sentido:

yr <- unique(os_critic$Year)
id <- unique(os_critic$Confirmed_Individual)

mcps <- data.frame(matrix(NA,nrow = length(id), ncol = length(yr)))
rownames(mcps) <- id
colnames(mcps) <- yr

for (i in 1:length(id)){
  ind <- os_critic[which(os_critic$Confirmed_Individual %in% id[i]), ]
  for (t in 1:length(yr)){
    ind_year <- ind[which(ind$Year %in% yr[t]),1]
    if(nrow(ind_year) < 10) next
    mcpind <- mcp(ind_year)
    mcps[which(rownames(mcps) %in% id[i]),which(colnames(mcps) %in% yr[t])] <- round(gArea(mcpind)/1000000,3) # in km2
  }
}

mcps

# O cogemos todo el radiotracking, o juntamos todos los años

#### COGIENDO SOLO RADIOTRACKING ###

os_critic <- os[which(os$Female_cubs_critic == 1), ]
os_critic <- os_critic[which(os_critic$Age_class %in% c("Adult", "Subadult") & 
                               os_critic$Confirmed_Individual != "Indetermined" & 
                               os_critic$Sex %in% c("F") &
                               os_critic$Method == "Radiotracking"
                               ), c(2,6) ]
os_critic <- spTransform(os_critic, CRS('+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'))

# resumen num. de observaciones
sum_radio <- os_critic@data %>% group_by(Confirmed_Individual, Year) %>%
  summarise(n())
sum_radio

id <- unique(os_critic$Confirmed_Individual)

mcps_radiotracking <- data.frame(matrix(NA,nrow = length(id), ncol = 1))
rownames(mcps_radiotracking) <- id
colnames(mcps_radiotracking) <- "area"

for (i in 1:length(id)){
  ind <- os_critic[which(os_critic$Confirmed_Individual %in% id[i]),1]
  mcpind <- mcp(ind)
  mcps_radiotracking[which(rownames(mcps_radiotracking) %in% id[i]),1] <- round(gArea(mcpind)/1000000,3) # in km2
  }

mcps_radiotracking 
# Al menos estos datos tienen sentido


#### JUNTANDO TODOS LOS AÑOS ###

os_critic <- os[which(os$Female_cubs_critic == 1), ]
os_critic <- os_critic[which(os_critic$Age_class %in% c("Adult", "Subadult") & 
                               os_critic$Confirmed_Individual != "Indetermined" & 
                               os_critic$Sex %in% c("F")), c(2,6) ]
os_critic <- spTransform(os_critic, CRS('+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'))

# resumen num. de observaciones
sum_allyears <- os_critic@data %>% group_by(Confirmed_Individual) %>%
  summarise(n())
sum_allyears

id <- unique(os_critic$Confirmed_Individual)

mcps_allyears <- data.frame(matrix(NA,nrow = length(id), ncol = 1))
rownames(mcps_allyears) <- id
colnames(mcps_allyears) <- "area"

for (i in 1:length(id)){
  ind <- os_critic[which(os_critic$Confirmed_Individual %in% id[i]),1]
  if(nrow(ind) < 10) next # al menos 10 obs para calcular home range
  mcpind <- mcp(ind)
  mcps_allyears[which(rownames(mcps_allyears) %in% id[i]),1] <- round(gArea(mcpind)/1000000,3) # in km2
}

mcps_allyears 

# COMPARANDO, mcps y numero de posiciones con la que está calculado

mcps_radiotracking
sum_radio

mcps_allyears
sum_allyears

