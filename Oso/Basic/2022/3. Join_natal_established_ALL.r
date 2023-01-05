
## -------------------------------------------------
##    Join positions of mother to natal positions cubs
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
os <- openxlsx::read.xlsx('Seguiment_Ossos_Pirineus_1996_2021_taula_final.xlsx')

# To check 
#checks <- os[sample(nrow(os), 60), ]
#write.csv(checks, "check.csv")

os <- os[-which(is.na(os$X)), ] # Remove NA
os <- os[-which(os$Coordinates.precision == 3 | os$Coordinates.precision == 2), ] # Remove coordinates precision 2 and 3

# Load radiotracking data
setwd("D:/MargSalas/Oso/Datos/GPS")
os_gps <- read.csv("Radiotracking_ossos_1996_2020_taula_final.csv", header = TRUE, row.names = NULL)
os_gps <- os_gps[,-1]
os_gps$Bear_name[which(os_gps$Bear_name == "Melba")] <- "Mellba" # Change name to fit with 2021

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
#write.csv(os@data, file = "Data_os_96_21.csv")


# Load basemap
map1 <- readOGR(dsn = "D:/MargSalas/Oso/Datos/GIS/Countries", layer = "clip_pyros2")

# Load info database
setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
info <- readxl::read_xlsx("D:/MargSalas/Oso/Datos/Tablas_finales/2022/Info_individuals_2021.xlsx", sheet = 1)
info <- info[,c(4:8,10)]
colnames(info)[6] <- "Year_death"

## ---- Summary identified individuals ----
osdat <- os@data

# Nº obs per individual, and identify those with > 9 obs
os_nobs <- data.frame(osdat %>% group_by(Confirmed_Individual) %>% summarise(n()))
id_more10 <- os_nobs$Confirmed_Individual[os_nobs$n.. > 9]
setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
save(id_more10, file = "id_more10.RData") # Save for future plots

# Nº obs per year
os_sum <- data.frame(osdat %>% group_by(Confirmed_Individual, Year) %>% summarise(n()))
colnames(os_sum)[3] <- "n_obs"
os_sum <- spread(os_sum, Year, n_obs)
os_sum[is.na(os_sum)] <- 0


## ---- Create new dataset adding the positions of the mother the year of birth ----

os <- os[-which(os$Confirmed_Individual == "Indetermined"), ]  # Remove indetermined

os$Confirmed_Individual_cub <- os$Confirmed_Individual # New column to indicate the cub name with the positions of the mother
os$Age_class_cub <- os$Age_class # New column to indicate the age class of the cub
os$Sex_cub <- os$Sex # New column to indicate the sex of the cub

# In the new dataset, we will use the column Confirmed_Individual_cub and Age_class_cub, 
# that contains the information of all mothers and cubs in the first year of birth
# This will be use only when looking at individual data, is it is repeated information (positions mother = positions cubs)

id <- unique(os$Confirmed_Individual)

for (i in 1:length(id)){
  
  #1. Identify the locations of the mother in the first year of the individual
  mother_id <- info$Mother[info$ID %in% id[i]] # Get ID of the mother
  if(is.na(mother_id)) next
  year_birth <- info$Year_birth[info$ID %in% id[i]] # Get year of birth
  locs_mother1 <- os[which(os$Confirmed_Individual_cub %in% mother_id &  # Select locs of mother the year of birth
                            os$Year %in% year_birth), ]                   # It is important to select Confirmed_Individual_cub, otherwise it selects the ones already joined in the previous individual i
  if(nrow(locs_mother1) == 0) next # If there are not locations of the mother the 1st year, pass onto the next individual
  locs_mother2 <- locs_mother1[which(locs_mother1$With_cubs_estimated == "<6month" | # Select only positions of the mother with 1st year cubs
                                       locs_mother1$With_cubs_estimated == "1"), ]
  
  #2. Assign new information to the locations
  locs_mother2$Confirmed_Individual_cub <- id[i] # Assign the id of the cub to the observation
  
  for (j in 1:nrow(locs_mother2)){              # Assign age class of the cub to the observation
    if(locs_mother2$With_cubs_estimated[j] == "<6month") { 
      locs_mother2$Age_class_cub[j] <- "Cub0" } else if (locs_mother2$With_cubs_estimated[j] == "1") {
        locs_mother2$Age_class_cub[j] <- "Cub1" }
    }
  
  locs_mother2$Sex_cub <- info$Sex[info$ID %in% id[i]]
    
  #3. Join this rows to the dataset
  os <- rbind(os,locs_mother2)
  }

# Save database to check and calculate natal and established centroids
os@data <- os@data[,c(23,25,24,5:21,1:4)]


setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
#write.csv(os@data, "Natal_established_coordinates.csv")


# Save different GIS layer per individual to explore in arcgis (including positions of the mother their first year)

setwd("D:/MargSalas/Oso/Datos/GIS/2022/Identified_individuals")
for (i in 1:length(id)){
  os_id <- os[which(os$Confirmed_Individual_cub == id[i]), ] # Select individual
  writeOGR(os_id, "D:/MargSalas/Oso/Datos/GIS/2022/Identified_individuals", paste("ID_", id[i], sep = ""), driver = "ESRI Shapefile")
}


## -------------------------------------------------
##        Join to database with all cub obsevations
## ------------------------------------------------- 

# Load monitoring data
setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
os <- openxlsx::read.xlsx('Seguiment_Ossos_Pirineus_1996_2021_taula_final_cubLocations.xlsx')

# To check 
#checks <- os[sample(nrow(os), 60), ]
#write.csv(checks, "check.csv")

os <- os[-which(is.na(os$X)), ] # Remove NA
os <- os[-which(os$Coordinates.precision == 3 | os$Coordinates.precision == 2), ] # Remove coordinates precision 2 and 3

# Load radiotracking data
setwd("D:/MargSalas/Oso/Datos/GPS")
os_gps <- read.csv("Radiotracking_ossos_1996_2020_taula_final.csv", header = TRUE, row.names = NULL)
os_gps <- os_gps[,-1]
os_gps$Bear_name[which(os_gps$Bear_name == "Melba")] <- "Mellba" # Change name to fit with 2021

## ---- JOIN SPATIAL DATA FROM BOTH DATASETS ----

# Keep column names common to both
colnames(os)
names_os <- c("ID_obs","Confirmed_Individual", "Sex", "Age","Age_class","Year","Date",
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
os_gps$ID_obs <- NA
os_gps <- os_gps[ ,c(22,1:6,9:12,8,13:19,21,7,20)] # Change column order as in monitoring db


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
#write.csv(os@data, file = "Data_os_96_21_cubLocations.csv")

# Load basemap
map1 <- readOGR(dsn = "D:/MargSalas/Oso/Datos/GIS/Countries", layer = "clip_pyros2")

# Load info database
setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
info <- readxl::read_xlsx("D:/MargSalas/Oso/Datos/Tablas_finales/2022/Info_individuals_2021.xlsx", sheet = 1)
info <- info[,c(4:8,10)]
colnames(info)[6] <- "Year_death"


## ---- Create new dataset adding the positions of the mother the year of birth ----

os <- os[-which(os$Confirmed_Individual == "Indetermined"), ]  # Remove indetermined

os$Confirmed_Individual_cub <- os$Confirmed_Individual # New column to indicate the cub name with the positions of the mother
os$Age_class_cub <- os$Age_class # New column to indicate the age class of the cub
os$Sex_cub <- os$Sex # New column to indicate the sex of the cub

# In the new dataset, we will use the column Confirmed_Individual_cub and Age_class_cub, 
# that contains the information of all mothers and cubs in the first year of birth
# This will be use only when looking at individual data, is it is repeated information (positions mother = positions cubs)

os <- os[-which(os$Confirmed_Individual == "Moonbots/Bonabe"), ]
id <- unique(os$Confirmed_Individual)

for (i in 1:length(id)){
  
  #1. Identify the locations of the mother in the first year of the individual
  mother_id <- info$Mother[info$ID %in% id[i]] # Get ID of the mother
  if(is.na(mother_id)) next
  year_birth <- info$Year_birth[info$ID %in% id[i]] # Get year of birth
  locs_mother1 <- os[which(os$Confirmed_Individual_cub %in% mother_id &  # Select locs of mother the year of birth
                             os$Year %in% year_birth), ]                   # It is important to select Confirmed_Individual_cub, otherwise it selects the ones already joined in the previous individual i
  if(nrow(locs_mother1) == 0) next # If there are not locations of the mother the 1st year, pass onto the next individual
  locs_mother2 <- locs_mother1[which(locs_mother1$With_cubs_estimated == "<6month" | # Select only positions of the mother with 1st year cubs
                                       locs_mother1$With_cubs_estimated == "1"), ]
  
  #2. Assign new information to the locations
  locs_mother2$Confirmed_Individual_cub <- id[i] # Assign the id of the cub to the observation
  
  for (j in 1:nrow(locs_mother2)){              # Assign age class of the cub to the observation
    if(locs_mother2$With_cubs_estimated[j] == "<6month") { 
      locs_mother2$Age_class_cub[j] <- "Cub0" } else if (locs_mother2$With_cubs_estimated[j] == "1") {
        locs_mother2$Age_class_cub[j] <- "Cub1" }
  }
  
  locs_mother2$Sex_cub <- info$Sex[info$ID %in% id[i]]
  
  #3. Join this rows to the dataset
  os <- rbind(os,locs_mother2)
}

# Save database to check and calculate natal and established centroids
os <- os[,c(1,24,26,25,6:22,2:5)]

setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
openxlsx::write.xlsx(os@data, 'Natal_established_coordinates_cubLocations.xlsx')

# Load the one where Elena added the new age classes

setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
os <- readxl::read_xlsx("D:/MargSalas/Oso/Datos/Tablas_finales/2022/Natal_established_coordinates_cubLocations_EPS.xlsx", sheet = 1)
id <- unique(os$Confirmed_Individual_cub)
coordinates(os) <- os[,c("x_long","y_lat")] # Spatial object
os@proj4string <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")


# Save different GIS layer per individual to explore in arcgis (including positions of the mother their first year)

setwd("D:/MargSalas/Oso/Datos/GIS/2022/Identified_individuals")
for (i in 1:length(id)){
  os_id <- os[which(os$Confirmed_Individual_cub == id[i]), ] # Select individual
  writeOGR(os_id, "D:/MargSalas/Oso/Datos/GIS/2022/Identified_individuals", paste("ID_", id[i], sep = ""), driver = "ESRI Shapefile")
}


