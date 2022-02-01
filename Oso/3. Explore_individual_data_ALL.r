
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
setwd("D:/Oso/Datos/Tablas_finales")
os <- read.csv("Seguiment_Ossos_Pirineus_1996_2020_taula_final.csv", header = TRUE, row.names = NULL)
os <- os[,-1] 

# To check 
#checks <- os[sample(nrow(os), 60), ]
#write.csv(checks, "check.csv")

os <- os[-which(is.na(os$X)), ] # Remove NA
os <- os[-which(os$Coordinates.precision == 3), ] # Remove coordinates precision 3

# Load radiotracking data
setwd("D:/Oso/Datos/GPS")
os_gps <- read.csv("Radiotracking_ossos_1996_2020_taula_final.csv", header = TRUE, row.names = NULL)
os_gps <- os_gps[,-1]

## ---- JOIN SPATIAL DATA FROM BOTH DATASETS ----

# Keep column names common to both
colnames(os)
names_os <- c("Confirmed_Individual", "Sex", "Age","Age_class","Year","Date_register",
         "X","Y","x_long","y_lat","Coordinate_system","Country","Region","With_cubs_estimated","N_cubs_estimated",
         "ID_cub1","ID_cub2","ID_cub3", "Method", "Obs_type", "Remarks")
unique(os$Method)

os <- os[,colnames(os) %in% names_os]

colnames(os_gps)
names_os_gps <- c("Bear_name","Sex","Age","Age_class","Year","Date_GMT","Tracking_system","Coordinate_system",
                  "X","Y", "x_long","y_lat","Country","Region","With_cubs_estimated", "N_cubs_estimated", "ID_cub1",
                  "ID_cub2", "ID_cub3", "Remarks")
os_gps <- os_gps[,colnames(os_gps) %in% names_os_gps]
colnames(os_gps)[c(1,6,7)] <- c("Confirmed_Individual", "Date_register", "Obs_type") # Change column names as in monitoring db
os_gps$Method <- "Radiotracking"
os_gps <- os_gps[ ,c(1:6,9:12,8,13:19,21,7,20)] # Change column order as in monitoring db

# Join
os_all <- rbind(os,os_gps)

# Arrange by date in case we plot it by colors
os_all$Date_register_formated <- os_all$Date_register 
os_all$Date_register_formated <- as.Date(os_all$Date_register_formated, format = "%d/%m/%Y")

os_all <- arrange(os_all, Date_register_formated, Year) 

# Set coordinates
coordinates(os_all) <- os_all[,c("x_long","y_lat")] # Spatial object
os_all@proj4string <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

os <- os_all # Call it os to make it fit with code script 3.1

# Load basemap
map1 <- readOGR(dsn = "D:/Oso/Datos/GIS/Countries", layer = "clip_pyros2")

# Load info database
setwd("D:/Oso/Datos")
info <- read.csv("Info_individuals.csv", header = TRUE, row.names = NULL, sep = ";")
info <- info[,c(4:9)]


## ---- Summary identified individuals ----
osdat <- os@data

# Nº obs per individual, and identify those with > 10 obs
os_nobs <- data.frame(osdat %>% group_by(Confirmed_Individual) %>% summarise(n()))
id_more10 <- os_nobs$Confirmed_Individual[os_nobs$n.. > 9]
setwd("D:/Oso/Datos")
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


setwd("D:/Oso/Datos/Tablas_finales")
write.csv(os@data, "Natal_established_coordinates.csv")


# Save different GIS layer per individual to explore in arcgis (including positions of the mother their first year)

setwd("D:/Oso/Datos/GIS/Identified_individuals")
for (i in 1:length(id)){
  os_id <- os[which(os$Confirmed_Individual_cub == id[i]), ] # Select individual
  writeOGR(os_id, "D:/Oso/Datos/GIS/Identified_individuals", paste("ID_", id[i], sep = ""), driver = "ESRI Shapefile")
}


## ---- Plot observations identified individuals ----

#### MALES ####

os_males <- os[which(os$Sex_cub == "M"), ]
id_males <- unique(os_males$Confirmed_Individual_cub)
id_males <- id_males[id_males %in% id_more10] # ID males with more than 10 observations (all)

# Visual 1 (Three categories: Cub, subadult, adult) 

setwd("D:/Oso/Datos/Plots")

pdf("ID_Males1.pdf")  

par(mfrow = c(3,3),
    oma = c(2,2,4,1),
    mar = c(0,0.5,1,0))

for (i in 1:9) {
  os_id <- os_males[which(os_males$Confirmed_Individual_cub == id_males[i]), ] # Select individual
  os_id_adult <- os_id[which(os_id$Age_class_cub == "Adult"), ]
  os_id_subadult <- os_id[which(os_id$Age_class_cub == "Subadult"), ]
  os_id_cub <- os_id[which(os_id$Age_class_cub == "Cub0" | os_id$Age_class_cub == "Cub1" | os_id$Age_class_cub == "Cub2"), ]
  
  plot(map1, col = "lightgrey", border = "grey")
  mtext(id_males[i], side = 3, line = 0.5, cex = 1)
  mtext(paste("(",info$Year_birth[info$ID %in% id_males[i]], "-", info$Year_death[info$ID %in% id_males[i]], ")"), side = 3, line = -1, cex = 0.8)
  
  points(os_id_adult, col = adjustcolor("purple", alpha.f = 1), pch = 19, cex = 0.8)
  points(os_id_subadult, col = adjustcolor("orange", alpha.f = 0.8), pch = 19, cex = 0.8)
  points(os_id_cub, col = adjustcolor("darkgreen", alpha.f = 1), lwd = 1.8, pch = 8, cex = 0.8)
  
  mtext("Males", side = 3, line = 1.5, cex = 1.2, outer = TRUE)
  
}

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',legend = c("Cub", "Subadult", "Adult"), col = c("darkgreen", "Orange", "Purple"), fill = c("darkgreen", "orange", "purple"), xpd = TRUE, horiz = TRUE, cex = 1.5, seg.len=1, bty = 'n')

par(mfrow = c(3,3),
    oma = c(2,2,4,1),
    mar = c(0,0.5,1,0))

for (i in 10:18) {
  os_id <- os_males[which(os_males$Confirmed_Individual_cub == id_males[i]), ] # Select individual
  os_id_adult <- os_id[which(os_id$Age_class_cub == "Adult"), ]
  os_id_subadult <- os_id[which(os_id$Age_class_cub == "Subadult"), ]
  os_id_cub <- os_id[which(os_id$Age_class_cub == "Cub0" | os_id$Age_class_cub == "Cub1" | os_id$Age_class_cub == "Cub2"), ]
  
  plot(map1, col = "lightgrey", border = "grey")
  mtext(id_males[i], side = 3, line = 0.5, cex = 1)
  mtext(paste("(",info$Year_birth[info$ID %in% id_males[i]], "-", info$Year_death[info$ID %in% id_males[i]], ")"), side = 3, line = -1, cex = 0.8)
  
  points(os_id_adult, col = adjustcolor("purple", alpha.f = 1), pch = 19, cex = 0.8)
  points(os_id_subadult, col = adjustcolor("orange", alpha.f = 0.8), pch = 19, cex = 0.8)
  points(os_id_cub, col = adjustcolor("darkgreen", alpha.f = 1), lwd = 1.8, pch = 8, cex = 0.8)
  
  #mtext("Males", side = 3, line = 1.5, cex = 1.2, outer = TRUE)
  
}

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',legend = c("Cub", "Subadult", "Adult"), col = c("darkgreen", "Orange", "Purple"), fill = c("darkgreen", "orange", "purple"), xpd = TRUE, horiz = TRUE, cex = 1.5, seg.len=1, bty = 'n')

par(mfrow = c(3,3),
    oma = c(2,2,4,1),
    mar = c(0,0.5,1,0))

for (i in 19:21) {
  os_id <- os_males[which(os_males$Confirmed_Individual_cub == id_males[i]), ] # Select individual
  os_id_adult <- os_id[which(os_id$Age_class_cub == "Adult"), ]
  os_id_subadult <- os_id[which(os_id$Age_class_cub == "Subadult"), ]
  os_id_cub <- os_id[which(os_id$Age_class_cub == "Cub0" | os_id$Age_class_cub == "Cub1" | os_id$Age_class_cub == "Cub2"), ]
  
  plot(map1, col = "lightgrey", border = "grey")
  mtext(id_males[i], side = 3, line = 0.5, cex = 1)
  mtext(paste("(",info$Year_birth[info$ID %in% id_males[i]], "-", info$Year_death[info$ID %in% id_males[i]], ")"), side = 3, line = -1, cex = 0.8)
  
  points(os_id_adult, col = adjustcolor("purple", alpha.f = 1), pch = 19, cex = 0.8)
  points(os_id_subadult, col = adjustcolor("orange", alpha.f = 0.8), pch = 19, cex = 0.8)
  points(os_id_cub, col = adjustcolor("darkgreen", alpha.f = 1), lwd = 1.8, pch = 8, cex = 0.8)
  
  #mtext("Males", side = 3, line = 1.5, cex = 1.2, outer = TRUE)
  
}

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend(-0.5,0.35,legend = c("Cub", "Subadult", "Adult"), col = c("darkgreen", "Orange", "Purple"), fill = c("darkgreen", "orange", "purple"), xpd = TRUE, horiz = TRUE, cex = 1.5, seg.len=1, bty = 'n')

dev.off()


# Visual 2: Three categories (Cub, Subadult, Adult) and gradient of colours within

color_cub <- c("lightyellow1", "yellow", "gold", "orange") # Create three color palettes
col_subadult <- c("darkorange1", "orangered","firebrick1", "red1")
color_adult <- c("red3", "red4", "darkred", "black")


pdf("ID_Males2.pdf")  

par(mfrow = c(3,3),
    oma = c(2,2,4,1),
    mar = c(0,0.5,1,0))

for (i in 1:length(id_males)) {
  os_id <- os_males[which(os_males$Confirmed_Individual_cub == id_males[i]), ] # Select individual
  os_id_adult <- os_id[which(os_id$Age_class_cub == "Adult"), ]
  os_id_subadult <- os_id[which(os_id$Age_class_cub == "Subadult"), ]
  os_id_cub <- os_id[which(os_id$Age_class_cub == "Cub0" | os_id$Age_class == "Cub1" | os_id$Age_class == "Cub2"), ]
  
  plot(map1, col = "lightgrey", border = "grey")
  mtext(id_males[i], side = 3, line = 0.5, cex = 1)
  mtext(paste("(",info$Year_birth[info$ID %in% id_males[i]], "-", info$Year_death[info$ID %in% id_males[i]], ")"), side = 3, line = -1, cex = 0.8)

  points(os_id_adult, col = color_adult, pch = 19, cex = 0.8)
  points(os_id_subadult, col = col_subadult, pch = 19, cex = 0.8)
  points(os_id_cub, col = color_cub, pch = 19, cex = 0.8)
  
  mtext("Males", side = 3, line = 1.5, cex = 1.2, outer = TRUE)
}

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
#legend('center',legend = c("Cub", "Subadult", "Adult"), col = c("darkgreen", "Orange", "Purple"), fill = c("Green", "Orange", "Purple"), xpd = TRUE, horiz = TRUE, cex = 1.5, seg.len=1, bty = 'n')

dev.off()


#### FEMALES ####

os_females <- os[which(os$Sex_cub == "F"), ]
id_females <- unique(os_females$Confirmed_Individual_cub)
id_females <- id_females[id_females %in% id_more10] # ID females with more than 10 observations (all)

pdf("ID_Females1.pdf")  

par(mfrow = c(3,3),
    oma = c(2,2,4,1),
    mar = c(0,0.5,1,0))

for (i in 1:9) {
  os_id <- os_females[which(os_females$Confirmed_Individual_cub == id_females[i]), ] # Select individual
  os_id_adult <- os_id[which(os_id$Age_class_cub == "Adult"), ]
  os_id_subadult <- os_id[which(os_id$Age_class_cub == "Subadult"), ]
  os_id_cub <- os_id[which(os_id$Age_class_cub == "Cub0" | os_id$Age_class_cub == "Cub1" | os_id$Age_class_cub == "Cub2"), ]
  
  plot(map1, col = "lightgrey", border = "grey")
  mtext(id_females[i], side = 3, line = 0.5, cex = 1)
  mtext(paste("(",info$Year_birth[info$ID %in% id_females[i]], "-", info$Year_death[info$ID %in% id_females[i]], ")"), side = 3, line = -1, cex = 0.8)
  
  points(os_id_adult, col = adjustcolor("purple", alpha.f = 1), pch = 19, cex = 0.8)
  points(os_id_subadult, col = adjustcolor("orange", alpha.f = 0.8), pch = 19, cex = 0.8)
  points(os_id_cub, col = adjustcolor("darkgreen", alpha.f = 1), lwd = 1.8, pch = 8, cex = 0.8)
  
  mtext("Females", side = 3, line = 1.5, cex = 1.2, outer = TRUE)
  
}

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',legend = c("Cub", "Subadult", "Adult"), col = c("darkgreen", "orange", "purple"), fill = c("darkgreen", "orange", "purple"), xpd = TRUE, horiz = TRUE, cex = 1.5, seg.len=1, bty = 'n')

par(mfrow = c(3,3),
    oma = c(2,2,4,1),
    mar = c(0,0.5,1,0))

for (i in 10:18) {
  os_id <- os_females[which(os_females$Confirmed_Individual_cub == id_females[i]), ] # Select individual
  os_id_adult <- os_id[which(os_id$Age_class_cub == "Adult"), ]
  os_id_subadult <- os_id[which(os_id$Age_class_cub == "Subadult"), ]
  os_id_cub <- os_id[which(os_id$Age_class_cub == "Cub0" | os_id$Age_class_cub == "Cub1" | os_id$Age_class_cub == "Cub2"), ]
  
  plot(map1, col = "lightgrey", border = "grey")
  mtext(id_females[i], side = 3, line = 0.5, cex = 1)
  mtext(paste("(",info$Year_birth[info$ID %in% id_females[i]], "-", info$Year_death[info$ID %in% id_females[i]], ")"), side = 3, line = -1, cex = 0.8)
  
  points(os_id_adult, col = adjustcolor("purple", alpha.f = 1), pch = 19, cex = 0.8)
  points(os_id_subadult, col = adjustcolor("orange", alpha.f = 0.8), pch = 19, cex = 0.8)
  points(os_id_cub, col = adjustcolor("darkgreen", alpha.f = 1), lwd = 1.8, pch = 8, cex = 0.8)
  
  #mtext("Females", side = 3, line = 1.5, cex = 1.2, outer = TRUE)
  
}

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',legend = c("Cub", "Subadult", "Adult"), col = c("darkgreen", "Orange", "Purple"), fill = c("darkgreen", "orange", "purple"), xpd = TRUE, horiz = TRUE, cex = 1.5, seg.len=1, bty = 'n')

par(mfrow = c(3,3),
    oma = c(2,2,4,1),
    mar = c(0,0.5,1,0))

for (i in 19:24) {
  os_id <- os_females[which(os_females$Confirmed_Individual_cub == id_females[i]), ] # Select individual
  os_id_adult <- os_id[which(os_id$Age_class_cub == "Adult"), ]
  os_id_subadult <- os_id[which(os_id$Age_class_cub == "Subadult"), ]
  os_id_cub <- os_id[which(os_id$Age_class_cub == "Cub0" | os_id$Age_class_cub == "Cub1" | os_id$Age_class_cub == "Cub2"), ]
  
  plot(map1, col = "lightgrey", border = "grey")
  mtext(id_females[i], side = 3, line = 0.5, cex = 1)
  mtext(paste("(",info$Year_birth[info$ID %in% id_females[i]], "-", info$Year_death[info$ID %in% id_females[i]], ")"), side = 3, line = -1, cex = 0.8)
  
  points(os_id_adult, col = adjustcolor("purple", alpha.f = 1), pch = 19, cex = 0.8)
  points(os_id_subadult, col = adjustcolor("orange", alpha.f = 0.8), pch = 19, cex = 0.8)
  points(os_id_cub, col = adjustcolor("darkgreen", alpha.f = 1), lwd = 1.8, pch = 8, cex = 0.8)
  
  #mtext("Females", side = 3, line = 1.5, cex = 1.2, outer = TRUE)
  
}

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend(-0.5,-0.3,legend = c("Cub", "Subadult", "Adult"), col = c("darkgreen", "Orange", "Purple"), fill = c("darkgreen", "orange", "purple"), xpd = TRUE, horiz = TRUE, cex = 1.5, seg.len=1, bty = 'n')

dev.off()
