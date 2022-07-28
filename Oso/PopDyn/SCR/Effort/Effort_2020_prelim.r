## -------------------------------------------------
##                 EFFORT 2020 HAIR
## ------------------------------------------------- 

rm(list = ls())


library(dplyr)
library(sp)
library(rgdal)

# Load basemap
map1 <- readOGR(dsn = "D:/MargSalas/Oso/Datos/GIS/Countries", layer = "clip_pyros2")

CRS_utm <- CRS("+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") # WGS 84 utm zone 31N
map1 <- spTransform(map1, CRS_utm)

## ---- Trap Data File ----

setwd("D:/MargSalas/Oso/Datos/Effort_raw")
ef <- read.csv("Revisions trampes 2020 V2.csv", sep = ";")
Encoding(ef$tipus_tr) <- "UTF-8" 

# For the traps that were located in january, I consider as if they where put in the first revision date
ef$data_pos[which(ef$data_pos == "01/01/2020")] <- ef$Data_rev1[which(ef$data_pos == "01/01/2020")]

# For the moment, keep only Codi_tr to identify. 
ef <- ef[,c(1,4,5,7,10,11,17:34)] # Keep relevant columns
ef$Year <- 2020

dates_effort <- ef[ ,c(4,7:24)]

ef <- ef[ ,c(1:3,25,5:6)]

# Colnames
colnames(ef)[1] <- "TrapID"
colnames(ef)[2] <- "Trap_type"
colnames(ef)[3] <- "Bait"


# Divide occasions

# Range of dates: min and max
dates_effort$data_pos <- as.Date(dates_effort$data_pos, format = "%d/%m/%Y") # MIN from posicionamiento trampas
min_date <- min(dates_effort$data_pos)

date <- unique(c(dates_effort$Data_rev1,dates_effort$Data_rev2,dates_effort$Data_rev3,dates_effort$Data_rev4,dates_effort$Data_rev5,dates_effort$Data_rev5,dates_effort$Data_rev6,dates_effort$Data_rev7,dates_effort$Data_rev8,dates_effort$Data_rev9))
date <- as.Date(date, format = "%d/%m/%Y")
date <- sort(date)
date <- date[complete.cases(date)]
max_date <- max(date)

# Period of time when the traps were open
all_days <- seq(min_date,max_date, by = "day")

e <- matrix(NA, nrow = nrow(ef), ncol = length(all_days))
row.names(e) <- ef$TrapID
colnames(e) <- as.character(all_days)

for (i in 1:nrow(e)){
  e[i,which(colnames(e) == dates_effort$data_pos[i])] <- 1 # Location of traps
  revisions <- as.Date(unlist(dates_effort[i,c(2,4,6,8,10,12,14,16,18)]), format = "%d/%m/%Y")
  
  e[i,which(colnames(e) %in% as.character(revisions))] <- 1
}

#  FILL 1 FROM FIRST TO LAST

known.state <- function(ch){
    state <- ch
    for (i in 1:dim(ch)[1]){
      n1 <- min(which(ch[i,]==1))
      n2 <- max(which(ch[i,]==1))
      state[i,n1:n2] <- 1
    }
    return(state)
  }
e1 <- known.state(e)


# HERE! 

# -> DIVIDE IN 15 DAYS SAMPLING OCCASIONS
# -> CREATE COVARIATE EFFORT WITH NUMBER OF VISITS


## ---- Encounter Data File ----

# Load monitoring data
setwd("D:/MargSalas/Oso/Datos/Tablas_finales")
os <- read.csv("Seguiment_Ossos_Pirineus_1996_2020_taula_final.csv", header = TRUE, row.names = NULL)
os <- os[which(os$Year == 2020 & os$Probable_Individual != "Indetermined"),-1]
os <- os[which(os$Method == "Sampling_station" & os$Obs_type == "Hair"), ]

coordinates(os) <- os[,c("x_long","y_lat")] # Spatial object
os@proj4string <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # Set coord system (WGS lat_long)


os <- spTransform(os, CRS_utm) #Transform to WGS utm

# Effort as spatial object to do spatial overlay
coordinates(ef) <- ef[,c("coord_x","coord_y")] # Spatial object
ef@proj4string <- CRS("+proj=utm +zone=31 +ellps=GRS80 +units=m +no_defs") # I THINK IS ETRS89 31N, BUT NOT SURE
ef <- spTransform(ef, CRS_utm)

#Plot
ext <- bbox(ef)
plot(map1, col = "lightgrey", border = "grey", 
     xlim = c((ext[1,1]),(ext[1,2])), ylim = c((ext[2,1]),(ext[2,2])))
points(ef, pch = 4)
points(os, pch = 18, col = adjustcolor("red", alpha.f = 0.4))

# Overlay: It is like the over function but I need to do it with character because it doesnt work
# This works when traps are in the exactly same location than observations

# Which traps have bear observations (and which obervations)
withWhich <-list()
for(i in 1:length(ef)){
  
  withWhich[[i]] <-  which(as.character(coordinates(os)[,1]) %in% as.character(coordinates(ef)[i,1]) & 
                           as.character(coordinates(os)[,2]) %in% as.character(coordinates(ef)[i,2]))
  
}

# Assign trap name to observations
os$TrapID <- NA
os$TrapID <- as.character(os$TrapID)
for (i in 1:length(withWhich)){
  if(length(withWhich[[i]]) == 0) next
  os$TrapID[withWhich[[i]]]<- ef[i,]$TrapID
}

length(which(!is.na(os$TrapID)))

# Only 19 out of 64 have a trap assigned...what do I do with the rest?
dat <- os@data


       