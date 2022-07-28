


## -------------------------------------------------
##          Explore individual dispersal 
## ------------------------------------------------- 

# Add to the individual layers the information on the natal territory (locations of the mum the 1st year)
# or the preliminary established territories (as locations that overlay)


rm(list = ls())

library(rgdal)
library(tidyverse)

# Load layers 

# Load monitoring data
setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
os <- read.csv("Natal_established_coordinates.csv", header = TRUE, row.names = NULL)
os <- os[,-1] 

# Set coordinates
coordinates(os) <- os[,c("x_long","y_lat")] # Spatial object
os@proj4string <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# ID with more than 10 0bs
setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2022")
load("id_more10.RData")

# Load mcp
setwd("D:/MargSalas/Oso/Datos/GIS/2022")
mcp_natal <- readRDS("mcp_natal.rds")
mcp_est1 <- readRDS("mcp_est1.rds")
mcp_est2 <- readRDS("mcp_est2.rds")

names_mcp_nat <- lapply(mcp_natal, `[[`, 1) # Take only the first argument of the list to know the name
names_mcp_est1 <- lapply(mcp_est1, `[[`, 1)
names_mcp_est2 <- lapply(mcp_est2, `[[`, 1)


# Save different GIS layer per individual to explore in arcgis (including positions of the mother their first year)

id <- unique(os$Confirmed_Individual)

id_loc <- list()

for (i in 1:length(id)){
  
  os_id <- os[which(os$Confirmed_Individual_cub == id[i]), ] # Select individual
  os_id@data$period <- "transit" 
  
  # 1. Identify locations belonging to NATAL mcp
  
  mcp_nat_id <- which(names_mcp_nat %in% id[i]) # MCP natal
  ifelse(length(mcp_nat_id) > 0, x <- mcp_natal[[which(names_mcp_nat %in% id[i])]], x <- 0) # store in x
  
  if (length(x) > 1) { # If there is a natal territory 
  os_id <- cbind(os_id,over(os_id,x[[2]])) # Select points that fall within natal mcp
  
  colnames(os_id@data)[which(colnames(os_id@data) %in% "id")] <- "natal_loc"
  
  os_id@data$period <- ifelse(is.na(os_id@data$natal_loc), os_id@data$period <- "transit", os_id@data$period <- "natal")
  
  os_id <- os_id[,-which(colnames(os_id@data) %in% "area")] }
  
  # 2. Identify locations belonging to ESTABLISHED mcp (last or 2 last years)
  
  mcp_est2_id <- which(names_mcp_est2 %in% id[i]) # MCP natal
  mcp_est1_id <- which(names_mcp_est1 %in% id[i]) # MCP natal
  
  if(length(mcp_est2_id) > 0 ) {
    
    z <- mcp_est2[[which(names_mcp_est2 %in% id[i])]]
    os_id <- cbind(os_id,over(os_id,z[[2]])) # Select points that fall within established mcp
    colnames(os_id@data)[which(colnames(os_id@data) %in% "id")] <- "est_loc"
    
    for (j in 1:nrow(os_id)){
      if(os_id@data$period[j] == "transit" & !is.na(os_id@data$est_loc[j])) {
        os_id@data$period[j] <- "est2"
      } }

    os_id <- os_id[,-which(colnames(os_id@data) %in% "area")]
    
  } else if (length(mcp_est2_id) < 0 & length(mcp_est1_id) > 0) {
    z <- mcp_est1[[which(names_mcp_est1 %in% id[i])]]
    os_id <- cbind(os_id,over(os_id,z[[2]])) # Select points that fall within established mcp
    colnames(os_id@data)[which(colnames(os_id@data) %in% "id")] <- "est_loc"
    
    for (k in 1:nrow(os_id)){
      if(os_id@data$period[k] == "transit" & !is.na(os_id@data$est_loc[k])) {
        os_id@data$period[k] <- "est1"
      } }
    
    os_id <- os_id[,-which(colnames(os_id@data) %in% "area")]
  }
  #writeOGR(os_id, "D:/MargSalas/Oso/Datos/GIS/2022/Identified_individuals_periods", paste("ID_", id[i], sep = ""), driver = "ESRI Shapefile")
  id_loc[[i]] <- os_id

}

## ---- Exploration of number of locations per period and year ----

summary_id <- list()
info <- list()

for (i in 1:length(id_loc)){
  
  dat_id <- id_loc[[i]]@data
  
  info[[1]] <- unique(id_loc[[i]]@data$Confirmed_Individual_cub)
  
  info[[2]] <- dat_id %>%        # Summary of number of obs per period and year, and if it reproduced
    group_by(period,Year, With_cubs_estimated) %>%
    summarise(n())
  
  info[[3]] <- dat_id %>%       # Summary of n obs per type of obs in transit period
    filter(period == "transit") %>%
    group_by(Year, Obs_type) %>%
    summarise(n())
  
  summary_id[[i]] <- info
  
}

summary_id[which(sapply(summary_id, `[[`, 1) == "Cachou")]

unique(dat_id$N_cubs_estimated)



  