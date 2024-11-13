


## -------------------------------------------------------
## Explore female dispersal patterns to argument discussion
## --------------------------------------------------------

rm(list = ls())

## -------------------------------------------------
##             1. Looking at AC locations    
## ------------------------------------------------- 

library(nimbleSCR)
library(sf)
library(rgdal)

setwd("D:/MargSalas/Oso/Datos/GIS/Countries")
Xbuf <- readOGR("Buffer_statespace.shp")
Xbuf2 <- readOGR("Buffer_8500_traps.shp")

nuc <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/nuc.shp")
per <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/per.shp")
sa <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/clip_pyros2_WGS84_31N_all.shp")


setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721")
load("zObserved_yAgeDeaths.RData")

zdatAGE[is.na(zdatAGE)] <- 0

# Load sxy estimated by the model and unscale

Tt <- 5
M.aug <- 300

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams_FINAL")
load("myResults_3-3.1_sxy.RData")

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Data/Systematic_FINAL_1721") # Load original habitat coordinates to unscale
load("habcoord.RData")

sampmat2 <- do.call(rbind, nimOutputSXY) 
s.which <- grep('sxy', colnames(sampmat2)) # ASP: index columns all sxy (sampmat matrix)
sampmat2_sxy <- sampmat2[, s.which]
dim(sampmat2_sxy)
mean_sampmat <- colMeans(sampmat2_sxy) # I will plot the mean location over iterations
names(mean_sampmat)

sxy <- array(NA, c(300, 2, 5))
for(t in 1:Tt){
  s.which.year <- grep(paste(t,"]", sep = ""), colnames(sampmat2_sxy)) # ASP: index columns all sxy (sampmat matrix)
  sxy[,,t] <- matrix(mean_sampmat[s.which.year] , M.aug, 2) 
}

dimnames(sxy)[[2]] <- c('x','y') 
sxy.uns <- scaleCoordsToHabitatGrid(coordsData = sxy,## this are your sxy
                                    coordsHabitatGridCenter = G,# this is your unscaled habitat (as you used when scaling the habitat/detector to the habitat. G?
                                    scaleToGrid = FALSE)$coordsDataScaled

# Check observed individuals that model places outside
yearnames <- c("2017", "2018", "2019", "2020", "2021")

sa_cropped <- st_crop(sa, st_bbox(per))
plot(st_geometry(sa_cropped))

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Results_section/Plots")
pdf("SI_AC_females.pdf",7,2)
par(mfrow = c(2,3),
    mar = c(0,0,1.5,0),
    oma = c(1,1,1,1))

for (t in 1:5){
  plot(st_geometry(sa_cropped), col = "beige", border = "white")
  plot(st_geometry(per), col = "#a6dba0", border = "#a6dba0", add = TRUE)
  plot(st_geometry(nuc), col = "#9970ab", border = "#9970ab", add = TRUE)
  mtext(yearnames[t], side = 3, line = 0)
  #plot(Xbuf, main = yearnames[t])
  #plot(Xbuf2, add = TRUE)
  p <- sxy.uns[which(zdatAGE[,t] == 1 & sex == 0),,t]
  sp <- SpatialPoints(p, proj4string=CRS(proj4string(Xbuf2)))
  points(sp, pch= 19, lwd = 0.7, col = "black")
  
  over(sp,Xbuf2)
}
dev.off()

# Identify which ones are outside
setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Spread")
load("ages.RData") # To know the actual names and compare with dispersal maps below
rownames(zdatAGE) <- c(rownames(age.cat),seq(62:300))

p.outside <- list()
o <- list()
w <- list()

for (t in 1:5){
  
  which.alive <- which(zdatAGE[,t] == 1 & sex == 0)
  p <- data.frame(sxy.uns[which.alive,,t])
  sp <- st_as_sf(p, coords = c(1,2), crs = "+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") 
  #sp <- SpatialPoints(p, proj4string=CRS(proj4string(Xbuf2)))
  inter <- lengths(st_intersects(sp, nuc)) > 0 # FALSE is no intersection point not present in any polygon)
  which.out <- which(inter == FALSE) # Index to identify
  w$age <- ageMatAug[which.alive[which.out],t]
  w$sex <- sex[which.alive[which.out]]
  w$ind <- names(which.alive)[which.out]
  o[[t]] <- w # Store info
  p.outside[[t]] <- sp[which.out,] # Store points
  
} 

## -------------------------------------------------
## 2. Looking at all observations (natal-established)   
## ------------------------------------------------- 

library(ggplot2)
#library(ggmap)
library(gridExtra)
library(rgdal)
#library(ggpubr)
library(dplyr)
#library(geosphere)
library("viridis")
#library(RColorBrewer)
library(sf)
#library(mapview)
#
#devtools::install_github('oswaldosantos/ggsn')


# Load monitoring data
setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2021")
os <- read.csv("Natal_established_coordinates_cubLocations.csv", header = TRUE, row.names = NULL, sep = ";")
os <- os[,-1]

# Set coordinates
coordinates(os) <- os[,c("x_long","y_lat")] # Spatial object
os@proj4string <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
os <- spTransform(os, CRS("+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))

# ID with more than 10 0bs
setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2021")
load("id_more10.RData")

# Load mcp
setwd("D:/MargSalas/Oso/Datos/GIS/2021")
mcp_natal <- readRDS("mcp_natal.rds")
mcp_est2 <- readRDS("mcp_est2.rds")

names_mcp_nat <- lapply(mcp_natal, `[[`, 1) # Take only the first argument of the list to know the name
names_mcp_est <- lapply(mcp_est2, `[[`, 1)


# Load info database
setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2021")
info <- readxl::read_xlsx("D:/MargSalas/Oso/Datos/Tablas_finales/2021/Info_individuals_2021.xlsx", sheet = 1)
info <- info[,c(4:8,10)]
colnames(info)[6] <- "Year_death"

# Load dispersal distances
setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2021")
d <- read.csv("disp_distance.csv", header = TRUE, row.names = NULL, sep = ",")
d <- d[,-1]
# Funtions
#'%!in%' <- function(x,y)!('%in%'(x,y))

# Load nucleus and periphery
sa <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/clip_pyros2_WGS84_31N_all.shp")
nuc <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/nuc.shp")
per <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/per.shp")

## -------------------------------------------------
##                    FEMALES
## ------------------------------------------------- 

os_females <- os[which(os$Sex_cub == "F"), ]
id_females <- unique(os_females$Confirmed_Individual_cub)
id_females <- id_females[id_females %in% id_more10] # ID females with more than 10 observations (all)
id_females <- id_females[id_females %in% d$ID] # Females for which we could estimate a dispersal distance

# Unify age classes to plot
unique(os_females$Age_class_cub)
os_females$Age_class_cub[which(os_females$Age_class_cub %in% c("Cub0","Cub1", "Cub2"))] <- "Cub"

# This is to put it anonymous
info$Year_birth <- as.integer(info$Year_birth)
info$Year_death[which(info$Year_death == "Alive")] <- 2021
info$Year_death <- as.integer(info$Year_death)

#pdf("ID_females2.pdf")  
f <- list()
natal <- est <- mean_loc <- list()

#for (xxx in 1:4) { # WIth this f*shit it doesn't work, but with apply it does so....
f <-  lapply(1:length(id_females), function(xxx){ 
  os_id <- os_females[which(os_females$Confirmed_Individual_cub == id_females[xxx]), ] # Select individual
  df <- os_id@data
  df$x_coord <- coordinates(os_id)[,1]
  df$y_coord <- coordinates(os_id)[,2]
  
  # ---- Plot average location by year ----
  
  # 1. Average x and y
  mean_loc[[xxx]] <- df %>%
    group_by(Year, Age_class_cub ) %>%
    summarise(
      count = n(),
      meanx = mean(x_coord, na.rm = TRUE),
      meany = mean(y_coord, na.rm = TRUE))
  na <- which(is.na(mean_loc[[xxx]]$Year)) 
  if(length(na) > 0) { # If there is NA in the data, remove them
    mean_loc[[xxx]][-which(is.na(mean_loc$Year)), ]
  }
  
  # 2. Plot

    # x[[2]] <- spTransform(x[[2]], CRS("+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")) # mcp
    
    # natal <- SpatialPoints(matrix(c(d[d$ID %in% id_females[xxx],c(2)], d[d$ID %in% id_females[xxx],c(3)]), nrow = 1, ncol = 2), 
    #                        proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
    # natal <- spTransform(natal, CRS("+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
    # natal <- as.data.frame(coordinates(natal))
    
    natal[[xxx]] <- SpatialPoints(matrix(c(d[d$ID %in% id_females[xxx],c(2)], d[d$ID %in% id_females[xxx],c(3)]), nrow = 1, ncol = 2), 
                           proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
    natal[[xxx]] <- spTransform(natal[[xxx]], CRS("+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
    natal[[xxx]] <- as.data.frame(coordinates(natal[[xxx]]))
    
    if(!is.na(d[d$ID %in% id_females[xxx],c(4)])) { # If there are more exact established coordinates, take est1
    
    est[[xxx]] <- SpatialPoints(matrix(c(d[d$ID %in% id_females[xxx],c(4)], d[d$ID %in% id_females[xxx],c(5)]), nrow = 1, ncol = 2), 
                           proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
    est[[xxx]] <- spTransform(est[[xxx]], CRS("+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
    est[[xxx]] <- as.data.frame(coordinates(est[[xxx]]))
    
    } else {  
      
      est[[xxx]] <- SpatialPoints(matrix(c(d[d$ID %in% id_females[xxx],c(6)], d[d$ID %in% id_females[xxx],c(7)]), nrow = 1, ncol = 2),  # Otherwise take est2
                                          proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
      est[[xxx]] <- spTransform(est[[xxx]], CRS("+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
      est[[xxx]] <- as.data.frame(coordinates(est[[xxx]]))
      
    }
    
    # Calculate years of study for individual (for anonymous ids)
    
    study_years <- info$Year_death[info$ID %in% id_females[xxx]] - info$Year_birth[info$ID %in% id_females[xxx]]
    
    
    #f[[xxx]] <-
    p <- ggplot() + 
      geom_sf(data = sa, colour = "white", fill = "beige") +
      geom_sf(data = per, colour = "#a6dba0", fill = "#a6dba0") +
      geom_sf(data = nuc, colour = "#9970ab", fill = "#9970ab") +
      coord_sf(xlim = c(st_bbox(per)[1], st_bbox(per)[3]),  ylim = c(st_bbox(per)[2] , st_bbox(per)[4]), expand = FALSE) +
      geom_point(data = mean_loc[[xxx]], aes(x = meanx , y = meany, colour = "red"), 
               size = 1.5) +
      geom_segment(aes(x = natal[[xxx]][,1], y = natal[[xxx]][,2], 
                       xend = est[[xxx]][,1], yend = est[[xxx]][,2]),
                   arrow = arrow(length = unit(0.3, "cm")),
                   lwd = 0.7) +
      theme(panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "none" ) + 
      # ggtitle(paste(id_females[xxx]," (",info$Year_birth[info$ID %in% id_females[xxx]], "-", info$Year_death[info$ID %in% id_females[xxx]], ")", sep = ""))
      ggtitle(paste("ID",xxx," (",study_years, " years)", sep = ""))
      
    
    p
})



  ##################################################################################################################################
    

  
#setwd("D:/MargSalas/Oso/Datos/Plots/2021/Detailed_ID2/females")
setwd("D:/MargSalas/Oso/OPSCR_project/Results/Results_section/Plots")

  pdf("SI_disp_females.pdf",
      width = 7, height = 7)
  
  grid.arrange(grobs = f, ncol = 3)

  dev.off()

