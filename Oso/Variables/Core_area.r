
## -------------------------------------------------
##                 CORE AREA
## ------------------------------------------------- 

# Identify the core zone of the population to create the variable "distance to
# core area" to explain spatial patterns in brown bear density

rm(list = ls())

library(tidyverse)
library(sf)
library(rgdal)
library(mapview)
library(lubridate)
library(adehabitatHR)
library(sp)
library(raster)
library(ggplot2)

## ---- Visual exploration of reintroduction patterns ----

setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2021")
os <- read.csv("Data_os_96_21_cubLocations.csv", header = TRUE, row.names = NULL)  %>% 
  filter(Confirmed_Individual != "Indetermined") %>%
  st_as_sf(coords = c("x_long","y_lat"), crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

os <- os[,-which(colnames(os) %in% c("ID_obs"))] # So that the script runs, because this column was added later and it doesn't fit otherwise

# Mellba y Ziva reintroducidas en 1996

R_fem96_97 <- os %>%
  filter(Year %in% c(1996,1997) & Confirmed_Individual %in% c("Mellba", "Ziva"))
data.frame(R_fem96_97 %>% group_by(Confirmed_Individual, Obs_type) %>% summarise(n()))


# Hembras del 1996-2005 (periodo hasta siguiente reintroducción en 2006):
  # Ziva + Cannelle + Caramelles -> Caramellita (matrilineage Mellba)

fem98_05 <- os %>%
  filter(Year %in% c(1998:2005) & Sex == "F")
data.frame(fem98_05 %>% group_by(Confirmed_Individual, Obs_type) %>% summarise(n()))


# Hembras reintroducidas en 2006: Palouma, Franscka, Hvala, Sarousse

R_fem06 <- os %>%
  filter(Year %in% c(2006) & Confirmed_Individual %in% c("Palouma", "Franska", "Hvala", "Sarousse"))
data.frame(R_fem06 %>% group_by(Confirmed_Individual, Obs_type) %>% summarise(n()))

R_fem06_dead <- os %>%
  filter(Year %in% c(2006) & Confirmed_Individual %in% c("Palouma", "Franska"))

# Hembras del 2006 - 2017 (periodo hasta siguiente reintroducción)
   # Seguimiento sistemático empieza en 2010, nuestro estudio en 2017

fem06_17 <- os %>%
  filter(Year %in% c(2006:2017) & Sex == "F")
data.frame(fem06_17 %>% group_by(Confirmed_Individual, Obs_type) %>% summarise(n()))

# Hembras reintroducidas en 2018: Sorita y Claverina

R_fem18 <- os %>%
  filter(Year %in% c(2018) & Confirmed_Individual %in% c("Sorita", "Claverina"))
data.frame(R_fem18 %>% group_by(Confirmed_Individual, Obs_type) %>% summarise(n()))

mapview(os, cex = 4) + 
  mapview(R_fem96_97, col.regions = "red", cex = 4) +  # -> Se detecta el núcleo principal
  mapview(fem98_05, col.regions = "green", cex = 4) +  # -> Posiciones en el mismo núcleo. Se detecta el núcleo oriental con Cannelle
  mapview(R_fem06, col.regions = "pink", cex = 4)  +   # -> Se amplia el núcleo principal (sobre todo por Palouma). Nuevo nucleo intermedio por Franska (que muere sin descendencia en 2007)
  mapview(fem06_17, col.regions = "orange") + # -> Se amplia el núcleo principal: Posiciones extra de Sarousse y Franska
  mapview(R_fem18, col.regions = "yellow", cex = 4)    # -> Nucleo oriental

# Save plots

map <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/clip_pyros2.shp")

setwd("D:/MargSalas/Oso/Datos/Plots/2022/Dist_core")
pdf("R_fem96_97.pdf")
plot(st_geometry(map), col = "grey")
mtext("Introduced Mellba and Ziva (1996)", side = 3, line = -1.5, cex = 2)
plot(os, pch = 19, col = "darkblue", add = TRUE)
plot(R_fem96_97, pch = 19, col = "red", add = TRUE)
dev.off()

setwd("D:/MargSalas/Oso/Datos/Plots/2022/Dist_core")
pdf("fem98_05.pdf")
plot(st_geometry(map), col = "grey")
mtext("Females 98-05", side = 3, line = -1.5, cex = 2)
plot(os, pch = 19, col = "darkblue", add = TRUE)
plot(fem98_05, pch = 19, col = "green", add = TRUE)
dev.off()

setwd("D:/MargSalas/Oso/Datos/Plots/2022/Dist_core")
pdf("R_fem06.pdf")
plot(st_geometry(map), col = "grey")
mtext("Females introduced 2006", side = 3, line = -1.5, cex = 2)
mtext("Palouma, Franscka, Hvala, Sarousse", side = 3, line = -3, cex = 1)
plot(os, pch = 19, col = "darkblue", add = TRUE)
plot(R_fem06, pch = 19, col = "magenta", add = TRUE)
plot(R_fem06_dead, pch = 19, col = "deeppink4", add = TRUE)
dev.off()

setwd("D:/MargSalas/Oso/Datos/Plots/2022/Dist_core")
pdf("fem06_17.pdf")
plot(st_geometry(map), col = "grey")
mtext("Females 06-17", side = 3, line = -1.5, cex = 2)
plot(os, pch = 19, col = "darkblue", add = TRUE)
plot(fem06_17, pch = 19, col = "orange", add = TRUE)
dev.off()

setwd("D:/MargSalas/Oso/Datos/Plots/2022/Dist_core")
pdf("R_fem18.pdf")
plot(st_geometry(map), col = "grey")
mtext("Introduced Sorita and Claverina (2018)", side = 3, line = -1.5, cex = 2)
plot(os, pch = 19, col = "darkblue", add = TRUE)
plot(R_fem18, pch = 19, col = "yellow", add = TRUE)
dev.off()


# Core area: Posiciones de hembras reproductoras (eliminar Franska, Palouma, Sarousse)
# -> Tiene poco sentido incluir las localizaciones de hembras reintroducidas no reproductoras? 
#    Sobre todo Franska y Palouma que murieron el mismo año
# -> Quitando estas, todas son hijas de Hvala y Caramelles (hija de Mellba)
# -> Posible delimitación: Distribución de hembras reintroducidas reproductivas,el año de su primera reproducción (asumiendo que están asentadas)
#    * Núcleo central: Distribición Ziva, Mellba y Hvala (tiene sentido porque antes de 2010 no hay seguimiento sistemático)
#    * Núcleo oriental: Distribución Sorita (Claverina no se ha reproducido)

# Pregunta: Reproducciones detectaadas endemicas??? Cannelle y Cannellito solo?

# Subset locations to delineate potential core area

core <- os %>%
  filter(Year %in% c(1997,1998, 2007,2008, 2019, 2020) 
         & Confirmed_Individual %in% c("Mellba", "Ziva", "Hvala", "Sorita") 
         & With_cubs_estimated %in% c("<6month", "1", "2"))


mapview(os) + 
  mapview(core, col.regions = "red") +
  mapview(fem98_05, col.regions = "green")

setwd("D:/MargSalas/Oso/Datos/Plots/2022/Dist_core")
pdf("core.pdf")
plot(st_geometry(map), col = "grey")
mtext("Core Area (locations)", side = 3, line = -1.5, cex = 2)
plot(os, pch = 19, col = "darkblue", add = TRUE)
plot(core, pch = 19, col = "red", add = TRUE)
dev.off()

## ---- Calculate kernel utilization distribution and MCP ----

setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2021")
os_spatial <- read.csv("Data_os_96_21_cubLocations.csv", header = TRUE, row.names = NULL)  %>% 
  filter(Confirmed_Individual != "Indetermined") 

os_spatial <- os_spatial[,-which(colnames(os_spatial) %in% c("ID_obs"))]

core_spatial <- os_spatial %>%
  filter(Year %in% c(1997,1998, 2007,2008, 2019, 2020) 
         & Confirmed_Individual %in% c("Mellba", "Ziva", "Hvala", "Sorita") 
         & With_cubs_estimated %in% c("<6month", "1", "2"))

coordinates(core_spatial) <- core_spatial[,c("x_long", "y_lat")]
proj4string(core_spatial) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# 1. Calculate kernel for each individual and join

core_kernel <- core_spatial[,2]

kref <- kernelUD(core_kernel, h="href")
kernel.poly <- getverticeshr(kref, percent = 95) 

core_area_kernelpol <- aggregate(kernel.poly, dissolve = TRUE)

mapview(core, col.regions = "red") + 
  mapview(core_area_kernelpol, col.regions = "pink") 
  + mapview(fem98_05, col.regions = "green") #) With this we see that cannelle overlaps also with occidental

setwd("D:/MargSalas/Oso/Datos/Plots/2022/Dist_core")
pdf("core_fullkernel95.pdf")
plot(st_geometry(map), col = "grey")
mtext("Full kernel 95%", side = 3, line = -1.5, cex = 2)
plot(os, pch = 19, col = "darkblue", add = TRUE)
plot(core, pch = 19, col = "red", add = TRUE)
plot(core_area_kernelpol, pch = 19, col = adjustcolor("yellow", alpha.f = 0.6), add = TRUE)
dev.off()

# Estimate saturation curve per individual
estimateHomeRange <- function(locations, n_locations, percent) {
  # Select subset of locations
  subset_locations <- locations[1:n_locations, ]
  # Perform KDE home range estimation
  hr <- kernelUD(subset_locations, h = "href")
  # Calculate home range size
  hr_size <- kernel.area(hr, percent = percent) 
  
  # Return home range size
  return(hr_size)
}
ind <- unique(core_kernel$Confirmed_Individual)

par(mfrow = c(2,2))
for(i in 1:length(ind)){
  core_kernel_ind <- core_kernel[which(core_kernel$Confirmed_Individual %in% ind[i]),]
  saturation_data <-  data.frame(
    Num_Locations = 50:nrow(core_kernel_ind),
    Home_Range_Size = NA)
  for(l in 1:nrow(saturation_data)){
    saturation_data[l,2] <- estimateHomeRange(locations = core_kernel_ind, n_locations = saturation_data[l,1], percent = 95)
  }
  plot(saturation_data$Home_Range_Size ~ saturation_data$Num_Locations,
       xlab = "Number of Locations", ylab = "Home Range Size", main = ind[i])
}



# 2. Calculate kernel for each population nucleus
  
core_spatial$Nucleus <- ifelse(core_spatial$Confirmed_Individual %in% c("Mellba", "Ziva", "Hvala"), 1, 2)

###### Count the number of locations per nucleus for methods
nrow(core_spatial[which(core_spatial$Nucleus == 1),])
nrow(core_spatial[which(core_spatial$Nucleus == 2),])
#######################################################

core_kernel2 <- core_spatial[,24]

kref2 <- kernelUD(core_kernel2, h="href")
kernel.poly2 <- getverticeshr(kref2, percent = 95) 

core_area_kernelpol2 <- aggregate(kernel.poly2, dissolve = TRUE)

setwd("D:/MargSalas/Oso/Datos/Plots/2022/Dist_core")
pdf("core_kernel95_n1n2.pdf")
plot(st_geometry(map), col = "grey")
mtext("Kernel 95% per nucleus", side = 3, line = -1.5, cex = 2)
plot(os, pch = 19, col = "darkblue", add = TRUE)
plot(core, pch = 19, col = "red", add = TRUE)
plot(core_area_kernelpol2, pch = 19, col = adjustcolor("yellow", alpha.f = 0.6), add = TRUE)
dev.off()

# Check saturation curve for home range for both nucleus
# Estimate saturation curve per individual

ind <- unique(core_kernel2$Nucleus)
nuc <- c("central", "occidental")
par(mfrow = c(1,2))
for(i in 1:length(ind)){
  core_kernel_ind <- core_kernel2[which(core_kernel2$Nucleus %in% ind[i]),]
  saturation_data <-  data.frame(
    Num_Locations = 50:nrow(core_kernel_ind),
    Home_Range_Size = NA)
  for(l in 1:nrow(saturation_data)){
    saturation_data[l,2] <- estimateHomeRange(locations = core_kernel_ind, n_locations = saturation_data[l,1], percent = 95)
  }
  plot(saturation_data$Home_Range_Size ~ saturation_data$Num_Locations,
       xlab = "Number of Locations", ylab = "Home Range Size", main = nuc[i])
}


# Try also to limit the % to 90 in the occidental core, because it is only one female
# and is over-represented?

core_kernel_n1 <- core_spatial[core_spatial$Nucleus == 1,24]
core_kernel_n2 <- core_spatial[core_spatial$Nucleus == 2,24]

kref <- kernelUD(core_kernel_n1, h="href")
kernel.poly <- getverticeshr(kref, percent = 95) 
core_area_kernelpol_n1 <- aggregate(kernel.poly, dissolve = TRUE)

kref <- kernelUD(core_kernel_n2, h="href")
kernel.poly <- getverticeshr(kref, percent = 80) 
core_area_kernelpol_n2 <- aggregate(kernel.poly, dissolve = TRUE)

setwd("D:/MargSalas/Oso/Datos/Plots/2022/Dist_core")
pdf("core_kernel_95n1_80n2_SI.pdf")
plot(st_geometry(map), col = "grey")
mtext("Kernel 95% (central), 80% (occidental)", side = 3, line = -1.5, cex = 2)
#plot(os, pch = 19, col = "darkblue", add = TRUE)
plot(core, pch = 19, col = "red", add = TRUE)
plot(core_area_kernelpol_n1, pch = 19, col = adjustcolor("yellow", alpha.f = 0.6), add = TRUE)
plot(core_area_kernelpol_n2, pch = 19, col = adjustcolor("yellow", alpha.f = 0.6), add = TRUE)
dev.off()

########## Plot for methods

sa <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/clip_pyros2_WGS84_31N_all.shp")
eur <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/esp_fr_2.shp") %>%
  st_transform(dpts, crs = crs(sa))
esp <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/ESP_adm/ESP_adm0.shp") %>%
  st_transform(dpts, crs = crs(sa))
fr <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/FRA_adm/FRA_adm0.shp") %>%
  st_transform(dpts, crs = crs(sa))

core <- core %>% st_transform(32631)
core_area_kernelpol_n1 <- spTransform(core_area_kernelpol_n1,CRS("+proj=utm +zone=31 +datum=WGS84 +units=m +no_defs +type=crs"))
core_area_kernelpol_n2 <- spTransform(core_area_kernelpol_n2,CRS("+proj=utm +zone=31 +datum=WGS84 +units=m +no_defs +type=crs"))

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Results_section/Plots")
pdf("SI_core_kernel_95n1_80n2.pdf",7,7)

plot(st_geometry(sa), col = "beige", border = "beige",  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000)) # To add andorra with same color (not present in eur)
plot(st_geometry(eur), col = "beige", border = "beige",  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000), add = TRUE)
plot(st_geometry(esp), col = "beige", border = "white",  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000), add = TRUE)
plot(st_geometry(fr), col = "beige", border = "white",  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000), add = TRUE)
plot(core, pch = 19, col = "brown", add = TRUE)
plot(core_area_kernelpol_n1, pch = 19, col = adjustcolor("yellow", alpha.f = 0.7), add = TRUE)
plot(core_area_kernelpol_n2, pch = 19, col = adjustcolor("yellow", alpha.f = 0.6), add = TRUE)

dev.off()


mapview(core, col.regions = "red") + 
  mapview(core_area_kernelpol, col.regions = "pink") +
  mapview(core_area_kernelpol2, col.regions = "green") +
  mapview(core_area_kernelpol_n1, col.regions = "violet") +
  mapview(core_area_kernelpol_n2, col.regions = "blue") +
  mapview(fem98_05, col.regions = "green") #) With this we see that cannelle overlaps also with occidental

core_kernel_n1n2 <- bind(core_area_kernelpol_n1, core_area_kernelpol_n2)
core_kernel_n1n2 <- aggregate(core_kernel_n1n2, dissolve = TRUE)
plot(core_kernel_n1n2) ## BETTER TO CALCULATE A KERNEL PER POPULATION NUCLEUS


# 3. Sum of individuals MCP total MCP
core_mcp <- core_spatial[,2]

mcp.poly <- mcp(core_mcp, percent = 95)
core_area_mcppol <- aggregate(mcp.poly, dissolve = TRUE)

mapview(core, col.regions = "red") +
  mapview(core_area_kernelpol2, col.regions = "green") +
  mapview(core_area_mcppol, col.regions = "yellow") ## KERNEL LOOKS BETTER AND IT TAKES INTO ACCOUNT THE DISTRIBUTION OF THE DATA  
 

# 4. Kernel per population nucleus + sample GPS locations to avoid over-representation

# After comment of Pierre-Yves Quennette.
# Especially for Sorita, in order to be able to estimate Kernel 95% for BOTH nucleus (more fancy)

f <- core_spatial@data[which(core_spatial$Obs_type == "GPS"),]

# I need to take their locations directly from the Radiotracking data, to get the fixes
setwd("D:/MargSalas/Oso/Datos/GPS")

os_spatial_rt <- read.csv("Radiotracking_ossos_1996_2020_taula_final.csv", header = TRUE, row.names = NULL)  

core_spatial_rt <- os_spatial_rt %>%
  filter(Year %in% c(1997,1998, 2007,2008, 2019, 2020) 
         & Bear_name %in% c("Mellba", "Ziva", "Hvala", "Sorita") 
         & With_cubs_estimated %in% c("<6month", "1", "2"))

hv <- core_spatial_rt[core_spatial_rt$Bear_name == "Hvala",] # 1 year (2007). 8 locations per day (or even more), GPS
so <- core_spatial_rt[core_spatial_rt$Bear_name == "Sorita",] # 2 years (2019 - 2020). 8 locations per day (or even more), GPS
mb <- core_spatial@data[core_spatial@data$Confirmed_Individual == "Mellba",] # 1 year 1997. 1 location per day or less, only VHF (NOT in radiotracking df)
zv <- core_spatial_rt[core_spatial_rt$Bear_name == "Ziva",] # 2 years (1997 - 1998). 1 location per DAY, only VHF

nrow(hv) + nrow(so) + nrow(mb) + nrow(zv) # Similar to df with all locations core_spatial (core_spatial has more, which makes sense)
# For Hvala and Sorita, sample dataset to get one location per day. I can do it in the original dataset, as I have the date there

#core_spatial@data$Date <- as.Date(core_spatial@data$Date, format = "%d/%m/%Y")
core_spatial_thinned <- core_spatial@data %>% group_by(Confirmed_Individual,Date_register_formated) %>% sample_n(1) %>% arrange(Confirmed_Individual,Year,Date_register_formated)

coordinates(core_spatial_thinned) <- core_spatial_thinned[,c("x_long", "y_lat")]
proj4string(core_spatial_thinned) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

core_spatial_thinned$Nucleus <- ifelse(core_spatial_thinned$Confirmed_Individual %in% c("Mellba", "Ziva", "Hvala"), 1, 2)
core_kernel3 <- core_spatial_thinned[,24]
kref3 <- kernelUD(core_kernel3, h="href")
kernel.poly3 <- getverticeshr(kref3, percent = 95)
core_area_kernelpol3 <- aggregate(kernel.poly3, dissolve = TRUE)

setwd("D:/MargSalas/Oso/Datos/Plots/2022/Dist_core")
pdf("core_thinned_kernel95_n1n2.pdf")
plot(st_geometry(map), col = "grey")
mtext("Kernel 95% per nucleus thinned", side = 3, line = -1.5, cex = 2)
plot(os, pch = 19, col = "darkblue", add = TRUE)
plot(core, pch = 19, col = "red", add = TRUE)
plot(core_area_kernelpol3, pch = 19, col = adjustcolor("yellow", alpha.f = 0.6), add = TRUE)

dev.off()


kernel.poly4 <- getverticeshr(kref3, percent = 80)
core_area_kernelpol4 <- aggregate(kernel.poly4, dissolve = TRUE)

setwd("D:/MargSalas/Oso/Datos/Plots/2022/Dist_core")
pdf("core_thinned_kernel80_n1n2.pdf")
plot(st_geometry(map), col = "grey")
mtext("Kernel 80% per nucleus thinned", side = 3, line = -1.5, cex = 2)
plot(os, pch = 19, col = "darkblue", add = TRUE)
plot(core, pch = 19, col = "red", add = TRUE)
plot(core_area_kernelpol4, pch = 19, col = adjustcolor("yellow", alpha.f = 0.6), add = TRUE)

dev.off()



## ---- Raster distance to core area ----
## ----  1. 95% central, 80% western ----
## FROM KERNEL
# Create empty raster with extent study area rasterize polygon core area
setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/EU_DEM")
dem <- raster("dem_WGS84_31N_Clip.tif")
rs <- dem 
rs.aggregate <- aggregate(rs, fact=8) # 200m resolution, similar to forest
res(rs.aggregate)
values(rs) <- NA # Raster to calculate 

# Rasterize polygon core area
core_area_kernelpol <- spTransform(core_kernel_n1n2, crs(rs) )
core_area_kernelpol <- as(core_area_kernelpol, "SpatialPolygonsDataFrame")

setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/DistCore")
writeOGR(core_area_kernelpol, "corepol", dsn = "D:/MargSalas/Oso/Datos/GIS/Variables/Europe/DistCore/corepol", driver = "ESRI Shapefile")

r <- rasterize(core_area_kernelpol,rs)
plot(r)
# To save as raster with 1 = core, 0 = not core
#values(r)[is.na(values(r))] <- 0
#setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/DistCore")
#writeRaster(r, filename = 'core', format = 'GTiff')

# Calculate distance
dist_corekernel <- distance(r)

setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/DistCore")
writeRaster(dist_corekernel, filename = 'dist_corekernel', format = 'GTiff')
plot(dist_corekernel)
plot(st_geometry(map), add = TRUE)

# Calculate the logaritmic distance (To give less importance to far locations)
setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/DistCore")
distcore <- raster("dist_corekernel.tif")
plot(distcore)

distcore@data@values <- distcore[] # The values were not stored, don't know why

logDistcore <- distcore # Create layer to store the log
logDistcore@data@values <- log(logDistcore@data@values) # Transform and store (2 ways)
logDistcore[] <- log(distcore[])

setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/DistCore")
writeRaster(logDistcore, filename = 'dist_corekernel_log', format = 'GTiff')

plot(logDistcore)

## ----  2. 80% central, 80% western, thinned ----

# Rasterize polygon core area
core_area_kernelpol_thinned <- spTransform(core_area_kernelpol4, crs(rs) )
core_area_kernelpol_thinned <- as(core_area_kernelpol_thinned, "SpatialPolygonsDataFrame")

r <- rasterize(core_area_kernelpol_thinned,rs)
plot(r)
#values(r)[is.na(values(r))] <- 0
dist_corekernel_thinned <- distance(r)
setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/DistCore")
writeRaster(dist_corekernel_thinned, filename = 'dist_corekernel_thinned', format = 'GTiff')

## ---- -> Check correlation between 1 and 2 ----

setwd("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/DistCore")
distcore <- raster("dist_corekernel.tif")
dist_corekernel_thinned <- raster("dist_corekernel_thinned.tif")
d <- stack(distcore, dist_corekernel_thinned)
cor1 <- cor(sampleRandom(d, size= ncell(distcore), method = "spearman") )

