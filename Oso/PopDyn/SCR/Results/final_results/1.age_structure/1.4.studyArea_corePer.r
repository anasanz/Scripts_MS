
## ---- Study area and location of core/periphery ----

rm(list = ls())

library(sf)
library(terra)
library(raster)

# Load europe
sa <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/clip_pyros2_WGS84_31N_all.shp")
eur <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/esp_fr_2.shp") %>%
  st_transform(dpts, crs = crs(sa))


esp <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/ESP_adm/ESP_adm0.shp") %>%
  st_transform(dpts, crs = crs(sa))
fr <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/FRA_adm/FRA_adm0.shp") %>%
  st_transform(dpts, crs = crs(sa))


nuc <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/nuc.shp")
per <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/per.shp")

# Plot
setwd("D:/MargSalas/Oso/OPSCR_project/Results/Results_section/Plots")
pdf("1.4.sa_coreper.pdf",7,7)

plot(st_geometry(sa), col = "beige", border = "beige",  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000)) # To add andorra with same color (not present in eur)
plot(st_geometry(eur), col = "beige", border = "beige",  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000), add = TRUE)
plot(st_geometry(esp), col = "beige", border = "white",  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000), add = TRUE)
plot(st_geometry(fr), col = "beige", border = "white",  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000), add = TRUE)
plot(st_geometry(nuc), col = "#9970ab", border = "#9970ab", add = TRUE)
plot(st_geometry(per), col = "#a6dba0", border = "#a6dba0", add = TRUE)
plot(st_geometry(esp), border = adjustcolor("white", alpha.f = 0.3),  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000), add = TRUE)
plot(st_geometry(fr), border = adjustcolor("white", alpha.f = 0.3),  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000), add = TRUE)
segments(x0 = 400000, x1 = 400000 + 50000,
         y0 = 4650000, y1 = 4650000, col = grey (0.3), lwd = 2)
text(400000 + 50000/2, 4665000, labels = "50 km")

dev.off()

st_bbox(sa)

# To plot with all the department borders:
#plot(st_geometry(sa), col = "beige", border = "white",  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000)) # To add andorra with same color (not present in eur)
#plot(st_geometry(eur), col = "beige", border = "white",  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000), add = TRUE)
#plot(st_geometry(nuc), col = "#9970ab", border = "#9970ab", add = TRUE)
#plot(st_geometry(per), col = "#a6dba0", border = "#a6dba0", add = TRUE)
#plot(st_geometry(eur), border = adjustcolor("white", alpha.f = 0.3),  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000), add = TRUE)


## ---- Size study area for methods ----

Xbuf2 <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/Buffer_8500_traps_sxyObs.shp") # This sampling buffer includes AC of observed individuals a bit outside the trapping array
Xbuf2$area <- st_area(Xbuf2) #m2
Xbuf2$area/1000000 # km2
Xbuf2$area/10000 # ha

## ---- Size core/periphery for methods ----

nuc <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/nuc.shp")
per <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/per.shp")

nuc$area <- st_area(nuc) #m2
per$area <- st_area(per) #m2

nuc$area/1000000 # km2
per$area/1000000 # km2

nuc$area/1000000 + per$area/1000000
