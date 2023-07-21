
## ---- Study area and location of core/periphery ----

rm(list = ls())

library(sf)


# Load europe
sa <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/clip_pyros2_WGS84_31N_all.shp")
eur <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/esp_fr_2.shp") %>%
  st_transform(dpts, crs = crs(sa))

nuc <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/nuc.shp")
per <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/per.shp")

# Plot
setwd("D:/MargSalas/Oso/OPSCR_project/Results/Results_section/Plots")
pdf("1.4.sa_coreper.pdf",7,7)

plot(st_geometry(sa), col = "beige", border = "white",  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000)) # To add andorra with same color (not present in eur)
plot(st_geometry(eur), col = "beige", border = "white",  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000), add = TRUE)
plot(st_geometry(nuc), col = "#9970ab", border = "#9970ab", add = TRUE)
plot(st_geometry(per), col = "#a6dba0", border = "#a6dba0", add = TRUE)
plot(st_geometry(eur), border = adjustcolor("white", alpha.f = 0.3),  xlim = c(st_bbox(sa)[1]-10000, st_bbox(sa)[3]+10000),  ylim = c(st_bbox(sa)[2] - 50000 , st_bbox(sa)[4] + 50000), add = TRUE)

dev.off()