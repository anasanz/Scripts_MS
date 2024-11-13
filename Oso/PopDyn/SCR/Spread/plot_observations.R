


# Map to show observations and increase in spread

library(sf)
library(tmap)
library(ggspatial)
library(ggplot2)


# Load europe
sa <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/clip_pyros2_WGS84_31N_all.shp")
eur <- st_read("D:/MargSalas/Oso/Datos/GIS/Countries/esp_fr_2.shp") %>%
  st_transform(dpts, crs = st_crs(sa))
nuc <- st_read("D:/MargSalas/Oso/Datos/GIS/Variables/Europe/DistCore/corepol/corepol.shp")
obs <- st_read("D:/MargSalas/Oso/Datos/GIS/2022/Seguiment_GIS_layer/Seguiment_Ossos_Pirineus_1996_2021_coordinates_final.shp") %>%
  st_transform(crs = st_crs(sa))
obs <- obs[which(obs$Obs_typ %in% c("Hair", "Photo", "Photo/video", "Video")), ]

year <- c(2017, 2018, 2019, 2020, 2021)

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Spread")

par(mfrow = c(3,2))

for (i in 1:length(year)){
  obs_year <- obs[which(obs$Year %in% year[i]), ]
  pdf(file = paste("obs_", year[i], ".pdf", sep = ""), 7, 7)
  ggplot(data = sa) +
    geom_sf() + 
    ggtitle(year[i]) +
    annotation_scale(location = "bl", width_hint = 0.2) +
    geom_sf(data = obs_year, colour = "yellow4") +
    geom_sf(data = nuc, colour = "beige", fill = adjustcolor("yellow2", alpha.f = 0.3))
  dev.off()
}

p1 <- ggplot(data = sa) +
  geom_sf() + 
  ggtitle(2017) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  geom_sf(data = obs[which(obs$Year %in% c(2017)), ], colour = "yellow4") +
  geom_sf(data = nuc, colour = "beige", fill = adjustcolor("yellow2", alpha.f = 0.3)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(2, 0, 0, 0), "pt"))

p2 <- ggplot(data = sa) +
  geom_sf() + 
  ggtitle(2018) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  geom_sf(data = obs[which(obs$Year %in% c(2018)), ], colour = "yellow4") +
  geom_sf(data = nuc, colour = "beige", fill = adjustcolor("yellow2", alpha.f = 0.3)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(2, 0, 0, 0), "pt"))

p3 <- ggplot(data = sa) +
  geom_sf() + 
  ggtitle(2019) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  geom_sf(data = obs[which(obs$Year %in% c(2019)), ], colour = "yellow4") +
  geom_sf(data = nuc, colour = "beige", fill = adjustcolor("yellow2", alpha.f = 0.3)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(2, 0, 0, 0), "pt"))


p4 <- ggplot(data = sa) +
  geom_sf() + 
  ggtitle(2020) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  geom_sf(data = obs[which(obs$Year %in% c(2020)), ], colour = "yellow4") +
  geom_sf(data = nuc, colour = "beige", fill = adjustcolor("yellow2", alpha.f = 0.3)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(2, 0, 0, 0), "pt"))


p5 <- ggplot(data = sa) +
  geom_sf() + 
  ggtitle(2021) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  geom_sf(data = obs[which(obs$Year %in% c(2021)), ], colour = "yellow4") +
  geom_sf(data = nuc, colour = "beige", fill = adjustcolor("yellow2", alpha.f = 0.3)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(2, 0, 0, 0), "pt"))

library("gridExtra")
grid.arrange(p1, p2, p3, p4, p5, 
             ncol = 2, nrow = 3)


## ---- Number of independent offspring ----

setwd("D:/MargSalas/Scripts_MS/Oso/PopDyn/SCR/Spread")
load("ages.RData")
age <- age.cat-1
apply(age, c(2), function(x){length(which(x == 3))})


