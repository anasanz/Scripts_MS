


## -------------------------------------------------
##          Explore individual dispersal 
## ------------------------------------------------- 

rm(list = ls())


library(ggplot2)
library(ggmap)
library(gridExtra)
library(rgdal)
library(ggpubr)
library(dplyr)
library(geosphere)
library("viridis")
library(RColorBrewer)



# Load monitoring data
setwd("D:/MargSalas/Oso/Datos/Tablas_finales")
os <- read.csv("Natal_established_coordinates.csv", header = TRUE, row.names = NULL)
os <- os[,-1] 

# Set coordinates
coordinates(os) <- os[,c("x_long","y_lat")] # Spatial object
os@proj4string <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# ID with more than 10 0bs
setwd("D:/MargSalas/Oso/Datos")
load("id_more10.RData")

# Load mcp
setwd("D:/MargSalas/Oso/Datos/GIS")
load("mcp_natal.RData")
load("mcp_est2.RData")

names_mcp_nat <- lapply(mcp_natal, `[[`, 1)
names_mcp_est <- lapply(mcp_est2, `[[`, 1)


# Load basemap
map1 <- readOGR(dsn = "D:/MargSalas/Oso/Datos/GIS/Countries", layer = "clip_pyros2")

## -------------------------------------------------
##                    MALES
## ------------------------------------------------- 

os_males <- os[which(os$Sex_cub == "M"), ]
id_males <- unique(os_males$Confirmed_Individual_cub)
id_males <- id_males[id_males %in% id_more10] # ID males with more than 10 observations (all)

# Visual 2: Three categories (Cub, Subadult, Adult) and gradient of colours within

color_cub <- c("lightyellow1", "yellow", "gold", "orange") # Create three color palettes
col_subadult <- c("darkorange1", "orangered","firebrick1", "red1")
color_adult <- c("red3", "red4", "darkred", "black")


#pdf("ID_Males2.pdf")  



for (i in 1:length(id_males)) {
  os_id <- os_males[which(os_males$Confirmed_Individual_cub == id_males[i]), ] # Select individual
  df <- os_id@data
  ext <- data.frame(bbox(os_id)) # Extent
  myMap <- get_map(location = bbox(os_id), source = "stamen", maptype = "terrain", crop = FALSE) # Basemap
  
  mcp_nat_id <- which(names_mcp_nat %in% id_males[i]) # MCP natal
  ifelse(mcp_nat_id > 0, x <- mcp_natal[[which(names_mcp_nat %in% id_males[i])]], x <- 0)
  
  
  # Plot location within study area
  p0 <- ggplot() + 
        geom_polygon(data = map1, aes(x = long, y = lat, group = group), colour = "black", fill = "lightgrey") +
        geom_rect(data = ext, aes(xmin=ext[1,1] , xmax=ext[1,2], ymin=ext[2,1], ymax=ext[2,2]), colour = "red", fill = "transparent", size = 1) +
        theme(panel.grid = element_blank(),
              axis.title=element_blank(),
              axis.text=element_blank(),
              axis.ticks=element_blank(),
              #plot.margin = unit(c(1,2,1,2), "cm")
              )    

  # Plot points my year and age class
  p1 <- ggmap(myMap)+
    geom_point(data = df, aes(x = x_long, y = y_lat, colour = factor(Year), shape = factor(Age_class_cub)), 
               size = 2) +
    scale_color_viridis(discrete = TRUE) +
    #scale_color_brewer(palette = "Dark2") +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          #plot.margin = unit(c(-1,-1,-1,-1), "cm")
          )
    

  # Plot average location by year
  
  # 1. Average x and y
  mean_loc <- df %>%
    group_by(Year) %>%
    summarise(
      count = n(),
      meanx = mean(x_long, na.rm = TRUE),
      meany = mean(y_lat, na.rm = TRUE))

  # 2. Plot

  p2 <- ggmap(myMap)+
    geom_point(data = mean_loc, aes(x = meanx , y = meany, colour = factor(Year)), 
               size = 2) +
    scale_color_viridis(discrete = TRUE) +
    #scale_color_brewer(palette = "Dark2") +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          #plot.margin = unit(c(-1,-1,-1,-1), "cm")
          ) +
    geom_polygon(data = x[[2]], aes(x = long, y = lat, group = group), colour = "black", fill = "transparent", size = 1)
    


  # Common plot:
  # TITLE
  title1 = text_grob(id_males[i], size = 20, face = "bold", hjust = c(5,0))
  # LEGEND: 
  #1. Create plot with legend
  ggp1_legend <- ggmap(myMap)+
                geom_point(data = df, aes(x = x_long, y = y_lat, colour = factor(Year), shape = factor(Age_class_cub)), 
                size = 2) +
                scale_color_viridis(name = "Year",discrete = TRUE) +
                #scale_color_brewer(palette = "Dark2") +
                theme(
                  legend.position = "top",
                  legend.justification=c(0.5,0.5),
                  legend.box = "vertical") + 
                  #legend.position=c(1,0)) +
                #scale_colour_discrete(name = "Year") +
                scale_shape_discrete(name = "Age class")
  #2. Create user-defined function, which extracts legends from ggplots
  extract_legend <- function(my_ggp) {
    step1 <- ggplot_gtable(ggplot_build(my_ggp))
    step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
    step3 <- step1$grobs[[step2]]
    return(step3)
  }
  
  #3. Apply user-defined function to extract legend
  shared_legend <- extract_legend(ggp1_legend)
  
  grid.arrange(p0,p1,p2, nrow = 2, ncol = 2,
               top = title1,
               shared_legend)
}

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
#legend('center',legend = c("Cub", "Subadult", "Adult"), col = c("darkgreen", "Orange", "Purple"), fill = c("Green", "Orange", "Purple"), xpd = TRUE, horiz = TRUE, cex = 1.5, seg.len=1, bty = 'n')


x <- c(1,1,1,2,2,2,3,3,3)
y <- c(1,2,3,1,2,3,1,2,3)
xy <- cbind(x,y)
S <- SpatialPoints(xy)
bbox(S)
