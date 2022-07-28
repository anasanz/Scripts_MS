


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

devtools::install_github('oswaldosantos/ggsn')



# Load monitoring data
setwd("D:/MargSalas/Oso/Datos/Tablas_finales/2021")
os <- read.csv("Natal_established_coordinates.csv", header = TRUE, row.names = NULL)
os <- os[,-1] 

# Set coordinates
coordinates(os) <- os[,c("x_long","y_lat")] # Spatial object
os@proj4string <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# ID with more than 10 0bs
setwd("D:/MargSalas/Oso/Datos")
load("id_more10.RData")

# Load mcp
setwd("D:/MargSalas/Oso/Datos/GIS/2021")
load("mcp_natal.RData")
load("mcp_est2.RData")

names_mcp_nat <- lapply(mcp_natal, `[[`, 1) # Take only the first argument of the list to know the name
names_mcp_est <- lapply(mcp_est2, `[[`, 1)


# Load basemap
map1 <- readOGR(dsn = "D:/MargSalas/Oso/Datos/GIS/Countries", layer = "clip_pyros2")

# Load info database
setwd("D:/MargSalas/Oso/Datos")
info <- read.csv("Info_individuals.csv", header = TRUE, row.names = NULL, sep = ";")
info <- info[,c(4:9)]

# Load dispersal distances
setwd("D:/MargSalas/Oso/Datos")
d <- read.csv("disp_distance.csv", header = TRUE, row.names = NULL, sep = ",")
d <- d[,-1]
# Funtions
'%!in%' <- function(x,y)!('%in%'(x,y))

extract_legend <- function(my_ggp) { # Extracts legends from ggplots
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}

## -------------------------------------------------
##                    MALES
## ------------------------------------------------- 

os_males <- os[which(os$Sex_cub == "M"), ]
id_males <- unique(os_males$Confirmed_Individual_cub)
id_males <- id_males[id_males %in% id_more10] # ID males with more than 10 observations (all)

# Unify age classes to plot
unique(os_males$Age_class_cub)
os_males$Age_class_cub[which(os_males$Age_class_cub %in% c("Cub0","Cub1", "Cub2"))] <- "Cub"

#pdf("ID_Males2.pdf")  


for (i in 1:length(id_males)) {
  os_id <- os_males[which(os_males$Confirmed_Individual_cub == id_males[i]), ] # Select individual
  df <- os_id@data
  ext <- data.frame(bbox(os_id)) # Extent
  myMap <- get_map(location = bbox(os_id), source = "stamen", maptype = "terrain", crop = FALSE) # Basemap
  
  mcp_nat_id <- which(names_mcp_nat %in% id_males[i]) # MCP natal
  
  ifelse(length(mcp_nat_id) > 0, x <- mcp_natal[[which(names_mcp_nat %in% id_males[i])]], x <- 0)
  
  
  # ---- Plot location within study area ----
  p0 <- ggplot() + 
        geom_polygon(data = map1, aes(x = long, y = lat, group = group), colour = "black", fill = "lightgrey", size = 0.3) +
        geom_rect(data = ext, aes(xmin=ext[1,1] , xmax=ext[1,2], ymin=ext[2,1], ymax=ext[2,2]), colour = "red", fill = "transparent", size = 0.5) +
        theme(panel.grid = element_blank(),
              axis.title=element_blank(),
              axis.text=element_blank(),
              axis.ticks=element_blank(),
              #plot.margin = unit(c(1,2,1,2), "cm")
              )    

  # ---- Plot points my year and age class ----
  p1 <- ggmap(myMap) +
    
    geom_point(data = df, aes(x = x_long, y = y_lat, colour = factor(Year), shape = factor(Age_class_cub)), 
               size = 1.5) +
    scale_color_viridis(discrete = TRUE) +
    #scale_color_brewer(palette = "Dark2") +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none"
          #plot.margin = unit(c(-1,-1,-1,-1), "cm")
          ) +
    ggtitle("All locations")
    

  # ---- Plot average location by year ----
  
  # 1. Average x and y
  mean_loc <- df %>%
    group_by(Year, Age_class_cub ) %>%
    summarise(
      count = n(),
      meanx = mean(x_long, na.rm = TRUE),
      meany = mean(y_lat, na.rm = TRUE))
  na <- which(is.na(mean_loc$Year)) 
  if(length(na) > 0) { # If there is NA in the data, remove them
    mean_loc[-which(is.na(mean_loc$Year)), ]
  }


  # 2. Plot(different with and without MCP)
  
  
  if (length(x) >= 1){ # If there is a mcp polygon
  
    if(id_males[i] %in% d$ID) { # And there is enough positions to create a natal-est arrow
      
          if(!is.na(d[d$ID %in% id_males[i],c(4)])) {
            
            p2 <- ggmap(myMap)+
              geom_point(data = mean_loc, aes(x = meanx , y = meany, colour = factor(Year), shape = factor(Age_class_cub)), 
                         size = 1.5) +
              scale_color_viridis(discrete = TRUE) +
              #scale_color_brewer(palette = "Dark2") +
              theme(panel.grid = element_blank(),
                    axis.title = element_blank(),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    legend.position = "none"
                    #plot.margin = unit(c(-1,-1,-1,-1), "cm")
              ) +
              geom_polygon(data = x[[2]], aes(x = long, y = lat, group = group), 
                           colour = "black", fill = "transparent", size = 0.5) +
              geom_segment(aes(x = d[d$ID %in% id_males[i],c(2)], y = d[d$ID %in% id_males[i],c(3)], 
                               xend = d[d$ID %in% id_males[i],c(4)], yend = d[d$ID %in% id_males[i],c(5)]),
                           arrow = arrow(length = unit(0.1, "cm")),
                           lwd = 0.2) +
              ggtitle("Average location")
            
          } else { p2 <- ggmap(myMap)+
                        geom_point(data = mean_loc, aes(x = meanx , y = meany, colour = factor(Year), shape = factor(Age_class_cub)), 
                                   size = 1.5) +
                        scale_color_viridis(discrete = TRUE) +
                        #scale_color_brewer(palette = "Dark2") +
                        theme(panel.grid = element_blank(),
                              axis.title = element_blank(),
                              axis.text = element_blank(),
                              axis.ticks = element_blank(),
                              legend.position = "none"
                              #plot.margin = unit(c(-1,-1,-1,-1), "cm")
                        ) +
                        geom_polygon(data = x[[2]], aes(x = long, y = lat, group = group), 
                                     colour = "black", fill = "transparent", size = 0.5) +
                        geom_segment(aes(x = d[d$ID %in% id_males[i],c(2)], y = d[d$ID %in% id_males[i],c(3)], 
                                         xend = d[d$ID %in% id_males[i],c(6)], yend = d[d$ID %in% id_males[i],c(7)]),
                                     arrow = arrow(length = unit(0.1, "cm")),
                                     lwd = 0.2) +
                        ggtitle("Average location") }} else {p2 <- ggmap(myMap)+
                                                            geom_point(data = mean_loc, aes(x = meanx , y = meany, colour = factor(Year), shape = factor(Age_class_cub)), 
                                                                           size = 1.5) +
                                                            scale_color_viridis(discrete = TRUE) +
                                                            #scale_color_brewer(palette = "Dark2") +
                                                            theme(panel.grid = element_blank(),
                                                                      axis.title = element_blank(),
                                                                      axis.text = element_blank(),
                                                                      axis.ticks = element_blank(),
                                                                      legend.position = "none"
                                                                      #plot.margin = unit(c(-1,-1,-1,-1), "cm")
                                                            ) +
                                                            ggtitle("Average location")}}
        
    
  # Common plot:
  # TITLE
  title1 = text_grob(paste(id_males[i]," (",info$Year_birth[info$ID %in% id_males[i]], "-", info$Year_death[info$ID %in% id_males[i]], ")", sep = ""), size = 15, face = "bold")

  # LEGEND: 
  #1. Create plot with legend
  ggp1_legend <- ggmap(myMap) +
                geom_point(data = df, aes(x = x_long, y = y_lat, colour = factor(Year), shape = factor(Age_class_cub)), 
                size = 2) +
                scale_color_viridis(name = "Year",discrete = TRUE) +
                scale_shape_discrete(name = "Age class") +
                labs(tag = "- MCP Natal") +
                #scale_color_brewer(palette = "Dark2") +
                theme(
                  plot.tag.position = c(0.1,1),
                  legend.position = "top",
                  legend.justification=c(0.5,0.5),
                  legend.box = "vertical")
                  
                  #plot.tag = element_text(size = 30) 
                  #legend.position=c(1,0)) +
                #scale_colour_discrete(name = "Year") +
 
  #2. Apply user-defined function to extract legend
  shared_legend <- extract_legend(ggp1_legend)
  
  setwd("D:/MargSalas/Oso/Datos/Plots/Detailed_ID2")
  pdf(paste(id_males[i], ".pdf", sep = ""),
      width = 8, height = 5)
  grid.arrange(p0,p1,p2, nrow = 2, ncol = 2,
               top = title1,
               shared_legend)
  dev.off()
}


