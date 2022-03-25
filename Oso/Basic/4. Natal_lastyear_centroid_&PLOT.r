## -------------------------------------------------
##      Rough centroid of natal territory 
##            and last year positions
## ------------------------------------------------- 

rm(list = ls())

library(geosphere)
library(adehabitatHR)
library(dplyr)
library(rgdal)

# Load monitoring data
setwd("D:/MargSalas/Oso/Datos/Tablas_finales")
os <- read.csv("Natal_established_coordinates.csv", header = TRUE, row.names = NULL)
os <- os[,-1] 

# Set coordinates
coordinates(os) <- os[,c("x_long","y_lat")] # Spatial object
os@proj4string <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# Load basemap
map1 <- readOGR(dsn = "D:/MargSalas/Oso/Datos/GIS/Countries", layer = "clip_pyros2")

# Load info database
setwd("D:/MargSalas/Oso/Datos")
info <- read.csv("Info_individuals.csv", header = TRUE, row.names = NULL, sep = ";")
info <- info[,c(4:9)]

## ---- Calculate centroids ----

id <- unique(os$Confirmed_Individual_cub)
dat <- data.frame(matrix(NA, nrow = length(id), ncol = 7))
colnames(dat) <- c("ID", "X_nat", "Y_nat", "X_est1", "Y_est1", "X_est2", "Y_est2")
dat$ID <- id

mcp_natal <- list()
o <- list()

# 1. Natal
for(i in 1:length(id)){
  os_id <- os[which(os$Confirmed_Individual_cub %in% id[i]), ] # Select individual
  pos_nat <- os_id[which(os_id$Age_class_cub == "Cub0" | os_id$Age_class_cub == "Cub1"| os_id$Age_class_cub == "Cub2"), ] # Select positions of 1st year
  if(nrow(pos_nat) < 5 ) next
  mx <- mean(coordinates(pos_nat)[,1])
  my <- mean(coordinates(pos_nat)[,2])
  dat[dat$ID %in% id[i], which(colnames(dat) %in% c("X_nat", "Y_nat"))] <- c(mx,my) #♣ Calculate centroid
  o[[1]] <- id[i]
  o[[2]] <- mcp(pos_nat)
  mcp_natal[[i]] <- o  # Store MCP to plot it
    }

# Check
par(mfrow = c(1,1))
os_id <- os[which(os$Confirmed_Individual_cub %in% id[i]), ] # Select individual
pos_nat <- os_id[which(os_id$Age_class_cub == "Cub0" | os_id$Age_class_cub == "Cub1"| os_id$Age_class_cub == "Cub2"), ] # Select positions of 1st year

plot(map1, col = "lightgrey", border = "grey")
points(pos_nat, col = adjustcolor("darkgreen", alpha.f = 0.5), pch = 19)
mx <- mean(coordinates(pos_nat)[,1])
my <- mean(coordinates(pos_nat)[,2])
points(mx, my, pch = 3)
plot(mcp_natal[[5]][[2]], add = TRUE)

# 2. Established

mcp_est1 <- list()
mcp_est2 <- list()
o <- list()

for(i in 1:length(id)){
  os_id <- os[which(os$Confirmed_Individual_cub %in% id[i]), ]
  pos_est <- os_id[which(os_id$Age_class_cub == "Subadult" | os_id$Age_class_cub == "Adult"), ]
  # Define the last and 2 last years of data
  last_year <- max(pos_est$Year, na.rm = TRUE) 
  last_2year <- c(last_year, last_year - 1)
  pos_est1 <- pos_est[which(pos_est$Year %in% last_year), ]
  pos_est2 <- pos_est[which(pos_est$Year %in% last_2year), ]
  
  if(nrow(pos_est2) < 5 ) next # First the positions in the 2 last years, as it is the limiting factor
  dat[dat$ID %in% id[i], which(colnames(dat) %in% c("X_est2", "Y_est2"))] <- centroid(coordinates(pos_est2))
  o[[1]] <- id[i]
  o[[2]] <- mcp(pos_est2)
  mcp_est2[[i]] <- o  # Store MCP to plot it
  
  if(nrow(pos_est1) < 5 ) next
  dat[dat$ID %in% id[i], which(colnames(dat) %in% c("X_est1", "Y_est1"))] <- centroid(coordinates(pos_est1))
  o[[1]] <- id[i]
  o[[2]] <- mcp(pos_est1)
  mcp_est1[[i]] <- o  # Store MCP to plot it
  
}

# Save MCP
setwd("D:/MargSalas/Oso/Datos/GIS")

save(mcp_natal, file = "mcp_natal.RData")
save(mcp_est1, file = "mcp_est1.RData")
save(mcp_est2, file = "mcp_est2.RData")

names_mcp <- lapply(mcp_natal, `[[`, 1)
which(names_mcp %in% "Caramelles")

## ---- Dispersal distances ----

d <- dat[which(!is.na(dat$X_nat) & !is.na(dat$X_est2)), ] #
d$Distance <- NA
for (i in 1:nrow(d)){
  if(!is.na(d[i,c(4)])) {
  dist <- distm(d[i,c(2,3)], d[i,c(4,5)], fun = distGeo) } else {dist <- distm(d[i,c(2,3)], d[i,c(6,7)], fun = distGeo)} # Calculate distance between the natal (Cub0 - Cub2) to the last year positions
  d$Distance[i] <- dist
}
d$Distance <- round(d$Distance,digits = 0)
info_sex <- info[,c(1,2)]
d <- left_join(d,info_sex)

setwd("D:/Oso/Datos")
write.csv(d, "disp_distance.csv")

# Calculate dispersal distance by sex

stderr <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

mean_dist <- d %>%
  group_by(Sex) %>%
  summarise(
    count = n(),
    mean = mean(Distance, na.rm = TRUE),
    sd = sd(Distance, na.rm = TRUE),
    se = stderr(Distance, na.rm = TRUE)
  )

# Plot

setwd("D:/Oso/Datos/Plots")
pdf("Disp_distances.pdf")

plot(mean_dist$mean, ylim = c(0,55000), xlim = c(0,3), axes = FALSE, xaxs="i", yaxs = "i",
     ylab = " ", xlab = " ", pch = 19)
axis(1, labels = c("Female", "Male"), at = c(1,2))
axis(2, pos = 0.5)
mtext("Dispersal distance (m)", side = 2, line = 0, cex = 1.2)
x <- c(1,2)
segments(x, mean_dist$mean - mean_dist$se * 2, 
         x, mean_dist$mean + mean_dist$se * 2, 
         lwd = 1.5)
arrows(x, mean_dist$mean - mean_dist$se * 2, 
       x, mean_dist$mean + mean_dist$se * 2, 
       lwd = 1.5, angle = 90, code = 3, length = 0.05)
dev.off()

# Anova

d2 <- data.frame(dist = d$Distance, sex = d$Sex)
m <- aov(dist ~ sex, data = d2)
summary(m)


## ---- Plot ----

# Plot all individuals, and add dispersal arrow for those for which I could calculate the centroids

#### MALES ####
setwd("D:/MargSalas/Oso/Datos")
load("id_more10.RData")

'%!in%' <- function(x,y)!('%in%'(x,y))

os_males <- os[which(os$Sex_cub == "M"), ]
id_males <- unique(os_males$Confirmed_Individual_cub)
id_males <- id_males[id_males %in% id_more10] # ID males with more than 10 observations (all)

# Three categories: Cub, subadult, adult

setwd("D:/Oso/Datos/Plots")

pdf("ID_Males_&Disp1.pdf",7,9)  

par(mfrow = c(4,3),
    oma = c(2,2,4,1),
    mar = c(0,0.5,1,0))

for (i in 1:12) {
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
  
  if(id_males[i] %!in% d$ID) next # If there is was enough information to calculate natal and established coordinates, next
  
  if(!is.na(d[d$ID %in% id_males[i],c(4)])) {
    arrows(x0 = d[d$ID %in% id_males[i],c(2)], y0 = d[d$ID %in% id_males[i],c(3)], 
           x1 = d[d$ID %in% id_males[i],c(4)], y1 = d[d$ID %in% id_males[i],c(5)],
           length = 0.1, col = "red", lwd = 1.7) 
  } else { arrows(x0 = d[d$ID %in% id_males[i],c(2)], y0 = d[d$ID %in% id_males[i],c(3)], 
                  x1 = d[d$ID %in% id_males[i],c(6)], y1 = d[d$ID %in% id_males[i],c(7)],
                  length = 0.1, col = "red", lwd = 1.7)
  }
  
  mtext("Males", side = 3, line = 1.5, cex = 1.2, outer = TRUE)
  
}

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',legend = c("Cub", "Subadult", "Adult"), col = c("darkgreen", "Orange", "Purple"), fill = c("darkgreen", "orange", "purple"), xpd = TRUE, horiz = TRUE, cex = 1.5, seg.len=1, bty = 'n')

dev.off()

pdf("ID_Males_&Disp2.pdf") 

par(mfrow = c(3,3),
    oma = c(2,2,4,1),
    mar = c(0,0.5,1,0))

for (i in 13:21) {
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
  
  if(id_males[i] %!in% d$ID) next # If there is was enough information to calculate natal and established coordinates, next
  
  if(!is.na(d[d$ID %in% id_males[i],c(4)])) {
    arrows(x0 = d[d$ID %in% id_males[i],c(2)], y0 = d[d$ID %in% id_males[i],c(3)], 
           x1 = d[d$ID %in% id_males[i],c(4)], y1 = d[d$ID %in% id_males[i],c(5)],
           length = 0.1, col = "red", lwd = 1.7) 
  } else { arrows(x0 = d[d$ID %in% id_males[i],c(2)], y0 = d[d$ID %in% id_males[i],c(3)], 
                  x1 = d[d$ID %in% id_males[i],c(6)], y1 = d[d$ID %in% id_males[i],c(7)],
                  length = 0.1, col = "red", lwd = 1.7)
  }
  #mtext("Males", side = 3, line = 1.5, cex = 1.2, outer = TRUE)
  
}

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',legend = c("Cub", "Subadult", "Adult"), col = c("darkgreen", "Orange", "Purple"), fill = c("darkgreen", "orange", "purple"), xpd = TRUE, horiz = TRUE, cex = 1.5, seg.len=1, bty = 'n')

dev.off()


#### FEMALES ####

os_females <- os[which(os$Sex_cub == "F"), ]
id_females <- unique(os_females$Confirmed_Individual_cub)
id_females <- id_females[id_females %in% id_more10] # ID females with more than 10 observations (all)

setwd("D:/Oso/Datos/Plots")

pdf("ID_Females_&Disp1.pdf",7,9)   

par(mfrow = c(4,3),
    oma = c(2,2,4,1),
    mar = c(0,0.5,1,0))

for (i in 1:12) {
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
  
  if(id_females[i] %!in% d$ID) next # If there is was enough information to calculate natal and established coordinates, next
  
  if(!is.na(d[d$ID %in% id_females[i],c(4)])) {
    arrows(x0 = d[d$ID %in% id_females[i],c(2)], y0 = d[d$ID %in% id_females[i],c(3)], 
           x1 = d[d$ID %in% id_females[i],c(4)], y1 = d[d$ID %in% id_females[i],c(5)],
           length = 0.1, col = "red", lwd = 1.7) 
  } else { arrows(x0 = d[d$ID %in% id_females[i],c(2)], y0 = d[d$ID %in% id_females[i],c(3)], 
                  x1 = d[d$ID %in% id_females[i],c(6)], y1 = d[d$ID %in% id_females[i],c(7)],
                  length = 0.1, col = "red", lwd = 1.7) }
  
  mtext("Females", side = 3, line = 1.5, cex = 1.2, outer = TRUE)
  
}

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',legend = c("Cub", "Subadult", "Adult"), col = c("darkgreen", "orange", "purple"), fill = c("darkgreen", "orange", "purple"), xpd = TRUE, horiz = TRUE, cex = 1.5, seg.len=1, bty = 'n')

dev.off()

setwd("D:/Oso/Datos/Plots")

pdf("ID_Females_&Disp2.pdf")  

par(mfrow = c(3,3),
    oma = c(2,2,4,1),
    mar = c(0,0.5,1,0))

for (i in 13:18) {
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
  
  if(id_females[i] %!in% d$ID) next # If there is was enough information to calculate natal and established coordinates, next
  
  if(!is.na(d[d$ID %in% id_females[i],c(4)])) {
    arrows(x0 = d[d$ID %in% id_females[i],c(2)], y0 = d[d$ID %in% id_females[i],c(3)], 
           x1 = d[d$ID %in% id_females[i],c(4)], y1 = d[d$ID %in% id_females[i],c(5)],
           length = 0.1, col = "red", lwd = 1.7) 
  } else { arrows(x0 = d[d$ID %in% id_females[i],c(2)], y0 = d[d$ID %in% id_females[i],c(3)], 
                  x1 = d[d$ID %in% id_females[i],c(6)], y1 = d[d$ID %in% id_females[i],c(7)],
                  length = 0.1, col = "red", lwd = 1.7) }
  
  #mtext("Females", side = 3, line = 1.5, cex = 1.2, outer = TRUE)
  
}


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend(-0.5,-0.3,legend = c("Cub", "Subadult", "Adult"), col = c("darkgreen", "Orange", "Purple"), fill = c("darkgreen", "orange", "purple"), xpd = TRUE, horiz = TRUE, cex = 1.5, seg.len=1, bty = 'n')

dev.off()
