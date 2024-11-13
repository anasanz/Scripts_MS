
## -------------------------------------------------
##    MODEL 7.4.3 BPA: CJS with temporal covariates 
##          Variation: Temporal + ID co-variate 
##        Covariate: Slope
## -------------------------------------------------


rm(list = ls())

library(dplyr)
library(tidyr)
library(rgdal)
library(raster)
library(rjags)
library(jagsUI)

## ---- DATA ----

# Load monitoring data

setwd("D:/MargSalas/Oso/Datos/Tablas_finales")
os <- read.csv("Data_os_96_20.csv", header = TRUE, row.names = NULL)
os <- os[,-1] 

os_id <- os[which(os$Confirmed_Individual != "Indetermined"), ] # Only identified
os_id$Method[which(os_id$Method == "sampling_station")] <- "Sampling_station"
os_id$Obs_type[which(os_id$Obs_type == "hair")] <- "Hair"

os_id <- os_id[which(os_id$Region == "Catalunya"), ]

coordinates(os_id) <- os_id[,c("x_long","y_lat")] # Spatial object
os_id@proj4string <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# Load and calculate slope

dem <- raster("D:/MargSalas/Oso/Datos/GIS/Variables/Catalunya/MDE30_rev2_ED50.tif")
dem_4326 <- projectRaster(dem, crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

slope <- terrain(dem_4326,'slope', unit = 'degrees', neighbors = 8, filename = 'slope')


# Extract coordinates

coord <- os_id[ ,c("x_long","y_lat")] 
cells <- cellFromXY(slope, coord) # 1. Tells the number of the cells where the coord. fall
slope_coords <- slope[cells] # 2. Gets the raster values of those cells

# Capture history

setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Functions")
source("capt_hist_bear.r")

d <- os_id@data

HairStCH <- capt_hist_bear(data = d, 
                           method = "Sampling_station", 
                           obs_type = c("Hair"))

ch <- data.frame(HairStCH$capt.hist[,c(1,3:13)])
rownames(ch) <- ch$Confirmed_Individual
ch <- ch[,-c(1)]
colnames(ch) <- c(2010:2020)

# --- Mean distance to paths per year and individual ----

# Data frame to fill

slo <- matrix(NA, nrow = nrow(ch), ncol = ncol(ch))
rownames(slo) <- rownames(ch)
colnames(slo) <- colnames(ch)

# Calculate distances

d <- cbind(d, slope_coords)

slope_sum <- d %>%
  group_by(Confirmed_Individual,Year) %>%
  summarise(
    mean = mean(slope_coords, na.rm = TRUE))

for (i in 1:nrow(slope_sum)){
  slo[which(rownames(slo) %in% slope_sum$Confirmed_Individual[i]),which(colnames(slo) %in% slope_sum$Year[i])] <- slope_sum$mean[i] 
}

# Select from 2010 to 2019 because it is a co-variate on survival, 
# which happens between years

slo <- slo[,c(1:10)]

# ---- Arrange data for JAGS ----

# Capture history
chb <- as.data.frame(HairStCH$capt.hist)
chb <- as.matrix(chb[,c(3:13)]) # Keep onlu from 2010 - 2020 (systematic sampling)

# Covariate distance to path

# There can't be NA in co-variate, so simulate values
mean_years <- apply(slo,2,mean, na.rm = TRUE)
sd_years <- apply(slo,2,sd, na.rm = TRUE)
number <- colSums(is.na(slo)) # How many values generate

for (i in 1:length(number)){
  col <- slo[,i] 
  col[is.na(col)] <- rnorm(number[i],mean_years[i],sd_years[i])
  slo[,i] <- col
}

# Standardize
slo_mean <- mean(slo)
slo_sd <- sd(slo)
slo_sc <- (slo - slo_mean) / slo_sd


# First observation
get.first <- function(x) min(which(x!=0)) # Important: Only initial values after first capture 
f_bear <- apply(chb, 1, get.first)

# JAGS data
setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Functions")
source("known.state.cjs.r")

data_bear <- list(y = chb, f = f_bear, nind = dim(chb)[1], n.occasions = dim(chb)[2], z = known.state.cjs(chb), x = slo_sc)

# Initial values

setwd("D:/MargSalas/Scripts_MS/Scripts_MS/Functions")
source("cjs.init.z.r")

inits <- function(){list(z = cjs.init.z(chb, f_bear), mu = rnorm(1), sigma = runif(1, 0, 5), beta = runif(1, -5, 5), mean.p = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("mean.phi", "mean.p", "phi.est", "sigma2", "beta")

# MCMC settings
ni <- 200000
nt <- 10
nb <- 10000
nc <- 3

# Call JAGS from R 
setwd("D:/MargSalas/Oso/SCR/Model")

outHairSt_slope_it<- jags(data_bear, inits, parameters, "cjs-cov-raneff[it].txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(outHairSt_slope_it, digits = 3)

# Summarize posteriors 

hist(outHairSt_slope_it$sims.list$mean.phi, nclass = 50, col = "gray", main = "", xlab = "Phi", las = 1, xlim = c(0,1))
abline(v = mean(outHairSt_slope_it$sims.list$mean.phi), col = "blue", lwd = 3)

hist(outHairSt_slope_it$sims.list$beta, nclass = 50, col = "gray", main = "", xlab = "Slope", las = 1, xlim = c(-12,1))
abline(v = mean(outHairSt_slope_it$sims.list$beta), col = "blue", lwd = 3)

# Save outputs
setwd("D:/MargSalas/Oso/SCR/Exp_analysis")
write.csv(outHairSt_slope_it$summary, file = "outHairSt_Cat_slope[it].csv")

