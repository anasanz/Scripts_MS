age_st <- function(nucleus, periphery, sampling_buffer){
  
  
  # 1. Total abundance in buffer area (sampling area)
  
  # Matrix to store abundance in the buffer each iteration and year
  NIn_trapBuf <- matrix(NA,nrow=dim(myResultsSXYZ$sims.list$z)[1],ncol=dim(myResultsSXYZ$sims.list$z)[3]) # nrow = iterations, ncol = year
  
  for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
    for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
      
      which.alive <- which(myResultsSXYZ$sims.list$z[ite,,t]==1) # Select only the individuals alive (z=1)
      
      which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
      
      sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(sampling_buffer))) # CONVERT SXY TO SPATIAL POINTS 
      
      which.In <- over(sp, sampling_buffer) # Check which ones are in the buffer
      
      NIn_trapBuf[ite,t] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
    }
  }
  
  # 2. Abundance of each age category in both polygons
  
  zones <- c(nucleus, periphery)
  
  ## CUBS
  ZZcubs <- myResultsSXYZ$sims.list$z
  ZZcubs[!myResultsSXYZ$sims.list$age.cat %in% c(1,2) ]  <- 3 # Considers all individuals that are not 1 and 2 as dead
  
  cub_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], length(zones))) # nrow = iterations, ncol = year, dim3 = nuc/periphery
  
  
  for(s in 1:length(zones)){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(ZZcubs[ite,,t]==1) # Select only the individuals alive (z=1)
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        cub_trapBuf[ite,t,s] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}
  
  
  ## SUBADULTS 
  
  ZZsub <- myResultsSXYZ$sims.list$z
  ZZsub[!myResultsSXYZ$sims.list$age.cat %in% c(3,4) ]  <- 3 # Considers all individuals that are not 3 and 4 (SUBADULTS) as dead
  
  subad_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], length(zones)))
  
  for(s in 1:length(zones)){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(ZZsub[ite,,t]==1) # Select only the individuals alive (z=1)
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        subad_trapBuf[ite,t,s] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}
  
  
  ## ADULTS
  
  ZZad <- myResultsSXYZ$sims.list$z
  ZZad[!myResultsSXYZ$sims.list$age.cat %in% c(5) ]  <- 3 # Considers all individuals that are not 3 and 4 (SUBADULTS) as dead
  
  ad_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], length(zones)))
  
  for(s in 1:length(zones)){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(ZZad[ite,,t]==1) # Select only the individuals alive (z=1)
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        ad_trapBuf[ite,t,s] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}
  
  ## ADULT FEMALES
  
  ZZadFEM <- myResultsSXYZ$sims.list$z
  ZZadFEM[!myResultsSXYZ$sims.list$age.cat %in% c(5) ]  <- 3 # Considers all individuals that are not 5 (ADULTS) as dead
  ZZadFEM[!myResultsSXYZ$sims.list$sex %in% c(0) ]  <- 3 # Considers all individuals that dont have a sex 0 (females) as dead
  
  adFEM_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], length(zones)))
  
  for(s in 1:length(zones)){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(ZZadFEM[ite,,t]==1) # Select only the individuals alive (z=1)
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        adFEM_trapBuf[ite,t,s] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}
  
  ## ADULT MALES
  
  ZZadMAL <- myResultsSXYZ$sims.list$z
  ZZadMAL[!myResultsSXYZ$sims.list$age.cat %in% c(5) ]  <- 3 # Considers all individuals that are not 5 (ADULTS) as dead
  ZZadMAL[!myResultsSXYZ$sims.list$sex %in% c(1) ]  <- 3 # Considers all individuals that dont have a sex 1 (males) as dead
  
  adMAL_trapBuf <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], dim(myResultsSXYZ$sims.list$z)[3], length(zones)))
  
  for(s in 1:2){
    for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
      for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
        
        which.alive <- which(ZZadMAL[ite,,t]==1) # Select only the individuals alive (z=1)
        
        which.aliveSXY <- myResultsSXYZ$sims.list$sxy[ite,which.alive,,t] # Retrieve the activity center for those individuals
        
        sp <- SpatialPoints(which.aliveSXY,proj4string=CRS(proj4string(zones[[s]]))) # CONVERT SXY TO SPATIAL POINTS 
        which.In <- over(sp, zones[[s]]) # Check which ones are in the buffer
        
        adMAL_trapBuf[ite,t,s] <- sum(which.In,na.rm = T) # The sum of the points in the buffer is the abundance that year and iteration. Store
      }}}
  
  
  # 3. Estimate abundance proportion of each age class 
  
  prop_years <- list()
  prop_years_zone <- list()
  
  prop <- array(NA,c(dim(myResultsSXYZ$sims.list$z)[1], 5, dim(myResultsSXYZ$sims.list$z)[3], length(zones))) # Dim1 = iterations, Dim2 = c(cub, subad, ad, adFEM, adMAL), Dim3 = Years, Dim4 = Core/Periphery
  
  for(s in 1:length(zones)){
    for(t in 1:dim(myResultsSXYZ$sims.list$z)[3]){
      for(ite in 1:dim(myResultsSXYZ$sims.list$z)[1]){
        prop[ite, 1, t, s] <- cub_trapBuf[ite,t,s]/NIn_trapBuf[ite,t]
        prop[ite, 2, t, s] <- subad_trapBuf[ite,t,s]/NIn_trapBuf[ite,t]
        prop[ite, 3, t, s] <- ad_trapBuf[ite,t,s]/NIn_trapBuf[ite,t]
        prop[ite, 4, t, s] <- adFEM_trapBuf[ite,t,s]/NIn_trapBuf[ite,t]
        prop[ite, 5, t, s] <- adMAL_trapBuf[ite,t,s]/NIn_trapBuf[ite,t]
      } } }
  
  return(prop)
  
}