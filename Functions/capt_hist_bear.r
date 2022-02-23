
## -------------------------------------------------
##     Generate capture histories for bear data
## ------------------------------------------------- 


capt_hist_bear <- function(data = data, method = method, obs_type = obs_type){
  require(dplyr)
  require(tidyr)
  
  # Select data
  data_method <- data[which(data$Method %in% method & data$Obs_type %in% obs_type), ] # Selected observations
  data_method <- data_method[,colnames(data_method) %in% c("Confirmed_Individual", "Year")] # Select relevant columns for ch
  
  data_method$Year <- as.numeric(data_method$Year) # Column occasion
  data_method$detect <- 1 # All are detections
  # Create capture history
  capt.hist.method <- data_method %>%
    # remove duplicates, which may occur when individuals are caught multiple times in an event
    # For example, your event may be a year and an individual may be caught multiple times in a year.
    distinct() %>%
    # spread out data. The fill = 0 adds rows for combinations of id and event where individuals were not observerd
    spread(Year, detect, fill = 0) %>% 
    # For every individual....
    group_by(Confirmed_Individual) 
  
  cat(nrow(capt.hist.method), "out of", length(unique(data$Confirmed_Individual)), "animals were identified by", obs_type, "in", method)
  return(list(capt.hist = capt.hist.method, Nident = nrow(capt.hist.method), Nsurv = sum(capt.hist.method$`2020`)))
}
