## -------------------------------------------------
##            Increase occassions for 2022
## ------------------------------------------------- 


## ---- Data simulation with estimates Mt 2021 ----


# Simulate data with the estimates of 2021 under the Mt model (best according to AIC)
# Check how they would improve by increasing occassions

set.seed(2022)
data.fn <- function(N = 104, p = c(0.173, 0.211, 0.038), T = 3){
  yfull <- yobs <- array(NA, dim = c(N, T) )
  p.vec <- array(NA, dim = T)
  for (j in 1:T){
    yfull[,j] <- rbinom(n = N, size = 1, prob = p[j])
  }
  ever.detected <- apply(yfull, 1, max)
  C <- sum(ever.detected)
  yobs <- yfull[ever.detected == 1,]
  cat(C, "out of", N, "animals present were detected.\n")
  capt.hist_sim <- as.data.frame(apply( yobs[ ,1:T ] , 1 , paste , collapse = "" )) # Create capture history for mark
  colnames(capt.hist_sim)[1] <- "ch"
  
  return(list(N = N, p = p, C = C, T = T, yfull = yfull, yobs = yobs, capt.hist_sim = capt.hist_sim))
}

data <- data.fn()
data$C
capt.hist_sim <- data$capt.hist_sim

p.time.shared<-list(formula=~time,share=TRUE) # Mt, p=c, varying over time 
Mt_sim<-mark(data=capt.hist_sim, model="Closed", model.parameters=list(p=p.time.shared), 
         brief=T)
Mt_sim$results$real
Mt_sim$results$derived

uncertainty_sim <- Mt_sim$results$derived$`N Population Size`$ucl - Mt_sim$results$derived$`N Population Size`$lcl

## ---- Add 1 sampling occasion, see if the estimates improve ----

set.seed(2022)
data_4oc <- data.fn (N = 104, p = c(0.173, 0.211, 0.038, 0.1), T = 4)
data_4oc$C
capt.hist_4oc <- data_4oc$capt.hist_sim

# Mt, p=c, varying over time 
p.time.shared<-list(formula=~time,share=TRUE)
Mt_sim_4oc<-mark(data=capt.hist_4oc, model="Closed", model.parameters=list(p=p.time.shared), 
             brief=T)
Mt_sim_4oc$results$real
Mt_sim_4oc$results$derived

uncertainty_sim_4oc <- Mt_sim_4oc$results$derived$`N Population Size`$ucl - Mt_sim_4oc$results$derived$`N Population Size`$lcl

## ---- Add 2 sampling occasion, see if the estimates improve ----

set.seed(2022)
data_5oc <- data.fn (N = 104, p = c(0.173, 0.211, 0.038, 0.1, 0.1), T = 5)
data_5oc$C
capt.hist_5oc <- data_5oc$capt.hist_sim

# Mt, p=c, varying over time 
p.time.shared<-list(formula=~time,share=TRUE)
Mt_sim_5oc<-mark(data=capt.hist_5oc, model="Closed", model.parameters=list(p=p.time.shared), 
                 brief=T)
Mt_sim_5oc$results$real
Mt_sim_5oc$results$derived

uncertainty_sim_5oc <- Mt_sim_5oc$results$derived$`N Population Size`$ucl - Mt_sim_5oc$results$derived$`N Population Size`$lcl

uncertainty_sim
uncertainty_sim_4oc
uncertainty_sim_5oc
