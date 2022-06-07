

## -------------------------------------------------
##                 Example oSCR
## ------------------------------------------------- 

library(oSCR) 
data(rbs_ecography)

rbs.edf
unique(rbs.edf$Session)
# One edf for all sessions (with session identifier) -> 4 sessions
# One tdf per session. 
# Different number of occasions (days) each session

## ---- Prepare format the data for analysis using the data2oscr() ----

rbs.data <- data2oscr(edf = rbs.edf,
                      #the EDF
                      sess.col = 1,
                      #session column
                      id.col = 2,
                      #individual column
                      occ.col = 3,
                      #occasion column
                      trap.col = 4, # BOARD
                      #trap column
                      sex.col = 5,
                      #sex column (optional)
                      tdf = list(rbs.tdf1, rbs.tdf2,
                                 #list of TDFs
                                 rbs.tdf3, rbs.tdf4),
                      K = c(7,5,6,4), # Different number of occasions per season
                      #occasion vector
                      ntraps = c(50,50,50,50),
                      #no. traps vector
                      trapcov.names = c("jday","jday2"), #covariate names 
                      tdf.sep = "/",  #char used separator
                      sex.nacode = "U")   #unknown sex code

# The one we need for analysis using oSCR is the scrFrame
rbs.data$traplocs
rbs.sf <- rbs.data$scrFrame

par(mfrow=c(2,2), mar=c(2,2,2,2), oma=c(0,0,0,0)) 
plot(rbs.sf,jit = 2)
plot(rbs.sf)
d <- rbs.sf$caphist[[1]]

## ---- Create state space ----

rbs.ss <- make.ssDF(scrFrame = rbs.sf, buffer = 4, res=0.5)
plot(rbs.ss, scrFrame = rbs.sf)

## ---- Model fitting ----


if(1 == 2){    #this takes a *long* time, see below for loading these models 
  m1 <- oSCR.fit(model = list(D~1,       #density
                          p0~b + jday + jday2,   #detection 
                          sig~1),               #space use 
  scrFrame = rbs.sf, ssDF = rbs.ss) 


m16 <- oSCR.fit(model = list(D~session, #density
                            p0~b + jday + jday2 + session + sex, #detection 
                            sig~session + sex),#space use 
                scrFrame = rbs.sf, ssDF = rbs.ss)
}

# Run of 16 candidate models, take very long to run! So included in the package

data(rbs_ecography_mods)

fl <- fitList.oSCR(list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16), drop.asu=T, rename=TRUE) #rename=T adds sesnible model names
ms <- modSel.oSCR(fl)

## ---- Predict ----

library(car) 
library(ggplot2) 
library(ggthemes)

#pick the model with the lowest AIC 
topmod <- fl[[which.min(sapply(fl,function(x) x$AIC))]]

## DETECTION

#make a dataframe of values for DETECTION predictions 
p.pred.df <- data.frame(jday = rep(seq(-1,7,length=100),2), #obs range 
                        jday2 = rep(seq(-1,7,length=100)^2,2),#obs range^2 
                        b=0,#pr(initial capture)
                        sex=rep(c(0,1),each=100)) #binary sex
#now predict 
p.preds <- get.real(model = topmod, type = "det", newdata = p.pred.df) 
p.preds$sex <- factor(p.preds$sex) 
levels(p.preds$sex) <- c("Female","Male") 
head(p.preds)

#now plot 
ggplot(p.preds, aes(x=jday*10, y=estimate, fill = sex)) + 
  geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.1, colour=NA) + 
  geom_line(size=1) + facet_grid(.~sex) + 
  theme_bw() + scale_color_fivethirtyeight() + 
  xlab("Days since Sept 1st")

## DENSITY

#make a dataframe of values for DENSITY predictions 
# Predictions of density will be on the ‘per-pixel’ scale
d.pred.df <- data.frame(session = rep(factor(1:4),2)) #session specific

#now predict 
d.preds <- get.real(model = topmod, type = "dens", newdata = d.pred.df, d.factor = 4) 
d.preds$sex <- factor(d.preds$sex) 
levels(d.preds$sex) <- c("Female","Male") 
head(d.preds)
