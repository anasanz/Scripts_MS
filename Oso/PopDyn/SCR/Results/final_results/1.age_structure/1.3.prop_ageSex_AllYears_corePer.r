
## ---- Proportion/Ageclass/sex all years (no bars) ----

rm(list = ls())

library(nimbleSCR)
library(nimble)
library(rgdal)

source("D:/MargSalas/Scripts_MS/Functions/plot.violins3.r")
source("D:/MargSalas/Scripts_MS/Functions/DoScale.r")
source("D:/MargSalas/Scripts_MS/Functions/PlotViolinsHoriz.r")

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams")
load("Nbuffer_newSize.RData") # Obtained from script: 3.3. Core_vS_Periphery_polygonBigCore_ageCatSt.r

prop <- list(ageSt_bear$prop[,1,,], ageSt_bear$prop[,2,,], ageSt_bear$prop[,3,,], ageSt_bear$prop[,4,,],
             ageSt_bear$prop[,5,,], ageSt_bear$prop[,6,,], ageSt_bear$prop[,7,,], ageSt_bear$prop[,8,,])

for(n in 1:length(prop)){
  if(all(complete.cases(prop[[n]])) == FALSE){ # delete rows with NA if there are
    delete <- unique(c(which(!complete.cases(prop[[n]][,,1])), which(!complete.cases(prop[[n]][,,2]))))
    prop[[n]] <- prop[[n]][-delete,,]
  }
}
prop_F <- prop[c(7,5,3,1)]
prop_M <- prop[c(8,6,4,2)]
colZn4 <- c("#9970ab", "#a6dba0")


setwd("D:/MargSalas/Oso/OPSCR_project/Results/Results_section/Plots")
pdf("1.3.Prop_AgeSex_ALLyears_corePer.pdf", 6, 3)

par(mfrow = c(1,1), 
    mar = c(3, 3, 3, 6),
    xpd=TRUE)

plot(-0.8, xlim=c(-0.8,1), ylim=c(0,(length(prop_F)+1)),
     type ="n", 
     yaxt="n", 
     xaxt="n", 
     xlab = " ", ylab = "", main = " ",
     cex.axis = 0.8, axes = FALSE)

axis(1, seq(-0.6,0.8, by = 0.2), labels = abs(seq(-0.6,0.8, by = 0.2)), 
     at = seq(-0.6,0.8, by = 0.2), las = 1, cex.axis = 1.2, pos = 0, lwd.ticks = 0.2, lwd = 0.2)


for(i in 1:length(prop_F)){
  
  
  plot.violinsHoriz(dat.list = list(apply(prop_F[[i]], c(1,3), mean)[,1]) ,x = i,at = i+0.2,
                    
                    add = T, violin.width = 0.18, horizontal = T,col="#9970ab", border.col = "#9970ab", alpha = 0.5)
  
  
  plot.violinsHoriz(dat.list =list(-apply(prop_M[[i]], c(1,3), mean)[,1]) ,x = i,at = i+0.2,
                    
                    add = T, violin.width = 0.18, horizontal = T,col="#9970ab",border.col = "#9970ab", alpha = 0.5)
  
  
  plot.violinsHoriz(dat.list =list(apply(prop_F[[i]], c(1,3), mean)[,2]) ,x = i,at = i-0.2,
                    
                    add = T, violin.width = 0.18, horizontal = T,col="#a6dba0",border.col = "#a6dba0", alpha = 0.5)
  
  
  plot.violinsHoriz(dat.list =list(-apply(prop_M[[i]], c(1,3), mean)[,2]) ,x = i,at = i-0.2,
                    
                    add = T, violin.width = 0.18, horizontal = T,col="#a6dba0",border.col = "#a6dba0", alpha = 0.5)
}

segments(x0=0,y0=0,x1=0,y1=5)

text(-0.4,5, "Males")
text(0.4,5, "Females")
mtext(c("All", "Adults", "Subadults", "Cubs"), side = 2, at = c(1,2,3,4), las = 1, line = -2)
#legend("topright", inset=c(-0.2,0), legend = c("Core", "Periphery"), fill = colZn4, border = NA)

dev.off()