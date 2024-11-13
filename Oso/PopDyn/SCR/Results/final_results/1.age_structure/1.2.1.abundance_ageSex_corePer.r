

## ---- Absolute number/Ageclass/year/sex ----

rm(list = ls())

library(nimbleSCR)
library(nimble)
library(rgdal)

source("D:/MargSalas/Scripts_MS/Functions/plot.violins3.r")
source("D:/MargSalas/Scripts_MS/Functions/DoScale.r")
source("D:/MargSalas/Scripts_MS/Functions/PlotViolinsHoriz.r")

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Models/3.openSCRdenscov_Age/2021/Cyril/3-3.1_allparams")
load("Nbuffer_newSize.RData") # Obtained from script: 3.3. Core_vS_Periphery_polygonBigCore_ageCatSt.r

colZn4 <- c("#9970ab", "#a6dba0")

setwd("D:/MargSalas/Oso/OPSCR_project/Results/Results_section/Plots")
pdf("1.2.Abundance_AgeSex_CorePer.pdf", 5, 5)

par(mfrow = c(4,2),
    mar = c(1,1,1.3,1),
    oma = c(2,3,2,2))

ageSt_list <- ageSt_bear[-c(1,2)]
nAgeSex <- length(ageSt_list)

zones <- dim(ageSt_list[[1]])[3]
offset <- seq(0,0.25,length.out = zones)
at <- c(1,2,3,4,5)

leg <- c(FALSE, TRUE, rep(FALSE, 6))
tit1 <- c("Females", "Males", rep(FALSE, 6))
tit2 <- c("Cubs", FALSE, "Subadults", FALSE, "Adults", FALSE, "All", FALSE)


scaleY_plots <- c(rep(max(ageSt_list[[1]], ageSt_list[[2]]),2), rep(max(ageSt_list[[3]], ageSt_list[[4]], na.rm = TRUE),2),
                  rep(max(ageSt_list[[5]], ageSt_list[[6]]),2), rep(max(ageSt_list[[7]], ageSt_list[[8]]),2))

for (n in 1:nAgeSex){
  
  if(all(complete.cases(ageSt_list[[n]])) == FALSE){ # delete rows with NA if there are
    delete <- unique(c(which(!complete.cases(ageSt_list[[n]][,,1])), which(!complete.cases(ageSt_list[[n]][,,2]))))
    ageSt_list[[n]] <- ageSt_list[[n]][-delete,,]
  }
  
  
  plot(1, ylim = c(0,scaleY_plots[n]+1), 
       xlim = c(0.5, at[5] + max(offset) + 0.5)  , 
       type ="n", 
       yaxt="n", 
       xaxt="n", 
       xlab = " ", ylab = "", main = " ",
       cex.axis = 0.8, axes = FALSE)
  
  axis(1, c(1:ncol(ageSt_list[[n]])), labels = c(2017:2021), 
       at = at + max(offset)/2, las = 1, cex.axis = 0.75, pos = 0, lwd.ticks = 0.2, lwd = 0.2 )
  axis(2, cex.axis = 0.75, pos = 0.75, lwd.ticks = 0.2, lwd = 0.2)
  
  
  m <- data.frame(matrix(NA, nrow = 5, ncol = 2))     
  for(s in 1:zones){
    for (i in 1:5){
      plot.violins3(list(ageSt_list[[n]][ ,i,s]),
                    x = i,
                    at = at[i]+offset[s],
                    violin.width = 0.05,
                    plot.ci = 0.85,
                    col = colZn4[s],
                    add = T,
                    alpha = 0.8,
                    scale.width = FALSE,
                    border.col = colZn4[s],
                    horizontal = FALSE)
      m[i,s] <- median(ageSt_list[[n]][ ,i,s]) 
    }
  }
  
  for(s in 1:zones){
    lo <- loess(m[,s]~ c(at+offset[s]))
    lines(predict(lo)~ c(at+offset[s]), col=colZn4[s], lwd=1.5)
  }
  #if(leg[n] == TRUE){
  #  legend("topright", inset = c(+0.1, -0.05), legend = c("Core", "Periphery"), fill = colZn4, border = NA, bty = "n", horiz = FALSE, cex = 0.9)
  #}
  if(tit1[n] != FALSE){
    mtext(tit1[n], side = 3, cex = 0.8, line = 1)
  }
  if(tit2[n] != FALSE){
    mtext(tit2[n], side = 2, cex = 0.75, line = 1.5)
  }
}

dev.off()

