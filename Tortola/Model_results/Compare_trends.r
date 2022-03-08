
rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)

setwd("D:/Otros/Tórtola/Results/Study2/Model_results")

t0 <- read.csv("trend0.csv")
t1 <- read.csv("trend1.csv")
t1.1 <- read.csv("trend1.1.csv")

t0$model <- "model0"
t1$model <- "model1"
t1.1$model <- "model1.1"

all <- list(t0,t1,t1.1)

for (i in 1:3){
  colnames(all[[i]])[3] <- unique(all[[i]]$model)
  all[[i]] <- all[[i]][,c(2,3)]
} 


all1 <- full_join(all[[1]],all[[2]])
all2 <- full_join(all1,all[[3]])

all2 <- t(all2)
colnames(all2) <- all2[1,]
all2 <- as.data.frame(all2[-1,])

cols <- colnames(all2)
all2_long <- gather(all2, key = "Site", value = "Beta")
all2_long$Model <- rep(rownames(all2),46)

setwd("D:/Otros/Tórtola/Results/Study2/Plots")
pdf("plot_model_transects.pdf",20,5)
    ggplot(all2_long, aes(fill=Model, y=Beta, x=Site)) + 
      geom_bar(position="dodge", stat="identity")
dev.off()


