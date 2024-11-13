
plot(distcoreMask, main = 1)
points(Xpoints[J.year[[1]], ], pch = 19, cex = 0.5)
# Plot capture locations
caps <- which(Y[1,,1] > 0) ## ASP: Get in which traps the ind was captured at year t
points(Xt[[1]][caps,], pch = 4, col = "red")

plot(distcoreMask, main = 2)
points(Xpoints[J.year[[2]], ], pch = 19, cex = 0.5)
# Plot capture locations
caps <- which(Y[1,,2] > 0) ## ASP: Get in which traps the ind was captured at year t
points(Xt[[2]][caps,], pch = 4, col = "red")

plot(distcoreMask, main = 3)
points(Xpoints[J.year[[3]], ], pch = 19, cex = 0.5)
# Plot capture locations
caps <- which(Y[1,,3] > 0) ## ASP: Get in which traps the ind was captured at year t
points(Xt[[3]][caps,], pch = 4, col = "red")

## INDIVIDUAL 4: Captured in year 1, NOT CAPTURED IN YEARs 2-3 (Where the model calculates a likelihood)

plot(distcoreMask, main = 1)
points(Xpoints[J.year[[1]], ], pch = 19, cex = 0.5)
# Plot capture locations
caps <- which(Y[4,,1] > 0) ## ASP: Get in which traps the ind was captured at year t
points(Xt[[1]][caps,], pch = 4, col = "red")

plot(distcoreMask, main = 2)
points(Xpoints[J.year[[2]], ], pch = 19, cex = 0.5)
# Plot capture locations
caps <- which(Y[4,,2] > 0) ## ASP: Get in which traps the ind was captured at year t
points(Xt[[2]][caps,], pch = 4, col = "red")

plot(distcoreMask, main = 3)
points(Xpoints[J.year[[3]], ], pch = 19, cex = 0.5)
# Plot capture locations
caps <- which(Y[4,,3] > 0) ## ASP: Get in which traps the ind was captured at year t
points(Xt[[3]][caps,], pch = 4, col = "red")


## INDIVIDUAL 3: Captured in year 1 and 3, NOT CAPTURED IN YEAR 2 (Where the model calculates a likelihood)

plot(distcoreMask, main = 1)
points(Xpoints[J.year[[1]], ], pch = 19, cex = 0.5)
# Plot capture locations
caps <- which(Y[3,,1] > 0) ## ASP: Get in which traps the ind was captured at year t
points(Xt[[1]][caps,], pch = 4, col = "red")

plot(distcoreMask, main = 2)
points(Xpoints[J.year[[2]], ], pch = 19, cex = 0.5)
# Plot capture locations
caps <- which(Y[3,,2] > 0) ## ASP: Get in which traps the ind was captured at year t
points(Xt[[2]][caps,], pch = 4, col = "red")

plot(distcoreMask, main = 3)
points(Xpoints[J.year[[3]], ], pch = 19, cex = 0.5)
# Plot capture locations
caps <- which(Y[3,,3] > 0) ## ASP: Get in which traps the ind was captured at year t
points(Xt[[3]][caps,], pch = 4, col = "red")

# CHECK FUNCTION GETLOCALOBJECTS
habitatMask = habitatMask
coords = Xt.sc[[1]]
dmax = 10
resizeFactor = 1

localTraps_test <- getLocalObjects(habitatMask, Xt.sc[[1]], resizeFactor = 1, dmax = 10)
localTraps_test$numLocalIndices

plot(habitatCoords[ ,2] ~ habitatCoords[ ,1], pch = 16, cex = 0.1)
points(coords[ ,2] ~ coords[ ,1], pch = 16, cex = 0.2, col = "red")
points(habitatCoords[766,2] ~ habitatCoords[766,1], col = "green")


