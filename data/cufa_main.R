

library('R.matlab')

setwd('/home/cm/Documents/Run/cufa/R')


#source("./center_mc.R")
source("./spa_mc.R")
source("./detcr_mc.R")
source("./knnsolo_mc.R")
##########################################################


xDat <- readMat("../data/amazon_SURF_L10.mat")
yDat <- readMat("../data/Caltech10_SURF_L10.mat")

X <- xDat$fts
Y <- yDat$fts
xL <- xDat$labels
yL <- yDat$labels
X <- t(X)
Y <- t(Y)


para <- list(rDim=50, k=5, lambda=0.00001, sup=FALSE)
# para <- list(rDim=30, k=5, lambda=0.00001, sup=TRUE, n = 50)

pq <- spa_mc(X, Y, para)


P <- pq$p
Q <- pq$q


X <- center_mc(X)
Y <- center_mc(Y)

#sx <- concat_mc(X)
#sy <- concat_mc(Y)

xx <- crossprod(P, X)
yy <- crossprod(Q, Y)

R <- knnsolo_mc(xx, xL, yy, yL, k = 1)


 

