

library('R.matlab')

# setwd('/home/cm/Documents/Run/kasa/IR')
setwd('M:/Run/kasa/IR')


source("./center_mc.R")
source("./knnsolo_mc.R")
source("./paradat_mc.R")
source("./getmx_mc.R")
# source("./spal_mc.R")
# source("./refv_mc.R")
source("./getmxl_mc.R")
source("./paraxy_mc.R")
source("./paray_mc.R")
source("./getmxy_mc.R")
source("./spa_mc.R")
source("./spak_mc.R")



#################################################################################
xDat <- readMat("../data/amazon_SURF_L10.mat")
yDat <- readMat("../data/Caltech10_SURF_L10.mat")

X <- xDat$fts
Y <- yDat$fts
xL <- xDat$labels
yL <- yDat$labels
X <- t(X)
Y <- t(Y)


# pdat <- paraxy_mc(X, Y, xL, yL, xSam+300, 100)
pdat <- paray_mc(X, Y, xL, yL, 0, 250)
X <- pdat$X
xL <- pdat$xL
nY <- pdat$nY
nL <- pdat$nL
Y <- pdat$Y
yL <- pdat$yL
l <- pdat$l

#################################################################################
# xDat <- readMat("../data/caltranGist.mat")
# 
# pdat <- paradat_mc(xDat, 500, 300)
# X <- pdat$X
# xL <- pdat$xL
# Y <- pdat$Y
# yL <- pdat$yL
# l <- pdat$l
#################################################################################

xSam <- dim(X)[2]
ySam <- dim(Y)[2]


para <- list(rDim=30, lambda=0.5, k = 3)
X <- center_mc(X)
Y <- center_mc(Y)

mx <- getmx_mc(X, para)
# mx <- getmxy_mc(X, X, para)
# mx <- getmxl_mc(X, para)

para <- mx$para
# xx <- mx$ux
# xU <- mx$u
# M <- mx$M
# Kx <- mx$Kx


t <- length(l)
alc <- 0
ind <- 0

nX <- X
nK <- mx$Kx
for (i in 1:t)
{
  
  tnd <- ind + 1
  n <- l[i]
  tY <- Y[, tnd:(ind+n)]
  tL <- yL[tnd:(ind+n)]
  
  # tY <- center_mc(tY)
  # yy <- spal_mc(X, tY, mx, para)
  
  # yy <- refv_mc(yy, tY, xU)
  
#   mx <- spa_mc(nX, nK, tY, mx, para)
  mx <- spak_mc(nX, nK, tY, mx, X, para)
  xx <- mx$xx
  yy <- mx$yy
  
  # Cxy <- spa_mc(X, tY, mx, para)
  
  R <- knnsolo_mc(xx, xL, yy, tL, k = 1)
  ac <- R$accu
  alc <- alc + ac
  
  ind <- ind + l[i]
  nX <- cbind(nX, tY)
  nK <- mx$nK
  
  print(i)
  rm(xx, yy, tY, tL)
  
}

accuracy <- alc / ySam




