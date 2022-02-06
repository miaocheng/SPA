# 
# This function implements the sequantial pattern alignment.
# Input: X - xDim * xSam
#        tY - temporal input
#        para - input parameters
# Output: yy - obtained results
# Coded by Miao Cheng
# Data: 2016-4-21
#######################################################################################

# library(MASS)

source("./contk_mc.R")


spa_mc <- function(X, Y, mx, para){
  
  r <- para$rDim
  
  xSam <- dim(X)[2]
  ySam <- dim(Y)[2]
  
  X <- center_mc(X)
  Y <- center_mc(Y)
  
  M <- mx$M
  xx <- mx$ux
  Kx <- mx$Kx
  
  K <- contk_mc(X, Y, "gaussian", para)
  
  # tmp <- Kx %*% M %*% K
  
  # tmx <- xx %*% t(xx)
  # tmx <- ginv(tmx) %*% xx
  
  # tmx <- t(xx)
  # tmx <- ginv(tmx)
  # yy <- tmx %*% tmp
  # 
  # return(yy)
  
  res <- svd(M)
  u <- res$u
  v <- res$v
  d <- res$d
  # d <- d[1:r]
  d <- d^.5
  s <- diag(d)
  
  G <- u %*% s
  xx <- t(G) %*% Kx
  yy <- t(G) %*% K
  res <- list(xx=xx, yy=yy)
  return(res)
  
  # Cxy <- tmp
  # return(res)
  
}


