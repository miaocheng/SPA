# 
# This function calculates the refined data.
# 
# Coded by Miao Cheng
# Data: 2016-4-6
################################################################################



refv_mc <- function(Y, tY, xU){
  
  yDim <- dim(Y)[1]
  ySam <- dim(Y)[2]
  
  tY <- as.matrix(tY)
  # xU <- as.matrix(xU)
  
  P <- Y %*% ginv(tY)
  P <- t(P)
  
  yy <- P %*% t(P) %*% tY
  
  # yS <- svd(Y, nu = min(yDim, ySam), nv = min(yDim, ySam))
  # u <- yS$u
  # yy <- u %*% t(u) %*% Y
  
  yy <- t(xU) %*% yy
  
  
  return(yy)
  
}

