# 
# This function implements the kasa algorithm.
# Input: X - xDim * xSam
#        xL - xSam * 1
#        Y - yDim * ySam
#        yL - ySam * 1
# Output: P - xDim * rDim
#         Q - yDim * rDim
# Coded by Miao Cheng
# Data: 2016-4-4
###########################################################################################

library(MASS)


source("./contk_mc.R")
source("./getmx_mc.R")


kasa_mc <- function(X, Y, para){
  
  lambda <- para$lambda
  r <- para$rDim
  
  xSam <- dim(X)[2]
  ySam <- dim(Y)[2]
  
  res <- getmx_mc(X, r, lambda)
  
  M <- res$M
  xx <- res$ux
  Kx <- res$Kx
  
  # Kx <- contk_mc(xx, xx, "cosine")
  K <- contk_mc(X, Y, "gaussian", para)
  
  # xx <- normx_mc(ux)
  tmp <- Kx %*% M %*% K
  tmx <- xx %*% t(xx)
  tmx <- ginv(tmx) %*% xx
  yy <- tmx %*% tmp
  
  
  # ut <- svd(tmp)
  # U <- ut$u
  # U <- t(U)
  # S <- ut$d
  # ss <- diag(S)
  # # ss <- ss[1:r]
  # ss <- ss^.5
  # D <- diag(ss)
  # 
  # yy <- U %*% D
  
  # yy <- normx_mc(yy)
  res <- list(xx=xx, yy=yy)
  return(res)
  
  
}




