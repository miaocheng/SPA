# 
# This function calculates the aligned matrix of original data set associated with
# localities preservation.
# Input: X - nDim * nSam
#        r - reduced dimensionality
#        lambda - regularization parameter
# Output: M - nDim * r
# Coded by Miao Cheng
# Data: 2016-4-18
#############################################################################################

library(MASS)

source("./contk_mc.R")
source("./getXP_mc.R")



getmxl_mc <- function(X, para){
  
  r <- para$rDim
  lambda <- para$lambda
  nSam <- dim(X)[2]
  
  X <- as.matrix(X)
  # xS <- svd(X, nu = r, nv = r)
  # u <- xS$u
  res <- getXP_mc(X, para)
  P <- res$P
  
  
  px <- crossprod(P, X)
  # ux <- normx_mc(ux)
  # 
  # xx <- crossprod(ux, ux)
  
  tmp <- colSums(X * X)
  sigma <- mean(tmp)
  para$sigma <- sigma
  
  # para$sigma <- 1
  Kx <- contk_mc(X, X, "gaussian", para)
  
  I <- diag(array(1, c(nSam, 1)))
  tmp <- Kx + lambda * I
  tmp <- ginv(tmp)
  # tt <- xx - lambda * tmp
  
  # M <- tmp %% tt %*% tmp
  # M <- tmp %% (xx) %*% tmp
  
  M <- tmp %*% t(px)
  M <- M %*% t(M)
  
  res <- list(M = M, ux = px, u = P, Kx = Kx, para = para)
  
  
  return(res)
  
  
}


