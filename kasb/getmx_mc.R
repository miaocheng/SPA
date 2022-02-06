# 
# This function calculates the aligned matrix of original data set.
# Input: X - nDim * nSam
#        r - reduced dimensionality
#        lambda - regularization parameter
# Output: M - nDim * r
# Coded by Miao Cheng
# Data: 2016-4-5
#############################################################################################

library(MASS)

source("./contk_mc.R")


getmx_mc <- function(X, para){
  
  r <- para$rDim
  lambda <- para$lambda
  nSam <- dim(X)[2]
  
  X <- as.matrix(X)
  svX <- svd(X, nu = r, nv = r)
  P <- svX$u
  
  px <- crossprod(P, X)
  # ux <- normx_mc(ux)
  # 
  # xx <- crossprod(ux, ux)
  
  tmp <- colSums(X * X)
  sigma <- mean(tmp) * 0.5
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
  # M <- Kx %*% M
  
  # svM <- svd(M)
  # # A <- M %*% t(M)
  # v <- svM$u
  # u <- svM$v
  # d <- svM$d
  # u <- P %*% u
  # svM <- list(u=u, v=v, d=d)
  
  svM <- svX
  u <- svM$u
  v <- svM$v
  d <- svM$d
  d <- d[1:r]
  svM <- list(u=u, v=v, d=d)
  
  svK <- svd(Kx)
  V <- svK$u
  d <- svK$d
  d <- d^.5
  D <- d
  
  V <- V[, 1:r]
  D <- D[1:r]
  
  svK <- list(V = V, D = D)
  
  # res <- list(svM = svM, M = M, ux = ux, u = u, Kx = Kx, para = para)
  res <- list(Kx = Kx, svM = svM, svK = svK, para = para)
  
  return(res)
  
  
}




