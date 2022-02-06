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
source("./lda_mc.R")



getmxy_mc <- function(X, Y, para){
  
  r <- para$rDim
  lambda <- para$lambda
  xSam <- dim(X)[2]
  ySam <- dim(Y)[2]
  
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  # output <- lda_mc(X, xL, r)
  # xU <- output$V
  
  xS <- svd(X, nu = r, nv = r)
  xU <- xS$u
  
  yS <- svd(Y, nu = r, nv = r)
  yU <- yS$u
  
  # res <- getXP_mc(X, para)
  # P <- res$P
  
  # P <- xU %*% t(xU) %*% yU
  P <- xU %*% t(xU)
  px <- crossprod(P, X)
  # R <- knnsolo_mc(xx, xL, yy, yL, k = 3)
  # ac <- R$accu
  
  # ux <- normx_mc(ux)
  # 
  # xx <- crossprod(ux, ux)
  
  Q <- yU %*% t(yU)
  py <- crossprod(Q, Y)
  pxy <- crossprod(px, py)
  
  tmp <- colSums(X * X)
  sigma <- mean(tmp) *0.5
  para$sigma <- sigma
  
  # para$sigma <- 1
  Kx <- contk_mc(X, X, "gaussian", para)
  Kxy <- contk_mc(X, Y, "gaussian", para)
  pKx <- ginv(Kx)
  pKxy <- ginv(Kxy)
  tKxy <- t(Kxy)
  # ptKxy <- ginv(tKxy)
  
  # tmp <- Kxy %*% ptKxy
  # tmp <- pxy + lambda * tmp
  # M <- pKx %*% tmp %*% pKxy
  
  # Ix <- diag(1, ySam, ySam)
  # Z <- array(0, c(xSam-ySam, ySam))
  # I <- rbind(Ix, Z)
  # tmp <- pxy + lambda * I
  
  tmp <- pxy - lambda * t(pKxy)
  M <- pKx %*% tmp %*% pKxy
  
  # M <- tmp %% tt %*% tmp
  # M <- tmp %% (xx) %*% tmp
  
  # M <- tmp %*% t(px)
  # M <- M %*% t(M)
  
  res <- list(M = M, ux = px, u = P, Kx = Kx, para = para)
  
  
  return(res)
  
  
}


