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
source("./nomx_mc.R")


spa_mc <- function(nX, nK, Y, mx, para){
  
  r <- para$rDim
  
  nSam <- dim(nX)[2]
  xDim <- dim(X)[1]
  xSam <- dim(X)[2]
  ySam <- dim(Y)[2]
  
  # X <- center_mc(X)
  Y <- center_mc(Y)
  
  ###### update G ######
  svM <- mx$svM
  # u <- svM$u
  u <- svM$u
  v <- svM$v
  d <- svM$d
  s <- diag(d)
  
  # tmp <- v[1:nSam, ]
  # G <- tmp %*% s
  G <- v %*% s
  rm(svM)
  
  ###### update svM ######
  I <- diag(1, xDim, xDim)
  tmp <- u %*% t(u)
  tmp <- (I - tmp) %*% Y
  res <- qr(tmp)
  Q <- qr.Q(res)
  R <- qr.R(res)
  
  tmp <- crossprod(u, Y)
  tmp <- as.matrix(cbind(s, tmp))
  Z <- matrix(0, nrow = ySam, ncol = r)
  ZR <- as.matrix(cbind(Z, R))
  F <- rbind(tmp, ZR)
  
  res <- svd(F)
  uu <- res$u
  vv <- res$v
  d <- res$d
  tmp <- cbind(u, Q)
  u <- tmp %*% uu
  Z <- matrix(0, nrow = xSam, ncol = ySam)     # 958 * 30???
  tmp <- cbind(v, Z)     # 
  Z <- matrix(0, nrow = ySam, ncol = r)
  I <- diag(1, ySam, ySam)
  OI <- cbind(Z, I)
  tmp <- rbind(tmp, OI)    
  v <- tmp %*% vv
  
  u <- u[, 1:r]
  v <- v[, 1:r]
  d <- d[1:r]
  svM <- list(u=u, v=v, d=d)
  rm(tmp, I, res, Q, R, Z, ZR, F, uu, vv, OI, u, v, d)     #####
  
  ##### calcualte K #####
  Kx <- mx$Kx
  K <- contk_mc(X, Y, "gaussian", para)
  Ky <- contk_mc(Y, Y, "gaussian", para)
  
  tK <- cbind(Kx, K)
  tmp <- cbind(t(K), Ky)
  kk <- rbind(tK, tmp)
  rm(tmp, Kx)
  
  ##### update comp of K #####
  svK <- mx$svK
  kv <- svK$V
  kd <- svK$D
  V <- kv
  D <- kd
  rm(svK)
  
  s <- D^-1
  s <- diag(s)
  Qxy <- s %*% t(V)
  Qxy <- Qxy %*% K
  tmp <- crossprod(Qxy, Qxy)
  tmp <- Ky - tmp
  res <- svd(tmp)
  u <- res$u
  d <- res$d
  d <- d^.5
  s <- diag(d)
  Qy <- u %*% s
  
  QD <- rbind(Qxy, Qy)
  # I <- diag(1, xSam+ySam, xSam+ySam)
  # I <- diag(1, xSam, xSam)
  I <- diag(1, ySam, ySam)
  Z <- matrix(0, ySam, xSam)
  OI <- cbind(Z, I)
  Z <- matrix(0, xSam, xSam + ySam)
  tmp <- rbind(Z, OI)
  # tmp <- I - QD
  tmp <- tmp %*% QD
  
  res <- qr(tmp)
  Q <- qr.Q(res)
  R <- qr.R(res)
  
  I <- diag(1, xSam, xSam)
  Z <- matrix(0, nrow = xSam, ncol = ySam)
  tmp <- cbind(I, Z)
  tmp <- tmp %*% QD
  rm(Qxy, QD, I, OI, Z, res, Qy)     #####
  
  s <- diag(D)
  AI <- cbind(s, tmp)
  Z <- matrix(0, nrow = ySam, ncol = xSam)
  ZR <- cbind(Z, R)
  tmp <- rbind(AI, ZR)
  res <- svd(tmp)
  uu <- res$u
  vv <- res$v
  d <- res$d
  
  # U <- Q %*% uu
  I <- diag(1, xSam+ySam, xSam)
  QU <- cbind(I, Q)
  U <- QU %*% uu
  
  Z <- matrix(0, xSam, ySam)
  tmp <- cbind(V, Z)
  Z <- matrix(0, nrow = ySam, ncol = xSam)
  I <- diag(1, ySam, ySam)
  OI <- cbind(Z, I)
  tmp <- rbind(tmp, OI)
  V <- tmp %*% vv
  D <- d
  
  svK <- list(Kx = kk, U = U, V = V, D = D)
  rm(AI, Z, ZR, OI, Q, R, uu, vv, d, U, V, D)     #####
  
  ##### calculate align #####
  kd <- kd^-2
  s <- diag(kd)
  T <- kv %*% s %*% t(kv)
  # T <- s %*% t(kv)
  G <- T %*% G
  
  # tmp <- Kx %*% M %*% K
  
  # tmx <- xx %*% t(xx)
  # tmx <- ginv(tmx) %*% xx
  
  # tmx <- t(xx)
  # tmx <- ginv(tmx)
  # yy <- tmx %*% tmp
  # 
  # return(yy)
  
  # res <- svd(M)
  # u <- res$u
  # v <- res$v
  # d <- res$d
  # # d <- d[1:r]
  # d <- d^.5
  # s <- diag(d)
  #
  # G <- u %*% s
  # G <- M
  
  xx <- t(G) %*% nK
  yy <- t(G) %*% K
  rm(G, T, kd, kv, s)
  ##### update K #####
  K <- contk_mc(Y, nX, "gaussian", para)
  nK <- rbind(nK, K)
  # nK <- t(nK)
  
  xx <- nomx_mc(xx)
  yy <- nomx_mc(yy)
  
  res <- list(xx=xx, yy=yy, nK = nK, Kx = kk, svM=svM, svK=svK)
  
  
  return(res)
  
  # Cxy <- tmp
  # return(res)
  
}


