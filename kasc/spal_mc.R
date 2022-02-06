# 
# This function implements the sequantial pattern alignment.
# Input: X - xDim * xSam
#        tY - temporal input
#        para - input parameters
# Output: yy - obtained results
# Coded by Miao Cheng
# Date: 2016-4-6
#######################################################################################

library(base)
library(rCUR)

source("./contk_mc.R")


#spal_mc <- function(X, Y, mx, para){
spal_mc <- function(nX, nK, Y, mx, para){
  
  
  rDim <-para$rDim
  oSam <- para$oSam
  
#   xSam <- dim(nX)[2]
  xDim <- dim(nX)[1]
  xSam <- dim(nX)[2]
  ySam <- dim(Y)[2]
#   n <- ySam * ratio
  
#   M <- mx$M
#   xx <- mx$ux
#   Kx <- mx$Kx
  
  
  ###### Select the ideal matching data ######
  K <- contk_mc(nX, Y, "gaussian", para)
  K <- as.matrix(K)
  res <- CUR(K, c=ySam, r=oSam, k=rDim, method = "ortho.top.scores")
#   rC <- getC(res)
#   rR <- getR(res)
#   ind <- rownames(rC)
  ind <- res@R.index
  ind <- as.numeric(ind)
#   selX <- subset(X, select=ind)
  sX <- nX[, ind]
  
  
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
  
  ov <- svM$ov
  oG <- ov %*% s
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
  Z <- matrix(0, nrow = ySam, ncol = rDim)
  ZR <- as.matrix(cbind(Z, R))
  F <- rbind(tmp, ZR)
  
  res <- svd(F)
  uu <- res$u
  vv <- res$v
  d <- res$d
  tmp <- cbind(u, Q)
  u <- tmp %*% uu
  
  Z <- matrix(0, nrow = xSam, ncol = ySam)     # 958 * 30???
  tmp <- cbind(v, Z)
  Z <- matrix(0, nrow = ySam, ncol = rDim)
  I <- diag(1, ySam, ySam)
  OI <- cbind(Z, I)
  tmp <- rbind(tmp, OI)    
  
  ###### ***** ######
  oZ <- matrix(0, nrow = oSam, ncol = ySam)
  tp <- cbind(ov, oZ)    #???
  ov <- tp %*% vv
  ###### ***** ######
  
  v <- tmp %*% vv
  u <- u[, 1:rDim]
  v <- v[, 1:rDim]
  d <- d[1:rDim]
  ov<- ov[, 1:rDim]
  svM <- list(u=u, v=v, d=d, ov=ov)
  
  rm(tmp, tp, I, res, Q, R, Z, ZR, F, uu, vv, OI, u, v, d)     #####
#   rm(tmp, I, res, Q, R, Z, ZR, F, uu, vv, u, v, d)
  
  ###### calcualte K ######
  Kx <- mx$Kx
#   K <- contk_mc(X, Y, "gaussian", para)
  Ky <- contk_mc(Y, Y, "gaussian", para)
  
  tK <- cbind(Kx, K)
  tmp <- cbind(t(K), Ky)
  kk <- rbind(tK, tmp)
  rm(tmp, Kx)
  
  ###### update comp of K ######
  svK <- mx$svK
  kv <- svK$V
  kd <- svK$D
  V <- kv
  D <- kd
  okv <- svK$okv
  rm(svK)
  
  ###### ***** ######
  sK <- K[ind, ]
  # V <- V[, 1:r]
  # D <- D[1:r]
  s <- D^-1
  s <- diag(s)
  Qxy <- s %*% t(okv)
  Qxy <- Qxy %*% sK     # ???
  
  tmp <- crossprod(Qxy, Qxy)
  tmp <- Ky - tmp
  res <- svd(tmp)
  u <- res$u
  d <- res$d
  d <- d^.5
  s <- diag(d)
  Qy <- u %*% s
  
  QD <- rbind(Qxy, Qy)
  len <- dim(QD)[1]
  II <- diag(1, len, len)
  len <- len - rDim
  I <- diag(1, rDim, rDim)
  Z <- matrix(0, len, rDim)
  I <- rbind(I, Z)
  tmp <- I %*% t(I)
  tmp <- II - tmp
  IQ <- tmp %*% QD
  
  res <- qr(IQ)
  Q <- qr.Q(res)
  R <- qr.R(res)
  
  tmp <- t(I) %*% QD
  rm(Qxy, Qy, QD, IQ, I, II, Z, res)     #####
  
  s <- diag(D)
  AI <- cbind(s, tmp)
  Z <- matrix(0, nrow = ySam, ncol = rDim)
  ZR <- cbind(Z, R)
  tmp <- rbind(AI, ZR)
  res <- svd(tmp, rDim, rDim)
  uu <- res$u
  vv <- res$v
  d <- res$d
  d <- d[1:rDim]
  
  # U <- Q %*% uu
  I <- diag(1, rDim, rDim)
  Z <- matrix(0, ySam, rDim)
  I <- rbind(I, Z)
  IQ <- cbind(I, Q)
  U <- IQ %*% uu
  U <- U[, 1:rDim]
  
  Z <- matrix(0, oSam, ySam)
  tmp <- cbind(okv, Z)
  okv <- tmp %*% vv
  
  rm(AI, Z, ZR, Q, R, uu, vv, d, U)
  ###### ***** ######
  
  # V <- V[, 1:r]
  # D <- D[1:r]
  
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
  len <- dim(QD)[1]
  II <- diag(1, len, len)
  len <- len - rDim
  I <- diag(1, rDim, rDim)
  Z <- matrix(0, len, rDim)
  I <- rbind(I, Z)
  tmp <- I %*% t(I)
  tmp <- II - tmp
  IQ <- tmp %*% QD
  
  # I <- diag(1, ySam, ySam)
  # Z <- matrix(0, ySam, xSam)
  # OI <- cbind(Z, I)
  # Z <- matrix(0, xSam, xSam + ySam)
  # tmp <- rbind(Z, OI)
  # # tmp <- I - QD
  # IQ <- tmp %*% QD
  
  res <- qr(IQ)
  Q <- qr.Q(res)
  R <- qr.R(res)
  
  # I <- diag(1, xSam, xSam)
  # Z <- matrix(0, nrow = xSam, ncol = ySam)
  # tmp <- cbind(I, Z)
  tmp <- t(I) %*% QD
  rm(Qxy, Qy, QD, IQ, I, II, Z, res)     #####
  
  s <- diag(D)
  AI <- cbind(s, tmp)
  Z <- matrix(0, nrow = ySam, ncol = rDim)
  ZR <- cbind(Z, R)
  tmp <- rbind(AI, ZR)
  res <- svd(tmp, rDim, rDim)
  uu <- res$u
  vv <- res$v
  d <- res$d
  d <- d[1:rDim]
  
  # U <- Q %*% uu
  I <- diag(1, rDim, rDim)
  Z <- matrix(0, ySam, rDim)
  I <- rbind(I, Z)
  IQ <- cbind(I, Q)
  U <- IQ %*% uu
  U <- U[, 1:rDim]
  
  Z <- matrix(0, xSam, ySam)
  tmp <- cbind(V, Z)
  Z <- matrix(0, nrow = ySam, ncol = rDim)
  I <- diag(1, ySam, ySam)
  OI <- cbind(Z, I)
  tmp <- rbind(tmp, OI)
  
  V <- tmp %*% vv
  V <- V[, 1:rDim]
  D <- d[1:rDim]
  
  svK <- list(Kx = kk, U = U, V = V, D = D, okv=okv)
  rm(AI, Z, ZR, OI, Q, R, uu, vv, d, U, V, D)     #####
#   rm(AI, Z, ZR, Q, R, uu, vv, d, U, V, D)     #####
  
  
  ###### calculate alignment ######
  kd <- kd^-2
  s <- diag(kd)
  T <- kv %*% s %*% t(kv)
  # T <- s %*% t(kv)
  G <- T %*% G
  xx <- t(G) %*% nK
  
  oT <- okv %*% s %*% t(okv)
  oG <- oT %*% oG
  yy <- t(oG) %*% sK
  rm(G, T, kd, kv, s)
  ###### update K ######
  K <- contk_mc(Y, nX, "gaussian", para)
  nK <- rbind(nK, K)
  # nK <- t(nK)
  
  
  res <- list(xx=xx, yy=yy, nK = nK, Kx = kk, svM=svM, svK=svK)
  
  
  return(res)
  
  
}


