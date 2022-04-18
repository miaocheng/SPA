# 
# This function implements the sequantial pattern alignment with given ranking.
# Input: X - xDim * xSam
#        tY - temporal input
#        para - input parameters
# Output: yy - obtained results
# Coded by Miao Cheng
# Data: 2016-4-21
###########################################################################################

# library(MASS)
library(rCUR)

source("./contk_mc.R")
source("./nomx_mc.R")


spar_mc <- function(nX, nK, Y, mx, oX, para){
  
  r <- para$rDim
  
  #   nSam <- dim(nX)[2]
  xDim <- dim(nX)[1]
  xSam <- dim(nX)[2]
  ySam <- dim(Y)[2]
  
  # X <- center_mc(X)
  # Y <- center_mc(Y)
  
  ###### Select the ideal matching data ######
  K <- contk_mc(nX, Y, "gaussian", para)
  K <- as.matrix(K)
  res <- CUR(K, c=ySam*0.6, r=xSam, k=r, method = "ortho.top.scores")
  #   rC <- getC(res)
  #   rR <- getR(res)
  #   ind <- rownames(rC)
  ind <- res@C.index
  ind <- as.numeric(ind)
  #   selX <- subset(X, select=ind)
  sY <- Y[, ind]
  sSam <- dim(sY)[2]
  
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
  
  tmp <- (I - tmp) %*% sY
  res <- qr(tmp)
  Q <- qr.Q(res)
  R <- qr.R(res)
  
  tmp <- crossprod(u, sY)
  tmp <- as.matrix(cbind(s, tmp))
  Z <- matrix(0, nrow = sSam, ncol = r)
  ZR <- as.matrix(cbind(Z, R))
  F <- rbind(tmp, ZR)
  
  res <- svd(F)
  uu <- res$u
  vv <- res$v
  d <- res$d
  tmp <- cbind(u, Q)
  u <- tmp %*% uu
  Z <- matrix(0, nrow = xSam, ncol = sSam)     # 958 * 30???
  tmp <- cbind(v, Z)     # 
  Z <- matrix(0, nrow = sSam, ncol = r)
  I <- diag(1, sSam, sSam)
  OI <- cbind(Z, I)
  tmp <- rbind(tmp, OI)    
  v <- tmp %*% vv
  
  u <- u[, 1:r]
  v <- v[, 1:r]
  d <- d[1:r]
  svM <- list(u=u, v=v, d=d)
  rm(tmp, I, res, Q, R, Z, ZR, F, uu, vv, OI, u, v, d)     #####
  
  ###### calcualte K ######
#   Kx <- mx$Kx
  K <- contk_mc(nX, sY, "gaussian", para)
  Ky <- contk_mc(sY, sY, "gaussian", para)
  
#   tK <- cbind(Kx, K)
#   tmp <- cbind(t(K), Ky)
#   kk <- rbind(tK, tmp)
#   rm(tmp, Kx)
  
  ###### update comp of K ######
  svK <- mx$svK
  kv <- svK$V
  kd <- svK$D
  V <- kv
  D <- kd
  rm(svK)
  
  # V <- V[, 1:r]
  # D <- D[1:r]
  
  s <- D^-1
  s <- diag(s)
  Qxy <- s %*% t(V)
  Qxy <- Qxy %*% K    # K !
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
  len <- len - r
  I <- diag(1, r, r)
  Z <- matrix(0, len, r)
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
  Z <- matrix(0, nrow = sSam, ncol = r)
  ZR <- cbind(Z, R)
  tmp <- rbind(AI, ZR)
  res <- svd(tmp, r, r)
  uu <- res$u
  vv <- res$v
  d <- res$d
  d <- d[1:r]
  
  # U <- Q %*% uu
  I <- diag(1, r, r)
  Z <- matrix(0, sSam, r)
  I <- rbind(I, Z)
  IQ <- cbind(I, Q)
  U <- IQ %*% uu
  U <- U[, 1:r]
  
  Z <- matrix(0, xSam, sSam)
  tmp <- cbind(V, Z)
  Z <- matrix(0, nrow = sSam, ncol = r)
  I <- diag(1, sSam, sSam)
  OI <- cbind(Z, I)
  tmp <- rbind(tmp, OI)
  
  V <- tmp %*% vv
  V <- V[, 1:r]
  D <- d[1:r]
  
  svK <- list(U = U, V = V, D = D)
#   svK <- list(Kx = kk, U = U, V = V, D = D)
  #   rm(AI, Z, ZR, OI, Q, R, uu, vv, d, QU, U, V, D)     #####
  rm(AI, Z, ZR, OI, Q, R, uu, vv, d, U, V, D)
  
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
  ###### update K ######
  K <- contk_mc(sY, oX, "gaussian", para)
  nK <- rbind(nK, K)
  # nK <- t(nK)
  
  # xx <- nomx_mc(xx)
  # yy <- nomx_mc(yy)
  
  #   res <- list(xx=xx, yy=yy, nK = nK, Kx = kk, svM=svM, svK=svK)
  res <- list(xx=xx, yy=yy, sY=sY, nK = nK, svM=svM, svK=svK)
  
  
  return(res)
  
  # Cxy <- tmp
  # return(res)
  
}


