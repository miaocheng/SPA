# 
# This function implements the projection of data matrix X.
# Input: X - nDim * nSam
#        para - parameters
# Output: P - Projection matrix
# Coded by Miao Cheng
# Data: 2016-4-18
#################################################################################

source("./eudist_mc.R")


getXP_mc <- function(X, para){
  
  r <- para$rDim
  k <- para$k
  
  nDim <- dim(X)[1]
  nSam <- dim(X)[2]
  
  D <- eudist_mc(X, X)
  for (ii in 1:nSam)
  {
    D[ii, ii] = 1e10
  }
  
  
  wets <- NULL
  W <- NULL
  for (i in 1:nSam)
  {
    tmp <- D[i, ]
    ind <- order(tmp)
    ind <- ind[1:k]
    
    tmp <- array(0, c(k, 3))
    tmp[, 1] <- i
    tmp[, 2] <- ind
    wets <- rbind(wets, tmp)
    
    tmp <- array(0, c(nSam, k))
    tmp[i, ] <- 1
    for (j in 1:k)
    {
      idx <- ind[j]
      tmp[idx, j] <- -1
    }
    W <- cbind(W, tmp)
    
    wets <- as.matrix(wets)
    W <- as.matrix(W)
    
  }
  
  L <- W %*% t(W)
  
  # S <- array(k, c(nSam, 1))
  # S <- S^.5
  # S <- diag(S)
  # val <- sqrt(k)
  # S <- diag(val, nSam, nSam)
  
  S <- diag(L)
  S <- sqrt(S)
  S <- diag(S)
  
  XS <- X %*% S
  
  SX <- svd(XS, nu = min(nDim, nSam), nv = min(nDim, nSam))
  u <- SX$u
  d <- SX$d
  # dd <- diag(d)
  dd <- d^-1
  dd <- diag(dd)
  
  pp <- u %*% dd
  PX <- crossprod(pp, X)
  T <- PX %*% L %*% t(PX)
  
  res <- eigen(T)
  V <- res$vectors
  D <- res$values
  
  n <- dim(V)[2]
  vv <- V[, n:1]
  # dd <- diag(D)
  dd <- D[n:1]
  dd <- diag(dd)
  
  # PX <- crossprod(vv, PX)
  vv <- vv[, 1:r]
  P <- pp %*% vv
  qres <- qr(P)
  P <- qr.Q(qres)
  PX <- crossprod(P, X)
  
  res <- list(PX=PX, P=P)
  return(res)
  
  
}





