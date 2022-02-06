#
# This function calculates the kernel matrix with given data sets.
# Input: X - nDim * xSam
#        Y - nDim * ySam
#        sigma - referred parameter
#        
# Output: output kernel matrix
# Coded by Miao Cheng
# Data: 2016-3-6
#####################################################################################


source("./eudist_mc.R")



contk_mc <- function(X, Y = X, ktype, para){
  
  # ktype <- para$ktype
  
  if (ktype == "gaussian")
  {
    sigma = para$sigma
  }
  
  xSam <- dim(X)[2]
  ySam <- dim(Y)[2]
  
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  D <- eudist_mc(X, Y)
  # #D2 <- as.data.frame(D)
  # #E <- apply(D2, 2, exp)
  # 
  #D <- as.matrix(D)
  
  if (ktype == "gaussian")
  {
    E <- D / (2*sigma^2)
    K <- exp(-E)
  }
  else if (ktype == "cosine")
  {
    XX <- crossprod(X, X)
    YY <- crossprod(Y, Y)
    #XY <- crossprod(X, Y)
    XY <- t(X) %*% Y
    
    # sX <- array(XX, c(xSam, ySam))
    # # sX <- t(sX)
    # sY <- array(YY, c(ySam, xSam))
    # sY <- t(sY)
    
    #XX <- array(diag(XX), c(xSam, ySam))
    dx <- diag(XX)
    XX <- repmat_mc(dx, ySam)
    # XX <- rep(dx, ySam)
    dy <- diag(YY)
    YY <- repmat_mc(dy, xSam)
    # YY <- rep(dy, xSam)
    YY <- t(YY)
    
    tmp <- sqrt(XX) * sqrt(YY)
    K <- abs(XY) / tmp
    
  }
  
  return(K)

}



repmat_mc <- function(v, times){

  
  M <- NULL
  for (i in 1:times)
  {
    M <- cbind(M, v)
  }
  
  #tmp <- t(tmp)
  return(M)
  
  
}

