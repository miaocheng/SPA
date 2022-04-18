#
# This function calculates the Euclidean distances between given inputs
# Input: X - nDim * xSam
#        Y - nDim * ySam
#        normsqure - FALSE
# Output: D - Obtained distance matrix
# 
# Coded by Miao Cheng
# Data: 2015-11-9
###################################################################################


eudist_mc <- function(X, Y = X, issqrt = FALSE)
{
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  xSam <- length(X[1,])
  ySam <- length(Y[1,])
  
  XX <- apply(X*X, 2, sum)
  YY <- apply(Y*Y, 2, sum)
#   XX <- as.matrix(XX)
#   YY <- as.matrix(YY)
  
  #XY <- crossprod(X, Y)
#   XY <- t(X) %*% Y
  XY <- inprod_mc(X, Y)
  
  sX <- array(XX, c(xSam, ySam))
  # sX <- t(sX)
  sY <- array(YY, c(ySam, xSam))
  sY <- t(sY)
  
  XY <- as.matrix(XY)
  sX <- as.matrix(sX)
  sY <- as.matrix(sY)
  
  res <- sX + sY - 2*XY
  
  if (issqrt == TRUE)
  {
    res <- sqrt(res)
  }
  
  # output <- list(D=res)
  # return(output)
  
  return(res)
  
}


inprod_mc <- function(X, Y)
{
  xSam <- dim(X)[2]
  ySam <- dim(Y)[2]
  res <- NULL
  
  for (i in 1:ySam)
  {
    tmp <- crossprod(X, Y[, i])
    tmp <- t(tmp)
    res <- cbind(res, tmp)
  }
  
  return(res)
}

