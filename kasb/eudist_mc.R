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
  #XY <- crossprod(X, Y)
  XY <- t(X) %*% Y
  
  sX <- array(XX, c(xSam, ySam))
  # sX <- t(sX)
  sY <- array(YY, c(ySam, xSam))
  sY <- t(sY)
  
  res <- sX + sY - 2*XY
  
  if (issqrt == TRUE)
  {
    res <- sqrt(res)
  }
  
  # output <- list(D=res)
  # return(output)
  
  return(res)
  
}

