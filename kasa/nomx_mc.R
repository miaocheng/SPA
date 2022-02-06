# 
# This function normalizes the given data set
# Input: X - nDim * nSam
# Output: xx
# Coded by Miao Cheng
# Data: 2016-4-26
#########################################################################################



nomx_mc <- function(X){
  
  
  nDim <- dim(X)[1]
  nSam <- dim(X)[2]
  
  tmp <- X * X
  s <- apply(tmp, 2, sum)
  
  ss <- s^-1
  ss <- diag(ss, nSam, nSam)
  E <- matrix(1, nSam, nDim)
  
  tmp <- ss %*% E
  tmp <- t(tmp)
  xx <- X * tmp
  
  return(xx)
  
  
}