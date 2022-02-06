
# This function implements the centering operation of data set. 
# Input: X - nDim * nSam
# Output: cX - centered data
# Coded by Miao Cheng
# Data: 2016-3-6
########################################################################################



center_mc <- function(X){
  
  X <- as.matrix(X)
  
  nDim <- dim(X)[1]
  nSam <- dim(X)[2]
  
  meanx <- rowMeans(X)
  cX <- X - repmat_mc(meanx, nSam)
  
  
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


