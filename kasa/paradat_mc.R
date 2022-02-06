# 
# This function implements the patition of given data into sequence data sets.
# Input: dat - nDim * nSam
#        L - data labels
#        m - amount of training data
#        n - amount of appended data for continuous alignment
# Output: X - original data
#         Y - sequence data
#         l - length of each data set
# Coded by Miao Cheng
# Data: 2016-4-6
#######################################################################################

source("./center_mc.R")


paradat_mc <- function(xDat, m, n){
  
  dat <- xDat$fea
  L <- xDat$labels
  
  dat <- as.matrix(dat)
  nSam <- dim(dat)[2]
  
  ind <- 1:nSam
  ind <- sample(ind, nSam)
  dat <- dat[, ind]
  L <- L[ind]
  
  dat <- center_mc(dat)
  
  X <- dat[, 1:m]
  xL <- L[1:m]
  Y <- dat[, (m+1):nSam]
  yL <- L[(m+1):nSam]
  
  t <- (nSam-m) %/% n
  l <- array(n, c(t, 1))
  rt <- nSam - m - t*n
  if (rt > 0)
  {
    l <- rbind(l, rt)
    # t <- t + 1
  }
  
  pdat <- list(X=X, xL=xL, Y=Y, yL=yL, l=l)
  return(pdat)
  
}


