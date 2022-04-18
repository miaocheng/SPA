# 
# This funciton partitions the given X and Y data sets.
# Input: X - input data
#        Y - input data
#        m - amount of first appended data
#        n - amount of appended data for continuous alignment
# Output: X - original data
#         Y - sequence data
#         l - length of each data set
# Coded by Miao Cheng
# Data: 2016-4-20
########################################################################################

source("./center_mc.R")


paray_mc <- function(X, Y, xL, yL, m, n){
  
  xSam <- dim(X)[2]
  ySam <- dim(Y)[2]
  
  
  ########### re-sort data set ##########
  
  ind <- 1:xSam
  ind <- sample(ind, xSam)
  X <- X[, ind]
  xL <- xL[ind]
  
  ind <- 1:ySam
  ind <- sample(ind, ySam)
  Y <- Y[, ind]
  yL <- yL[ind]
  
  ######################################
  
  nY <- Y[, 1:m]
  nL <- yL[1:m]
  tmy <- Y[, (m+1):ySam]
  tml <- yL[(m+1):ySam]
  Y <- tmy
  yL <- tml
  
  nSam <- dim(Y)[2]
  t <- nSam %/% n
  
  # s <- nSam %% n
  # 
  # if (s > 0)
  # {
  #   l <- array(n, c(t, 1))
  #   l <- cbind(l, s)
  # }
  # else if (s == 0)
  # {
  #   l <- array(n, c(t, 1))
  # }
  
  l <- array(n, c(t, 1))
  rt <- nSam - t*n
  if (rt > 0)
  {
    l <- rbind(l, rt)
    # t <- t + 1
  }
  
  pdat <- list(X=X, xL=xL, nY=nY, nL=nL, Y=Y, yL=yL, l=l)
  return(pdat)
  
  
}




