
# This function implements the K-nearest neighbor classifier.
# Input: X - xDim * xSam
#        L - xSam
#        Y - yDim * ySam
#        tL - ySam
#        k - number of nearest neighbors
# Output: accu - number of samples that are correctly recognized
#         accuracy - ratio of correct recognition
#         
# Coded by Miao Cheng
# Data: 2015-11-9


library(stats)


knnsolo_mc <- function(X, xL, Y, yL, k = 1)
{
  xSam <- length(X[1,])
  ySam <- length(Y[1,])
  
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  xL <- as.matrix(xL)
  yL <- as.matrix(yL)
  
  tdist <- array(0, c(xSam, 1))
  class <- array(0, c(k, 1))
  
  accu <- 0
  err <- 0
  
  for (i in 1:ySam)
  {
    for (j in 1:xSam)
    {
      tmp = Y[, i] - X[, j]
      tdist[j] = norm(crossprod(tmp, tmp))
    }
    dist <- as.matrix(sort(tdist))
    idx <- as.matrix(order(tdist))
    
    mdist <- dist[xSam] + 100
    cou <- 0
    kdist <- dist[1:k]
    rm(dist)
    
    for (j in 1:k)
    {
      
      idj = idx[j]
      tmp = xL[idj]
      class[j] = tmp
      
      if (yL[i] == tmp)
      {
        cou <- cou + 1
      }
      
    }
    uclass <- unique(class)
    cnum <- dim(uclass)[1]
    if (cou > k/2)
    {   result <- yL[i]   }
    else if (cnum == k)
    {   result <- class[1]   }
    else
    {
      ccou <- array(0, c(cnum, 1))
      
      for (j in 1:cnum)
      {
        ind <- which(class == uclass[j], arr.ind = TRUE)
        ccou[j] <- length(ind)
      }
      
      C <- max(ccou)
      I <- which(ccou == C, arr.ind = TRUE)
      if (ccou[I] == cou)
      {
        ind1 <- which(ccou == C, arr.ind = TRUE)
        
        for (m in 1:length(ind1))
        {
          ind2 <- which(class == uclass[ind1[m]])
          cdist <- kdist[ind2]
          if (cdist[1] < mdist)
          {
            mdist <- cdist[1]
            result <- class[ind2[1]]
          }
        }
      }
      else
      { result <- uclass[I] }
      
    }
    
    if (result == yL[i])
    {   accu <- accu + 1   }
    else
    {   err <- err + 1   }
    
  }
  
  accuracy <- accu / ySam
  
  output <- list(accu=accu, accuracy=accuracy)
  return(output)
  
}

