# 
# This function implements the sequantial pattern alignment.
# Input: X - xDim * xSam
#        tY - temporal input
#        para - input parameters
# Output: yy - obtained results
# Coded by Miao Cheng
# Data: 2016-4-6
#######################################################################################

# library(MASS)

source("./contk_mc.R")


spal_mc <- function(X, Y, mx, para){
  
  M <- mx$M
  xx <- mx$ux
  Kx <- mx$Kx
  
  K <- contk_mc(X, Y, "gaussian", para)
  
  # xx <- normx_mc(ux)
  tmp <- Kx %*% M %*% K
  
  # tmx <- xx %*% t(xx)
  # tmx <- ginv(tmx) %*% xx
  
  tmx <- t(xx)
  tmx <- ginv(tmx)
  
  yy <- tmx %*% tmp
  
  
  return(yy)
  
}


