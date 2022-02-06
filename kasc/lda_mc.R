
#This function implements the efficient version of linear discriminant analysis (LDA),
#in which high-dimensional data can be handled fast.
# Input:   X - nDim * nSam
#          L - Labels of training set
#          r - Reduced dimension
# Output:  rX - Obtained data with reduced dimension
#          ldaV - Projection vectors of LDA
#          wD - Eigenvalues of associated eigenvectors
#          
# Coded by Miao Cheng
# Data: 2015-11-5


# library(base)


lda_mc <- function(X, L, r = 0)
{
  nSam = length(X[1,])
  nDim = length(X[,1])
  X <- as.matrix(X)
  L <- as.matrix(L)
  meanx <- rowMeans(X)
  
  xL <- unique(L)
  nL <- length(xL)
  B <- NULL
  W <- NULL
  for (i in 1:nL)      #seq(from=1, to=nL, by=1))
  {
    # ml <- L %in% xL(i)
    # ml <- match(L, xL(i), nomatch=0)
    ml <- which(L==xL[i], arr.ind=TRUE)
    nC <- length(ml)
    cX <- subset(X, select = ml)
    
    # Constructing between and within scatter
    meanc <- rowMeans(cX)
    B <- cbind(B, meanc)
    wX <- cX - rep(meanc, nC)
    W <- cbind(W, wX)
  }
  B <- B - rep(meanx, nL)
  
  resB = svd(B)
  U <- resB$u
  S <- resB$d
  
  # Calculating the discriminant subspaces spanned by eigens of LDA
  # W <- crossprod(U, W)
  # ss <- solve(S)
  # ss <- diag(S)
  ss <- S^(-1)
  ss <- diag(ss)
  bV <- U %*% ss
  W <- crossprod(bV, W)
  
  resW <- svd(W)
  U <- resW$u
  S <- resW$d
  rW <- length(U[1,])
  wV <- subset(U, select = c(rW:1))
  # wD <- diag(S)
  wD <- rev(S)
  # wD <- diag(wD)
  
  V <- bV %*% wV
  
  if(r != 0 && r < rW)
  {
    V <- V[,1:r]
    wD <- wD[1:r]
  }
  
  qres <- qr(V)
  V <- qr.Q(qres)
  
  rX <- crossprod(V, X)
  
  output <- list(rX=rX, V=V, D=wD)
  
  return(output)
  
  
}





