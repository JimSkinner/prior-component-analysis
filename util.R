distanceMatrix <- function(X, X2=NA, max.dist=Inf, max.points=NA) {
  stopifnot(!(all(is.na(nrow(X)))))

  symmetric=FALSE
  if (missing(X2) || all(is.na(X2))) {
    symmetric=TRUE
    X2 = X
    d = nrow(X)
  } else {
    d = max(nrow(X), nrow(X2))
  }

  if (!is.infinite(max.dist)) {
    # Take a guess at the maximum number of points needed. Calculate the volume
    # of data-space neighboring a single point. The maximum number of ponts is
    #    #points x (#points x neighboring.volume)
    #  =
    #    #points x (expected #neighbours)
    if (is.na(max.points)) {
      d_s        = ncol(X) # Spatial dimensions
      #nBallVol   = pi^(0.5*d_s) * max.dist^d_s / gamma(0.5*d_s+1)
      #max.points = ceiling(nBallVol*d*d)
      nCubeVol   = (2*max.dist)^d_s
      max.points = ceiling(nCubeVol*d)*d*2
    }

    # Use fields.rdist.near to find all distances within max.dist, stored as triples
    Dtriples = fields.rdist.near(X, X2, delta=max.dist, max.points = max.points)

    if (symmetric) {
      upInd = Dtriples$ind[,1] <= Dtriples$ind[,2]
      D = sparseMatrix(i=Dtriples$ind[upInd,1], j=Dtriples$ind[upInd,2],
                       dims=c(d,d), x=Dtriples$ra[upInd], symmetric=TRUE,
                       index1=TRUE)
    } else {
      D = sparseMatrix(i=Dtriples$ind[,1], j=Dtriples$ind[,2],
                       dims=c(nrow(X), nrow(X2)), x=Dtriples$ra, index1=TRUE)
    }
  } else {
    if (symmetric) {
      D = Matrix(rdist(X))
    } else {
      D = Matrix(rdist(X, X2))
    }
  }
  return(D)
}
