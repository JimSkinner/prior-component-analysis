library(fields)
library(igraph)

exp.cov <- function(X, X2, beta) {
  D = Matrix(rdist(X, X2))
  D@x = exp(beta[1]) * exp(-(D@x^2)/(2*exp(beta[2])^2))
  return(D)
}

exp.cov.d <- function(X, X2, beta) {
  D = Matrix(rdist(X, X2))

  dK1 = exp.cov(X, X2, beta)
  dK2 = (D^2)*exp(-2*beta[2]) * dK1
  return(list(dK1, dK2))
}


cov.matern <- function(X, beta) {
  max.dist   = exp(beta[1])
  smoothness = exp(beta[2]) + 1
  dimension  = ncol(grid)
  D = distanceMatrix(X, max.dist=max.dist)
  scale.const = wendland.eval(0, n=dimension, k=smoothness, derivative=0)
  covx = wendland.eval(D@x, n=dimension, k=smoothness, derivative=0)/scale.const
  D@x = covx
  return(D)
}

cov.matern.wend = function(X, beta) {
  sigma0     = exp(beta[1])
  max.dist   = exp(beta[2])
  smoothness.wend = exp(beta[3]) + 1
  smoothness.mat  = exp(beta[4])
  return(sigma0 *
         cov.wendland(X, c(1, beta[2:3]), D) *
         cov.matern(X, c(1, beta[c(2,4)]), D))
}

cov.exp.wend = function(grid, beta) {
  sigma0     = exp(beta[1])
  max.dist   = exp(beta[2])
  smoothness.wend = exp(beta[3]) + 1
  smoothness.mat  = exp(beta[4])
  return(sigma0 *
         cov.wendland(X, c(1, beta[2:3]), D) *
         cov.matern(X, c(1, beta[c(2,4)]), D))
}

cov.adj = function(grid, beta) {
  sigma0     = exp(beta[1])
  coupling   = 1 / (1 + exp(-beta[2]))
  D = distanceMatrix(X, max.dist=max.dist)
  D@x = sigma0 * coupling
  diag(D) = sigma0
  return(D)
}

#MR.cov <- function(X, l) {
#  stopifnot(is.matrix(X))
#  stopifnot(length(l)==ncol(X))
#  if (all(l==0)) return(diag.spam(nrow(X)))
#
#  #R = cleanup(spind2spam(fields.rdist.near(X%*%diag(1/l), delta=1, max.points=0.05*nrow(X)^2)))
#
#  # mean.neighbour is volume of a 2-ball times the number of
#  # neighbours per unit volume
#  # TODO: Can work w/ squashed circle instead of max(l)
#  nBallVol = ceiling(pi^(0.5*ncol(X)) * max(l)^ncol(X) / gamma(0.5*ncol(X) + 1))
#
#  R = cleanup(spind2spam(fields.rdist.near(
#        X%*%diag(1/l, nrow=ncol(X), ncol=ncol(X)),
#        delta=1,
#        mean.neighbor=nBallVol*nrow(X))))
#        #mean.neighbor=ceiling(pi*nrow(X)*max(l)^2))))
#        #mean.neighbor=ceiling(2*pi*max(l)*nrow(X)))))
#
#  K = ((2+cos(2*pi*R))*(1-R)/3 + sin(2*pi*R)/(2*pi)) + diag.spam(nrow(X))
#  return(as.dgCMatrix.spam(K))
#}

MR.cov.1d <- function(X, l) {
  D   = distanceMatrix(X, max.dist=l)
  D@x = ((2+cos(2*pi*D@x))*(1-D@x)/3 + sin(2*pi*D@x)/(2*pi))
  return(D)
}

MR.cov.1d.d <- function(X, l) {
  D   = distanceMatrix(X, max.dist=l)
  D@x = cos(2*pi*D@x) - (2*pi*(1-D@x)*sin(2*pi*D@x) + 2 + cos(2*pi*D@x))/3
  return(D)
}

cov.sexp <- function(x1, x2=NA, beta=c(1, 1)) {
  exp(beta[1]) * Matrix(Exp.cov(x1, x2, theta=exp(beta[2]), p=2))
}

cov.sexp.d <- function(x1, x2=NA, beta=c(1, 1)) {
  beta[1] * Matrix(Exp.cov(x1, x2, theta=beta[2], p=2))
}

exp.MR.cov <- function(X, X2=NA, beta=c(), D=NA, max.dist=Inf, max.points=NA) {
  if (all(is.na(D))) {
    D = distanceMatrix(X, X2, max.dist=max.dist, max.points=max.points)
  }
  MRpart = ((2+cos(2*pi*D@x/max.dist))*(1-D@x/max.dist)/3 + sin(2*pi*D@x/max.dist)/(2*pi))
  expPart = exp(beta[1]) * exp(-(D@x^2)/(2*exp(beta[2])^2))
  D@x = MRpart * expPart
  return(D)
}

exp.MR.cov.d <- function(X, X2=NA, beta=c(), D=NA, max.dist=0.1, max.points=NA) {
  if (all(is.na(D))) {
    D   = distanceMatrix(X, max.dist=max.dist, max.points=max.points)
  }
  MRpart = ((2+cos(2*pi*D@x/max.dist))*(1-D@x/max.dist)/3 + sin(2*pi*D@x/max.dist)/(2*pi))

  expPart = exp(beta[1]) * exp(-(D@x^2)/(2*exp(beta[2])^2))

  dK1 = sparseMatrix(x=expPart * MRpart,
                     i=D@i, p=D@p, dims=D@Dim, symmetric=TRUE, index1=FALSE)
  dK2 = sparseMatrix(x=(D@x^2)*exp(-2*beta[2]) * expPart * MRpart,
                     i=D@i, p=D@p, dims=D@Dim, symmetric=TRUE, index1=FALSE)

  return(list(dK1, dK2))
}

matern.MR.cov <- function(X, beta, max.dist = 0.1) {
  D   = distanceMatrix(X, max.dist=max.dist)
  MRpart = ((2+cos(2*pi*D@x/max.dist))*(1-D@x/max.dist)/3 + sin(2*pi*D@x/max.dist)/(2*pi))
  maternPart = exp(beta[1])*Matern(D, range=exp(beta[2]), smoothness=exp(beta[3]))
  D@x = MRpart * maternPart
  return(D)
}

faims.cov <- function(X, X2=NA, beta=c(), D=NA, max.dist=0.1, max.points=NA) {
  MRD = D
  if (all(is.na(MRD))) {
    MRD = distanceMatrix(X, X2, max.dist=max.dist, max.points=max.points)
  }
  MRD@x = ((2+cos(2*pi*MRD@x/max.dist))*(1-MRD@x/max.dist)/3 + sin(2*pi*MRD@x/max.dist)/(2*pi))

  sigma0 = exp(beta[1])
  l1 = exp(beta[2])
  l2 = exp(beta[3])
  l3 = exp(beta[4])
  l4 = exp(beta[5])

  # TODO: Have I got x- and y- the right way around??
  # TODO: Respect taper support
  dSquared = (rdist(X[,1])/(l1 + 0.5*l2*outer(X[,1], X[,1], '+'))^2 +
              rdist(X[,2])/((l3 + 0.5*l4*outer(X[,2], X[,2], '+'))^2))
  expPart  = sigma0 * exp(-0.5*dSquared)

  M = MRD * expPart
  return(M)
}
