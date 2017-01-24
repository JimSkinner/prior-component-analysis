library(fields)
library(igraph)

cov.wendland <- function(X, beta, D=NA) {
  sigma0     = exp(beta[1])
  max.dist   = exp(beta[2])
  smoothness = exp(beta[3]) + 1
  dimension  = ncol(grid)
  if (is.na(D)) {D = distanceMatrix(X, max.dist=max.dist)}
  scale.const = wendland.eval(0, n=dimension, k=smoothness, derivative=0)
  covx = wendland.eval(D@x, n=dimension, k=smoothness, derivative=0)/scale.const
  D@x = covx
  return(sigma0 * D)
}

cov.matern <- function(X, beta, D=NA) {
  sigma0     = exp(beta[1])
  range      = exp(beta[2])
  smoothness = exp(beta[3])
  if (is.na(D)) {D = distanceMatrix(X, max.dist=max.dist)}
  covx = Matern(D@x, scale=1, range=range, smoothness=smoothness)
  D@x = covx
  return(sigma0 * D)

}

cov.exp.wend = function(X, beta, D=NA) {
  sigma0     = exp(beta[1])
  max.dist   = exp(beta[2])
  smoothness.wend = exp(beta[3]) + 1
  smoothness.mat  = exp(beta[4])
  return(sigma0 *
         cov.wendland(X, c(1, beta[2:3]), D) *
         cov.matern(X, c(1, beta[c(2,4)]), D))
}

cov.exp.wend = function(D, beta) {
  sigma0 = exp(beta[1])
  lenscale = exp(beta[2])
  smoothness.taper  = exp(beta[3]) + 1

  

  ## mean.neighbour is volume of a 2-ball times the number of
  ## neighbours per unit volume
  #delta = 0.01
  #nBallVol = (pi^(0.5*ncol(grid)) * delta^ncol(grid)) / gamma(0.5*ncol(grid) + 1)


  #D = fields.rdist.near(grid, delta=delta,
  #      mean.neighbor=ceiling(nBallVol*nrow(grid)))


  sigma0*as.dgCMatrix.spam(stationary.taper.cov(x1=grid,
    Covariance="Exponential", p=1, theta=lenscale,
    Taper="Wendland", Taper.args=list(theta=0.1, k=smoothness.taper, dimension=ncol(grid))))
  ##TODO: Make me a CORRELATION Matrix
}

cov.matern.wend = function(grid, beta) {
  sigma0 = exp(beta[1])
  lenscale = exp(beta[2])
  smoothness.matern = exp(beta[3])
  smoothness.taper  = exp(beta[4]) + 1
  sigma0*as.dgCMatrix.spam(stationary.taper.cov(x1=grid,
    Dist.args=NULL,
    Covariance="Matern", theta=lenscale, smoothness=smoothness.matern,
    Taper="Wendland", Taper.args=list(theta=0.1, k=smoothness.taper, dimension=ncol(grid))))
}

dim.faims = c(512, 51)
A = get.adjacency(make_lattice(dim.faims))
covar.fn.neighbours = function(gamma) {
  return(Diagonal(ncol(X)) + gamma*A)
}

