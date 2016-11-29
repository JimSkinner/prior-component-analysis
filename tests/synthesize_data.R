library(MASS)
library(fields)
library(Matrix)
library(spam)

synthesize_data_kern <- function(n, k, dim, lenScale=rep(1, length(dim)), noisesd=0, kern=Curry(Exp.cov, p=2)) {
  grid = make.surface.grid(lapply(dim, function(d_) seq(0, 1, length=d_)))
  K    = kern(grid%*%diag(1/sqrt(lenScale)))
  if (is.spam(K)) K = as.dgCMatrix.spam(K)
  data.synth = synthesize_data(n, k, K, noisesd)
  return(list(X    = data.synth$X,
              W    = data.synth$W,
              V    = data.synth$V,
              grid = grid,
              K    = K))
}

synthesize_data <- function(n, k, K, noisesd=0) {
  d   = ncol(K)
  diag(K) = diag(K) + 1e-5
  R_K = chol(K)
  W   = t(R_K) %*% matrix(rnorm(d*k), nrow=d, ncol=k)
  V   = mvrnorm(n=n, mu=rep(0, k), Sigma=diag(k))
  X   =  V %*% t(W)

  X.noisy = X + matrix(rnorm(n*d, sd=sqrt(noisesd)), nrow=n, ncol=d)
  return(list(X   = X.noisy,
              W   = W,
              V   = V))
}
