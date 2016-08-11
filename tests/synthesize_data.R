library(MASS)
library(fields)
library(Matrix)
library(spam)

synthesize_data_kern <- function(n, k, dim, lenScale=rep(1, length(dim)), sigSq=0, kern=Exp.cov) {
  grid = make.surface.grid(lapply(dim, function(d_) seq(0, 1, length=d_)))
  K    = kern(grid%*%diag(1/sqrt(lenScale)))
  if (is.spam(K)) K = as.dgCMatrix.spam(K)
  data.synth = synthesize_data(n, k, K, sigSq)
  return(list(X    = data.synth$X,
              W    = data.synth$W,
              V    = data.synth$V,
              grid = grid,
              K    = K))
}

synthesize_data <- function(n, k, K, sigSq=0) {
  d   = ncol(K)
  R_K = chol(K)
  W   = t(R_K) %*% matrix(rnorm(d*k), nrow=d, ncol=k)
  V   = mvrnorm(n=n, mu=rep(0, k), Sigma=diag(k))
  X   =  V %*% t(W) +
           matrix(rnorm(n*d, sd=sqrt(sigSq)), nrow=n, ncol=d)

  X.noisy = X + matrix(rnorm(n*d, sd=sqrt(sigSq)), nrow=n, ncol=d)
  return(list(X   = X.noisy,
              W   = W,
              V   = V))
}

real_data <- function() {

}
