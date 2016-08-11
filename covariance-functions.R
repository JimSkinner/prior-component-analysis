

MR.cov <- function(X, l) {
  if (all(l==0)) return(diag.spam(nrow(X)))
  R = cleanup(spind2spam(fields.rdist.near(X%*%diag(1/l), delta=1, max.points=0.05*nrow(X)^2)))
  K = sigma0*((2+cos(2*pi*R))*(1-R)/3 + sin(2*pi*R)/(2*pi)) + diag.spam(nrow(X))
  return(as.dgCMatrix.spam(K))
}
