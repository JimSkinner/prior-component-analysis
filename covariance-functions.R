MR.cov <- function(X, l) {
  stopifnot(is.matrix(X))
  stopifnot(length(l)==ncol(X))
  if (all(l==0)) return(diag.spam(nrow(X)))

  #R = cleanup(spind2spam(fields.rdist.near(X%*%diag(1/l), delta=1, max.points=0.05*nrow(X)^2)))

  # mean.neighbour is volume of a 2-ball times the number of
  # neighbours per unit volume
  # TODO: Can work w/ squashed circle instead of max(l)
  nBallVol = ceiling(pi^(0.5*ncol(X)) * max(l)^ncol(X) / gamma(0.5*ncol(X) + 1))

  R = cleanup(spind2spam(fields.rdist.near(
        X%*%diag(1/l, nrow=ncol(X), ncol=ncol(X)),
        delta=1,
        mean.neighbor=nBallVol*nrow(X))))
        #mean.neighbor=ceiling(pi*nrow(X)*max(l)^2))))
        #mean.neighbor=ceiling(2*pi*max(l)*nrow(X)))))

  K = ((2+cos(2*pi*R))*(1-R)/3 + sin(2*pi*R)/(2*pi)) + diag.spam(nrow(X))
  return(as.dgCMatrix.spam(K))
}

MR.cov.deriv <- function(X, l) {
  stopifnot(is.matrix(X))
  stopifnot(length(l)==ncol(X))

  n = nrow(X)

  nBallVol = min(pi^(0.5*ncol(X)) * max(l)^ncol(X) / gamma(0.5*ncol(X) + 1), 1)

  R = cleanup(spind2spam(fields.rdist.near(
        X%*%diag(1/l, nrow=ncol(X), ncol=ncol(X)),
        delta=1,
        mean.neighbor=ceiling(nBallVol)*nrow(X))))

  supp = spam2spind(R)$ind # Support of R

  derivs = lapply(1:length(l), function(i) NA)
  dKdl = list()
  for (i in 1:length(l)) {
    # Construct matrix of pairwise distances between values in column i
    Ri = spam(list(ind=supp,
                   val=((X[supp[,1],i] - X[supp[,2],i])/l[i])^2),
              nrow=n, ncol=n)

    # Eq. 25
    dKdl[[i]] = as.dgCMatrix.spam(spam(ncol=n, nrow=n, list(ind=supp,
       val=((4*(pi*(1-R[supp])*cos(pi*R[supp]) + sin(pi*R[supp])) *
             sin(pi*R[supp]) * Ri[supp] / (3 * R[supp] * l[i]))))))
  }
  return(dKdl)
}
