library(Matrix)
library(optimx)
library(ucminf)

source("util.R")

prca <- function(X, k, locations, covar.fn, covar.fn.d, beta0=c(),
                 trace=0, report_iter=10, max.dist=Inf,
                 maxit=10, tol=1e-2, ucminf.control=list()) {

  ## Define commonly used variables.
  #X = scale(X, scale=FALSE) # Centered data: nxd
  n = nrow(X)               # Number of samples
  d = ncol(X)               # Original dimensionality

  ## Perform sanity checks.
  stopifnot(ncol(X) >= k) # Cannot deal with complete/overcomplete case
  stopifnot(nrow(X) >= k) # TODO: Check if I can actually deal witht his and if I can change to an equality

  # Define a few defaults for the optimx routine in optimizing beta, and
  # overwrite them with any values specified by the user.
  overwrite = ucminf.control
  ucminf.control = list(trace=0, grtol=1e-3, xtol=1e-4, maxeval=20)
  ucminf.control[names(overwrite)] = overwrite

  ## Initialize W and sigma^2 from PPCA
  covar.svd = svd(X/sqrt(n), nu=0, nv=k)
  covar.eigval = covar.svd$d^2
  sigSq = mean(covar.eigval[-(1:k)])
  W     = covar.svd$v %*% diag(sqrt(covar.eigval[1:k] - sigSq), ncol=k, nrow=k)

  if (sigSq < 1e-10) {warning("The data provided lie close to a subspace of",
    "dimensionality equal to or lower than the k provided; prca may fail due",
    "to producing a degenerate probability model. Maybe pick a smaller k?")}

  if (trace>=1) {print(paste("Starting prca with", length(beta0), "hyperparameters"))}

  beta = beta0
  D    = distanceMatrix(locations, max.dist=max.dist)
  K    = covar.fn(locations, beta=beta, D=D, max.dist=max.dist)
  stopifnot(is(K, "Matrix"))

  lp   = prca.log_posterior(X, K, W, sigSq) # Current log posterior
  lps  = numeric(maxit) # Record of log posteriors for monitoring convergence

  converged = FALSE
  iteration = 0
  while (!converged) {
    ##################
    ## EM for sigma^2
    ##################

    ## Expectation Step
    WtW = crossprod(W)
    Minv = chol2inv(chol(WtW + sigSq*diag(k)))
    E_V1 = X %*% W %*% Minv
    E_V2 = lapply(1:n, function(i_) sigSq*Minv + tcrossprod(E_V1[i_,]))

    ## Maximization step for sigma^2
    E_V2sum = Reduce('+', E_V2)

    sigSq = (norm(X, 'F')^2 -
             2*sum(vapply(1:n, function(n_) E_V1[n_,] %*% t(W) %*% X[n_,], numeric(1))) +
             sum(vapply(1:d, function(d_) W[d_,] %*% E_V2sum %*% W[d_,], numeric(1))))/(n*d)

    ##################
    ## EM for W
    ##################

    ## Expectation Step
    WtW = crossprod(W)
    if (all(WtW==0)) {stop(paste("SPCA has failed due to numerical instability.",
      "Try dividing X by it's largest singular value to improve numerical",
      "stability"))}
    Minv = chol2inv(chol(WtW + sigSq*diag(k)))
    E_V1 = X %*% W %*% Minv
    E_V2 = lapply(1:n, function(i_) sigSq*Minv + tcrossprod(E_V1[i_,]))

    ## Maximization step for W
    xvsum = Reduce('+', lapply(1:n, function(i_) tcrossprod(X[i_,], E_V1[i_,])))
    vvsum.eig = eigen(Reduce('+', E_V2), symmetric=TRUE)
    vvsuminv.eig = list(values=rev(1/vvsum.eig$values),
                        vectors=vvsum.eig$vectors[,k:1])

    C.tilde  = (K %*% xvsum %*% vvsuminv.eig$vectors %*%
                diag(vvsuminv.eig$values, ncol=k, nrow=k))

    ## I have been unable to speed this up by doing Kc <- Cholesky(K) and
    ## calculating each W.tilde_i with a diagonal update to Kc. Using updown
    ## gives the correct answer but is MUCH slower and using update gives the
    ## wrong answer.
    W.tilde = vapply(1:k, function(i_) {
      Kc = Cholesky(K, Imult=sigSq*vvsuminv.eig$values[i_], LDL=FALSE, perm=TRUE)
      as.vector(Matrix::solve(Kc, C.tilde[,i_], system="A"))
    }, numeric(d))
    W = W.tilde %*% t(vvsuminv.eig$vectors)

    ##################
    ## EM for beta
    ##################
    if (length(beta0) > 0) { # There are HPs to tune
      ## Expectation Step
      WtW = crossprod(W)
      Minv = chol2inv(chol(WtW + sigSq*diag(k)))
      E_V1 = X %*% W %*% Minv
      E_V2 = lapply(1:n, function(i_) sigSq*Minv + tcrossprod(E_V1[i_,]))

      # The part of the expected log likelihood which varies with beta
      min.f <- function(beta_) {
        K_ = covar.fn(locations, beta=beta_, D=D, max.dist=max.dist)

        return(-prca.log_prior(K_, W))

        #tryCatch({
        #  Kc   = Cholesky(K_, LDL=FALSE, pivot=TRUE)
        #}, error=function(msg) browser())
        #perm = Kc@perm + 1

        #logDet    = as.numeric(Matrix::determinant(Kc, logarithm=TRUE)$modulus)*2
        #trWtKinvW = norm(Matrix::solve(Kc, W[perm,,drop=FALSE], system='L'), 'F')^2
        #return(k*logDet + trWtKinvW)

        #Ksub   = nearPD(covar.fn(locations.sub, beta_)) #TODO: Precondition??
        #logDet = log(sum(Ksub$eigenvalues))

        #R = chol(as.matrix(Ksub$mat), permut=TRUE)

        #trWtKinvW = norm(Matrix::solve(Kc, W.sub[perm,,drop=FALSE], system='L'), 'F')^2

      }

      min.f.d <- function(beta_) {
        K_  = covar.fn(locations, beta=beta_, D=D, max.dist=max.dist)
        dK_ = covar.fn.d(locations, beta=beta_, D=D, max.dist=max.dist)

        Kc  = Cholesky(K_, LDL=TRUE, pivot=TRUE)

        deriv = numeric(length(beta_))
        for (i in 1:length(beta_)) { # TODO: Neaten me
          a = k * sum(diag(solve(Kc, dK_[[i]], system="A")))

          b = solve(Kc, W, system="A")

          c = sum(vapply(1:k, function(k_) {
            as.numeric(b[,k_] %*% dK_[[i]] %*% b[,k_])
          }, numeric(1)))

          deriv[i] = a-c
        }
        return(0.5*deriv)
      }

      # TODO: Relative tolerance = relative change in LL from updating W? (Maybe lower-bounded by ~1e6)
      beta.opt = ucminf(par=beta, fn=min.f, gr=min.f.d, control=ucminf.control)
      beta     = beta.opt$par # ucminf
      K        = covar.fn(locations, beta=beta, D=D, max.dist=max.dist)
    }

    ####################################
    ## Convergence criteria & printouts if trace > 0
    ####################################
    lpNew = prca.log_posterior(X, K, W, sigSq)
    if (iteration >= maxit) { # Check for maximum iterations reached. If so, print.
      if (trace>0) {
        print(paste("Convergence criteria reached:", iteration, "iterations"))
      }
      converged=TRUE
    } else if (trace==1 && (iteration%%report_iter)==0) {
      print(paste("Iteration ", iteration, ": log likelihood = ",
                  round(lpNew, 4), " (increase=", round(lpNew-lp, 4),")", sep=''))
    }

    lp = lpNew
    iteration = iteration + 1
    lps[iteration] = lp
  }

  if (iteration < maxit) { # Trim unused lps entries
    lps = lps[1:iteration]
  }

  ## Remove nonidentifiability VW^T = (VR)(WR)^T by setting R=I
  W.svd = svd(W)
  W     = W.svd$u %*% diag(W.svd$d, nrow=k, ncol=k)
  E_V1  = E_V1 %*% W.svd$v

  # TODO
  # Identify directionality of each component by fixing sign of 1st element to be +ve
  #P = diag(sign(W[1,]), nrow=k, ncol=k)
  #W = W %*% P
  #E_V1 = E_V1 %*% P

  ll = prca.log_likelihood(X, W, sigSq)
  dof = d*k - 0.5*k*(k-1) + 3 + length(beta0) # Degrees of Freedom for PPCA + #HPs
  bic = -2*ll + dof*log(n)

  prcaObj = list(W     = W,
                 sigSq = sigSq,
                 mu    = attr(X, "scaled:center"),
                 V     = E_V1,
                 #Vvar  = E_V2,
                 ll    = ll,
                 lp    = lp,
                 lps   = lps,
                 bic   = bic,
                 beta  = beta)
  class(prcaObj) = "prca"
  return(prcaObj)
}

prca.log_likelihood <- function(X, W, sigSq) {
  # Careful! Data X should be centered
  d = nrow(W)
  k = ncol(W)
  n = nrow(X)

  W.sv  = svd(W, nu=0, nv=0)$d
  R   = chol(crossprod(W)/sigSq + diag(k))

  # This monstrosity is an efficient calculation of the log likelihood
  lla = -0.5*n*d*log(2*pi)
  llb = -0.5*n*((d-k)*log(sigSq) + sum(log(W.sv^2 + sigSq)))
  llc = -0.5*((norm(X, 'F')^2)/sigSq -
              (norm(forwardsolve(t(R), t(W)%*%t(X)), 'F')^2)/(sigSq^2))

  return(lla + llb + llc)
}

prca.log_prior <- function(K, W) {
  d = nrow(W)
  k = ncol(W)

  if (is(K, "sparseMatrix")) {
    Kc = Cholesky(K, LDL=TRUE, pivot=TRUE)

    # Fast calculation of the log determinant of K
    logDetK = as.numeric(-determinant(solve(Kc, system='D'), log=TRUE)$modulus)
    #logDetK2 = as.numeric(determinant(K, log=TRUE)$modulus) # TODO: Benchmark, and confirm equality

    # Fast calculation of Tr[W' K^{-1} W]
    KinvW     = solve(Kc, W, system="A")
    trWtKinvW = sum(vapply(1:k, function(k_) (W[,k_] %*% KinvW[,k_])[1,1],
                           numeric(1)))
  } else {
    stop("Still need to implement dealing w/ non-sparse matrices")
  }


  return(-0.5*( k*d*log(2*pi) +
                k*logDetK +
                trWtKinvW))

}

prca.log_posterior <- function(X, K, W, sigSq) {
  return(prca.log_likelihood(X, W, sigSq) + prca.log_prior(K, W));
}

predict.prca <- function(prcaObj, samples) {
  if (missing(samples)) {
    return(prcaObj$V)
  }

  stopifnot(is.matrix(samples))

  k     = ncol(prcaObj$W) # Latent dimensionality
  muMat = matrix(rep(prcaObj$mu, nrow(samples)), nrow=nrow(samples), byrow=TRUE)

  # Latent representation of new samples
  Mchol  = chol(crossprod(prcaObj$W) + prcaObj$sigSq*diag(k))
  latent = t(backsolve(Mchol, forwardsolve(t(Mchol), t((samples-muMat)%*%prcaObj$W))))

  # 'samples' projected onto PrCA model
  proj = tcrossprod(latent, prcaObj$W) + muMat
  return(proj)
}
