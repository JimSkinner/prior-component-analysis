library(Matrix)
library(ucminf)
library(numDeriv)

source("util.R")

# TODO: If passed a prca obj, continue inference.

prca <- function(X, k, locations, covar.fn, covar.fn.d=NA, beta0=c(),
                 trace=0, report_iter=10, max.dist=Inf,
                 maxit=10, tol=1e-2, ucminf.control=list()) {

  ## Define commonly used variables.
  Xc = scale(X, scale=FALSE) # Centered data: nxd
  mu = attr(Xc, "scaled:center")
  n  = nrow(X)               # Number of samples
  d  = ncol(X)               # Original dimensionality

  ## Perform sanity checks.
  stopifnot(ncol(X) > k) # Cannot deal with complete/overcomplete case
  stopifnot(nrow(X) >= k) # TODO: Check if I can deal with equality case

  if (length(ucminf.control)>0) {stop("ucminf.control not yet implemented")}
  ## Define a few defaults for ucminf in optimizing beta, and
  ## overwrite them with any values specified by the user.
  #overwrite = ucminf.control
  #ucminf.control = list(trace=0, grtol=1e-3, xtol=1e-4, maxeval=20)
  #ucminf.control[names(overwrite)] = overwrite

  ## Initialize W and sigma^2 from PPCA
  covar.svd = svd(Xc/sqrt(n), nu=0, nv=k)
  covar.eigval = covar.svd$d^2
  sigSq = sum(covar.eigval[-(1:k)])/(d-k)
  W     = covar.svd$v %*% diag(sqrt(covar.eigval[1:k] - sigSq), ncol=k, nrow=k)

  if (sigSq < 1e-10) {warning("The data provided lie close to a subspace of",
    "dimensionality equal to or lower than the k provided; prca may fail due",
    "to producing a degenerate probability model. Maybe pick a smaller k?")}

  if (trace>=1) {print(paste("Starting prca with", length(beta0), "hyperparameters"))}

  beta = beta0
  D    = distanceMatrix(locations, max.dist=max.dist)
  K    = covar.fn(locations, beta=beta, D=D, max.dist=max.dist)
  stopifnot(is(K, "Matrix"))

  lp   = prca.log_posterior(X, K, W, mu, sigSq) # Current log posterior
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
    E_V1 = Xc %*% W %*% Minv
    E_V2 = lapply(1:n, function(i_) sigSq*Minv + tcrossprod(E_V1[i_,]))

    ## Maximization step for sigma^2
    E_V2sum = Reduce('+', E_V2)

    sigSq = (norm(Xc, 'F')^2 -
             2*sum(vapply(1:n, function(n_) E_V1[n_,] %*% t(W) %*% Xc[n_,], numeric(1))) +
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
    E_V1 = Xc %*% W %*% Minv
    E_V2 = lapply(1:n, function(i_) sigSq*Minv + tcrossprod(E_V1[i_,]))

    ## Maximization step for W
    xvsum = Reduce('+', lapply(1:n, function(i_) tcrossprod(Xc[i_,], E_V1[i_,])))
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
      if (is(K, "sparseMatrix")) {
        Kc = Cholesky(K, Imult=sigSq*vvsuminv.eig$values[i_], LDL=FALSE, perm=TRUE)
        return(as.vector(Matrix::solve(Kc, C.tilde[,i_], system="A")))
      } else {
        KplusDiag = K
        diag(KplusDiag) = diag(KplusDiag) + sigSq*vvsuminv.eig$values[i_]
        Kc = base::chol(KplusDiag, pivot=TRUE)
        pivot = attr(Kc, "pivot")
        unpivot = order(pivot)
        return(backsolve(Kc, forwardsolve(t(Kc), C.tilde[pivot,i_]))[unpivot])
      }
    }, numeric(d))
    W = W.tilde %*% t(vvsuminv.eig$vectors)

    ##################
    ## EM for beta
    ##################
    if (length(beta0) > 0) { # There are HPs to tune
      ## Expectation Step
      WtW = crossprod(W)
      Minv = chol2inv(chol(WtW + sigSq*diag(k)))
      E_V1 = Xc %*% W %*% Minv
      E_V2 = lapply(1:n, function(i_) sigSq*Minv + tcrossprod(E_V1[i_,]))

      # The part of the expected log likelihood which varies with beta
      min.f <- function(beta_) {
        #f <- function(beta_) {
        #  K_ = covar.fn(locations, beta=beta_, D=D, max.dist=max.dist)
        #  return(-prca.log_prior(K_, W))
        #}
        #H = hessian(f, beta_)
        #print(kappa(H))
        #tryCatch({
        K_ = covar.fn(locations, beta=beta_, D=D, max.dist=max.dist)
        if (any(is.nan(K_@x))) {
          #print(paste(paste(beta_, collapse=', '), ": NaN"))
          return(NaN)
        }
        tryCatch({
          prca.log_prior(K_, W)
        }, error=function(msg) browser())

        #print(paste(paste(beta_, collapse=', '), ":", -prca.log_prior(K_, W)))
        return(-prca.log_prior(K_, W))
      }

      if (is.function(covar.fn.d)) {
        min.f.d <- function(beta_) {
          K_  = covar.fn(locations, beta=beta_, D=D, max.dist=max.dist)
          dK_ = covar.fn.d(locations, beta=beta_, D=D, max.dist=max.dist)

          deriv = numeric(length(beta_))
          if (is(K_, "sparseMatrix")) {
            K_c  = Cholesky(K_, LDL=TRUE, pivot=TRUE)
            for (i in 1:length(beta_)) {
              # Tr[ K^-1 dK/dBeta ]
              working1 = k * sum(diag(solve(K_c, dK_[[i]], system="A")))

              # K^-1 W
              working2 = solve(K_c, W, system="A")

              # Tr[ W^T K^-1 dK/dBeta K^-1 K^-1 W ]
              working3 = sum(vapply(1:k, function(k_) {
                as.numeric(working2[,k_] %*% dK_[[i]] %*% working2[,k_])
              }, numeric(1)))

              deriv[i] = 0.5*(working1 - working3)
            }
          } else { # Non-sparse K case
            K_c = base::chol(as.matrix(K_), pivot=TRUE)
            for (i in 1:length(beta_)) {
              # Meanings as above
              working1 = k * sum(diag(cholSolve(K_c, dK_[[i]])))
              working2 = cholSolve(K_c, W)

              working3 = sum(vapply(1:k, function(k_) {
                as.numeric(working2[,k_] %*% dK_[[i]] %*% working2[,k_])
              }, numeric(1)))

              deriv[i] = 0.5*(working1 - working3)
            }
          }
          #print(paste("*****", paste(beta_, collapse=', '), ":", paste(deriv, collapse=',')))
          return(deriv)
        }
        #gg = make.surface.grid(list(x=seq(beta[1]*2, beta[1]/2, length=10),
        #                            y=seq(beta[2]/2, beta[2]*2, length=10)))
        #zz = apply(gg, 1, min.f)
        #plot.surface(as.surface(gg, zz), type='I')

        #beta.opt = multimin(x=beta, f=min.f, df=min.f.d, method="bfgs",
        #                    step.size=0.2, prec=1e-2)
        #beta.opt = optimx(par=beta, fn=min.f, gr=min.f.d, method="Nelder-Mead", control=list(trace=0))

        if (iteration==0) {
          H0 = nearPD(hessian(min.f, beta))$mat
          invHess = solve(H0)[lower.tri(H0,diag=TRUE)]
          stepmax = 0.2
        }
        beta.opt = ucminf(par=beta, fn=min.f, gr=min.f.d, hessian=2,
                          control=list(invhessian.lt=invHess, stepmax=stepmax))
        invHess = beta.opt$invhessian.lt
        stepmax = beta.opt$info['stepmax']*1.1
        if (all(beta.opt$par == beta)) {stepmax = stepmax*0.75}
      } else {
        beta.opt = ucminf(par=beta, fn=min.f)
      }

      # TODO: Relative tolerance = relative change in LL from updating W? (Maybe lower-bounded by ~1e6)
      beta     = beta.opt$par # ucminf
      #beta     = as.numeric(coef(beta.opt)) # optimx
      #beta = beta.opt$x #multimin
      K        = covar.fn(locations, beta=beta, D=D, max.dist=max.dist)
    }

    ####################################
    ## Convergence criteria & printouts if trace > 0
    ####################################
    lpNew = prca.log_posterior(X, K, W, mu, sigSq)
    if (iteration >= maxit) { # Check for maximum iterations reached. If so, print.
      if (trace>0) {
        print(paste("Convergence criteria reached:", iteration, "iterations"))
      }
      converged=TRUE
    } else if (trace==1 && (iteration%%report_iter)==0) {
      print(paste("Iteration ", iteration, ": log likelihood = ",
                  round(lpNew, 4), " (increase=", round(lpNew-lp, 4),")", sep=''))
    } else if (trace==2 && (iteration%%report_iter)==0) {
      print(paste("Iteration ", iteration, ": log likelihood = ",
                  round(lpNew, 4), " (increase=", round(lpNew-lp, 4),")", sep=''))
      print(paste("    beta=[", paste(beta, collapse=', '), ']', sep=''))
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

  ll = prca.log_likelihood(X, W, mu, sigSq)
  dof = d*k - 0.5*k*(k-1) + 3 + length(beta0) # Degrees of Freedom for PPCA + #HPs
  bic = -2*ll + dof*log(n)
  bic2 = -2*ll + dof*log(n)

  #hess = hessian(function(z) {
  #  prca.log_posterior(X, K,
  #                     W=matrix(z[1:(d*k)], nrow=d, ncol=k),
  #                     mu=z[d*k+1],
  #                     sigSq=z[d*k+2])
  #}, x(W, mu, sigSq))
  #logDetHess = determinant(hess, logarith=TRUE)
  #logZ = (prca.log_posterior(X, K, W, mu, sigSq) +
  #        (0.5*(d*k+2))*log(2*pi) +
  #        logDetHess)

  prcaObj = list(X     = X,
                 W     = W,
                 sigSq = sigSq,
                 mu    = mu,
                 V     = E_V1,
                 #Vvar  = E_V2,
                 ll    = ll,
                 lp    = lp,
                 lps   = lps,
                 bic   = bic,
                 beta  = beta,
                 K     = K)
  class(prcaObj) = "prca"
  return(prcaObj)
}

prca.log_likelihood <- function(X, W, mu, sigSq) {
  if (is.vector(X)) {X = matrix(X, nrow=1)}

  Xc = sweep(X, 2, mu)

  d = nrow(W)
  k = ncol(W)
  n = nrow(X)

  W.sv = svd(W, nu=0, nv=0)$d
  R    = chol(crossprod(W) + sigSq*diag(k))

  const     = n*d*log(2*pi)
  nLogDetC  = n*((d-k)*log(sigSq) + sum(log(W.sv^2 + sigSq)))
  trXCinvXt = (norm(Xc, 'F')^2 -
               norm(forwardsolve(t(R), t(W))%*%t(Xc), 'F')^2)/sigSq
  return(-0.5*(const + nLogDetC + trXCinvXt))
}

prca.log_prior <- function(K, W) {
  d = nrow(W)
  k = ncol(W)


  if (is(K, "sparseMatrix")) {

    Kc = Cholesky(K, LDL=TRUE, pivot=TRUE)

    # Fast calculation of the log determinant of K
    #logDetK = as.numeric(-determinant(solve(Kc, system='D'), log=TRUE)$modulus)

    logDetK = as.numeric(determinant(K, logarithm=TRUE)$modulus)
    if (logDetK==Inf) {return(-Inf)}

    #logDetK2 = as.numeric(determinant(K, log=TRUE)$modulus) # TODO: Benchmark, and confirm equality

    # Fast calculation of Tr[W' K^{-1} W]
    KinvW     = solve(Kc, W, system="A")
    trWtKinvW = sum(vapply(1:k, function(k_) (W[,k_] %*% KinvW[,k_])[1,1],
                           numeric(1)))
  } else {
    # TODO: Poor numrical problems? Try a QRdecomp & solve() or an LDL
    # decomposition (KFAS package)
    Kc = base::chol(K, pivot=TRUE)
    pivot = attr(Kc, "pivot")

    logDetK   = 2*sum(log(diag(Kc)))
    trWtKinvW = norm(forwardsolve(t(Kc), W[pivot,]), 'F')^2

    #KinvW     = cholSolve(Kc, W)
    #trWtKinvW = sum(vapply(1:k, function(k_) (W[,k_] %*% KinvW[,k_])[1,1],
    #                       numeric(1)))
  }


  return(-0.5*( k*d*log(2*pi) +
                k*logDetK +
                trWtKinvW))

}

prca.log_posterior <- function(X, K, W, mu, sigSq) {
  return(prca.log_likelihood(X, W, mu, sigSq) + prca.log_prior(K, W));
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

prca.log_evidence <- function(X, K, W, mu, sigSq) {
  #' Compute the laplace approximation to the log evidence given the MAP
  #' parameters K, mu, sigSq as well as the prior covariance matrix K.
  #' Note that this is multiplied by an UN-KNOWN CONSTANT due to the flat
  #' priors over mu and sigSq. However, this unknown constant is always
  #' the same regardless of k and K, so this may be used to compute
  #' meaningful bayes factors between SPCA models.

  if (is(X, "prca")) {
    prcaObj = X
    X = prcaObj$X
    K = prcaObj$K
    W = prcaObj$W
    mu = prcaObj$mu
    sigSq = prcaObj$sigSq
  }

  n = nrow(X)
  d = ncol(X)
  k = ncol(W)

  # Centered X
  #Xc = sweep(X, 2, mu)

  ## recomute C^{-1}, which is used all over the place
  #R = Matrix::chol(crossprod(W) + sigSq*diag(k))
  #Cinv = Matrix(Diagonal(d)/sigSq - crossprod(forwardsolve(t(R), t(W)))/sigSq)

  #logDetH = 0

  ## Compute each of the blocks of H corresponding to each w_i, and the log
  ## determinant of this block to the comulative log determinant
  #HW = list()
  #for (i in 1:k) {
  #  Hw = Matrix(solve(K) +
  #              Cinv*as.numeric(W[,i]%*%Cinv%*%t(Xc)%*%Xc%*%Cinv%*%W[,i] -
  #                              n*W[,i]%*%Cinv%*%W[,i] + n) +
  #              Cinv%*%(
  #                t(Xc)%*%Xc%*%Cinv%*%W[,i,drop=F]%*%t(W[,i,drop=F]) +
  #                W[,i,drop=F]%*%t(W[,i,drop=F])%*%Cinv%*%t(Xc)%*%Xc +
  #                as.numeric(W[,i]%*%Cinv%*%W[,i] - 1)*t(Xc)%*%Xc -
  #                n*W[,i,drop=F]%*%t(W[,i,drop=F])
  #              )%*%Cinv
  #  , forceCheck=TRUE)
  #  HW[[i]] = Hw
  #  logDetH = logDetH + as.numeric(determinant(Hw, logarithm=TRUE)$modulus)
  #}

  ## Compute the mu block of H, add the log det to the cumulative total
  #Hmu = n*Cinv
  #HmuLogDet = as.numeric(determinant(n*Cinv, logarithm=TRUE)$modulus)
  #logDetH   = logDetH + HmuLogDet

  ## Compute the sigSq block of H & add log det to cumulative total
  #HsigSq       = sum(diag(Xc%*%Cinv%*%Cinv%*%Cinv%*%t(Xc))) -
  #                   0.5*n*sum(diag(Cinv%*%Cinv))
  #HsigSqLogDet = log(sum(diag(Xc%*%Cinv%*%Cinv%*%Cinv%*%t(Xc))) -
  #                   0.5*n*sum(diag(Cinv%*%Cinv)))
  #logDetH      = logDetH + HsigSqLogDet

  ## Laplace-approximated log evidence
  #logZ = (prca.log_posterior(X, K, W, mu, sigSq) +
  #        (0.5*(d*k+d+1))*log(2*pi) -
  #        0.5*logDetH)
  #return(logZ)

  require(numDeriv)
  H2 = -Matrix(hessian(function(z) {
    W_ = matrix(z[1:d*k], nrow=d, ncol=k)
    mu_ = z[(d*k+1):(d*k+d)]
    sigSq_ = z[d*k+d+1]
    return(prca.log_posterior(X, K, W_, mu_, sigSq_))
  }, c(W, mu, sigSq)))
  logZ2 = (prca.log_posterior(X, K, W, mu, sigSq) +
           (0.5*(d*k+d+1))*log(2*pi) -
           0.5*as.numeric(determinant(H2, logarithm=TRUE)$modulus))
  return(logZ2)

  # Hessian of the log posterior
  #require(numDeriv)
  #H = hessian(function(z) {
  #  prca.log_posterior(X, K,
  #                     W=matrix(z[1:(d*k)], nrow=d, ncol=k),
  #                     mu=z[(d*k+1):(d*k+d)],
  #                     sigSq=z[d*k+d+1])
  #}, c(W, mu, sigSq))
  #logDetHess = as.numeric(determinant(H, logarithm=TRUE)$modulus)
  #logZ = (prca.log_posterior(X, K, W, mu, sigSq) +
  #        (0.5*(d*k+d+1))*log(2*pi) -
  #        0.5*logDetHess)
}

prca.log_bayes_factor <- function(X, K1, W1, mu1, sigSq1, K2, W2, mu2, sigSq2) {
  only2argsspecified = (!missing(X) & !missing(K1) & missing(W1) & missing(mu1)
                        & missing(sigSq1) & missing(K2) & missing(W2)
                        & missing(mu2) & missing(sigSq2))

  if(only2argsspecified) {
    model1 = X
    model2 = K1

    X1     = model1$X
    K1     = model1$K
    W1     = model1$W
    mu1    = model1$mu
    sigSq1 = model1$sigSq

    X2     = model2$X
    K2     = model2$K
    W2     = model2$W
    mu2    = model2$mu
    sigSq2 = model2$sigSq
  }

  ev1 = prca.log_evidence(X1, K1, W1, mu1, sigSq1)
  ev2 = prca.log_evidence(X2, K2, W2, mu2, sigSq2)
  return(ev1-ev2)
}
