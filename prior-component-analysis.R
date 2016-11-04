library(functional)
library(Matrix)
library(optimx)

prca <- function(X, k, covar.fn, beta.init=c(), maxit=10, tol=1e-2, trace=0,
                 report_iter=10, warnDiag=TRUE) {
  stopifnot(ncol(X) >= k)

  X = scale(X, scale=FALSE) # Centered data: nxd
  n = nrow(X)               # number of samples
  d = ncol(X)               # original dimensionality

  svd.X = svd(X)
  W     = svd.X$v[,1:k,drop=FALSE] %*% diag(svd.X$d[1:k], ncol=k, nrow=k)
  #W     = svd.X$v[,1:k,drop=FALSE] # TODO: Which init is better?
  sigSq = mean(svd.X$d[-(1:k)]^2)

  if (sigSq < 1e-10) {warning("The data provided lie close to a subspace of",
    "dimensionality equal to or lower than the k provided; prca may fail due",
    "to producing a degenerate probability model. Maybe pick a smaller k?")}

  if (trace>=1) print(paste("Starting prca with", length(beta), "hyperparameters"))
  if (length(beta.init)==0) {
    # covar.fn has no hyperparameters
    K = covar.fn()
  } else {
    beta  = beta.init
    K     = covar.fn(beta)
  }

  # Test for ill conditioning. Need well conditioned matrix if we are HP tuning.
  if (length(beta.init)>0 && condest(K)$est > 10^4.5) {
    stop(paste("The covariance matrix constructed with the covariance function",
      "and starting parameters provided is ill-conditioned. The first",
      "iteration requires a well-conditioned covariance matrix."))
  }

  # TODO: Handle regular matrix case (cast to the correct Matrix)
  stopifnot(is(K, "Matrix"))
  stopifnot(isSymmetric(K))

  if (warnDiag & Matrix::isDiagonal(K)) warning(paste("The covariance matrix",
    "constructed from covar.fun with parameters beta.init is diagonal. This",
    "can sometimes cause beta optimisation to get stuck. If all inputs are",
    "truly independent, this may not be a good technique to use."))

  lp  = -Inf # Log likelihood
  lps = numeric(maxit)

  logPost = -Inf

  converged = FALSE
  iteration = 0
  while (!converged) {
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

    if (trace >=2) {
      newLogPost = prca.log_posterior(X, K, W, sigSq)
      cat("Iteration ", iteration, ", Updated sigma^2. Log Posterior=", round(newLogPost, 5), " (increase=", round(newLogPost-logPost, 5),")\n", sep='')
      logPost = newLogPost
    }

    ## Expectation Step 2
    WtW = crossprod(W)
    Minv = chol2inv(chol(WtW + sigSq*diag(k)))
    E_V1 = X %*% W %*% Minv
    E_V2 = lapply(1:n, function(i_) sigSq*Minv + tcrossprod(E_V1[i_,]))

    # Maximization step for W
    xvsum = Reduce('+', lapply(1:n, function(i_) tcrossprod(X[i_,], E_V1[i_,])))
    vvsum.eig = eigen(Reduce('+', E_V2), symmetric=TRUE)
    vvsuminv.eig = list(values=rev(1/vvsum.eig$values),
                        vectors=vvsum.eig$vectors[,k:1])
    C.tilde  = K %*% xvsum %*% vvsuminv.eig$vectors %*% diag(vvsuminv.eig$values, ncol=k, nrow=k)

    # TODO: Squash this bug:
    # Error in vapply(1:k, function(i_) { : values must be length 3258,
    # but FUN(X[[1]]) result is length 0
    tryCatch({

    W.tilde = vapply(1:k, function(i_) { # TODO: Can make this faster using R_K
      Matrix::solve(K + sigSq*vvsuminv.eig$values[i_]*Diagonal(d), C.tilde[,i_])@x
    }, numeric(d))

    }, error=function(msg) browser())

    W = W.tilde %*% t(vvsuminv.eig$vectors)

    if (trace >=2) {
      newLogPost = prca.log_posterior(X, K, W, sigSq)
      cat("Iteration ", iteration, ", Updated W. Log Posterior=", round(newLogPost, 5), " (increase=", round(newLogPost-logPost, 5),")\n", sep='')
      logPost = newLogPost
    }

    if (length(beta.init)!=0) { # covar.fn has hyperparameters to tune
      # The part of the expected log likelihood with varies with beta

      ## Expectation Step 3
      WtW = crossprod(W)
      Minv = chol2inv(chol(WtW + sigSq*diag(k)))
      E_V1 = X %*% W %*% Minv
      E_V2 = lapply(1:n, function(i_) sigSq*Minv + tcrossprod(E_V1[i_,]))

      min.f <- function(beta_) {
        K_ = covar.fn(beta_)

        success = tryCatch({
          K_chol = Matrix::chol(K_, pivot=FALSE, cache=FALSE) # Pivot?
          TRUE
        }, error = function(err) {
          FALSE
        })

        if (success) {
          return(2*k*sum(log(diag(K_chol)))
                 + norm(solve(t(K_chol), W), type='F')^2)
        } else {
          return(Inf)
        }
      }

      # TODO: Relative tolerance = relatice change in LL from updating W? (Maybe lower-bounded by ~1e6)
      beta.opt = suppressWarnings(optimx(par=beta, fn=min.f, method=c("Nelder-Mead"),
                                         control=list(trace=0, kkt=FALSE, reltol=1e-5,
                                                      starttests=FALSE, maxit=500)))
      beta     = as.numeric(coef(beta.opt)[1,])
      K        = covar.fn(beta)

      if (trace >=2) {
        newLogPost = prca.log_posterior(X, K, W, sigSq)
        cat("Iteration ", iteration, ", Updated beta. Log Posterior=", round(newLogPost, 5), " (increase=", round(newLogPost-logPost, 5),")\n", sep='')
        logPost = newLogPost
      }
    }

    lpNew  = prca.log_posterior(X, K, W, sigSq)

    ## Convergence criteria & printouts if trace > 0
    if (trace==1 && (iteration%%report_iter)==0) {
      print(paste("Iteration ", iteration, ": log likelihood = ",
                  round(lpNew, 4), " (increase=", round(lpNew-lp, 4),")", sep=''))
    }

    if ((lpNew - lp < tol) & (lpNew > lp)) {
      if (trace>0) {
        print(paste("Convergence criteria reached: delta log likelihood <", tol))
      }
      converged=TRUE
    } else if (iteration >= maxit) {
      if (trace>0) {
        print(paste("Convergence criteria reached:", iteration, "iterations"))
      }
      converged=TRUE
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

  dof = d*k - 0.5*k*(k-1) + 3 + length(beta.init) # Degrees of Freedom for PPCA + length(beta.init)
  bic = -2*lp + dof*log(n)

  prcaObj = list(W     = W,
                 sigSq = sigSq,
                 mu    = attr(X, "scaled:center"),
                 V     = E_V1,
                 Vvar  = E_V2,
                 lp    = lp,
                 lps   = lps,
                 bic   = bic,
                 beta  = beta)
  class(prcaObj) = "prca"
  return(prcaObj)
}

prca.log_likelihood <- function(X, W, sigSq) {
  d = nrow(W)
  k = ncol(W)
  n = nrow(X)

  W.sv  = svd(W, nu=0, nv=0)$d

  # This monstrosity is an efficient calculation of the log likelihood
  lla = -0.5*n*d*log(2*pi)
  llb = -0.5*n*((d-k)*log(sigSq) + sum(log(W.sv^2 + sigSq)))
  #llc = -0.5*(norm(X, 'F')^2 - sum(W.sv^2/(W.sv^2 + sigSq)))/(sigSq*n)

  R = chol(crossprod(W) + diag(k))
  llc2 = -0.5*((norm(X, 'F')^2)/sigSq -
               (norm(forwardsolve(t(R), t(W)%*%t(X)), 'F')^2)/(sigSq^2))

  #ll2 = sum(vapply(1:n, function(i_) dmvnorm(X[i_,], sigma=tcrossprod(W) + sigSq*diag(d), log=TRUE), numeric(1)))

  #print("*********************")
  #print("Log lik diff:")
  #print(llc - llc2)
  #print("*********************")

  return(lla + llb + llc2)
}

prca.log_prior <- function(K, W) {
  d = nrow(W)
  k = ncol(W)

  K_chol = Matrix::chol(K, pivot=FALSE)
  return(-0.5*( d*k*log(2*pi) +
                2*k*sum(log(diag(K_chol))) +
                norm(solve(t(K_chol), W), type='F')^2))
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
