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
    "to producing a degenerate probability model.")}

  if (trace==1) print(paste("Starting prca with", length(beta), "hyperparameters"))
  if (length(beta.init)==0) {
    # covar.fn has no hyperparameters
    K = covar.fn()
  } else {
    beta  = beta.init
    K     = covar.fn(beta)
  }

  # Test for ill conditioning
  if (condest(K)$est > 10^4.5) {
    stop("The covariance matrix constructed with the covariance function and starting parameters provided is ill-conditioned. The first iteration requires a well-conditioned covariance matrix.")
  }

  # TODO: Handle regular matrix case (cast to the correct Matrix)
  stopifnot(is(K, "Matrix"))
  stopifnot(isSymmetric(K))

  if (warnDiag & Matrix::isDiagonal(K)) warning(paste("The covariance matrix",
    "constructed from covar.fun with parameters beta.init is diagonal. This",
    "can sometimes cause beta optimisation to get stuck. If all inputs are",
    "truly independant, this may not be a good technique to use."))

  lp  = -Inf # Log likelihood
  lps = numeric(maxit)

  converged = FALSE
  iteration = 0
  while (!converged) {
    WtW = crossprod(W)

    ## Expectation Step
    Minv = chol2inv(chol(WtW + sigSq*diag(k)))
    E_V1 = X %*% W %*% Minv
    E_V2 = lapply(1:n, function(i_) sigSq*Minv + tcrossprod(E_V1[i_,]))

    ## Maximization step
    vvsuminv = chol2inv(chol(Reduce('+', E_V2)))
    xvsum    = Reduce('+', lapply(1:n, function(i_) tcrossprod(X[i_,], E_V1[i_,])))
    vvsuminv.eig = eigen(vvsuminv, symmetric=TRUE)
    C.tilde  = K %*% xvsum %*% vvsuminv %*% vvsuminv.eig$vectors

    # Update sigSq
    sigSq = (
      norm(X, 'F')^2 -
      sum(vapply(1:n, function(n_) {
        2*E_V1[n_,] %*% t(W) %*% X[n_,] -
        sum(vapply(1:k, function(i_) E_V2[[n_]][i_,] %*% WtW[i_,], numeric(1)))
      }, numeric(1)))
    )/(n*d)

    # Update W
    W.tilde = vapply(1:k, function(i_) { # TODO: Can make this faster using R_K
      Matrix::solve(K + sigSq*vvsuminv.eig$values[i_]*Diagonal(d), C.tilde[,i_])@x
    }, numeric(d))
    W = W.tilde %*% t(vvsuminv.eig$vectors)

    restricted.beta = FALSE
    if (length(beta.init)!=0) { # covar.fn has hyperparameters to tune
      min.f <- function(beta_) {
        K_ = covar.fn(beta_)
        K_cond = condest(K_)$est
        if (K_cond > 10^4.5) {return(Inf)} # Too ill conditioned to work with
        K_chol = tryCatch({
          Matrix::chol(K_, pivot=FALSE, cache=FALSE) # Pivot?
        }, error = function(err) browser()) # TODO: Keep the TryCatch just in case
        return(2*k*sum(log(diag(K_chol)))
               + norm(solve(t(K_chol), W), type='F')^2)
      }

      browser()
      beta.opt = suppressWarnings(optimx(par=beta, fn=min.f, method=c("Nelder-Mead"),
                                         itnmax=5, control=list(trace=0, kkt=FALSE,
                                                                starttests=FALSE)))
      beta     = beta.opt[,1:length(beta)]
      K        = covar.fn(beta)
    }

    lpNew  = prca.log_posterior(X, K, W, sigSq)

    ## Convergence criteria & printouts if trace > 0
    if (trace>1 && (iteration%%report_iter)==0) {
      print(paste("Iteration ", iteration, ": log likelihood = ",
                  round(lpNew, 4), " (delta=", round(lpNew-lp, 4),")", sep=''))
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

  if (restricted.beta) warning(paste("K became ill-conditioned when optimizing",
    "beta, so the value of beta has been restricted. It is likely that the",
    "optimum value of beta lies in a region that cannot be dealt with",
    "numerically."))

  ## Remove nonidentifiability VW^T = (VR)(WR)^T by setting R=I
  W.svd = svd(W)
  W     = W.svd$u %*% diag(W.svd$d, nrow=k, ncol=k)
  E_V1  = E_V1 %*% W.svd$v

  # Identify directionality of each component by fixing sign of 1st element to be +ve
  #P = diag(sign(W[1,]), nrow=k, ncol=k)
  #W = W %*% P
  #E_V1 = E_V1 %*% P

  dof = d*k - 0.5*k*(k-1) + 3 + length(beta.init) # Degrees of Freedom for PPCA
  bic = -2*lp + dof*log(n)

  return(list(W     = W,
              sigSq = sigSq,
              mu    = attr(X, "scaled:center"),
              V     = E_V1,
              Vvar  = E_V2,
              lp    = lp,
              lps   = lps,
              bic   = bic,
              beta  = beta))
}

prca.log_likelihood <- function(X, W, sigSq) {
  d = nrow(W)
  k = ncol(W)
  n = nrow(X)

  W.sv  = svd(W, nu=0, nv=0)$d

  # This monstrosity is an efficient calculation of the log likelihood
  lla = -0.5*n*d*log(2*pi)
  llb = -0.5*n*((d-k)*log(sigSq) + sum(log(W.sv^2 + sigSq)))
  llc = -0.5*(norm(X, 'F')^2 - sum(W.sv^2/(W.sv^2 + sigSq)))/(sigSq*n)

  return(lla + llb + llc)
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