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

  if (trace==1) print(paste("Starting prca with", length(beta), "hyperparameters"))
  if (length(beta.init)==0) {
    # covar.fn has no hyperparameters
    K = covar.fn()
  } else {
    beta  = beta.init
    K     = covar.fn(beta)
  }

  # if ( # TODO: Check type of K is Matrix

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
    tryCatch({
    Minv = chol2inv(chol(WtW + sigSq*diag(k)))
    E_V1 = X %*% W %*% Minv
    E_V2 = lapply(1:n, function(i_) sigSq*Minv + tcrossprod(E_V1[i_,]))
    }, error=function(msg) {browser()})

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
      solve(K + sigSq*vvsuminv.eig$values[i_]*Diagonal(d), C.tilde[,i_])@x
      #solve.spam(K + sigSq*vvsuminv.eig$values[i_]*diag.spam(d), C.tilde[,i_])
    }, numeric(d))
    W = W.tilde %*% t(vvsuminv.eig$vectors)

    restricted.beta = FALSE
    if (!all(is.na(beta.init))) {
      # covar.fn has hyperparameters to tune
      min.f <- function(beta_) {
        K_ = covar.fn(beta_)

        error = tryCatch({
          K_chol = Matrix::chol(K_, pivot=FALSE); FALSE
        }, error   = function(msg) TRUE, # If K is poorly conditioned, set `error'
           warning = function(msg) TRUE  # to TRUE, and return +Inf.
        )
        if (error) {
          restricted.beta = TRUE
          return(Inf) # Too poorly conditioned to perform reliable chol
        }
        return(2*k*sum(log(diag(K_chol)))
               + norm(solve(t(K_chol), W), type='F')^2)
      }

      beta.opt = suppressWarnings(optimx(par=beta, fn=min.f, method=c("Nelder-Mead"),
                                         itnmax=5, control=list(trace=0, kkt=FALSE,
                                                                starttests=FALSE)))
      beta     = beta.opt[,1:length(beta)]
      K        = covar.fn(beta)
    }

    lpNew  = ppca.log_posterior(X, K, W, sigSq, tune.beta)

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

  lp1 = ppca.log_posterior(X, K, W, sigSq, tune.beta)
  ll1  = ppca.log_likelihood(X, W, sigSq)

  ## Remove nonidentifiability VW^T = (VR)(WR)^T by setting R=I
  W.svd = svd(W)
  W     = W.svd$u %*% diag(W.svd$d, nrow=k, ncol=k)
  E_V1  = E_V1 %*% W.svd$v

  # Identify directionality of each component by fixing sign of 1st element to be +ve
  #P = diag(sign(W[1,]), nrow=k, ncol=k)
  #W = W %*% P
  #E_V1 = E_V1 %*% P

  return(list(W     = W,
              sigSq = sigSq,
              mu    = attr(X, "scaled:center"),
              V     = E_V1,
              Vvar  = E_V2,
              lp    = lp,
              lps   = lps,
              ll    = ppca.log_likelihood(X, W, sigSq),
              beta  = beta))
}

ppca.log_likelihood <- function(X, W, sigSq) {
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

ppca.log_prior <- function(K, W) {
  d = nrow(W)
  k = ncol(W)

  K_chol = Matrix::chol(K, pivot=FALSE)
  return(-0.5*( d*k*log(2*pi) +
                2*k*sum(log(diag(K_chol))) +
                norm(solve(t(K_chol), W), type='F')^2))
}

ppca.log_posterior <- function(X, K, W, sigSq, tune.beta) {
  return(ppca.log_likelihood(X, W, sigSq) +
         ifelse(tune.beta, ppca.log_prior(K, W), 0))
}

MR.cov <- function(X, l) {
  if (all(l==0)) return(diag.spam(nrow(X)))
  R = cleanup(spind2spam(fields.rdist.near(X%*%diag(1/l), delta=1, max.points=0.05*nrow(X)^2)))
  K = sigma0*((2+cos(2*pi*R))*(1-R)/3 + sin(2*pi*R)/(2*pi)) + diag.spam(nrow(X))
  return(as.dgCMatrix.spam(K))
}
