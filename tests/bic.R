library(fields)
library(spam)
library(igraph)
library(MASS)

source("../prior-component-analysis.R", chdir=TRUE)
source("synthesize_data.R")

prca.bic <- function(X, k, covar.fn, beta.init) {
  bic = numeric(length(k))
  for (ki in 1:length(k)) {
    DoF = ncol(X)*k[ki] - 0.5*k[ki]*(k[ki]-1) + 3 + length(beta.init) # Degrees of Freedom for PPCA

    out.prca  = prca(X, k[ki], covar.fn, beta.init, trace=0, tol=10,
                     report_iter=1, maxit=15, warnDiag=FALSE)
    bic[ki] = -2*out.prca$lp + DoF*log(nrow(X))
  }
  return(bic)
}

bic.test <- function(X, grid, dim, ks) {
  covar.fn.wend = function(lenscale) {
    if (any(lenscale==0)) return(Diagonal(ncol(X)))
    as.dgCMatrix.spam(wendland.cov(x1=grid%*%diag(1/lenscale)))
  }

  A = get.adjacency(make_lattice(c(dim[1], dim[2])))
  covar.fn.graph = function(gamma) {
    return(Diagonal(ncol(X)) + gamma*A)
  }

  covar.fn.ppca = function() {
    Diagonal(ncol(X))
  }

  bic.wend  = prca.bic(X, ks, covar.fn.wend,  c(0.05, 0.05))
  bic.graph = prca.bic(X, ks, covar.fn.graph, c(0.05))
  bic.ppca  = prca.bic(X, ks, covar.fn.ppca,  c())

  plot(1, type='n', xlim=range(ks), ylim=range(bic.wend, bic.graph, bic.ppca), xlab="k", ylab="bic")
  lines(ks, bic.wend,  col='red')
  lines(ks, bic.graph, col='green')
  lines(ks, bic.ppca,  col='blue')
  lines(c(k, k), c(0, 1), lty=2, col="black")
  legend("top",
         col=c("red", "green", "blue"),
         lty=1,
         legend=c("Wendland", "graph", "ppca"))
}

# Samples:     High
# Dimension:   Low
# Noise:       Low
n        = 100
k        = 7
dim      = c(10, 10)
sigSq    = 0.1
lenScale = c(0.1, 0.1)
synth    = synthesize_data_kern(n, k, dim, lenScale, sigSq)
X        = synth$X
grid     = synth$grid
bic.test(X, grid, dim, 1:9)

# Samples:     Low
# Dimension:   High
# Noise:       Low
n        = 20
dim      = c(100, 100)
synth    = synthesize_data_kern(n, k, dim, lenScale, sigSq)
X        = synth$X
grid     = synth$grid
bic.test(X, grid, dim, 1:9)

# Samples:     Low
# Dimension:   Low
# Noise:       High
dim      = c(10, 10)
sigSq    = 0.3
synth    = synthesize_data_kern(n, k, dim, lenScale, sigSq)
X        = synth$X
grid     = synth$grid
bic.test(X, grid, dim, 1:9)

# Real data
# Samples:     Low
# Dimension:   High
# Noise:       Low
source("~/Documents/uvFAIMS/HE/read_data.R")
data.df = read_data()
X       = scale(do.call(rbind, data.df$vecs)[,((52224/2)+1):52224], scale=FALSE)
n       = nrow(data.faims)
dim     = c(512, 51)
grid    = make.surface.grid(list(x=seq(0, 1, length=dim[1]),
                                 y=seq(0, 1, length=dim[2])))
bic.test(X, grid, dim, 1:9)
