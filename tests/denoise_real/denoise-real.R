library(functional)
library(fields)

source("../../prior-component-analysis.R", chdir=TRUE)
source("../../covariance-functions.R")
source("~/Documents/uvFAIMS/HE/read_data.R")
source("~/Documents/feature_learning/ppca/ppca-simple.R")

source("/home/jim/R/Library/MultiImage.R")

#' Load data
data.df = read_data()
X = scale(do.call(rbind, data.df$vecs)[,((52224/2)+1):52224], scale=FALSE)
#X = X/svd(X)$d[1]

dim.faims = c(512, 51)
grid = expand.grid(list(x=seq(0, 1, length=dim.faims[1]),
                        y=seq(0, 1, length=dim.faims[2])))

keep = apply(X, 2, sd) > 0.0077
X = X[,keep]
grid = grid[keep,]

#' Use BIC to pick latent dimensionality w/ PPCA (to make test most fair);
#' pick k=9.
bics = numeric(14)
for (k in 1:14) {
  out.ppca  = ppca(X, k)
  bics[k] = out.ppca$bic
}

plot(bics)

#' Add noise to data and de-noise. Use k=8 because BIC picked this
k = 9
nLevels = 6
noiseLevels = seq(0, 0.055*(max(X)-min(X)), length=nLevels)
reconstErrs = data.frame(pca=numeric(length(noiseLevels)),
                         ppca=numeric(length(noiseLevels)),
                         prca=numeric(length(noiseLevels)))

for (nli in length(noiseLevels)) {
  X.noisy = X + matrix(rnorm(nrow(X)*ncol(X), sd=noiseLevels[nli]), nrow=nrow(X), ncol=ncol(X))

  print(paste("svd", nli))
  out.svd  = svd(X.noisy)

  print(paste("ppca", nli))
  out.ppca  = ppca(X.noisy, k)

  print(paste("prca", nli))
  out.prca = prca(X.noisy, k, grid, exp.MR.cov, exp.MR.cov.d,
                  beta0=log(c(mean(X.noisy^2), 0.5)), maxit=8,
                  report_iter=1, max.dist=3, trace=1, ucminf.control=list(
                    trace=0, grtol=1e-2, xtol=1e-2, maxeval=4
                  ))

  rec.svd  = out.svd$u[,1:k] %*% diag(out.svd$d[1:k]) %*% t(out.svd$v)[1:k,]
  rec.ppca = tcrossprod(out.ppca$V, out.ppca$W)
  rec.prca = tcrossprod(out.prca$V, out.prca$W)

  #par(mfcol=c(1,3))
  #plot.surface(as.surface(grid, rec.svd[1,]), type='I', main="svd")
  #plot.surface(as.surface(grid, rec.ppca[1,]), type='I', main="ppca")
  #plot.surface(as.surface(grid, rec.prca[1,]), type='I', main="prca")

  reconstErrs[nli, c("pca", "ppca", "prca")] = c(
    norm(X - rec.svd, 'F'),
    norm(X - rec.ppca, 'F'),
    norm(X - rec.prca, 'F')
  )
}

#' Plot reconstruction error vs noise
plot(1, type='n', ylim=range(reconstErrs), xlim=range(noiseLevels), ylab="reconstruction error", xlab="added noise sd")
lines(noiseLevels, reconstErrs$pca, col='red')
lines(noiseLevels, reconstErrs$ppca, col='blue')
lines(noiseLevels, reconstErrs$prca, col='green')
legend('topleft',
       legend=c('pca', 'ppca', 'prca'),
       col=c('red', 'blue', 'green'),
       lty=1)
