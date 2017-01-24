library(fields)

source("../../prior-component-analysis.R", chdir=TRUE)
source("../../covariance-functions.R")

source("/home/jim/Documents/feature_learning/ppca/ppca-simple.R")
source("/home/jim/R/Library/MultiImage.R")

data(ozone2)

#' Missing data: Replace missing values by column means
X = ozone2$y
for (col in 1:ncol(X)) {
  ind = is.na(X[,col])
  X[ind,col] = mean(X[,col], na.rm=TRUE)
}

#' Plot a few samples
#+ plot-samples, fig.asp=0.22
locations = ozone2$lon.lat
locations = apply(locations, 2, function(col) (col-min(col))/(max(col)-min(col)))
multiImage(locations, X[1,], X[2,], X[3,], X[4,])

#' Look at BIC for SPCA with and for PPCA. We can see that SPCA picks k=3, PPCA picks k=9 and SPCA always has a higher BIC. This is actually expected, since both SPCA and PPCA will have (approximately) the same model complexity penalty, but PPCA will have a greater model fit.
source("../../prior-component-analysis.R", chdir=TRUE)
source("../../covariance-functions.R", chdir=TRUE)

source("/home/jim/Documents/feature_learning/ppca/ppca-simple.R")

ks = 1:12

bics.spca = numeric(length(ks))
bics.ppca = numeric(length(ks))
for (ki in seq_along(ks)) {
  out.spca = prca(X, ks[ki], locations, exp.MR.cov, exp.MR.cov.d,
                  beta0=log(c(0.1, 0.5)), maxit=15, max.dist=0.5,
                  trace=1, report_iter=15)
  bics.spca[ki] = out.spca$bic

  out.ppca = ppca(X, ks[ki])
  bics.ppca[ki] = out.ppca$bic
}

plot(ks, bics.spca, xlab='k', ylab='bic', type='b', main="BIC",
     ylim=range(bics.spca, bics.ppca))
points(ks, bics.ppca, col='green', type='b')
legend("top",
       legend=c("SPCA", "PPCA"),
       col=c("black", "green"),
       pch=1, lty=1)

#' Try looking at crossvalidated log likelihoods instead, as a better measure of model fit. Better! SPCA always has a greater holdout data likelihood. Interesting that they peak at the same place (k=7).
liks.spca = numeric(length(ks))
liks.ppca = numeric(length(ks))
for (ki in seq_along(ks)) {
  trainInd = sample(nrow(X), floor(nrow(X)*0.8))

  trainX = X[trainInd,]

  out.spca = prca(trainX, ks[ki], locations, exp.MR.cov, exp.MR.cov.d,
                  beta0=log(c(0.1, 0.5)), maxit=15, max.dist=0.5,
                  trace=1, report_iter=15)
  liks.spca[ki] = prca.log_likelihood(X[-trainInd,], out.spca$W, out.spca$sigSq)

  out.ppca = ppca(X, ks[ki])
  liks.ppca[ki] = prca.log_likelihood(X[-trainInd,], out.ppca$W, out.ppca$sigSq)
}

plot(ks, liks.spca, xlab='k', ylab='bic', type='b',
     main="Heldout log likelihoods", ylim=range(liks.spca, liks.ppca))
points(ks, liks.ppca, col='green', type='b')
legend("top",
       legend=c("SPCA", "PPCA"),
       col=c("black", "green"),
       pch=1, lty=1)

#' Rerun w/ k=7
out.spca = prca(X, 7, locations, exp.MR.cov, exp.MR.cov.d,
                beta0=log(c(0.1, 0.5)), maxit=10, max.dist=0.5,
                trace=1, report_iter=5)

#' We have learnt quite a long length scale; a good sign!
dist = matrix(seq(0, 1, length=100), ncol=1)
covar = exp.MR.cov(X=matrix(0), X2=dist, beta=out.spca$beta, max.dist=0.5)

plot(dist, covar, type='l', main='Learned covariance function')

##' Another PrCA/PPCA comparison; cross-validate and project test-samples onto model. Compare distributions of distances to original points (generalisation error)
#n        = nrow(X)
#err.wend = numeric(n)
#err.ppca = numeric(n)
#for (fold in 1:n) {
#  # PrCA
#  out.wend = prca(X[-fold,,drop=FALSE], 7, covar.fn, beta.init=c(1, log(0.3)),
#                  maxit=10)
#  rec.wend = predict(out.wend, X[fold,,drop=FALSE])
#  err.wend[fold] = sqrt(sum(rec.wend - X[fold,,drop=FALSE])^2)
#
#  # PPCA
#  out.ppca = prca(X[-fold,,drop=FALSE], 3, covar.fn.ppca, beta.init=c(),
#                  maxit=10, warnDiag=FALSE)
#  rec.ppca = predict(out.ppca, X[fold,,drop=FALSE])
#  err.ppca[fold] = sqrt(sum(rec.ppca - X[fold,,drop=FALSE])^2)
#}
#
## Histogram Colored (blue and red)
#hist(err.wend, 100, col=rgb(1,0,0,0.5), xlim=c(0, max(err.wend, err.ppca)),
#     ylim=c(0,45), xlab="Reconstruction error")
#hist(err.ppca, 100, col=rgb(0,0,1,0.5), add=T)
#box()
#
#boxplot(err.wend, err.ppca, names=c("PrCA", "PPCA"))
