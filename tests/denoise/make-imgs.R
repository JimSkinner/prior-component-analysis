library(fields)
# TODO: Deal with margins

# Plot some samples
pdf("img/samples.pdf", width=4, height=4)
par(mfcol=c(2,2), mar=c(1,1,1,1))
plot.surface(as.surface(grid, X[1,]), type='I', xlab='', ylab='', xaxt='n', yaxt='n')
plot.surface(as.surface(grid, X[2,]), type='I', xlab='', ylab='', xaxt='n', yaxt='n')
plot.surface(as.surface(grid, X[3,]), type='I', xlab='', ylab='', xaxt='n', yaxt='n')
plot.surface(as.surface(grid, X[4,]), type='I', xlab='', ylab='', xaxt='n', yaxt='n')
dev.off()

#' Look at the noisiest setting. Plot some samples.
pdf("img/noisy-samples.pdf", width=4, height=4)
par(mfcol=c(2,2), mar=c(1,1,1,1))
plot.surface(as.surface(grid, X.noisy[1,]), type='I', xlab='', ylab='', xaxt='n', yaxt='n')
plot.surface(as.surface(grid, X.noisy[2,]), type='I', xlab='', ylab='', xaxt='n', yaxt='n')
plot.surface(as.surface(grid, X.noisy[3,]), type='I', xlab='', ylab='', xaxt='n', yaxt='n')
plot.surface(as.surface(grid, X.noisy[4,]), type='I', xlab='', ylab='', xaxt='n', yaxt='n')
dev.off()

##' Compare real PCs vs ones reconstructed from PCA & SPPCA
pdf("img/PCs.pdf", width=4, height=4)
par(mfcol=c(2,2), mar=c(1,1,1,1))
plot.surface(as.surface(grid, synth$W[,1]), type='I', xlab='', ylab='', xaxt='n', yaxt='n')
plot.surface(as.surface(grid, out.svd$v[,1]), type='I', xlab='', ylab='', xaxt='n', yaxt='n')
plot.surface(as.surface(grid, out.ppca$W[,1]), type='I', xlab='', ylab='', xaxt='n', yaxt='n')
plot.surface(as.surface(grid, out.prca2$W[,1]), type='I', xlab='', ylab='', xaxt='n', yaxt='n')
dev.off()

#' Plot reconsturction error
pdf("img/rec-error.pdf", width=4, height=4)
plot(1, type='n', ylim=range(reconstErrs), xlim=range(noiseLevels),
     ylab="reconstruction error", xlab="added noise sd")
lines(noiseLevels, reconstErrs$pca,  col='red', type='b')
lines(noiseLevels, reconstErrs$ppca, col='green', type='b')
lines(noiseLevels, reconstErrs$prca2, col='purple', type='b')
legend('topleft',
       legend=c('pca', 'ppca', 'spca'),
       col=c('red', 'green', 'purple'),
       lty=1)
dev.off()

#' Plot reconsturction error including crap result
pdf("img/rec-error-poor.pdf", width=4, height=4)
plot(1, type='n', ylim=range(reconstErrs), xlim=range(noiseLevels),
     ylab="reconstruction error", xlab="added noise sd")
lines(noiseLevels, reconstErrs$pca,   col='red', type='b')
lines(noiseLevels, reconstErrs$ppca,  col='green', type='b')
lines(noiseLevels, reconstErrs$prca2, col='purple', type='b')
lines(noiseLevels, reconstErrs$prca,  col='blue', type='b')
legend('topleft',
       legend=c('pca', 'ppca', 'spca', 'spca-wendland'),
       col=c('red', 'green', 'purple', 'blue'),
       lty=1)
dev.off()
