library(functional)
library(fields)

lenScale = 0.03
covar.real = function(x) Exp.cov(0, x/sqrt(lenScale), p=2)
covar.wend = function(x, supp) c(wendland.cov(0, x/supp))

x = seq(0, 20, length=1000)
min.f = function(supp) {
  vapply(supp, function(s) {
    mean((covar.real(x) - covar.wend(x, s))^2)
  }, numeric(1))
}

optObj = optimize(min.f, c(0.1, 10))
supp.optimum = optObj$minimum

par(mfcol=c(1,2))

curve(min.f, 0.1, 10)
abline(v=supp.optimum, col='blue', lty=3)

curve(covar.real, 0, 0.5, xlab="d", ylab="covariance", col="blue")
lines(seq(0, 0.5, length=1000), covar.wend(seq(0, 0.5, length=1000), supp.optimum), col='red')
legend("topright",
       col=c("blue", "red"),
       lty=1,
       legend=c("SE(d)", "Wendland(d/0.4572787)"))
