library(sbde) # using for data
library(ggplot2)

tau.pred <- unique(c(seq(0.00001, 0.0001, by=0.00001), 
                          seq(0.0001, 0.001, by=0.0001), 
                          seq(0.001, 0.01, by=0.001), 
                          seq(0.01, 0.99, by=0.01),
                          seq(0.99, 0.999, by= 0.001), 
                          seq(0.999, 0.9999, by=0.0001),
                          seq(0.9999, 0.99999, by=0.00001)))


#####################################################################
# POSITIVE REALS ONLY

# Generate sample data using inverse transform sampling

# GPD quantile and density functions
Q0 <- function(p,nu) {  if(nu==0) { -1*log(1-p) } else{ ((1-p)^(-1*nu) -1)/nu}}
f0 <- function(x, nu) {
  if(nu==0) {val <- exp(-x) } else{ val <-(1 + nu*x)^(-1*(1/nu + 1))}
  if(nu<0) val <- (x <=(-1/nu))*val
  return(val)
}

# Note that sbde only works for fixed support and heavy tails (positive tail parameter)
# generate with tail parameter of 0.2
set.seed(10)
nu. <- 0.2; sigma. <- 1
n <- 1000
y <- Q0(runif(n), nu=nu.)*sigma.
hist(y)

trueQ.pred <-  Q0(tau.pred, nu=nu.)*sigma.
truef <- function(x) {f0(x, nu=nu.)/sigma.}

# Estimate density function
set.seed(20)
fit.gp <- sbde(y, nsamp = 2e3, thin = 10,  nknots = 11, fbase="gpd", incr=0.01, 
               prox.range = qbeta(c(.001,.999),6,4), blocking="single3")

summary(fit.gp, more=TRUE)
fit.gp <- update(fit.gp, nmc=1e3, append=FALSE)
summary(fit.gp, more=TRUE) # converge diagnostics seem to have stabilized

coef.gp <- coef(fit.gp, burn.perc=0,nmc=1000) 
coef.gp$parametric

# Density estimates
par(mfrow=c(1,1))
preds.gp <- predict(fit.gp, burn.perc=0,nmc=1000)
hist(y,freq=FALSE,breaks=40)
matplot(preds.gp$y, preds.gp$fest, type="l", col=1, lty=c(1,2,2),lwd=2, add=T)
curve(truef, add=T, from=0, to=150, col="blue")

# Semi-parametric quantiles
quant.gp <- quantile.sbde(tau.pred, fit.gp, burn.perc = 0.5, nmc=1000)
quant.gp <- data.frame(p=tau.pred, trueQ=trueQ.pred, meanQ=apply(quant.gp$Qsamp,1,mean), 
                         quant.gp$Qest, CIwidth=quant.gp$Qest[,"Up95%"] - quant.gp$Qest[,"Lo95%"])
tail(quant.gp,10) 

# Semi-parametric model fit diagnostics
quantlev.gp <- quantlev.sbde(fit.gp, burn.perc = 0.0, nmc=500)
quantlev.gp <- data.frame(y=fit.gp$y, meanQ=apply(quantlev.gp$qlsamp,1,mean), 
                       quantlev.gp$qlest, CIwidth=quantlev.gp$qlest[,"Up95%"] - quantlev.gp$qlest[,"Lo95%"])

ggplot() + geom_qq(aes(sample=quantlev.gp$meanQ), distribution=stats::qunif) + 
  ylab("actual") + geom_abline(intercept=0, slope=1)

############################################
# FULL REAL LINE

# Create 1 sample dataset and fit lgp to seed covariance matrix in remaining iterations
set.seed(10)
nu. <- 8; sigma. <- 1; gam0. <- 4      # Mode 2, heavy tailed @ right
y <- c(rt(500, df=nu.)*sigma. + gam0., rnorm(500,0,1))  # Mode 1 normal/exponential tail @ left
n <- length(y)
truef <- function(x) {0.5*dt((x-gam0.), df=nu.) + 0.5*dnorm(x,0,1)}
trueF <- function(x) {0.5*pt((x-gam0.), df=nu.) + 0.5*pnorm(x,0,1)}
trueQ.pred <-  sapply(tau.pred, function(prob) {
  uniroot(function(x,p=prob) {trueF(x) - p}, lower=-40, upper=40, tol=2.220446e-16)$root})

set.seed(30)
fit.gp <- sbde(y, nsamp = 1e3, thin = 1,  nknots = 11, fbase="t", incr=0.01, 
               prox.range = qbeta(c(.001,.999),6,4), blocking="single3")
# Note this takes much longer to run than a fbase = "gpd" or logistic because it
# calls the qt (quantile-of-t) function, which is quite slow

summary(fit.gp, more=TRUE) # needs updating - not yet converged
  