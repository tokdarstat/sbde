require(sbde)
require(evd)

n <- 1000
ix <- (runif(n) < 0.5)
y <- rep(NA, n)
y[ix] <- rnorm(sum(ix), -4, 1)
y[!ix] <- 4 + 2 * rt(sum(!ix), df = 3.0)

o <- sbde(y, nknots = 11, nsamp = 1000, thin = 30) 
summary(o, more = TRUE)
print(coef(o))
pp <- predict(o)

par(mfrow = c(1,1))
hist(y, 200, freq = FALSE, col = "gray", border = "white")
for(i in 1:3) lines(pp$y, pp$fest[,i], lty = 1 + (i != 2))
