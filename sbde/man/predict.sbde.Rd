\name{predict.sbde}
\Rdversion{0.0-1}
\alias{predict.sbde}
\title{Posterior predictive Summary for Quantile based Density Estimation}
\description{Extract posterior predictive density estimate for \code{\link{sbde}}}
\usage{ 
\method{predict}{sbde}(object, burn.perc = 0.5, nmc = 200, yRange = range(object$y), yLength = 401, ...)
} 
\arguments{
    \item{object}{a fitted model of the class 'sbde'.}
    \item{burn.perc}{a positive fraction indicating what fraction of the saved draws are to be discarded as burn-in}
    \item{nmc}{integer giving the number of samples, post burn-in, to be used in Monte Carlo averaging}
    \item{yRange}{Range of values over which posterior predictive density is to be evaluated.}
    \item{yLength}{Number of grid points spanning yRange for posterior predictive density evaluation.}
    \item{...}{no additional parameters are used.}
}
\value{
Returns a list with three items:
\item{y}{vector giving the grid over which the posterior predictive density is evaluated.}
\item{fsamp}{a matrix with \code{yLength} many rows and \code{nmc} many columns. Each column corresponds to a draw of the response density from the posterior predictive.}
\item{fest}{summary of the posterior predictive density given by point-wise median, 2.5th and 97.5th percentiles.}
}

\seealso{\code{\link{sbde}} and \code{\link{summary.sbde}}.}
\examples{ 
# Plasma data analysis

data(plasma)
Y <- plasma$BetaPlasma
Y <- Y + 0.1 * rnorm(length(Y)) ## remove atomicity

# model fitting with 50 posterior samples from 100 iterations (thin = 2)
fit.sbde <- sbde(Y, 50, 2)
pred <- predict(fit.sbde)
hist(Y, freq = FALSE, col = "gray", border = "white", ylim = c(0, max(pred$fest)))
lines(pred$y, pred$fest[,2])
lines(pred$y, pred$fest[,1], lty = 2)
lines(pred$y, pred$fest[,3], lty = 2)
}
\keyword{programming}
