\name{coef.qde}
\Rdversion{1.0-5}
\alias{coef.qde}
\title{Coefficient Extraction from qde Model Fit}
\description{Post process MCMC output from \code{\link{qde}} to create summaries of the quantile estimates.}
\usage{ 
\method{coef}{qde}(object, burn.perc = 0.5, nmc = 200, reduce = TRUE, ...)
} 
\arguments{
\item{object}{a fitted model of the class \code{qde}.}
\item{burn.perc}{a positive fraction indicating what fraction of the saved draws are to be discarded as burn-in}
\item{nmc}{integer giving the number of samples, post burn-in, to be used in Monte Carlo averaging}
\item{reduce}{logical indicating if the tail-expanded grid of tau values is to be reduced to the regular increment grid}
\item{...}{not currently implemented}
}
\value{
Extracts posterior draws of intercept and slope functions from saved draws of model parameters. A plot may be obtained if \code{plot = TRUE}. A list is returned invisibly with three fields.

\item{beta.samp}{a matrix with \code{nmc} many columns and \code{length(tau.grid)} many rows.}
\item{beta.est}{a 3-column matrix of median, 2.5th and 97.5th percentiles of the posterior distribution of \eqn{\beta_0}{{beta}[0]}}
\item{parametric}{a matrix with 3 columns giving the estimate (posterior median) and the lower and upper end points of the 95\% posterior credible interval for \eqn{\gamma_0}{{gamma}[0]}, \eqn{\sigma}{sigma}, and, \eqn{\nu}{nu}.}
}

\seealso{\code{\link{qde}} and \code{\link{summary.qde}} for model fitting under qrjoint. Also see \code{\link{getBands}} for plotting credible bands for coefficients.}
\examples{
## Plasma data analysis
data(plasma)
Y <- plasma$BetaPlasma
Y <- Y + 0.1 * rnorm(length(Y)) ## remove atomicity

# model fitting with 50 posterior samples from 100 iterations (thin = 2)
fit.qde <- qde(Y, 50, 2)
betas <- coef(fit.qde)
signif(betas$parametric, 3)
}
\keyword{programming}
