\name{sbde}
\Rdversion{0.0-1}
\alias{sbde}
\alias{update.sbde}
\title{Bayesian Semiparametric Density Estimation}
\description{Provides a semiparametric estimation of the density function of independent univariate data.
}
\usage{
sbde(y, nsamp = 1e3, thin = 10, cens = rep(0,length(y)),
    wt = rep(1,length(y)), incr = 0.01, par = "prior", nknots = 6,
    hyper = list(sig = c(.1,.1), lam = c(6,4), kap = c(0.1,0.1,1)),
    prox.range = c(.2,.95), acpt.target = 0.15, ref.size = 3,
    blocking = "single", temp = 1, expo = 2, blocks.mu, blocks.S,
    fix.nu = FALSE, fbase = c("t", "t+", "gpd", "unif"), verbose = TRUE)


\method{update}{sbde}(object, nadd, append = TRUE, ...)
}
\arguments{
\item{y}{numeric vector of response data.}
\item{nsamp}{number of posterior samples to be saved; defaults to 1000.}
\item{thin}{thinning rate for the Markov chain sampler -- one posterior sample is saved per \code{thin} iterations. Defaults to 10. The Markov chain sampler runs for a total of \code{nsamp * thin} many iterations.}
\item{cens}{censoring status of response. Must be a vector of 0s and 1s of length same as length(y), with 1 indicating right censoring, and, 0 indicating no censoring. Defaults to all zeros.}
\item{wt}{weights attached to the observation units, expected to be non-negative numbers, and defaults to a vector of ones.}
\item{incr}{tau grid increment. Defaults to 0.01.}
\item{par}{character string indicating how the sampler is to be initialized. Only two options are currently supported: "prior" to initialize at a random draw from the prior; "RQ" to initialize at a model space approximation of the estimates from \code{\link{rq}}.}
\item{nknots}{number of knots to be used for low rank approximation of the Gaussian process priors. Defaults to 6.}
\item{hyper}{hyperparameters of the prior distribution. Must be a list with some of all of the following fields: \code{sig}: a two vector giving the parameters of the inverse-gamma distribution on sigma-square that is used when \code{shrink}=TRUE, \code{lam}: a two vector giving the parameters of the beta distribution on proximity = \eqn{\exp(-0.01* \lambda^2)}{exp(-0.01 * lambda^2)}, and \code{kap}:  a vector to be coerced into a \code{3 * nkap} matrix, with \code{nkap} being the number of components in the mixture of gamma prior on \code{kappa}, and each column of the matrix gives the shape, rate and mixing weight of a component.}
\item{prox.range}{for specifying the range of length-scale parameter of the Gaussian process prior.}
\item{acpt.target}{target acceptance rate of the adaptive Metropolis sampler; defaults to 0.15}
\item{ref.size}{adaptation rate of the adaptive Metropolis sampler. The proposal density is updated once every \code{ref.size} iterations. Could be a single number or a vector of length same as the number of blocks.}
\item{blocking}{type of blocking to be applied. Either a character string specifying one to be chosen from the supplied menu (see Details), or a list giving user specified blocks. In the latter case, each element of the list is a logical vector of length equal to the total number of model parameters, which equals \code{nknots+3} indicating which model parameters belong to the block.}
\item{temp}{temperature of the log-likelihood function. The log-likelihood function is raised to the power of \code{temp}. Defaults to 1.}
\item{expo}{the exponent to be used in the covariance kernel of the Gaussian process priors. Defaults to 2, giving the standard squared-exponential covariance kernel.}
\item{blocks.mu}{initial block specific means in the form of  a list. If left unspecified then will be automatically generated as a list of vectors of zeros of appropriate lengths matching the corresponding block sizes.}
\item{blocks.S}{initial block specific covariance matrices in the form of a list. If left unspecified then will be automatically generated as a list of identity matrices of appropriate dimensions matching the corresponding block sizes. When \code{blocking} is chosen as one of the menu items of the form "std*", known prior covariance information and estimated variance matrices from \code{\link{rq}} are used.}
\item{fix.nu}{either the logical FALSE indicating that nu should be learned, or a positive real number giving the fixed value of nu, which is then excluded from MCMC updates}
\item{fbase}{either "t" (default) or "t+" (for half-t distributions on the positive real lines) or "gpd" (for generalized pareto distributions with location zero and parametrized by nu = 1 / shape) or  "unif" to indicate what base distribution is to be used.}
\item{verbose}{logical indicating whether MCMC progress should be printed, defaults to TRUE}
\item{object}{a fitted model of the class 'qde'.}
\item{nadd}{number of additional MCMC samples.}
\item{append}{logical indicating whether new samples should be appended to old ones. If FALSE then old samples are discarded.}
\item{...}{no additional arguments are allowed}
}
\value{
\code{sbde(y, ...)} returns a `sbde' class object to be used by \code{\link{coef}} and \code{\link{summary}}.

Returned object is a list containing the following variables.

\item{par}{latest draw of the parameter vector}
\item{y}{response vector}
\item{cens}{censoring status vector}
\item{wt}{vector of observation weights}
\item{hyper}{completed list of hyper-parameters}
\item{dim}{model dimension vector of the form c(n, p, length of tau grid, position of \eqn{\tau_0}{tau0} on the grid, nknots, length of lambda grid, nkap, total number of MCMC iterations, thin, nsamp)}
\item{gridmats}{details of covariance matrix factors etc, intended for internal use.}
\item{tau.g}{the tau grid}
\item{muV}{list of means for parameter blocks}
\item{SV}{list of covariance matrices for parameter blocks}
\item{blocks}{list of blocks}
\item{blocks.size}{vector of block lengths}
\item{dmcmcpar}{numeric vector containing details of adaptive MCMC runs, equals c(temp, decay rate of adaptation, vector of target acceptance rates for the blocks, vector of increment scales used in adaptation). Intended strictly for internal use.}
\item{imcmcpar}{numeric vector containing details of adaptive MCMC runs, equals c(number of parameter blocks, ref.size, indicator on whether details are to be printed during MCMC progress, rate of details printing, a vector of counters needed for printing). Intended strictly for internal use.}
\item{parsamp}{a long vector containing the parameter draws. Could be coerced into a matrix of dim \code{npar * nsamp}. Intended primarily for use by \code{\link{summary}} and \code{\link{coef}}.}
\item{acptsamp}{a long vector containing rates of acceptance statistics for parameter blocks. Could be coerced into a matrix of dim \code{nblocks * nsamp}. Not very informative, because thinning times and adaptation times may not be exactly synced.}
\item{lpsamp}{vector of log posterior values for the saved MCMC draws.}
\item{fbase.choice}{integer 1 for "t" and 2 for "unif" base.}
\item{prox}{vector of proximity (exp(-0.01*lambda^2)) grid values}
\item{reg.ix}{positions of the regular tau grid on the expanded tail-appended grid}
\item{runtime}{run time of the MCMC}
}

\details{
SOMETHING about the model.

In running the MCMC, the following menu choices are available for blocking the parameter vector. For this special case \code{p = ncol(X) = 0}, some of the menu choices are actually the same, in particular, "std0" is same as "single", "std1" is same as "single2", and, "std2" is same as "single3".

\code{"single"}: a single block containing all parameters

\code{"single2"}: one block containing all parameters and an additional block containing only (gamma[0], gamma, sigma, nu)

\code{"single3"}: like \code{"single2"}, but the second block is split into two further blocks, one with \eqn{(\gamma_0, \gamma)}{(gamma[0], gamma)}, the other with \eqn{(\sigma, \nu)}{(sigma, nu)}

\code{"std0"}: Same as \code{"single"}.

\code{"std1"}: Same as \code{"single2"}.

\code{"std2"}: Same as \code{"single3"}.

\code{"std3"}: total 3 blocks. First block is \eqn{W_{*0}}{W*[0]}, last two are \eqn{(\gamma_0, \gamma)}{(gamma[0], gamma)} and \eqn{(\sigma, \nu)}{(sigma, nu)}

\code{"std4"}: total 3 blocks. First block is \eqn{(W_{*0}, \gamma_0)}{(W*[0], gamma[0])}, last two are \eqn{(\gamma_0, \gamma)}{(gamma[0], gamma)} and \eqn{(\sigma, \nu)}{(sigma, nu)}

\code{"std5"}: total 4 blocks. First three are same as \code{"std4"} and one additional block containing all parameters.
}
\references{Some references.}
\seealso{\code{\link{summary.sbde}}, \code{\link{coef.sbde}} and \code{\link{predict.sbde}}.}
\examples{
## Plasma data analysis

data(plasma)
}
\keyword{programming}

