\name{getBands}
\Rdversion{0.1-1}
\alias{getBands}
\title{Posterior Credible Bands}
\description{Calculate and display credible bands of a vector of parameters from a sample of draws. Most suitable when the vector represents a discretized version of a function.} 
\usage{ 
getBands(b, col = 2, lwd = 1, plot = TRUE, add = FALSE, 
x = seq(0,1,len=nrow(b)), remove.edges = TRUE, ...)
} 
\arguments{
\item{b}{a matrix of sampled draws of a vector, columns giving samples and rows giving elements of the vector}
\item{col}{color of the median line and 95\% bands, usual color codes could be used}
\item{lwd}{line width for the median line}
\item{plot}{logical indicating whether plots are to be drawn, default is \code{TRUE}}
\item{add}{logical indicating whether plot is to be added to existing plot, default is \code{FALSE}}
\item{x}{indexing the parameter coordinates. When \code{b} represents a (discretized) function, \code{x} can be taken as the function argument values. Needed when plot is to be created. Default sets it to a uniform grid of points over [0,1].}
\item{remove.edges}{logical indicating whether the first and last entries of \code{b} are to be removed from plotting. This is helpful in qrjoint plots, where the two extremes could be \code{Inf} or nearly \code{Inf}.}
\item{...}{limited number of plotting parameters}
}
\value{
returns median, 2.5th and 97.5th percentiles as a 3-column matrix.
}

\seealso{\code{\link{coef.qrjoint}}}

\examples{ 
## toy example

x <- 2*pi*seq(0,1,.01)
fsamp <- replicate(100, rnorm(1,0,0.1) + rnorm(1,1,0.2) * cos(x))
getBands(fsamp)
getBands(fsamp, x = x, col = 3, remove.edges = FALSE, xlab = "x", ylab = "f", bty = "n")
getBands(fsamp, x = x, col = "darkgreen", remove.edges = FALSE, xlab = "x", ylab = "f")
}



\keyword{programming}