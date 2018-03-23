\name{waic}
\Rdversion{0.1-1}
\alias{waic}
\title{Watanabe Information Criterion}
\description{Calculates two versions of the Watanabe information criteria from MCMC draws.} 
\usage{ 
waic(logliks, print = TRUE)
} 
\arguments{
\item{logliks}{a matrix of observation level log-likelihood values, the columns are MCMC iterations and the rows are observations in the data}
\item{print}{logical whether to print the results}
}
\value{
Returns the two version of the WAIC
}

\seealso{\code{\link{summary.qrjoint}}}

\references{Gelman, A., Hwang, J., and Vehtari, A. (2014). Understanding predictive information criterion for Bayesian models. \emph{Stat Comput}, 24, 997-1016.}

\examples{ 
# Plasma data analysis

# recoding variables
data(plasma)
plasma$Sex <- as.factor(plasma$Sex)
plasma$SmokStat <- as.factor(plasma$SmokStat)
plasma$VitUse <- 3 - plasma$VitUse
plasma$VitUse <- as.factor(plasma$VitUse)

# creating predictors and response (beta carotene concentration in the plasma)
X <- model.matrix(BetaPlasma ~ Age + Sex + SmokStat + Quetelet + VitUse + Calories + 
        Fat + Fiber + Alcohol + Cholesterol + BetaDiet, data = plasma)[,-1]
Y <- plasma$BetaPlasma

# model fitting with 50 posterior samples from 100 iterations (thin = 2)
fit.qrj <- qrjoint(X, Y, 50, 2)
sm <- summary(fit.qrj, plot = FALSE)

# the call to summary already shows the waic for the fitted model, it also returns 
# the observation level log-likelihood vales. To calculate waic from last 20 draws 
# we can use:

ll <- sm$ll
waic(ll[,31:50])
}



\keyword{programming}