\name{plasma}
\alias{plasma}
\docType{data}
\title{
Plasma Concentration of Beta-Carotene and Retinol}
\description{
Plasma concentration levels of beta-carotene and retinol along with dietary intake and drug usage measurements for 315 patients who had an elective surgical procedure during a three-year period to biopsy or remove a lesion of the lung, colon, breast, skin, ovary or uterus that was found to be non-cancerous.
}
\usage{data(plasma)}
\format{
  A data frame with 315 observations on the following 14 variables.
  \describe{
    \item{\code{Age}}{age (years)}
    \item{\code{Sex}}{sex (\code{1}=Male, \code{2}=Female)}
    \item{\code{SmokStat}}{smoking status (\code{1}=Never, \code{2}=Former, \code{3}=Current)}
    \item{\code{Quetelet}}{Quetlet index, aka, BMI (weight / height^2)}
    \item{\code{VitUse}}{vitamin use (\code{0}=No, \code{1}=Yes, not often, \code{2}=Yes, fairly often)}
    \item{\code{Calories}}{number of calories consumed per day}
    \item{\code{Fat}}{grams of fat consumed per day}
    \item{\code{Fiber}}{grams of fiber consumed per day}
    \item{\code{Alcohol}}{number of alcoholic drinks consumed per week}
    \item{\code{Cholesterol}}{cholesterol consumed (mg per day)}
    \item{\code{BetaDiet}}{dietary beta-carotene consumed (mcg per day)}
    \item{\code{RetDiet}}{dietary retinol consumed (mcg per day) }
    \item{\code{BetaPlasma}}{plasma beta-carotene concentration (ng/ml)}
    \item{\code{RetPlasma}}{plasma retinol concentration (ng/ml)}
  }
}

\details{
Dietary intakes are self-reported. Results from analyzing this data should be used with caution!
}
\source{Statlib database}
\references{
Nierenberg, D. W., T. A. Stukel, J. A. Baron, B. J. Dain, and E. R. Greenberg (1989). Determinants of plasma levels of beta-carotene and retinol. \emph{American Journal of Epidemiology}, 130(3), 511--521.
}
\examples{
data(plasma)
}
\keyword{datasets}
