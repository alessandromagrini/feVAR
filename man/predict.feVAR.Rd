\name{predict.feVAR}
\alias{predict.feVAR}
\title{Predict method for class feVAR}
\description{Perform predictions based on a fixed effects vector autoregressive model.}
\usage{
\method{predict}{feVAR}(object, newdata=NULL, n.ahead=0, unit.id=NULL, level=0.95, ...)
}
\arguments{
  \item{object}{Object of class \code{feVAR}.}
  \item{newdata}{Object of class \code{data.frame} on which predictions should be based. Ignored if \code{n.ahead} is equal to 0. If \code{NULL} (the default) and \code{n.ahead} is equal to 0, predictions are based on data used for parameter estimation (in-sample predictions).}
  \item{n.ahead}{Integer non-negative number indicating the number of steps ahead of predictions. If equal to 0, in-sample predictions are returned (see argument \code{newdata}).}
  \item{unit.id}{Name or identificative number of the units of observation for which predictions should be made. Ignored in case of a single unit of observation. If \code{NULL} (the default), predictions are performed for all units of observation.}
  \item{level}{Numerical value between 0 and 1 indicating the level of prediction intervals. Default is 0.95.}
  \item{...}{Further parameters to be passed to the generic method \code{predict}.}
}
\value{A named list of data frames with one component for each endogenous variable. Each data frame has several columns including the unit name (only in case of panel data), the time of observation, point predictions, lower and upper bound of prediction intervals.}
%\details{}
%\note{}
%\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
%\references{}
\seealso{\link{feVAR}.}
\examples{
\dontrun{
data(agrisus2020)
varNames <- colnames(agrisus2020)[4:15]
m_agr <- feVAR(varNames, unit="Country", time="Year", exogenous="GDP",
  data=agrisus2020, box.cox=0, ndiff=1, nlags=1)

predict(mod_agr)             ## in-sample prediction
predict(mod_agr, n.ahead=3)  ## 3 steps ahead
}
}
