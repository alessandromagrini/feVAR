\name{predict.feVAR}
\alias{predict.feVAR}
\title{Predict method for class feVAR}
\description{Perform in-sample prediction or h-step ahead forecasts based on an estimated fixed effects vector autoregressive model.}
\usage{
\method{predict}{feVAR}(object, n.ahead=0, unit.id=NULL, newdata=NULL, level=0.95, subset=NULL, ...)
}
\arguments{
  \item{object}{Object of class \code{feVAR}.}
  \item{n.ahead}{The number of steps ahead of predictions. If equal to 0 (the default), in-sample predictions are returned. If an invalid value is provided, it is automatically replaced by the nearest valid one.}
  \item{unit.id}{Name or numerical identifier of the units of observation for which predictions must be made. Ignored in case of a single unit of observation. If \code{NULL} (the default), predictions are performed for all units of observation.}
  \item{newdata}{Object of class \code{data.frame} where to look at future values of exogenous variables. Ignored if \code{n.ahead} is equal to 0. If \code{NULL} (the default) and \code{n.ahead} is greater than 0, exogenous variables are assumed to keep the latest observed value.}
  \item{level}{Confidence level of prediction intervals. Default is 0.95. If outside the interval \code{[0,1]}, it is set to 0.95.}
  \item{subset}{For internal use only.}
  \item{...}{Further parameters to be passed to the generic method \code{predict}.}
}
\value{An object of class \code{predict.feVAR}, which is a list with the following four components:
\itemize{
  \item{\code{predicted}: }{a named list of objects of class \code{data.frame}, one for each endogenous variable. Each data frame includes one record for each observation, and the following variables: unit name (only in case of panel data), time of observation, point predictions, lower and upper bound of prediction intervals;}
  \item{\code{observed}: }{an object of class \code{data.frame} containing the observed data;}
  \item{\code{n.ahead}: }{the number of steps ahead of the predictions;}
  \item{\code{call}: }{component \code{call} of the object of class \code{feVAR} passed to argument \code{object}. See 'Value' of the function \code{\link{feVAR}}.}
  }
}
%\details{}
\note{The class \code{predict.feVAR} has its own \code{plot} method. See \code{\link{plot.predict.feVAR}} for details.}
%\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
%\references{}
\seealso{\code{\link{feVAR}}; \code{\link{plot.predict.feVAR}}.}
\examples{
\dontrun{
data(agrisus2020)
x_agr <- colnames(agrisus2020)[4:15]
m_agr <- feVAR(var.names=x_agr, unit="Country", time="Year", exogenous="GDP",
  data=agrisus2020, box.cox=0, ndiff=1, nlags=1)

predict(mod_agr)             ## in-sample predictions
predict(mod_agr, n.ahead=3)  ## 3 steps ahead forecasts
}
}
