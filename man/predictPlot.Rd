\name{predictPlot}
\alias{predictPlot}
\title{Graphics of predictions}
\description{Graphic of in-sample or out-of-sample predictions for a fixed effects vector autoregressive models.}
\usage{predictPlot(model, unit.id, var.names=NULL, n.ahead=0, newdata=NULL,
  level=0.95, xlim=NULL, ylim=NULL, add.grid=TRUE, cex.axis=NULL, cex.points=0.6,
  xlab=NULL, ylab=NULL, titles=NULL, fit.col="dodgerblue", out.col="red",
  obs.col="grey40", new.col="grey40", fit.lty=1, out.lty=1, obs.lty=1, new.lty=2,
interval.col="grey70", las=NULL, mar=c(3.5,3.5,2,2), mgp=c(2.3,0.8,0), ...)}
\arguments{
  \item{model}{An object of class \code{feVAR}.}
  \item{unit.id}{Name or identificative number of one unit of observation for which predictions should be displayed. Ignored in case of a single unit of observation.}
  \item{var.names}{Charcater vector indicating the name of endogenous variables for which prediction should be displayed.}
  \item{n.ahead}{Integer non-negative number indicating the number of steps ahead of predictions. If equal to 0, in-sample predictions are returned.}
  \item{newdata}{Object of class \code{data.frame} containing new observations to be added to the graphics. If \code{NULL} (the default), no further observations are added. Ignored if \code{n.ahead} is equal to 0.}
  \item{level}{Numeric value between 0 and 1 indicating the level of confidence intervals. If equal to 0, confidence intervals are not provided. Default is 0.95.}
  \item{xlim}{Numeric vector of length 2 indicating the limits of the x-axis.}
  \item{ylim}{Numeric vector of length 2 indicating the limits of the y-axis.}
  \item{add.grid}{Logical value indicating whether the grid should be added. Default is \code{TRUE}.}
  \item{cex.axis}{Numeric vector of length 2 indicating the expansion factor for x-axis and y-axis labels, respectively. Default is (1,1).}
  \item{cex.points}{Numeric value indicating the expansion factor for points. If 0, points are suppressed. Default is 0.6.}
  \item{xlab}{Label for the x-axis.}
  \item{ylab}{Label for the y-axis.}
  \item{titles}{Character vector containing the title for each graphic. If \code{NULL} (the default), the names of endogenous variables are used as titles.}
  \item{fit.col}{Color for the time series of fitted values. Default is black.}
  \item{out.col}{Color for the time series of predictions. Default is red.}
  \item{obs.col}{Color for the time series of observed data. Default is black.}
  \item{new.col}{Color for the time series of new data. Default is black.}
  \item{fit.lty}{Line type for the time series of fitted values. Default is 1.}
  \item{out.lty}{Line type for the time series of predictions. Default is 1.}
  \item{obs.lty}{Line type for the time series of observed data. Default is 1.}
  \item{new.lty}{Line type for the time series of new data. Default is 1.}
  \item{interval.col}{Color of prediction intervals. Default is \code{grey40}.}
  \item{las}{Numeric value indicating the orientation of the tick mark labels. Possible values are: parallel to the axis (0), horizontal (1), perpendicular to the axis (2), vertical (3). If \code{NULL} (the default), it is set automatically based on the length of axis labels.}
  \item{mar}{Numeric vector of length 4 indicating the margin sizes in the following order: bottom, left, top, right. Default is c(3.5,3.5,2,2).}
  \item{mgp}{Numeric vector of length 3 indicating the location of axis labels, tick mark labels and tick marks, respectively, relative to the plot window. Default is c(2.3,0.8,0).}
  \item{...}{Further parameters to be passed to the generic method \code{plot}.}
}
%\value{}
%\details{}
%\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
%\references{}
\seealso{\link{feVAR}; \link{predict.feVAR}.}
\examples{
\dontrun{
data(agrisus2020)
varNames <- colnames(agrisus2020)[4:15]
m_agr <- feVAR(varNames, unit="Country", time="Year", exogenous="GDP",
  data=agrisus2020, box.cox=0, ndiff=1, nlags=1)

# predictions for the first four variables on unit 'Italy'
predictPlot(m_agr, unit.id="Italy", var.names=varNames[1:4])  ## in-sample
predictPlot(m_agr, unit.id=15, var.names=1:4)                 ## same as before
predictPlot(m_agr, unit.id=15, var.names=1:4, n.ahead=3)      ## 3 steps ahead
}
}