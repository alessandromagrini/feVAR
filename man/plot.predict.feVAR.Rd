\name{plot.predict.feVAR}
\alias{plot.predict.feVAR}
\title{Plot method for class predict.feVAR}
\description{Display in-sample predictions or h-step ahead forecasts based on an estimated fixed effects vector autoregressive model.}
\usage{
\method{plot}{predict.feVAR}(x, var.names=NULL, unit.id, newdata=NULL, start=NULL, ylim=NULL, add.grid=TRUE,
  las=0, cex.axis=c(1,1), cex.lab=c(1,1), xlab=NULL, ylab=NULL, titles=NULL,
  obs.col="grey40", fit.col="dodgerblue", out.col="red", new.col="grey40",
  interval.col="grey70", obs.lty=1, fit.lty=5, out.lty=1, new.lty=2, obs.cex=0.4,
  new.cex=0.4, mfrow=NULL, mar=c(3.5,3.5,2,2), mgp=c(2.3,0.8,0), ...)}
\arguments{
  \item{x}{Object of class \code{predict.feVAR}.}
  \item{var.names}{Vector indicating the names or the numerical identifiers of the endogenous variables for which predictions must be displayed.}
  \item{unit.id}{Name or numerical identifier of a single cross-sectional unit for which predictions must be displayed. Ignored if there is only one cross-sectional unit.}
  \item{newdata}{Object of class \code{data.frame} containing new observations to be added to the graphics. If \code{NULL} (the default), no further observations are added besides the existing ones. Ignored if \code{x} contains in-sample predictions.}
  \item{start}{Time point from which the graphic must start. If \code{NULL} (the default), the whole observed time series is displayed.}
  \item{ylim}{Numeric vector of length 2 indicating the limits of the y-axis.}
  \item{add.grid}{Logical value indicating whether the grid must be added. Default is \code{TRUE}.}
  \item{las}{Orientation of the tick mark labels. Possible values are: parallel to the axis (0), horizontal (1), perpendicular to the axis (2), vertical (3). If \code{NULL} (the default), it is set automatically based on the length of axis labels.}
  \item{cex.axis}{Numeric vector of length 2 indicating the expansion factor for x-axis and y-axis labels, respectively. Default is (1,1).}
  \item{cex.lab}{Numeric vector of length 2 indicating the expansion factor for x and y labels, respectively. Default is (1,1).}
  \item{xlab}{Label for the x-axis.}
  \item{ylab}{Label for the y-axis.}
  \item{titles}{Character vector containing the title for each graphic. If \code{NULL} (the default), the names of endogenous variables are used as titles.}
  \item{obs.col}{Color for the time series of observed data. Default is black.}
  \item{fit.col}{Color for the time series of fitted values. Default is black.}
  \item{out.col}{Color for the time series of predictions. Default is red.}
  \item{new.col}{Color for the time series of new data. Default is black.}
  \item{interval.col}{Color of prediction intervals. Default is \code{grey40}.}
  \item{obs.lty}{Line type for the time series of observed data. Default is 1.}  
  \item{fit.lty}{Line type for the time series of fitted values. Default is 1.}
  \item{out.lty}{Line type for the time series of predictions. Default is 1.}
  \item{new.lty}{Line type for the time series of new data. Default is 1.}
  \item{obs.cex}{Expansion factor for points in the time series of observed data. If 0, points are suppressed. Default is 0.4.}
  \item{new.cex}{Expansion factor for points in the time series of new data. If 0, points are suppressed. Default is 0.4.}
  \item{mfrow}{Numerical vector of the form \code{c(nr,nc)}: graphics are displayed in an \code{nr}-by-\code{nc} array on the device by rows. If \code{NULL} (the default), it is determined automatically.}
  \item{mar}{Numeric vector of length 4 indicating the margin sizes in the following order: bottom, left, top, right. Default is c(3.5,3.5,2,2).}
  \item{mgp}{Numeric vector of length 3 indicating the location of axis labels, tick mark labels and tick marks, respectively, relative to the plot window. Default is c(2.3,0.8,0).}
  \item{...}{Further parameters to be passed to the generic method \code{plot}.}
}
%\value{}
%\details{}
%\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
%\references{}
\seealso{\code{\link{feVAR}}; \code{\link{predict.feVAR}}.}
\examples{
\dontrun{
data(agrisus2020)
x_agr <- colnames(agrisus2020)[4:15]
m_agr <- feVAR(var.names=x_agr, unit="Country", time="Year", exogenous="GDP",
  data=agrisus2020, box.cox=0, ndiff=1, nlags=1)

# predictions for the first four variables on unit 'Italy'
pr1 <- predict(m_agr, unit.id="Italy", var.names=varNames[1:4])  ## in-sample predictions
plot(pr1)
pr2 <- predict(m_agr, unit.id=15, var.names=1:4)                 ## equivalent
plot(pr2)
pr3 <- predict(m_agr, unit.id=15, var.names=1:4, n.ahead=3)      ## 3 steps ahead forecasts
plot(pr3)
}
}