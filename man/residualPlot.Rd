\name{residualPlot}
\alias{residualPlot}
\title{Graphics of residual diagnostics}
\description{Display graphics of residual diagnostics for a fixed effect vector autoregressive model.}
\usage{residualPlot(model, type="ts", var.names=NULL, max.nlags=NULL, signif=0.05,
  acf.return=FALSE, ylim=NULL, cex.points=0.6, add.grid=TRUE,  xlab=NULL, ylab=NULL,
  titles=NULL, las=0, mar=c(3.5,3.5,2,2), mgp=c(2.3,0.8,0), ...)}
\arguments{
  \item{model}{An object of class \code{feVAR}.}
  \item{type}{Type of diagnostic, which should be one among 'ts' (1): time series plot of residuals, 'acf' (2): autocorrelograms of residuals, 'qq' (3): normal quantile plot of residuals, and 'fitVSres' (4): fitted values versus residuals.}
  \item{var.names}{Charcater vector indicating the name of endogenous variables for which prediction should be displayed.}
  \item{max.nlags}{Integer non-negative number indicating the number of lags to be displayed in the ACF. If \code{NULL} (the default), it is automatically determined. Ignored if \code{type} is not equal to 'acf'.}
  \item{signif}{Numeric value between 0 and 1 indicating the significance level of the autocorrelation test. Ignored if \code{type} is not equal to 'acf'. Default is 0.05.}
  \item{acf.return}{Logical value indicating whether autocorrelations should be returned. Ignored if \code{type} is not equal to 'acf'. Default is \code{FALSE}.}
  \item{ylim}{Numeric vector of length 2 indicating the limits of the y-axis.}
  \item{cex.points}{Numeric value indicating the expansion factor for points. If equal to 0, points are suppressed. Ignored if \code{type} is not equal to 'qq' or 'fitVSres'. Default is 0.6.}
  \item{add.grid}{Logical value indicating whether the grid should be added. Ignored if \code{type} is equal to 'acf'. Default is \code{TRUE}.}
  \item{xlab}{Label for the x-axis.}
  \item{ylab}{Label for the y-axis.}
  \item{titles}{Character vector containing the title for each graphic. If \code{NULL} (the default), the names of endogenous variables are used as titles.}
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
data(agrisus2020)
varNames <- colnames(agrisus2020)[4:15]
m_agr <- feVAR(varNames, unit="Country", time="Year", exogenous="GDP",
  data=agrisus2020, box.cox=0, ndiff=1, nlags=1)

\dontrun{ 
residualPlot(m_agr, type="ts", cex.main=1.1)        ## time series plot
residualPlot(m_agr, type="acf", cex.main=1.1)       ## autocorrelogram
residualPlot(m_agr, type="qq", cex.main=1.1)        ## normal quantile plot
residualPlot(m_agr, type="fitVSres", cex.main=1.1)  ## fitted versus residuals
}
}