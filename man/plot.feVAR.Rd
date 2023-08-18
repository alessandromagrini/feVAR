\name{plot.feVAR}
\alias{plot.feVAR}
\title{Plot method for class feVAR}
\description{Display residual diagnostics for an estimated fixed effect vector autoregressive model.}
\usage{\method{plot}{feVAR}(x, type=c("ts","acf","acf2","fitVSres","qqnorm"), var.names=NULL, max.nlags=NULL,
  signif=0.05, ylim=NULL, cex.points=0.6, add.grid=TRUE, xlab=NULL, ylab=NULL,
  titles=NULL, las=0, mfrow=NULL, mar=c(3.5,3.5,2,2), mgp=c(2.3,0.8,0), ...)}
\arguments{
  \item{x}{Object of class \code{feVAR}.}
  \item{type}{Type of diagnostic, which should be one among 'ts' (1): time series plot of residuals, 'acf' (2): autocorrelograms of residuals, 'acf2' (3): autocorrelograms of squared residuals, 'fitVSres' (4): fitted values versus residuals, and 'qqnorm' (5): normal quantile plot of residuals.}
  \item{var.names}{Character vector indicating the name of endogenous variables for which the diagnostic must be displayed.}
  \item{max.nlags}{The number of lags to be displayed in the autocorrelograms. If \code{NULL} (the default), it is determined automatically. If an invalid value is provided, it is automatically replaced by the nearest valid one. Ignored if \code{type} is not equal to 'acf' or 'acf2'.}
  \item{signif}{The significance level of the test for no autocorrelation in the autocorrelograms. Default is 0.05. If outside the interval \code{[0,1]}, it is set to 0.05. Ignored if \code{type} is not equal to 'acf' or 'acf2'.}
  \item{ylim}{Numerical vector of length 2 indicating the limits of the y-axis.}
  \item{cex.points}{Expansion factor for points. If equal to 0, points are suppressed. Ignored if \code{type} is not equal to 'qqnorm' or 'fitVSres'. Default is 0.6.}
  \item{add.grid}{Logical value indicating whether the grid must be added. Ignored if \code{type} is equal to 'acf' or 'acf2'. Default is \code{TRUE}.}
  \item{xlab}{Label for the x-axis.}
  \item{ylab}{Label for the y-axis.}
  \item{titles}{Character vector containing the title for the graphics. If the vector has names, titles are assigned based on them, otherwise they are assigned in the same order of the endogenous variables. If \code{NULL} (the default), the names of endogenous variables are used as titles.}
  \item{las}{Orientation of the tick mark labels. Possible values are: parallel to the axis (0), horizontal (1), perpendicular to the axis (2), vertical (3). If \code{NULL} (the default), it is set automatically based on the length of axis labels.}
  \item{mfrow}{Numerical vector of the form \code{c(nr,nc)}: graphics are displayed in an \code{nr}-by-\code{nc} array on the device by rows. If \code{NULL} (the default), it is determined automatically.}
  \item{mar}{Numerical vector of length 4 indicating the margin sizes in the following order: bottom, left, top, right. Default is c(3.5,3.5,2,2).}
  \item{mgp}{Numerical vector of length 3 indicating the location of axis labels, tick mark labels and tick marks, respectively, relative to the plot window. Default is c(2.3,0.8,0).}
  \item{...}{Further parameters to be passed to the generic method \code{plot}.}
}
%\value{}
\details{
In case of panel data, the autocorrelation of residuals at lag \eqn{k} is computed as:
\deqn{\rho_k=\frac{\sum_{i=1}^n\sum_{t=1}^{T_i-k}e_{i,t}~e_{i,t+k}}{\sum_{i=1}^n\sum_{t=1}^{T_i}e_{i,t}^2}}
where \eqn{e_{i,t}} is the residual for unit \eqn{i} at time \eqn{t}, \eqn{~T_i} is the number of time points for unit \eqn{i}, and \eqn{n} is the number of units.
The autocorrelation of squared residuals is computed analogously by replacing \eqn{e_{i,t}} with \eqn{e_{i,t}^2} for each \eqn{i} and \eqn{t}.
}
\note{
In order to meet the hypothesis of serially uncorrelated errors, autocorrelograms of residuals must show no significant spikes. Note that the violation of the hypothesis of serially uncorrelated errors leads to inconsistent estimates.

In order to meet the hypothesis of homoschedastic errors, autocorrelograms of squared residuals must show no significant spikes, and the graphics of fitted values versus residuals must show no patterns. Note that the violation of the hypothesis of homoschedasticity leads to inefficient but still consistent estimates. The Box-Cox transformation may solve heteroschedasticity (see argument \code{box.cox} of the function \code{\link{feVAR}}).

In order to meet the hypothesis of normal errors, the normal quantile plot of residuals must look approximatively linear. Note that the violation of the hypothesis of normal errors does not affect neither consistency nor efficiency of the estimates, but significance tests and confidence intervals are no more exact and are valid only asymptotically.
}
%\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
%\references{}
\seealso{\code{\link{feVAR}}.}
\examples{
data(agrisus2020)
x_agr <- colnames(agrisus2020)[4:15]
m_agr <- feVAR(var.names=x_agr, unit="Country", time="Year", exogenous="GDP",
  data=agrisus2020, box.cox=0, ndiff=1, nlags=1)

# time series plots of residuals
plot(m_agr, type="ts", cex.main=1.1)

# autocorrelograms of residuals
plot(m_agr, type="acf", cex.main=1.1)

# autocorrelograms of squared residuals
plot(m_agr, type="acf2", cex.main=1.1)

# fitted versus residuals
plot(m_agr, type="fitVSres", cex.main=1.1)

# normal quantile plots of residuals
plot(m_agr, type="qqnorm", cex.main=1.1)
}