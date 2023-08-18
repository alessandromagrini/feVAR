\name{autocorTest}
\alias{autocorTest}
\title{Ljung-Box test for a feVAR model}
\description{Ljung-Box autocorrelation test on residuals for an estimated fixed effects vector autoregressive model.}
\usage{autocorTest(model, max.lag=NULL)}
\arguments{
  \item{model}{Object of class \code{feVAR}.}
  \item{max.lag}{Maximum lag at which to perform the test, which must be no lower than 1 and no greater than \code{T-5}, where \code{T} is the minimum number of time points across all units of observations. If \code{NULL} (the default), it is set to \code{sqrt(T)}. If an invalid value is provided, it is automatically replaced by the nearest valid one.}
}
\value{An object of class \code{autocorTest.feVAR}, which is a list with two components: \code{residual}, containing the results of the test on residuals, and \code{sq.residual}, containing the results of the test on squared residuals.
Each of these two components is a list with the following two components:
\itemize{
  \item{\code{statistic}: }{matrix containing the test statistic for each lag (by row) and each endogenous variable (by column);}
  \item{\code{p.value} :}{matrix containing the p-value for each lag (by row) and each endogenous variable (by column).}
  }
}
\note{The class \code{autocorTest.feVAR} has its own \code{print} method that reports the p-value of the test at lag \code{max.lag} for each endogenous variable.}
\details{
The Ljung-Box test (Ljung & Box, 1978) checks the hypothesis of no autocorrelation up to a specific lag.
The test is based on the fact that, under the hypothesis of no autocorrelation up to lag \eqn{k}, it holds:
\deqn{Q_k=T(T+2)\sum_{j=1}^k\frac{\hat{\rho}_j^2}{T-k}\approx\chi^2_k}
where \eqn{\hat{\rho}_j} is the sample autocorrelation at lag \eqn{j}, \eqn{~T} is the total number of time points across all units of observation, and \eqn{\chi^2_k} is a Chi-squared random variable with \eqn{k} degrees of freedom.
See 'Details' of the function \code{\link{plot.feVAR}} for details on how to compute the autocorrelation of residuals and of squared residuals in case of panel data.

Performing the test on the residuals provides an indication about the hypothesis of serially uncorrelated errors, while performing the test on the squared residuals provides an indication about the hypothesis of homoschedastic errors.
}
%\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
\references{
%Box G.E.P., Pierce D.A. (1970). Distribution of residual autocorrelations in autoregressive integrated moving average time series models. \emph{Journal of the American Statistical Association}, 65(332): 1509-1526. DOI: 10.1080/01621459.1970.10481180
%
Ljung G.M., Box G.E.P. (1978). On a measure of a lack of fit in time series models. \emph{Biometrika}, 65(2): 297-303. DOI: 10.1093/biomet/65.2.297
}
\seealso{\code{\link{feVAR}}.}
\examples{
data(agrisus2020)
x_agr <- colnames(agrisus2020)[4:15]
m_agr <- feVAR(var.names=x_agr, unit="Country", time="Year", exogenous="GDP",
  data=agrisus2020, box.cox=0, ndiff=1, nlags=1)

m_agr_bpt <- autocorTest(m_agr)
m_agr_bpt                      ## p-value at lag 5
m_agr_bpt$residual$p.value     ## p-values at each lag (test on residuals)
m_agr_bpt$sq.residual$p.value  ## p-values at each lag (test on squared residuals)
}
