\name{residBPtest}
\alias{residBPtest}
\title{Box-Pierce test for a feVAR model}
\description{Box-Pierce test on residuals for an estimated fixed effects vector autoregressive model.}
\usage{residBPtest(model, max.lag=NULL)}
\arguments{
  \item{model}{Object of class \code{feVAR}.}
  \item{max.lag}{Maximum lag at which to perform the test, which must be no lower than 1 and no greater than \code{T-5}, where \code{T} is the minimum time length across all units of observations. If \code{NULL} (the default), it is set to \code{sqrt(T)}. If an invalid value is provided, it is automatically replaced by the nearest valid one.}
}
\value{An object of class \code{residBPtest.feVAR}, which is a list with two components: \code{residual}, containing the results of the test on residuals, and \code{sq.residual}, containing the results of the test on squared residuals.
Each of these two components is a list with the following two components:
\itemize{
  \item{\code{statistic}: }{matrix containing the test statistic for each lag (by row) and each endogenous variable (by column);}
  \item{\code{p.value} :}{matrix containing the p-value for each lag (by row) and each endogenous variable (by column).}
  }
}
\note{The class \code{residBPtest.feVAR} has its own \code{print} method that reports, for each endogenous variable, the minimum lag with significant autocorrelation in the residuals and in the squared residuals.}
\details{The Box-Pierce test (Box & Pierce, 1970) checks the hypothesis of no autocorrelation at specific lags. Performing the test on the residuals provides an indication about the hypothesis of serially uncorrelated errors, while performing the test on the squared residuals provides an indication about the hypothesis of homoschedastic errors.}
%\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
\references{
Box G.E.P., Pierce D.A. (1970). Distribution of residual autocorrelations in autoregressive integrated moving average time series models. \emph{Journal of the American Statistical Association}, 65(332): 1509-1526. DOI: 10.1080/01621459.1970.10481180
}
\seealso{\code{\link{feVAR}}.}
\examples{
\dontrun{
data(agrisus2020)
x_agr <- colnames(agrisus2020)[4:15]
m_agr <- feVAR(var.names=x_agr, unit="Country", time="Year", exogenous="GDP",
  data=agrisus2020, box.cox=0, ndiff=1, nlags=1)

# perform the box-pierce test on residuals
m_agr_bpt <- residBPtest(m_agr)
m_agr_bpt  ## summary of the test

# p-value of the test at each lag
m_agr_bpt$residual$p.value     ## on residuals
m_agr_bpt$sq.residual$p.value  ## on squared residuals
}
}
