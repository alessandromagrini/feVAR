\name{VMA}
\alias{VMA}
\title{Vector moving average for a feVAR model}
\description{Compute the coefficients of the vector moving average representation for an estimated fixed effects vector autoregressive model.}
\usage{VMA(model, n.ahead=10)}
\arguments{
  \item{model}{Object of class \code{feVAR}.}
  \item{n.ahead}{The number of steps ahead, which must be no lower than 1. Default is 10. If an invalid value is provided, it is automatically replaced by the nearest valid one.}
  }
\value{Array with three dimensions that correspond, respectively, to the number of steps ahead, to the endogenous variable receiving the impulse, and to the endogenous variables originating the impulse.}
%\note{}
%\details{}
%\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
%\references{}
\seealso{\code{\link{feVAR}}, \code{\link{IRF}}.}
\examples{
\dontrun{
data(agrisus2020)
x_agr <- colnames(agrisus2020)[4:15]
m_agr <- feVAR(var.names=x_agr, unit="Country", time="Year", exogenous="GDP",
  data=agrisus2020, box.cox=0, ndiff=1, nlags=1)

m_agr_vma <- VMA(m_agr, n.ahead=10)
m_agr_vma["1",,]  ## VMA coefficient at lag 1
m_agr_vma["10",,]  ## VMA coefficient at lag 10
}
}
