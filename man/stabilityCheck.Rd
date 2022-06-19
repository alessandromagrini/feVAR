\name{stabilityCheck}
\alias{stabilityCheck}
\title{Stability check for a fixed effect vector autoregressive model}
\description{Check whether an estimated fixed effect vector autoregressive model is stationary (stable).}
\usage{stabilityCheck(model)}
\arguments{
  \item{model}{Object of class \code{feVAR}.}
  }
\value{Numerical vector including the modulus of each eigenvector of the companion matrix.}
%\note{}
\details{A vector autoregressive model is stationary (stable) if all eigenvalues of the companion matrix have modulus less than 1.}
%\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
%\references{}
\seealso{\link{feVAR}.}
\examples{
data(agrisus2020)
varNames <- colnames(agrisus2020)[4:15]
m_agr <- feVAR(varNames, unit="Country", time="Year", exogenous="GDP",
  data=agrisus2020, box.cox=0, ndiff=1, nlags=1)

stabilityCheck(m_agr)  ## all moduli are lower than 1 -> stable
}
