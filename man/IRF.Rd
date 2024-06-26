\name{IRF}
\alias{IRF}
\title{Impulse response functions for a feVAR model}
\description{Compute impulse response functions (IRFs) for an estimated fixed effects vector autoregressive model.}
\usage{IRF(model, n.ahead=10, cumulative=FALSE, orthogonal=FALSE, order=NULL, nboot=100,
  level=0.95)}
\arguments{
  \item{model}{Object of class \code{feVAR}.}
  \item{n.ahead}{The number of steps ahead, which must be no lower than 1. Default is 10. If an invalid value is provided, it is automatically replaced by the nearest valid one.}
  \item{cumulative}{Logical value indicating whether cumulative IRFs must be computed. Default is \code{FALSE}.}
  \item{orthogonal}{Logical value indicating whether orthogonal IRFs must be computed. Default is \code{FALSE}.}
  \item{order}{Vector indicating the order of orthogonalization. It is possible to provide either the names or the numerical identifiers of endogenous variables. If \code{NULL} (the default), the original order of endogenous variables is used. Ignored if \code{orthogonal} is \code{FALSE}.}
  \item{nboot}{The number of bootstrap resamples to compute confidence intervals. Default is 100. If not numeric or 0, confidence intervals are not computed.}
  \item{level}{The confidence level. Default is 0.95. If outside the interval \code{[0,1]}, it is set to 0.95.}
  }
\value{Object of class \code{IRF.feVAR}, which is a list with the following three components:
\itemize{
  \item{\code{irf}, }{including the computed IRFs,} 
  \item{\code{irf_sx}, }{including the lower bound of confidence intervals for IRFs,}
  \item{\code{irf_dx}, }{including the upper bound of confidence intervals for IRFs.}
}
Each of these three components is an array with three dimensions that correspond, respectively, to the number of steps ahead, to the endogenous variable receiving the impulse, and to the endogenous variables originating the impulse.}
\note{The class \code{IRF.feVAR} has its own \code{plot} method. See \code{\link{plot.IRF.feVAR}} for details.}
\details{
Impulse response functions (IRFs) at \eqn{k} steps ahead, say \eqn{\Phi_k=\frac{\partial\bm{y}_{i,t+k}}{\partial\bm{\varepsilon}_{i,t}}}, are computed based on the following recursive formula:
\deqn{\Phi_k=\sum_{j=1}^{\min(k,p)}\Phi_{k-j}\bm{B}_j}
where \eqn{\Phi_0} is a lower triangular \eqn{m \times m} matrix with 1s in the main diagonal,
\eqn{\bm{B}_j} is the \eqn{m \times m} matrix of coefficients at lag \eqn{j}, and \eqn{p} is the lag order.
Matrix \eqn{\Phi_0} is set equal to the \eqn{m \times m} identity matrix in order to obtain simple (non-orthogonal) IRFs, otherwise it is obtained through Cholesky decomposition.
}
%\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
%\references{}
\seealso{\code{\link{feVAR}}, \code{\link{plot.IRF.feVAR}}.}
\examples{
\dontrun{
data(agrisus2020)
x_agr <- colnames(agrisus2020)[4:15]
m_agr <- feVAR(var.names=x_agr, unit="Country", time="Year", exogenous="GDP",
  data=agrisus2020, box.cox=0, ndiff=1, nlags=1)

m_agr_irf <- IRF(m_agr, nboot=500)        ## compute IRFs (500 bootstrap resamples
                                          ##   for more accurate confidence intervals)
m_agr_irf$irf[,"TFP_2015",]               ## IRFs generated by 'TFP_2015'
m_agr_irf$irf[,"NetCapital_GVA",]         ## IRFs received by 'NetCapital_GVA'
m_agr_irf$irf[,"TFP_2015","Poverty_rur"]  ## IRFs from 'TFP_2015' to 'Poverty_rur'
}
}
