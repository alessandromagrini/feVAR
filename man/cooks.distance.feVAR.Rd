\name{cooks.distance.feVAR}
\alias{cooks.distance.feVAR}
\title{Cook's distance method for class feVAR}
\description{Compute Cook's distances for an estimated ficed effects vector autoregressive model.}
\usage{
\method{cooks.distance}{feVAR}(model, ...)
}
\arguments{
  \item{model}{Object of class \code{feVAR}.}
  \item{...}{Further parameters to be passed to the generic method \code{cooks.distance}.}
}
\value{An object of class \code{data.frame} including the Cook's distance of each observation (by row) in each equation (by column), plus the maximal one across all equations (last column).}
\details{The Cook's distance (Cook, 1979) measures the leverage of each observation, i.e. how much its deletion impacts on the estimated regression hyperplane.}
%\note{}
%\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
\references{
Cook R.D. (1979). Influential observations in linear regression. \emph{Journal of the American Statistical Association}, 74(365): 169-174. DOI: 10.2307/2286747
}
\seealso{\code{\link{feVAR}}.}
\examples{
\dontrun{
data(agrisus2020)
x_agr <- colnames(agrisus2020)[4:15]
m_agr <- feVAR(var.names=x_agr, unit="Country", time="Year", exogenous="GDP",
  data=agrisus2020, box.cox=0, ndiff=1, nlags=1)

m_agr_cook <- cooks.distance(m_agr)  ## cook's distances in each equation
plot(m_agr_cook$'(maximum)')         ## display maximal distance per observation
}
}
