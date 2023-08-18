\name{LAG}
\alias{LAG}
\title{Generate lagged values}
\description{Generate lagged values of a quantitative variable.}
\usage{LAG(x, lag, unit=NULL, del=0)}
\arguments{
  \item{x}{Numerical vector representing temporally ordered data.}
  \item{lag}{The lag to generate. If 0, \code{x} is returned with no changes. If an invalid value is provided, it is automatically replaced by the nearest valid one.}
  \item{unit}{Numerical or character vector of the same length of \code{x} containing the identifiers of the units of observation. If \code{NULL} (the default), a single unit of observation is assumed.}
  \item{del}{The number of initial observations to be set to missing. Default is 0. Note that the first \code{lag} observations are inevitably lost and therefore set to missing.}
}
\value{Numerical vector of the same length of \code{x}.}
%\note{}
%\details{}
%\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
%\references{}
%\seealso{\code{\link{feVAR}}.}
\examples{
# single unit of observation
x <- rnorm(10)  ## simulate a white noise on 10 time points
LAG(x, 3)       ## lag 3

# panel data
x <- rnorm(30)
u <- rep(1:3, each=10)  ## 3 units, each with 10 time points
LAG(x, 3, unit=u)       ## lag 3
}
