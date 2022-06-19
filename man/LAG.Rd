\name{LAG}
\alias{LAG}
\title{Matrix of lagged values}
\description{Generate the matrix of lagged values for a quantitative variable.}
\usage{LAG(x, p, unit=NULL, cut=0)}
\arguments{
  \item{x}{Numerical vector representing temporally ordered data.}
  \item{p}{Numerical value indicating the number of lags. If \code{p} is not integer, it is rounded to the nearest integer. If \code{p} is lower than 1, \code{x} is returned with no changes.}
  \item{unit}{Character containing the name of the units of observation. If \code{NULL} (the default), a single unit of observation is assumed.}
  \item{cut}{Numerical value indicating the number of initial time points to be set to missing. Default is 0.}
}
\value{A matrix with \code{length(x)} rows and \code{p} columns provided that \code{p} is greater or equal to 1. Otherwise, the vector \code{x} with no changes.}
%\note{}
%\details{}
%\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
%\references{}
%\seealso{\link{feVAR}.}
\examples{
# single unit of observation
x <- rnorm(10)  ## simulate a white noise on 10 time points
LAG(x, 3)       ## 3 lags

# panel data
x <- rnorm(30)
u <- rep(1:3, each=10)  ## 3 units, each with 10 time points
LAG(x, 3, unit=u)       ## 3 lags
}
