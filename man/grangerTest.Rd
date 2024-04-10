\name{grangerTest}
\alias{grangerTest}
\title{Granger causality test for a feVAR model}
\description{Granger causality test for an estimated fixed effects vector autoregressive model.}
\usage{grangerTest(model, from=NULL, to=NULL)}
\arguments{
  \item{model}{Object of class \code{feVAR}.}
  \item{from}{Numerical vector including the name of the endogenous variables considered as predictors. If \code{NULL} (the default), all endogenous variables are considered as predictors.}
  \item{to}{Numerical vector including the name of the endogenous variables considered as predictands. If \code{NULL} (the default), all endogenous variables are considered as predictands.}
}
\value{An object of class \code{grangerTest.feVAR}, which is a list with the following components:
\itemize{
  \item{\code{from}: }{endogenous variables considered as predictors in the tests;}
  \item{\code{to}: }{endogenous variables considered as predictands in the tests;}
  \item{\code{F.test}: }{a list of objects of class \code{anova}, one for each predictand variable, including the result of the test;}
  \item{\code{p.value} :}{a vector including the p-value of the test for each predictand variable.}
  }
}
\note{When a variable is both predictand and predictors, it is excluded from the predictor set.
For instance, when both arguments \code{from} and \code{to} are \code{NULL}, each endogenous variable is tested for granger causality given all other endogenous variables.

The class \code{grangerTest.feVAR} has its own \code{print} method that reports the p-value of the test for each predictand variable.}
%\details{}
%\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
%\references{}
\seealso{\code{\link{feVAR}}.}
\examples{
data(agrisus2020)
x_agr <- colnames(agrisus2020)[4:15]
m_agr <- feVAR(var.names=x_agr, unit="Country", time="Year", exogenous="GDP",
  data=agrisus2020, box.cox=0, ndiff=1, nlags=1)

# each endogenous variable is tested for granger causality
#   given all other endogenous variables
grangerTest(m_agr)
}
