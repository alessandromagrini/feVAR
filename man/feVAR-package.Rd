\name{gammadlm-package}
\alias{gammadlm-package}
\docType{package}
\title{Fixed Effects Vector Autoregressive Models}
\description{Estimation and analysis of fixed effects vector autoregressive models.}
\details{
  \tabular{ll}{
  Package: \tab feVAR\cr
  Type: \tab Package\cr
  Version: \tab 0.1.1\cr
  Date: \tab 2022-06-19\cr
  License: \tab GPL-2
  }
Fixed effects Vector Autoregressive Models are Vector Autoregressive Models (VAR) for panel data, firstly introduced by Holtz-Eakin et al. (1988).
Let \eqn{X_1,\ldots,X_m} be \eqn{m} endogenous variables and \eqn{Z_1,\ldots,Z_q} be \eqn{q} exogenous variables.
Also, let \eqn{\bm{y}_{i,t}} and \eqn{\bm{z}_{i,t}} be, respectively, the vector of endogenous and exogenous variables observed on unit \eqn{i} at time \eqn{t}.
A fixed effect Vector Autoregressive Model (feVAR) of order \eqn{p\in\mathrm{N}_0} is defined as:
\deqn{\bm{y}_{i,t} = \bm{\alpha}_i+\sum_{j=1}^p\bm{B}_j\bm{y}_{i,t-j}+\bm{G}\bm{x}_{i,t}+\bm{\varepsilon}_{i,t}}
where:
\itemize{
\item{\eqn{\bm{\alpha}_i} is a \eqn{p}-dimensional vector including the intercepts of the endogenous variables for unit \eqn{i}};
\item{\eqn{\bm{B}_j}} is a \eqn{p\times p} matrix including the regression coefficients among the endogenous variables at time lag \eqn{j};
\item{\eqn{\bm{G}}} is a \eqn{q\times q} matrix including the regression coefficients of exogenous variables on endogenous variables;
\item{\eqn{\bm{\varepsilon}_{i,t}} is a \eqn{p}-dimensional vector of random errors for unit \eqn{i} at time \eqn{t} such that \eqn{\mathrm{E}(\bm{\varepsilon}_{i,t})=0} and \eqn{\mathrm{E}(\bm{\varepsilon}_{i,t}\cdot\bm{\varepsilon}_{k,t})=0} \eqn{\forall i,k,t}.}
}
The order \eqn{p} represents the number of lags up to which endogenous variables are supposed to influence each other, which can be determined based on information criteria.

The main functions of the package are:
\itemize{
  \item{\link{unirootTest}, }{to check stationarity of the time series;}
  \item{\link{feVAR}, }{to estimate a fixed effect autoregressive model through bias-corrected ordinary least squares (Dhaene & Jochmans, 2016);}
  \item{\link{residualPlot}, }{to display residual diagnostics.}
  \item{\link{stabilityCheck}, }{to check stationarity (stability) of an estimated model;}
  \item{\link{predict.feVAR} and \link{predictPlot}, }{to perform and display predictions.}
}
}
\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
\references{
G. Dhaene, and K. Jochmans (2016). Bias-corrected estimation of panel vector autoregressions. \emph{Economics Letters}, 145: 98-103. DOI: 10.1016/j.econlet.2016.06.010

D. Holtz-Eakin, W. Newey, and H. S. Rosen (1988). Estimating vector autoregressions with panel data. \emph{Econometrica}, 56(6), 1371-1395. DOI: 10.2307/1913103
}
