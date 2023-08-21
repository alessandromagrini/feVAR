\name{feVAR-package}
\alias{feVAR-package}
\docType{package}
\title{Fixed Effects Vector Autoregressive Models}
\description{Estimation and analysis of fixed effects vector autoregressive models. The package implements the bootstrap bias correction for ordinary least squares estimation proposed by Everaert & Pozzi (2007) and the Expectation-Maximization (EM) algorithm (Dempster et al., 1977) to perform estimation in case of missing values. Several numerical and graphical tools to analyse the estimated models are available. Both panel, even unbalanced, and non-panel data are allowed.}
\details{
  \tabular{ll}{
  Package: \tab feVAR\cr
  Type: \tab Package\cr
  Version: \tab 0.1.1\cr
  Date: \tab 2023-08-18\cr
  License: \tab GPL-2
  }
Fixed effects vector autoregressive models (feVARs) are an extension of vector autoregressive models to panel data, firstly introduced by Holtz-Eakin et al. (1988).
The main problem of feVARs is that ordinary least squares (OLS) estimation is biased even with a large number of units (Nickell, 1981), and, although the bias diminishes as the number of time points gets larger, simulations show significant bias when the time series has moderate size (e.g. 30 time points, see Judson & Owen, 1999).
For this reason, generalized methods of moments (GMM) estimators like difference GMM (Arellano & Bond, 1991) and system GMM (Blundell & Bond, 1998) are often employed in place of OLS estimation.
Unfortunately, unbiasedness of GMM estimators strictly depends on the characteristics of the multivariate process underlying the data (weak instrument problem), and may be computationally unfeasible in short panels due to singularity.
In this case, a simple and effective solution is represented by the bootstrap bias correction for OLS estimation proposed by Everaert & Pozzi (2007).

The main functions of the package are:
\itemize{
  \item{\code{\link{unirootTest}}: }{to check stationarity of the time series;}
  \item{\code{\link{feVAR}}: }{to estimate a feVAR model;}
  \item{\code{\link{plot.feVAR}}: }{to display residual diagnostics for an estimated feVAR model;}
  \item{\code{\link{IRF}} and \code{\link{plot.IRF.feVAR}}: }{to compute and display impulse response functions for an estimated feVAR model;}
  \item{\code{\link{predict.feVAR}} and \code{\link{plot.predict.feVAR}}, }{to perform and display in-sample predictions or h-step ahead forecasts based on an estimated feVAR model;}
  \item{\code{\link{forecastError}}, }{to estimate the h-step ahead forecast error for an estimated feVAR model.}
}
}
\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
\references{
Arellano M., Bond S. (1991). Some tests of specification for panel data: Monte Carlo evidence and an application to employment equations. \emph{Review of Economic Studies}, 58(2): 277. DOI: 10.2307/2297968

Blundell R., Bond S. (1998). Initial conditions and moment restrictions in dynamic panel data models. \emph{Journal of Econometrics}, 87(1): 115-143. DOI: 10.1016/S0304-4076(98)00009-8

Dempster A.P., Laird N.M., Rubin D.B. (1977). Maximum likelihood from incomplete data via the EM algorithm. \emph{Journal of the Royal Statistical Society}, Series B, 39(1): 1-38. 

Everaert G., Pozzi L. (2007). Bootstrap-based bias correction for dynamic panels. \emph{Journal of Economic Dynamics & Control}, 31: 1160-1184. DOI: 10.1016/j.jedc.2006.04.006

Holtz-Eakin D., Newey W., Rosen H.S. (1988). Estimating vector autoregressions with panel data. \emph{Econometrica}, 56(6), 1371-1395. DOI: 10.2307/1913103

Judson R.A., Owen A.L. (1999). Estimating dynamic panel data models: a guide for macroeconomists. \emph{Economics Letters}, 65: 9-15. DOI: 10.1016/S0165-1765(99)00130-5

Nickell S. (1981). Biases in dynamic models with fixed effects. \emph{Econometrica}, 49: 1417-1426.
}
