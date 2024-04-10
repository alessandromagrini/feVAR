\name{feVAR-package}
\alias{feVAR-package}
\docType{package}
\title{Fixed Effects Vector Autoregressive Models}
\description{Estimation and analysis of fixed effects vector autoregressive models. The package implements the half panel jackknife bias correction (Dhaene & Jochmans, 2015) to obtain unbiased estimates, and the Expectation-Maximization (EM) algorithm (Dempster et al., 1977) to perform estimation in case of missing values. Several numerical and graphical tools to analyse the estimated models are available. Both panel, even unbalanced, and non-panel data are allowed.}
\details{
  \tabular{ll}{
  Package: \tab feVAR\cr
  Type: \tab Package\cr
  Version: \tab 0.1.1\cr
  Date: \tab 2024-04-11\cr
  License: \tab GPL-2
  }
Fixed effects vector autoregressive models (feVARs) are an extension of vector autoregressive models to panel data, firstly introduced by Holtz-Eakin et al. (1988).
Parameters of a feVAR model can be estimated using, for each equation, the Least Squares Dummy Variable (LSDV) estimator, i.e., ordinary least squares with dummy variables to denote cross-sectional units.
Parameters of a feVAR model can be estimated using, for each equation, the Least Squares Dummy Variable (LSDV) estimator, i.e., ordinary least squares with dummy variables to denote cross-sectional units.
Unfortunately, due to the presence of autoregressive terms, LSDV estimates have a finite sample bias of magnitude O(1/T) (Hahn & Kuersteiner, 2002), where T is the number of time points, which may be relevant even when T is moderately large (Judson & Owen, 1999).
The Half Panel Jackknife (HPJ) bias correction (Dhaene & Jochmans, 2015) can be used to reduce the bias of the LSDV estimator to O(1/T^2).
Compared to Generalized Method of Moments (GMM) estimators that require a large number of cross-sectional units to be efficient, the HPJ bias correction performs well in general, because it reduces the bias without inflating the asymptotic variance.

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
Dempster A.P., Laird N.M., Rubin D.B. (1977). Maximum likelihood from incomplete data via the EM algorithm. \emph{Journal of the Royal Statistical Society}, Series B, 39(1): 1-38. 

Dhaene G., Jochmans K. (2015). Split-panel jackknife estimation of fixed-effect models. \emph{The Review of Economic Studies}, 82(3): 991-1030. DOI: 10.1093/restud/rdv007

Hahn J., Kuersteiner G. (2002). Asymptotically unbiased inference for a dynamic panel model with fixed effects when both n and T are large. \emph{Econometrica}, 70(4): 1639-1657. DOI: 10.1111/1468-0262.00344

Holtz-Eakin D., Newey W., Rosen H.S. (1988). Estimating vector autoregressions with panel data. \emph{Econometrica}, 56(6), 1371-1395. DOI: 10.2307/1913103

Judson R.A., Owen A.L. (1999). Estimating dynamic panel data models: a guide for macroeconomists. \emph{Economics Letters}, 65: 9-15. DOI: 10.1016/S0165-1765(99)00130-5

Nickell S. (1981). Biases in dynamic models with fixed effects. \emph{Econometrica}, 49: 1417-1426.

Sigmund M., Ferstl R. (2021). Panel vector autoregression in R with the package panelvar. \emph{The Quarterly Review of Economics and Finance}, 80: 693-720. DOI: 10.1016/j.qref.2019.01.001.
}
