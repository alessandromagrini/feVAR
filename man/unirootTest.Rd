\name{unirootTest}
\alias{unirootTest}
\title{Unit root tests for panel data}
\description{Perform augmented Dickey-Fuller (ADF) and Kwiatkowski-Phillips-Schmidt-Shin (KPSS) tests for unit roots on a multivariate time series, which may have a panel structure.}
\usage{unirootTest(var.names, unit=NULL, time=NULL, data, box.cox=1, ndiff=0, max.nlags=NULL)}
\arguments{
  \item{var.names}{Character vector including the name of the variables to undergo the tests.}
  \item{unit}{Character containing the name of the variable that identifies the units of observation. If \code{NULL} (the default), a single unit of observation is assumed.}
  \item{time}{Character containing the name of the time variable, which must be in numeric or date format. If \code{NULL} (the default), data are assumed to be temporally ordered.}
  \item{data}{Object of class \code{data.frame} containing the variables in \code{var.names}, \code{unit} and \code{time}.}
  \item{box.cox}{Named vector including the \eqn{\lambda} values of the Box-Cox transformation (Box & Cox, 1964) for variables in \code{var.names}. If \code{box.cox} has no names and length greater than one, the same ordering as in \code{var.names} is assumed. If \code{box.cox} has no names and length equal to one, the same \eqn{\lambda} value is used for all variables in \code{var.names}. Valid \eqn{\lambda} values are in the interval \code{[-2,2]}. All invalid values are automatically replaced with the nearest valid ones. Default is 1, meaning no transformation for all variables. See 'Details'.}
  \item{ndiff}{Named vector including the number of differences for variables in \code{var.names}. If \code{ndiff} has no names and length greater than one, the same ordering as in \code{var.names} is assumed. If \code{ndiff} has no names and length equal to one, the same number of differences is used for all variables in \code{var.names}. The number of differences cannot exceed \code{T-5}, where \code{T} is the minimum number of time points across all units of observations. All invalid values are automatically replaced with the nearest valid ones. Default is 0, meaning no differencing for all variables.}
  \item{max.nlags}{The maximum lag to consider in the tests, which must be no greater than \code{T-5}, where \code{T} is the minimum number of time points across all units of observations. If \code{NULL} (the default), it is set to \code{sqrt(T)}. If an invalid value is provided, it is automatically replaced by the nearest valid one. See 'Details'.}
}
\value{An object of class \code{unirootTest}, which is a list including a list for each variable in \code{var.names}, each with the following components:
\itemize{
\item{\code{statistic}: }{test statistic for each test};
\item{\code{lag.selected}: }{lag selected for each test};
\item{\code{p.value}: }{p-value for each test};
\item{\code{box.cox}: }{offset and lambda values of the Box-Cox transformation for each variable subjected to the tests};
\item{\code{ndiff}: }{order of differencing for each variable subjected to the tests}.
}
}
\note{The class \code{unirootTest} has its own \code{print} method that reports the p-values of the tests for each endogenous variable.}
\details{Only quantitative variables can undergo the tests.
In order to keep statistical power as high as possible in presence of missing values, they are replaced by the most recent observed value (Ryan & Giles, 1998).

The Box-Cox transformation for variable is \eqn{Y} defined as:
\deqn{f(y;\lambda,\eta)=I(\lambda\neq 0)\cdot\frac{(y+\eta)^\lambda-1}{\lambda}+I(\lambda=0)\cdot\log(y+\eta)}
where \eqn{y} is a specific value of \eqn{Y}, \eqn{~I(\cdot)} is the indicator function, and \eqn{\eta} is an offset equal to zero if \eqn{Y} contains no negative values, otherwise equal to \eqn{-min(Y)}. If \eqn{\lambda=0} (logarithmic transformation), then \eqn{\eta} is further increased by half the lowest non-zero value of \eqn{Y+min(Y)}.
The value of \eqn{\eta} is computed automatically and cannot be specified by the user.
If \eqn{\lambda=1}, no transformation is applied.

In case of panel data, the ADF test is performed by combining the individual p-values according to the method by Demetrescu \emph{et al.} (2006), while the KPSS test is performed by combining the individual test statistics as proposed by Hadri (2000).

The tests account for a deterministic linear trend, unless the order of differencing is greater than 0. The lag length at which to perform the tests is selected through AIC-based backward elimination starting from the lag length passed to argument \code{max.lag}.

The null hypothesis of the ADF test (Dickey & Fuller, 1981) is the presence of unit roots, while the null hypothesis of the KPSS test (Kwiatkowski \emph{et al.}, 1992) is the absence of unit roots.
Therefore, p-value higher than 0.05 for the ADF test or p-value lower than 0.05 for the KPSS test suggest the presence of unit roots and the need of further differencing.}
%\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
\references{
Demetrescu M., Hassler U., Tarcolea A. (2006). Combining significance of correlated statistics with application to panel data. \emph{Oxford Bulletin of Economics and Statistics}, 68(5), 647-663. DOI: 10.1111/j.1468-0084.2006.00181.x
  
Dickey D.A., Fuller W.A. (1981). Likelihood ratio statistics for autoregressive time series with a unit root. \emph{Econometrica}, 49(4): 1057-1072. DOI: 10.2307/1912517

Hadri K. (2000). Testing for stationarity in heterogeneous panel data. \emph{Econometrics Journal}, 3: 148-161. DOI: 10.1111/1368-423X.00043 

Kwiatkowski D., Phillips P.C.B., Schmidt P., Shin Y. (1992). Testing the null hypothesis of stationarity against the alternative of a unit root. \emph{Journal of Econometrics}, 54(1-3): 159-178. DOI: 10.1016/0304-4076(92)90104-Y

Ryan K.F., Giles D.E.A (1998). Testing for unit roots in economic time-series with missing observations. In: T. B. Fomby & R. C. Hill (eds.) \emph{Messy Data (Advances in Econometrics, Vol. 13)}, pages 203-242. DOI: 10.1108/S0731-9053(1999)0000013010
}
\seealso{\code{\link{feVAR}}.}
\examples{
data(agrisus2020)
x_agr <- names(agrisus2020)[4:15]

# tests on variables in level
#  -> most variables are nonstationary
unirootTest(x_agr, unit="Country", time="Year", data=agrisus2020)

# tests on variables in logarithmic differences
#  -> stationarity confirmed for all variables excepting 'Manager_ratio'
unirootTest(x_agr, unit="Country", time="Year", data=agrisus2020,
  box.cox=0, ndiff=1)
}