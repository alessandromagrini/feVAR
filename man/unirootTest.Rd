\name{unirootTest}
\alias{unirootTest}
\title{Unit root and stationarity tests for panel data}
\description{Augmented Dickey-Fuller (ADF) and Kwiatkowski-Phillips-Schmidt-Shin (KPSS) tests for a multivariate time series, that may have a panel structure.}
\usage{unirootTest(var.names, unit=NULL, time=NULL, data, box.cox=1, ndiff=0, max.nlags=NULL)}
\arguments{
  \item{var.names}{Character vector including the name of the variables to be differenced.}
  \item{unit}{Character containing the name of the variable that identifies the units of observation. If \code{NULL} (the default), a single unit of observation is assumed.}
  \item{time}{Character containing the name of the time variable, which must be in numeric or date format. If \code{NULL} (the default), data are assumed to be temporally ordered.}
  \item{data}{Object of class \code{data.frame} containing the variables in \code{var.names}, \code{unit} and \code{time}.}
  \item{box.cox}{Named vector with non-negative real values indicating the parameters of the Box-Cox transformation (Box & Cox, 1964) for variables in \code{x.names}. If \code{box.cox} has no names and length greater than one, the same ordering as in \code{x.names} is assumed. If \code{box.cox} has no names and length equal to one, the same parameter is used for all variables in \code{x.names}. Value 0 of the parameter equates to the logarithmic transformation, while value 1 means no transformation. Default is 1 for all variables in \code{x.names}.}
  \item{ndiff}{Named vector with non-negative integer values indicating the number of differences for variables in \code{x.names}. If \code{ndiff} has no names and length greater than one, the same ordering as in \code{x.names} is assumed. If \code{ndiff} has no names and length equal to one, the same number of differences is used for all variables in \code{x.names}. Value 0 means no differencing. Default is 0 for all variables in \code{x.names}.}
  \item{max.nlags}{Non-negative integer value representing the maximum lag length at which to perform the tests (see 'Details'). If \code{NULL} (the default), it is taken as the squared root of the length of the time series.
%\code{trunc((length(x)-1)^(1/3))}
.}
}
\value{One list for each variable in \code{var.names}, each with three components:
\itemize{
\item{\code{statistic}: }{test statistic for each test};
\item{\code{lag.selected}: }{lag length selected for each test};
\item{\code{p.value}: }{p-value of each test, which is a single value if \code{unit} is \code{NULL}, otherwise one value for each unit of observation plus another one indicating the combined p-value};
\item{\code{box.cox}: }{parameter of the Box-Cox transformation for each variable subjected to the tests};
\item{\code{ndiff}: }{order of differencing for each variable subjected to the tests}.
}
}
\note{The first order difference of logarithmic values (\code{box.cox}=0 and \code{ndiff}=1) approximates the relative changes with respect to the previous time point.

If a variable contains negative values, the Box-Cox transformation will be not applied and a warning is returned.

If the number of differencing exceeds \code{n-5}, where \code{n} is the sample size, differencing will be not applied and a warning is returned.
}
\details{The variables subjected to the tests must be quantitative.
Missing values internal to the time series are imputed through linear interpolation, otherwise they are deleted out.

The lag length at which to perform the tests is selected through AIC-based backward elimination starting from the lag length specified in argument \code{max.lag}.

The null hypothesis of the ADF test (Dickey & Fuller, 1981) is the presence of unit roots, while the null hypothesis of the KPSS test (Kwiatkowski \emph{et al.}, 1992) is the absence of unit roots.
Therefore, p-value higher than 0.05 for the ADF test or p-value lower than 0.05 for the KPSS test suggest the presence of unit roots and the need of further differencing.

In case of panel data, p-values are combined according to the method by Demetrescu \emph{et al.} (2006).}
%\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
\references{
M. Demetrescu, U. Hassler, and A. Tarcolea (2006). Combining significance of correlated statistics with application to panel data. \emph{Oxford Bulletin of Economics and Statistics}, 68(5), 647-663. DOI: 10.1111/j.1468-0084.2006.00181.x
  
D. A. Dickey, and W. A. Fuller (1981). Likelihood ratio statistics for autoregressive time series with a unit root. \emph{Econometrica}, 49(4): 1057-1072. DOI: 10.2307/1912517

%K. Hadri, and R. Larsson (2005). Testing for stationarity in heterogeneous panel data where the time dimension is finite. \emph{The Econometrics Journal}, 8(1): 55-69.

D. Kwiatkowski, P. C. B. Phillips, P. Schmidt, and Y. Shin (1992). Testing the null hypothesis of stationarity against the alternative of a unit root. \emph{Journal of Econometrics}, 54(1-3): 159-178. DOI: 10.1016/0304-4076(92)90104-Y
}
\examples{
data(agrisus2020)
varNames <- names(agrisus2020)[4:15]

# tests on variables in level
unirootTest(varNames, unit="Country", time="Year", data=agrisus2020)

# tests on variables in logarithmic differences
unirootTest(varNames, unit="Country", time="Year", data=agrisus2020,
  box.cox=0, ndiff=1)
}