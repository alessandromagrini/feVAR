\name{feVAR}
\alias{feVAR}
\title{Estimation of a feVAR model}
\description{Estimation of a fixed effects vector autoregressive model through bias-corrected ordinary least squares (Dhaene & Jochmans, 2016). If there are missing values in the endogenous variables, estimation is performed through the EM algorithm. Both panel, even unbalanced, and non-panel data are allowed.}
\usage{feVAR(var.names, unit=NULL, time=NULL, exogenous=NULL, data, max.nlags=NULL,
  nlags=NULL, trend=c("none","global","unit"), ic=c("bic","aic","aicc","hqic"),
  box.cox=1, ndiff=0, max.ndiff=2, auto.restrict=FALSE, imputation=TRUE,
  em.tol=1e-4, em.maxiter=100, quiet=FALSE)}
\arguments{
  \item{var.names}{Character vector of length 2 or greater including the names of the endogenous variables, which must be quantitative variables.}
  \item{unit}{Character containing the name of the variable that identifies the units of observation. If \code{NULL} (the default), a single unit of observation is assumed.}
  \item{time}{Character containing the name of the time variable, which must be in numeric or date format. If \code{NULL} (the default), data are assumed to be temporally ordered.}
  \item{exogenous}{Character vector including the names of the exogenous variables, which can be either quantitative or qualitative variables. If \code{NULL} (the default), no exogenous variable is included in the model.}
  \item{data}{Object of class \code{data.frame} containing the variables in \code{var.names}, \code{exogenous}, \code{unit} and \code{time}.}
  \item{max.nlags}{The maximum lag order \code{p} of the model, which must be no greater than \code{T-5}, where \code{T} is the minimum number of time points across all units of observations. If an invalid value is provided, it is automatically replaced by the nearest valid one. If \code{NULL} (the default) and \code{nlags} is also \code{NULL}, it is set to \code{sqrt(T)}. If \code{NULL} and \code{nlags} is not \code{NULL}, the best lag order is selected automatically based on information criteria. Ignored if there are missing values for the endogenous variables. See 'Details'.}
  \item{nlags}{The lag order \code{p} of the model. If \code{NULL} (the default), the best lag order is selected based on information criteria (see argument \code{max.nlags}), unless there are missing values for the endogenous variables. In this case, it is set equal to \code{max.nlags}. See 'Details'.}
  \item{trend}{The type of trend, which may be one among 'none' (no trend, i.e. \eqn{\bm{\delta}_i=\bm{0}}), 'global' (same linear trend across all units of observation, i.e. \eqn{\bm{\delta}_i=\bm{\delta}~\forall i}), and 'unit' (unit-specific linear trend, i.e. unconstrained \eqn{\bm{\delta}_i}). If 'unit' and there is only one unit of observation, it is set to 'global'. See 'Details'.}
  \item{ic}{The information criterion to be used in the automated selction of the lag order \code{p} of the model, which should be one among 'bic' (Bayesian information criterion), 'aic' (Akaike information criterion), 'aicc' (corrected Akaike information criterion), and 'hqic' (Hannan & Quinn information criterion). Ignored if \code{nlags} is not \code{NULL}. Default is 'bic'.}
  \item{box.cox}{Named vector including the \eqn{\lambda} values of the Box-Cox transformation (Box & Cox, 1964) for variables in \code{var.names} and in \code{exogenous}. If \code{box.cox} has no names and length greater than one, the same ordering as in \code{c(var.names,exogenous)} is assumed. If \code{box.cox} has no names and length equal to one, the same \eqn{\lambda} value is used for all variables in \code{var.names} and in \code{exogenous}. Valid \eqn{\lambda} values are in the interval \code{[-2,2]}. All invalid values are automatically replaced with the nearest valid ones. Default is 1, meaning no transformation for all variables. If \code{NULL}, \eqn{\lambda} values are determined automatically for each variable through maximum likelihood estimation. See the function \code{\link{unirootTest}} for more details.}
  \item{ndiff}{Named vector including the number of differences for variables in \code{var.names} and in \code{exogenous}. If \code{ndiff} has no names and length greater than one, the same ordering as in \code{c(var.names,exogenous)} is assumed. If \code{ndiff} has no names and length equal to one, the same number of differences is used for all variables in \code{var.names} and in \code{exogenous}. The number of differences cannot exceed \code{T-5}, where \code{T} is the minimum number of time points across all units of observations. All invalid values are automatically replaced with the nearest valid ones. Default is 0, meaning no differencing for all variables. If \code{NULL}, the order of differencing is determined automatically for each variable based on sequential unit root tests. See the function \code{\link{unirootTest}} for more details.}
  \item{max.ndiff}{Maximum order of differencing. Default is 2. Ignored if \code{ndiff} is not \code{NULL}. See 'Note'.}
  \item{auto.restrict}{Logical value indicating whether automated parameter restriction based on information criteria must be performed. Default is \code{FALSE}.}
  \item{imputation}{Logical value indicating whether imputation of missing values must be performed. Default is \code{TRUE}.}
  \item{em.tol}{Tolerance of the EM algorithm, which must be a non-negative value. If negative, it is set to 1e-4 (the default).}
  \item{em.maxiter}{Maximum number of iterations for the EM algorithm, which must be an integer value no lower than 1. If not integer, it is rounded to the nearest integer. Default is 100.}
  \item{quiet}{Logical value indicating whether prompt messages must be suppressed. Default is \code{FALSE}.}
  }
\value{An object of class \code{feVAR} including the following components:
\itemize{
  \item{\code{equations}: }{list of objects of class \code{lm}, one for each endogenous variable, including the estimated equations;}
  \item{\code{call}: }{list with the following components:
  \itemize{
  \item{\code{var.names}, }{including the names of the endogenous variables,}
  \item{\code{unit}, }{including the name of the variable that identifies the units of observation,}
  \item{\code{time}, }{including the name of the time variable,}
  \item{\code{exogenous}, }{including the names of the exogenous variables,}
  \item{\code{trend}, }{including the type of trend (see argument \code{trend}),}
  \item{\code{box.cox}, }{including the parameters of the Box-Cox transformation for each variable (see argument \code{box.cox}),}
  \item{\code{ndiff}, }{including the order of differencing for each variable (see argument \code{ndiff}).}
  }}
  \item{\code{intercepts}: }{matrix containing the estimated intercepts for each unit (by row) and endogenous variable (by column);}
  \item{\code{Beta}: }{point estimation of regression coefficients among endogenous variables (matrix \eqn{B});}
  \item{\code{Sigma}: }{estimated covariance matrix of random errors;}
  \item{\code{data.orig}: }{original data (before Box-Cox transformation and differencing);}
  \item{\code{data.used}: }{data used for model estimation (after Box-Cox transformation and differencing);}
  \item{\code{companion}: }{the companion matrix of the model;}
  \item{\code{eigen.module}: }{module of eigenvalues of the companion matrix. If all of them are lower than one, then the vector autoregressive model is stable.}
  }
}
\details{
Fixed effects vector autoregressive models are an extension of vector autoregressive models to panel data, firstly introduced by Holtz-Eakin et al. (1988).
Let \eqn{Y_1,\ldots,Y_m} be \eqn{m} quantitative endogenous variables and \eqn{X_1,\ldots,X_q} be \eqn{q}, either quantitative or qualitative, exogenous variables.
Also, let \eqn{\bm{y}_{i,t}} be an \eqn{m}-dimensional vector including the values of the endogenous variables observed on unit \eqn{i} at time \eqn{t}, and \eqn{\bm{x}_{i,t}} be a \eqn{q}-dimensional vector including the values of the exogenous variables observed on unit \eqn{i} at time \eqn{t}.
A fixed effect Vector Autoregressive Model (feVAR) of order \eqn{p\in\mathrm{N}_0} is defined as:
\deqn{\bm{y}_{i,t} = \bm{\alpha}_i+\bm{\delta}_i t+\sum_{j=1}^p\bm{B}_j\bm{y}_{i,t-j}+\bm{G}\bm{x}_{i,t}+\bm{\varepsilon}_{i,t}}
where:
\itemize{
\item{\eqn{\bm{\alpha}_i} and \eqn{\bm{\delta}_i} are \eqn{p}-dimensional vectors of parameters specific for unit \eqn{i};}
\item{\eqn{\bm{B}_j} is a \eqn{m\times m} matrix of parameters common to all units, \eqn{j=1,\ldots,p};}
\item{\eqn{\bm{G}} is a \eqn{q\times q} matrix of parameters common to all units;}
\item{\eqn{\bm{\varepsilon}_{i,t}} is a \eqn{p}-dimensional vector of random errors for unit \eqn{i} at time \eqn{t} such that:
\deqn{\mathrm{E}(\bm{\varepsilon}_{i,t})=0~\forall i,t~~~~~~\mathrm{E}(\bm{\varepsilon}_{i,t}\bm{\varepsilon}_{i,t})=\Sigma~\forall i,t~~~~~~\mathrm{E}(\bm{\varepsilon}_{i,s}\bm{\varepsilon}_{i,t})=0~\forall i, s\neq t}}
}
A feVAR model can be viewed as a system of regression equations, one for each endogenous variable, where the regressors include the lags up to \eqn{p} of all endogenous variables and the exogenous variables.
In particular:
\itemize{
\item{parameter \eqn{\bm{\alpha}_i} contains the intercepts for unit \eqn{i} across all equations;}
\item{parameter \eqn{\bm{\delta}_i} contains the linear trend coefficients for unit \eqn{i} across all equations;}
\item{parameter \eqn{\bm{B}_j} contains, by column, the coefficients at lag \eqn{j} of the endogenous variables across all equations;}
\item{parameter \eqn{\bm{G}} contains, by column, the coefficients of the exogenous variables across all equations.}
}

The function \code{feVAR} allows both panel, even unbalanced, and non-panel data (single unit of observation). The bias-corrected ordinary least squares (Dhaene & Jochmans, 2016) is applied only in the former case, while ordinary least squares is used in the latter case.

In case of missing values in the endogenous variables, provided that \code{imputation} is set to \code{TRUE}, estimation is performed through the Expectation-Maximization (EM) algorithm (Dempster et al., 1977) as follows:
\enumerate{
  \item{\emph{Initialization}: }{missing values for each endogenous variable are filled with the sample mean;}
  \item{\emph{Expectation step}: }{Box-Cox transformation and/or differencing are eventually applied and the feVAR model is estimated;}
  \item{\emph{Maximization step}: }{in-sample predictions are used to fill missing values internal to the time series, and backward h-steps ahead forecasts are used to fill missing values at the beginning of the time series;}
  \item{Repeat from step 2) until convergence of the likelihood.}
  }
During this procedure, the lag order of the feVAR model is held constant without possibility to select it automatically based on information criteria. However, it is possible to apply automated parameter restriction (if \code{auto.restrict} is set to \code{TRUE}) after convergence of the algorithm. Estimation resulting from the EM algorithm has the property to maximize the expected value of the likelihood with respect to the distribution of missing data. Currently, the EM algorithm only handles missing values for the endogenous variables, while missing values in the exogenous variables are imputed with the most recent value before estimating the feVAR model.

S3 methods are available for class \code{feVAR} include:
\code{summary}, \code{coef}, \code{confint}, \code{nobs}, \code{residuals}, \code{fitted.values}, \code{logLik}, \code{extractAIC}, \code{AIC}, \code{BIC}, \code{\link{cooks.distance.feVAR}}, \code{\link{predict.feVAR}}.
}
\note{Second-order stationarity (expected value and autocorrelation function independent of time) of all time series is a basic assumption of the model, that is guaranteed if no time series contains unit roots. Before calling the function \code{\link{feVAR}}, the user is strongly recommended to check the absence of unit roots in each time series through the function \code{\link{unirootTest}}.

A popular combination of Box-Cox transformation and differencing is the first order difference of logarithmic values (\code{box.cox=0} and \code{ndiff=1}), which approximates the relative changes with respect to the previous time point:
\deqn{\Delta\log y_t\equiv\log y_t-\log y_{t-1}\approx \frac{y_t}{y_{t-1}}}
%In this case, the feVAR model is formulated as:
%\deqn{\Delta\log\bm{y}_{i,t} = \bm{\delta}_i+\sum_{j=1}^p\bm{B}_j ~\Delta\log\bm{y}_{i,t-j}+\bm{G}~\Delta\log\bm{x}_{i,t}+\bm{\varepsilon}_{i,t}}
%\deqn{\mathrm{E}(\bm{\varepsilon}_{i,t})=0~\forall i,t~~~~~~\mathrm{E}(\bm{\varepsilon}_{i,t}\bm{\varepsilon}_{i,t})=\Sigma~\forall i,t~~~~~~\mathrm{E}(\bm{\varepsilon}_{i,s}\bm{\varepsilon}_{i,t})=0~\forall i, s\neq t}

When the endogenous variable in the l.h.s. of an equation is differenced one time, the original intercepts (parameters denoted by letter \eqn{\alpha}) are canceled out and the new ones become the unit-specific linear trend coefficients (parameters denoted by letter \eqn{\delta}).
If the variable is differenced twice or more, unit-specific linear trend coefficients are canceled out, too, and the model will contain no intercepts.
However, in both cases the original intercepts (parameters denoted by letter \eqn{\alpha}) are recovered using OLS properties and made available in the component \code{intercepts} of the object of class \code{feVAR} (see 'Value').
}
%\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
\references{
%C. W. J. Granger, and P. Newbold (1974). Spurious regressions in econometrics. \emph{Journal of Econometrics}, 2(2): 111-120.
%
H. Akaike (1974). A new look at the statistical model identification. \emph{IEEE Transactions on Automatic Control}, 19(6): 716-723. DOI: 10.1109/TAC.1974.1100705

Cavanaugh, J.E. (1997), Unifying the derivations of the Akaike and corrected Akaike information criteria. \emph{Statistics & Probability Letters}, 31(2): 201-208. DOI: 10.1016/s0167-7152(96)00128-9

Dempster A.P., Laird N.M., Rubin D.B. (1977). Maximum likelihood from incomplete data via the EM algorithm. \emph{Journal of the Royal Statistical Society}, Series B, 39(1): 1-38. 

Dhaene G., Jochmans K. (2016). Bias-corrected estimation of panel vector autoregressions. \emph{Economics Letters}, 145: 98-103. DOI: 10.1016/j.econlet.2016.06.010

Hannan E.J., Quinn B.G. (1979). The determination of the order of an autoregression. \emph{Journal of the Royal Statistical Society}, Series B, 41: 190-195.

Holtz-Eakin D., Newey W., Rosen H.S. (1988). Estimating vector autoregressions with panel data. \emph{Econometrica}, 56(6), 1371-1395. DOI: 10.2307/1913103

Schwarz G.E. (1978). Estimating the dimension of a model. \emph{Annals of Statistics}, 6(2): 461-464. DOI: 10.1214/aos/1176344136
}
\seealso{\code{\link{unirootTest}}; \code{\link{autocorTest}}; \code{\link{plot.feVAR}}; \code{\link{cooks.distance.feVAR}}; \code{\link{IRF}}; \code{\link{predict.feVAR}}; \code{\link{forecastError}}.}
\examples{
##  EXAMPLE WITH A SINGLE UNIT OF OBSERVATION  ##

data(djstock)
x_dj <- colnames(djstock)[-1]

# fit a model with p=1 on data in logarithmic differences ('box.cox'=0 and 'ndiff'=1)
m_dj <- feVAR(var.names=x_dj, time="Date", data=djstock,
  box.cox=0, ndiff=1, nlags=1)
summary(m_dj)  ## summary of estimation

\dontrun{
# automated selection of the lag order (max p=3)
m_dj_auto <- feVAR(var.names=x_dj, time="Date", data=djstock,
  box.cox=0, ndiff=1, ic="bic", max.nlags=3)
summary(m_dj_auto)  ## summary of estimation

# automated parameter restriction
m_dj_auto_r <- feVAR(var.names=x_dj, time="Date", data=djstock,
  box.cox=0, ndiff=1, ic="bic", max.nlags=3, auto.restrict=T)
summary(m_dj_auto_r)  ## summary of estimation
}

##  EXAMPLE WITH PANEL DATA  ##

data(agrisus2020)
x_agr <- colnames(agrisus2020)[4:15]

# fit a model with p=1 on data in logarithmic differences ('box.cox'=0 and 'ndiff'=1)
m_agr <- feVAR(var.names=x_agr, unit="Country", time="Year", exogenous="GDP",
  data=agrisus2020, box.cox=0, ndiff=1, nlags=1)
summary(m_agr)  ## summary of estimation

\dontrun{
# automated selection of the lag order (max p=3)
m_agr_auto <- feVAR(var.names=x_agr, unit="Country", time="Year", exogenous="GDP",
  data=agrisus2020, box.cox=0, ndiff=1, ic="bic", max.nlags=3)
summary(m_agr_auto)  ## summary of estimation

# automated parameter restriction
m_agr_auto_r <- feVAR(var.names=x_agr, unit="Country", time="Year", exogenous="GDP",
  data=agrisus2020, box.cox=0, ndiff=1, ic="bic", max.nlags=3, auto.restrict=TRUE)
summary(m_agr_auto_r)  ## summary of estimation
}
}