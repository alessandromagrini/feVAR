\name{feVAR}
\alias{feVAR}
\title{Estimation of a fixed effect vector autoregressive model}
\description{Estimation of a fixed effect vector autoregressive model through bias-corrected ordinary least squares (Dhaene & Jochmans, 2016). Both panel and non-panel data are allowed.}
\usage{feVAR(var.names, unit=NULL, time=NULL, exogenous=NULL, data, max.nlags=NULL, nlags=NULL,
  add.intercept=TRUE, box.cox=1, ndiff=0, local.adapt=FALSE, ic="bic", quiet=FALSE)}
\arguments{
  \item{var.names}{Character vector of length 2 or greater including the names of the endogenous variables, which must be quantitative variables.}
  \item{unit}{Character containing the name of the variable that identifies the units of observation. If \code{NULL} (the default), a single unit of observation is assumed.}
  \item{time}{Character containing the name of the time variable, which must be in numeric or date format. If \code{NULL} (the default), data are assumed to be temporally ordered.}
  \item{exogenous}{Character vector including the names of the exogenous variables, which can be either quantitative or qualitative variables. If \code{NULL} (the default), no exogenous variable is included in the model.}
  \item{data}{Object of class \code{data.frame} containing the variables in \code{var.names}, \code{exogenous}, \code{unit} and \code{time}.}
  \item{max.nlags}{Non-negative value indicating the maximum lag order \code{p} of the model. If \code{NULL} (the default) and argument \code{nlags} is not \code{NULL}, the maximum lag order is determined automatically and the best one is selected based on information criteria.}
  \item{nlags}{The lag order \code{p} of the model. If \code{NULL} (the default), the best lag order is selected based on information criteria (see argument \code{max.nlags}).}
  \item{add.intercept}{Logical value indicating whether the intercept should be included in the model. Default is \code{TRUE}. See 'Note'.}
  \item{box.cox}{Named vector with non-negative real values indicating the parameters of the Box-Cox transformation (Box & Cox, 1964) for variables in \code{var.names} and in \code{exogenous}. If \code{box.cox} has no names and length greater than one, the same ordering as in \code{c(var.names,exogenous)} is assumed. If \code{box.cox} has no names and length equal to one, the same parameter is used for all variables in \code{var.names} and in \code{exogenous}. Value 0 of the parameter equates to the logarithmic transformation, while value 1 means no transformation. Default is 1 for all variables in \code{var.names} and in \code{exogenous}. See 'Note'.}
  \item{ndiff}{Named vector with non-negative integer values indicating the number of differences for variables in \code{var.names} and in \code{exogenous}. If \code{ndiff} has no names and length greater than one, the same ordering \code{c(var.names,exogenous)} is assumed. If \code{ndiff} has no names and length equal to one, the same number of differences is used for all variables in \code{var.names} and in \code{exogenous}. Value 0 means no differencing. Default is 0 for all variables in \code{var.names} and in \code{exogenous}. See 'Note'.}
  \item{local.adapt}{Logical value indicating whether the selection of the lag order should be performed locally (i.e., for each variable in each equation). Default is \code{FALSE}.}  
  \item{ic}{The information criterion to be used in the automated selction of the lag order \code{p} of the model, which should be one among 'bic' (Bayesian information criterion), 'hqic' (Hannan & Quinn information criterion), and 'aic' (Akaike information criterion). Ignored if \code{nlags} is not \code{NULL}. Default is 'bic'.}
  \item{quiet}{Logical value indicating whether prompt messages should be suppressed. Default is \code{FALSE}.}
  }
\value{An object of class \code{feVAR} including the following components:
\itemize{
  \item{\code{models}: }{list of objects of class \code{lm}, including ordinary least squares estimation of each equation;}
  \item{\code{call}: }{list including the main arguments passed to the function \code{feVAR};}
  \item{\code{nlags}: }{the lag order \code{p} of the model;}
  \item{\code{intercepts}: }{matrix including the estimated intercepts for each unit (by row) and endogenous variable (by column);}
  \item{\code{Beta}: }{point estimation of regression coefficients among endogenous variables (matrix \eqn{B});}
  \item{\code{Gamma}: }{point estimation of regression coefficients of exogenous variables on endogenous variables (matrix \eqn{G});}  
  \item{\code{Sigma}: }{estimated covariance matrix of random errors;}
  \item{\code{ic}: }{value of information criteria for the estimated model;}
  \item{\code{data}: }{original data (before Box-Cox transformation and differencing);}
  \item{\code{data.used}: }{data used for parameter estimation (after Box-Cox transformation and differencing).}
  }
}
\details{
S3 methods available for class \code{feVAR} include:

- \code{summary} provides a summary of parameter estimation;

- \code{coef} returns point estimation of parameters for each equation;

- \code{residuals} returns the residuals;

- \code{fitted.values} returns the fitted values;

- \code{confint} returns confidence intervals of parameters for each equation;

- \link{predict.feVAR} performs predictions.

Also, the function \code{stabilityCheck} can be used to check whether an estimated vector autoregressive model is stationary (stable).
}
\note{
Second-order stationarity (expected value and autocorrelation function independent of time) of all time series is a basic assumption of the model, that is guaranteed if no time series contains unit roots. Before calling the function \link{feVAR}, the user is strongly recommended to check the absence of unit roots in each time series through the function \link{unirootTest}.
Stationarity of an estimated model (stability) can be check through the function \link{stabilityCheck}.

The first order difference of logarithmic values (\code{box.cox}=0 and \code{ndiff}=1) approximates the relative changes with respect to the previous time point.

If a variable contains negative values, the Box-Cox transformation will be not applied and a warning is returned.

If the number of differencing exceeds \code{n-5}, where \code{n} is the sample size, differencing will be not applied and a warning is returned.

If an endogenous variable is not differenced, the intercepts capture the effect of time-invariant factors besides exogenous variables. Otherwise, the intercept represents the drift (i.e., the coefficient of a linear deterministic trend). To avoid the drift, argument \code{add.intercept} should be set to \code{FALSE} when differencing is applied to all endogenous variables.
}
%\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
\references{
%C. W. J. Granger, and P. Newbold (1974). Spurious regressions in econometrics. \emph{Journal of Econometrics}, 2(2): 111-120.
%
H. Akaike (1974). A new look at the statistical model identification. \emph{IEEE Transactions on Automatic Control}, 19(6): 716-723. DOI: 10.1109/TAC.1974.1100705

G. Dhaene, and K. Jochmans (2016). Bias-corrected estimation of panel vector autoregressions. \emph{Economics Letters}, 145: 98-103. DOI: 10.1016/j.econlet.2016.06.010

E. J. Hannan, and B. G. Quinn (1979). The determination of the order of an autoregression. \emph{Journal of the Royal Statistical Society}, Series B, 41: 190-195.

G. E. Schwarz (1978). Estimating the dimension of a model. \emph{Annals of Statistics}, 6(2): 461-464. DOI: 10.1214/aos/1176344136
}
\seealso{\link{unirootTest}; \link{stabilityCheck}; \link{predict.feVAR}; \link{predictPlot}; \link{residualPlot}.
}
\examples{
##  EXAMPLE WITH A SINGLE UNIT OF OBSERVATION  ##

data(djstock)
varNames <- colnames(djstock)[-1]

# fit a model with p=1 on data in logarithmic differences ('box.cox'=0 and 'ndiff'=1)
m_dj <- feVAR(varNames, time="Date", data=djstock,
  box.cox=0, ndiff=1, nlags=1)
summary(m_dj)  ## summary of parameter estimation

# automated selection of the lag order
m_dj_auto <- feVAR(varNames, time="Date", data=djstock,
  box.cox=0, ndiff=1, ic="bic", max.nlags=5)
summary(m_dj_auto)  ## summary of parameter estimation


##  EXAMPLE WITH PANEL DATA  ##

data(agrisus2020)
varNames <- colnames(agrisus2020)[4:15]

# fit a model with p=1 on data in logarithmic differences ('box.cox'=0 and 'ndiff'=1)
m_agr <- feVAR(varNames, unit="Country", time="Year", exogenous="GDP",
  data=agrisus2020, box.cox=0, ndiff=1, nlags=1)
summary(m_agr)  ## summary of parameter estimation

# automated selection of the lag order
m_agr_auto <- feVAR(varNames, unit="Country", time="Year", exogenous="GDP",
  data=agrisus2020, box.cox=0, ndiff=1, ic="bic", max.nlags=4)
summary(m_agr_auto)  ## summary of parameter estimation
}