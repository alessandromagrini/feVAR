\name{forecastError}
\alias{forecastError}
\title{Forecast error for a feVAR model}
\description{Estimate h-step ahead forecast error for an estimated fixed effects vector autoregressive model.}
\usage{forecastError(model, n.ahead=1, quiet=FALSE)}
\arguments{
  \item{model}{Object of class \code{feVAR}.}
  \item{n.ahead}{The number of steps ahead, which must be no lower than 1 (the default). If an invalid value is provided, it is automatically replaced by the nearest valid one.}
  \item{quiet}{Logical value indicating whether prompt messages must be suppressed. Default is \code{FALSE}.}
  }
\value{A list with one component for each step ahead. Each component of the list is a matrix including the value of forecast error metrics (by column) for each endogenous variable (by row). Let \eqn{\hat{y}_t(h)} be the (pseudo) prediction at \eqn{h} steps ahead performed based on the information up to time \eqn{t}.
Available metrics include:
\itemize{
  \item{'rmse', }{the root mean squared error (RMSE): \deqn{RMSE(h)=\sqrt{\frac{1}{T-h}\sum_{t=1}^{T-h}(y_{t+h}-\hat{y}_t(h))^2}}}
  \item{'mae', }{the mean absolute error (MAE): \deqn{MAE(h)=\frac{1}{T-h}\sum_{t=1}^{T-h}|y_{t+h}-\hat{y}_t(h)|}}
  \item{'rmsse', }{the root mean squared scaled error (RMSSE), equal to the RMSE divided by the RMSE of the naive predictor: \deqn{RMSSE(h)=\sqrt{\frac{\frac{1}{T-h}\sum_{t=1}^{T-h}(y_{t+h}-\hat{y}_t(h))^2}{\frac{1}{T-h}\sum_{t=1}^{T-h}(y_{t+h}-y_t)^2}}}}
  \item{'mase', }{the mean absolute scaled error (MAE), equal to the MAE divided by the MAE of the naive predictor: \deqn{MASE(h)=\frac{\frac{1}{T-h}\sum_{t=1}^{T-h}|y_{t+h}-\hat{y}_t(h)|}{\frac{1}{T-h}\sum_{t=1}^{T-h}|y_{t+h}-y_t|}}}
  \item{'mape', }{the mean absolute percentage error (MAPE), computed only for strictly positive endogenous variables: \deqn{MAPE(h)=\frac{100}{T-h}\sum_{t=1}^{T-h}\frac{|y_{t+h}-\hat{y}_t(h)|}{y_{t+h}}}}
  }
}
%\note{}
%\details{}
%\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
%\references{}
\seealso{\code{\link{feVAR}}.}
\examples{
\dontrun{
data(agrisus2020)
x_agr <- colnames(agrisus2020)[4:15]
m_agr <- feVAR(var.names=x_agr, unit="Country", time="Year", exogenous="GDP",
  data=agrisus2020, box.cox=0, ndiff=1, nlags=1)

forecastError(m_agr, n.ahead=3)  ## forecast error at 3 steps ahead
}
}
