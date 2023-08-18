# feVAR
__Estimation and analysis of fixed effects vector autoregressive models__

`feVAR` is an R package for estimation and analysis of fixed effects Vector AutoRegressive (feVAR) models (see, e.g., Holtz-Eakin et al., 1988). The package implements the bias-corrected ordinary least squares estimation proposed by Dhaene & Jochmans (2016) and the Expectation-Maximization (EM) algorithm (Dempster et al., 1977) to perform estimation in case of missing values. Several numerical and graphical tools to analyse the estimated models are available. Both panel, even unbalanced, and non-panel data are allowed.

References:
- Holtz-Eakin D., Newey W., Rosen H.S. (1988). Estimating vector autoregressions with panel data. _Econometrica_, 56(6), 1371-1395. DOI: <a href="https://doi.org/10.2307/1913103">10.2307/1913103</a>
- Dhaene G., Jochmans K. (2016). Bias-corrected estimation of panel vector autoregressions. _Economics Letters_, 145: 98-103. DOI: <a href="https://doi.org/10.1016/j.econlet.2016.06.010">10.1016/j.econlet.2016.06.010</a>
- Dempster A.P., Laird N.M., Rubin D.B. (1977). Maximum likelihood from incomplete data via the EM algorithm. _Journal of the Royal Statistical Society_, Series B, 39(1): 1-38. 

R (The R Project for Statistical Computing) needs to be installed on your system in order
to use the `feVAR` package. R can be downloaded from https://www.r-project.org/.

To install the `feVAR` package, open the console of R and type:
```
install.packages("devtools")  ## do not run if package 'devtools' is already installed
library(devtools)
install_github("alessandromagrini/feVAR")
```
More details on the methodology can be found in the manual of the package.

Some examples of use of the package are reported below.

For any request or feedback, please write to <alessandro.magrini@unifi.it> (Alessandro Magrini)
_________________________________________________________________

Load the `feVAR` package
```
library(feVAR)
```
Load data on EU agricultural sustainability indicators
```
data(agrisus2020)
```
Unit root tests to check stationarity of the time series
```
# names of endogenous variables
x_agr <- colnames(agrisus2020)[4:15]

# tests on variables in level
#  -> most variables are nonstationary
unirootTest(x_agr, unit="Country", time="Year", data=agrisus2020)

# tests on variables in logarithmic differences
#  -> stationarity confirmed for all variables excepting 'Manager_ratio'
unirootTest(x_agr, unit="Country", time="Year", data=agrisus2020, box.cox=0, ndiff=1)
```
Estimation of a feVAR model
```
# fit a feVAR model with 1 lag on data in logarithmic differences ('box.cox'=0 and 'ndiff'=1)
m_agr <- feVAR(var.names=x_agr, unit="Country", time="Year", exogenous="GDP",
  data=agrisus2020, box.cox=0, ndiff=1, nlags=1)
summary(m_agr)  ## summary of estimation

# automated selection of the lag order (max lag order: 3)
m_agr_auto <- feVAR(var.names=x_agr, unit="Country", time="Year", exogenous="GDP",
  data=agrisus2020, box.cox=0, ndiff=1, ic="bic", max.nlags=3)
summary(m_agr_auto)  ## summary of estimation

# automated parameter restriction
m_agr_auto_r <- feVAR(var.names=x_agr, unit="Country", time="Year", exogenous="GDP",
  data=agrisus2020, box.cox=0, ndiff=1, ic="bic", max.nlags=3, auto.restrict=TRUE)
summary(m_agr_auto_r)  ## summary of estimation
```
Residual diagnostics
```
# time series plot of residuals
plot(m_agr_auto_r, type="ts")
# autocorrelogram of residuals (to check serial uncorrelation)
plot(m_agr_auto_r, type="acf")
# autocorrelogram of squared residuals (to check linearity and homoschedasticity)
plot(m_agr_auto_r, type="acf2")
# fitted versus residuals (to check homoschedasticity)
plot(m_agr_auto_r, type="fitVSres", cex.main=1.1)
# normal quantile plot of residuals
plot(m_agr_auto_r, type="qq", cex.main=1.1)
```
In-sample predictions and h-step ahead forecasts
```
# predictions for the first four variables on unit 'Italy'
pr1 <- predict(m_agr_auto_r, unit.id="Italy", var.names=varNames[1:4])  ## in-sample predictions
plot(pr1)
pr2 <- predict(m_agr_auto_r, unit.id=15, var.names=1:4)                 ## equivalent
plot(pr2)
pr3 <- predict(m_agr_auto_r, unit.id=15, var.names=1:4, n.ahead=3)      ## 3 steps ahead forecasts
plot(pr3)
```
Impulse response functions
```
m_agr_irf <- IRF(m_agr_auto_r, nboot=500)  ## compute IRFs (500 bootstrap resamples)
m_agr_irf$irf[,"TFP_2015",]                ## IRFs generated by 'TFP_2015'
m_agr_irf$irf[,"NetCapital_GVA",]          ## IRFs received by 'NetCapital_GVA'
m_agr_irf$irf[,"TFP_2015","Poverty_rur"]   ## IRFs from 'TFP_2015' to 'Poverty_rur'
```
