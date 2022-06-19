# feVAR
__Estimation and analysis of fixed effects vector autoregressive models__

`feVAR` is an R package implementing bias-corrected least squares estimation for fixed effects vector autoregressive models (Dhaene & Jochmans, 2016), together with several post-estimation functionalities. Both panel and non-panel data are allowed.

G. Dhaene, and K. Jochmans (2016). Bias-corrected estimation of panel vector autoregressions. _Economics Letters_, 145: 98-103. DOI: <a href="https://doi.org/10.1016/j.econlet.2016.06.010">10.1016/j.econlet.2016.06.010</a>

R (The R Project for Statistical Computing) needs to be installed on your system in order
to use the `feVAR` package. R can be downloaded from https://www.r-project.org/.

To install the `feVAR` package, open the console of R and type:
```
install.packages("devtools")  ## do not run if package 'devtools' is already installed
library(devtools)
install_github("alessandromagrini/feVAR")
```

For any request or feedback, please write to <alessandro.magrini@unifi.it> (Alessandro Magrini)

Below, you find some examples of use of the package.
_________________________________________________________________

Load the `feVAR` package
```
library(feVAR)
```
Load data on EU agricultural sustainability indicators
```
data(agrisus2020)
```
Estimation of a fixed effects vector autoregressive model:
```
# names of endogenous variables
varNames <- colnames(agrisus2020)[4:15]

# fit a model with p=1 on data in logarithmic differences ('box.cox'=0 and 'ndiff'=1)
#   exogenous variable: GDP (optional)
m_agr <- feVAR(varNames, unit="Country", time="Year", exogenous="GDP",
  data=agrisus2020, box.cox=0, ndiff=1, nlags=1)

# automated selection of the lag order based on BIC
m_agr_auto <- feVAR(varNames, unit="Country", time="Year", exogenous="GDP",
  data=agrisus2020, box.cox=0, ndiff=1, ic="bic", max.nlags=4)
```
Summary of parameter estimation:
```
summary(m_agr)
```
Stability check:
```
eig0 <- stabilityCheck(m_agr)
sum(eig0>=1)  ## all eigenvalues are lower than 1 -> stable
```
Graphical diagnostics:
```
# time series plot of residuals
residualPlot(m_agr, type="ts", cex.main=1.1)
# autocorrelogram of residuals
residualPlot(m_agr, type="acf", cex.main=1.1)
# normal quantile plot of residuals
residualPlot(m_agr, type="qq", cex.main=1.1)
# fitted versus residuals
residualPlot(m_agr, type="fitVSres", cex.main=1.1)
```
Graphical functionalities for prediction:
```
# predictions for the first four variables on unit 'Italy'
predictPlot(m_agr, unit.id="Italy", var.names=varNames[1:4]) ## in-sample
predictPlot(m_agr, unit.id=15, var.names=1:4)                ## same as before
predictPlot(m_agr, unit.id=15, var.names=1:4, n.ahead=3)     ## 3 steps ahead
```
