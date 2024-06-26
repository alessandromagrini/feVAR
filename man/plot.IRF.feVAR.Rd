\name{plot.IRF.feVAR}
\alias{plot.IRF.feVAR}
\title{Plot method for class IRF.feVAR}
\description{Display impulse response functions (IRFs) for an estimated fixed effects vector autoregressive model.}
\usage{
\method{plot}{IRF.feVAR}(x, from=NULL, to=NULL, n.ahead=NULL, labels=NULL, ylim=NULL, add.grid=TRUE,
  add.titles=TRUE, cex.points=0, mfrow=NULL, mar=c(3.5,3.5,2,2),
  mgp=c(2.3,0.8,0), ...)
}  
\arguments{
  \item{x}{Object of class \code{IRF.feVAR}.}
  \item{from}{Vector containing the names or the numerical identifiers of the endogenous variables originating the IRFs. If \code{NULL} (the default), all the endogenous variables are considered.}
  \item{to}{Character vector containing the names of the endogenous variables receiving the IRFs. Numerical identifiers can be provided instead of names. If \code{NULL} (the default), all the endogenous variables are considered.}
  \item{n.ahead}{The number of steps ahead, which must be no lower than 1 and no greater than the number of steps ahead in \code{x} (the default). If an invalid value is provided, it is automatically replaced by the nearest valid one.}
  \item{labels}{Character vector containing the labels for the endogenous variables. If the vector has names, labels are assigned based on them, otherwise they are assigned in the same order of the endogenous variables. If \code{NULL} (the default), the names of endogenous variables are used as labels.}
  \item{ylim}{Numerical vector of length 2 indicating the limits of the y-axis.}
  \item{add.grid}{Logical value indicating whether the grid must be added. Default is \code{TRUE}.}
  \item{add.titles}{Logical value indicating whether a title must be added to each graphic. Default is \code{TRUE}.}
  \item{cex.points}{Expansion factor for points. If 0, points are suppressed. Default is 0.}
  \item{mfrow}{Numerical vector of the form \code{c(nr,nc)}: graphics are displayed in an \code{nr}-by-\code{nc} array on the device by rows. If \code{NULL} (the default), it is determined automatically.}
  \item{mar}{Numerical vector of length 4 indicating the margin sizes in the following order: bottom, left, top, right. Default is c(3.5,3.5,2,2).}
  \item{mgp}{Numerical vector of length 3 indicating the location of axis labels, tick mark labels and tick marks, respectively, relative to the plot window. Default is c(2.3,0.8,0).}
  \item{...}{Further parameters to be passed to the generic method \code{plot}.}
}
%\value{}
%\details{}
%\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
%\references{}
\seealso{\code{\link{feVAR}}; \code{\link{IRF}}.}
\examples{
\dontrun{
data(agrisus2020)
x_agr <- colnames(agrisus2020)[4:15]
m_agr <- feVAR(var.names=x_agr, unit="Country", time="Year", exogenous="GDP",
  data=agrisus2020, box.cox=0, ndiff=1, nlags=1)
m_agr_irf <- IRF(m_agr, nboot=500)    ## compute IRFs

# display IRFs generated by 'TFP_2015'
plot(m_agr_irf, from="TFP_2015")
plot(m_agr_irf, from=1)  ## equivalent

# display IRFs received by 'NetCapital_GVA'
plot(m_agr_irf, to="NetCapital_GVA")
plot(m_agr_irf, from=2)  ## equivalent

# display IRFs generated by 'TFP_2015' and 'NetCapital_GVA' and
#  received by 'Poverty_rur', 'Organic_p' and 'GHG_UAA'
plot(m_agr_irf, from=c("TFP_2015","NetCapital_GVA"),
  to=c("Poverty_rur", "Organic_p" and "GHG_UAA"))
plot(m_agr_irf, from=c(1,2), to=c(8,10,11))  ## equivalent
}
}