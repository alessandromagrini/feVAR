\name{agrisus2020}
\alias{agrisus2020}
\docType{data}
\title{EU agricultural sustainability data}
\description{Data on several indicators covering the economic, social and environmental dimensions of agricultural sustainability for 27 EU countries plus United Kingdom in the period 2004-2020.}
\usage{data(agrisus2020)}
\format{
  Object of class \code{data.frame} with a total of 476 observations on the following 16 variables:
  \describe{
    \item{\code{Country}}{Country name.}
    \item{\code{Country_code}}{Country code.}
    \item{\code{Year}}{Time of measurement (year).}
    \item{\code{TFP_2015}}{Total Factor Productivity (TFP) index of agriculture (2015=100). Source: USDA, Economic Research Service (ERS).}
    \item{\code{NetCapital_GVA}}{Net capital stocks in agriculture (2015 US dollars) to gross value added of agriculture (2015 US dollars). Source: Faostat.}
    \item{\code{Manager_ratio}}{Ratio young/elderly for farm managers (number of managers with less than 35 years by 100 managers with 55 years and over). Source: CMEF.}
    \item{\code{FactorIncome_paid_2010}}{Real income of agricultural factors per paid annual work unit (index 2010=100). Source: Eurostat.}
    \item{\code{EntrIncome_unpaid_2010}}{Net entrepreneurial income of agriculture per unpaid annual work unit (index 2010=100). Source: Eurostat.}
    \item{\code{Income_rur}}{Median equivalised net income in rural areas (power purchasing standard). Source: Eurostat.}
    \item{\code{Unempl_rur}}{At-risk-of-poverty rate in rural areas (\%). Source: Eurostat.}
    \item{\code{Poverty_rur}}{Unemployment rate in rural areas (\%). Source: Eurostat.}
    \item{\code{RenewProd_UAA}}{Production of renewable energy from agriculture (share of utilized agricultural area, \%). Source: CMEF.}
    \item{\code{Organic_p}}{Area under organic cultivation (share of utilized agricultural area, \%). Source: Faostat.}
    \item{\code{GHG_UAA}}{Greenhouse gas emissions due to agriculture (million CO2 equivalent grams per hectare of utilized agricultural area). Source: Faostat.}
    \item{\code{GNB_UAA}}{Gross nitrogen balance (tonnes of nutrient per hectare of utilized agricultural area). Source: Eurostat.}
    \item{\code{GDP}}{Gross domestic product (million 2015 international US dollars). Source: World Bank.}
  }
}
\note{This dataset was employed in Magrini & Giambona (2022) and in Magrini (2022).}
\references{
European Commission (2022). Eurostat database. https://ec.europa.eu/eurostat/data/database

European Commission (2020). Common Monitoring and Evaluation Framework (CMEF) for the CAP 2014-2020. https://agridata.ec.europa.eu/extensions/DataPortal/cmef_indicators.html

Food and Agriculture Organization (2022). Faostat statistical database. https://www.fao.org/faostat/en/#home

Magrini A., Giambona F. (2022). A composite indicator to assess sustainability of agriculture in European Union countries. \emph{Social Indicators Research}, 163: 1003-10362. DOI: 10.1007/s11205-022-02925-6

Magrini A. (2022). Assessment of agricultural sustainability in European Union countries: A group-based multivariate trajectory approach. \emph{Advances in Statistical Analysis}, 106: 673-703. DOI: 10.1007/s10182-022-00437-9

United States Department of Agriculture (2022). International agricultural productivity. USDA, Economic Research Service (ERS), Washington, DC, US. https://www.ers.usda.gov/data-products/international-agricultural-productivity
}