\name{veneto_weather}
\docType{data}
\alias{veneto_weather}

\title{Weekly Aggregated Climate and Meteorological Data in Veneto, Italy}

\description{
The dataset contains weekly averages of various meteorological variables for week 28 of 2021, collected from 75 monitoring stations across the Veneto region in Italy.
}

\format{
An \code{sf} (simple features) object with 75 observations and 7 variables:
\describe{
  \item{ET0}{Evapotranspiration levels (mm) (aggregated weekly).}
  \item{tmax}{Maximum temperature (°C) (weekly max).}
  \item{tmin}{Minimum temperature (°C) (weekly min).}
  \item{hmax}{Maximum humidity (\%) (weekly max).}
  \item{hmin}{Minimum humidity (\%) (weekly min).}
  \item{log_prec}{Log-transformed precipitation values (aggregated weekly).}
  \item{geometry}{Simple feature column containing point locations of monitoring stations (sfc_POINT). The coordinates are in meters (Gauss Boaga - EPSG: 3003)}
}
}

\source{The raw data can be downloaded from the Environmental Protection Agency of Veneto Region
(ARPA Veneto) website.}

\details{
The evapotranspiration levels were estimated by ARPA Veneto according to the Hargreaves model.
}

\usage{
data(veneto_weather)
}

\examples{
data(veneto_weather)
plot(veneto_weather["tmax"])
}

\keyword{datasets}



