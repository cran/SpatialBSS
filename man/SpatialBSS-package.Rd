\name{SpatialBSS-package}
\alias{SpatialBSS-package}
\docType{package}
\title{
Blind Source Separation for Multivariate Spatial Data
}

\description{
Blind source separation for multivariate spatial data based on simultaneous/joint diagonalization of local covariance matrices. This package is an implementation of the methods described in Nordhausen, Oja, Filzmoser and Reimann (2015) <doi:10.1007/s11004-014-9559-5>, Bachoc, Genton, Nordhausen, Ruiz-Gazen and Virta (2020) <doi:10.1093/biomet/asz079> and
Muehlmann, Bachoc and Nordhausen (2022) <doi:10.1016/j.spasta.2021.100574> as well as some related methods.
}


\details{
\tabular{ll}{
Package: \tab SpatialBSS\cr
Type: \tab Package\cr
Version: \tab 0.16-0 \cr
Date: \tab 2025-03-26\cr
License: \tab GPL (>= 2)\cr
}

This package provides functions to solve the Blind Source Separation problem for multivariate spatial data. These methods are designed to work with random fields that are observed on irregular locations. Moreover, the random field is assumed to show weak second order stationarity. The main functions of this package are:
\itemize{
  \item \code{\link{sbss}}: This function derives a set of local scatter matrices that are based on spatial kernel functions, where the spatial kernel functions can be chosen. Then this set of local covariance matrices as well as the sample covariance matrix are simultaneously/jointly diagonalized. Local covariance matrices as well as local difference matrices are implemented.
  \item \code{\link{sbss_asymp}}, \code{\link{sbss_boot}}: These functions test for white noise components in the estimated latent field estimated by the \code{\link{sbss}} function based on asymptotic results or bootstrap inference principles.
  \item \code{\link{snss_sd}}, \code{\link{snss_jd}}, and \code{\link{snss_sjd}}: These functions estimate the latent random field assuming a spatial non-stationary source separation model. This is done by splitting the domain into a number of sub-domains and diagonalizing the corresponding covariance and/or local covariance matrices for each sub-domain.
  \item \code{\link{robsbss}}: Uses robust estimates of local covariance matrices to solve the SBSS problem.
}


Joint diagonalization is computed with the \code{\link[JADE]{frjd}} (fast real joint diagonalization) algorithm from the package \code{\link[JADE:JADE-package]{JADE}}.

The random field can be either a pair of numeric matrices giving the coordinates and field values or an object of class \code{\link[sp]{SpatialPointsDataFrame}} or \code{\link[sf]{sf}}.

}

\author{
Christoph Muehlmann, Mika Sipila, Claudia Cappello, Sandra De Iaco, Klaus Nordhausen, Sara Taskinen, Joni Virta

Maintainer: 
Klaus Nordhausen \email{klausnordhausenR@gmail.com}
}

\references{
Bachoc, F., Genton, M. G, Nordhausen, K., Ruiz-Gazen, A. and Virta, J. (2020), \emph{Spatial Blind Source Separation}, Biometrika, 107, 627-646, \doi{10.1093/biomet/asz079}.

Muehlmann, C., Bachoc, F. and Nordhausen, K. (2022), \emph{Blind Source Separation for Non-Stationary Random Fields}, Spatial Statistics, 47, 100574, \doi{10.1016/j.spasta.2021.100574}.

Muehlmann, C., Bachoc, F., Nordhausen, K. and Yi, M. (2024), \emph{Test of the Latent Dimension of a Spatial Blind Source Separation Model}, Statistica Sinica, 34, 837-865, \doi{10.5705/ss.202021.0326}.

Muehlmann, C., Filzmoser, P. and Nordhausen, K. (2024), \emph{Spatial Blind Source Separation in the Presence of a Drift}, Austrian Journal of Statistics, 53, 48-68, \doi{10.17713/ajs.v53i2.1668}.

Nordhausen, K., Oja, H., Filzmoser, P. and Reimann, C. (2015), \emph{Blind Source Separation for Spatial Compositional Data}, Mathematical Geosciences 47, 753-770, \doi{10.1007/s11004-014-9559-5}.

Sipila, M., Muehlmann, C. Nordhausen, K. and Taskinen, S. (2024). \emph{Robust second order stationary spatial blind source separation using generalized sign matrices}, Spatial Statistics, 59, 100803, \doi{10.1016/j.spasta.2023.100803}.
}

\keyword{ package }
