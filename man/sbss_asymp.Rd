\name{sbss_asymp}
\alias{sbss_asymp}
\alias{sbss_asymp.default}
\alias{sbss_asymp.SpatialPointsDataFrame}
\alias{sbss_asymp.sf}

\title{
Asymptotic Test for the White Noise Dimension in a Spatial Blind Source Separation Model
}

\description{
\code{sbss_asymp} uses asymptotic theory for the spatial blind source separation (SBSS) methodology to test if the last \eqn{p - q} entries of the latent random field are white noise assuming that the \eqn{p}-variate observed random field follows a SBSS model.
}

\usage{
sbss_asymp(x, ...)

\method{sbss_asymp}{default}(x, coords, q, kernel_parameters, 
           kernel_list = NULL, \dots)
\method{sbss_asymp}{SpatialPointsDataFrame}(x, \dots)
\method{sbss_asymp}{sf}(x, \dots)
}

\arguments{
  \item{x}{
either a numeric matrix of dimension \code{c(n, p)} where the \code{p} columns correspond to the entries of the random field and the \code{n} rows are the observations, an object of class \code{\link[sp]{SpatialPointsDataFrame}} or an object of class \code{\link[sf]{sf}}. 
}
  \item{coords}{
a numeric matrix of dimension \code{c(n,2)} where each row represents the coordinates of a point in the spatial domain. Only needed if \code{x} is a matrix and the argument \code{kernel_list} is \code{NULL}.
}
  \item{q}{
an integer between \code{0} and \code{p - 1} specifying the number of hypothetical signal components (null hypothesis) in the latent random field.
}
  \item{kernel_parameters}{
a numeric vector that gives the parameters for the ring kernel function. At least length of two, see details.
}
  \item{kernel_list}{
a list of spatial kernel matrices with dimension \code{c(n,n)}, see details. Usually computed by the function \code{\link{spatial_kernel_matrix}}.
}
  \item{\dots}{
further arguments for the fast real joint diagonalization algorithm that jointly diagonalizes the local covariance matrices. See details and \code{\link[JADE]{frjd}}.
}
}

\details{
This function uses the SBSS methodology in conjunction with local covariance matrices based on ring kernel functions to estimate the \eqn{p}-variate latent random field \eqn{s = x^{wh} w}, where \eqn{x^{wh}} is the whitened version of the data and \eqn{w} is the estimated unmixing matrix. The considered (adapted) local covariance matrices write as
\deqn{ LCov^* = 1/(n F^{1/2}_n) \sum_{i,j} I(r_i < d_{i,j} \le r_o) (x(s_i)-\bar{x}) (x(s_j)-\bar{x})'}
with
\deqn{ F_n = 1 / n \sum_{i,j} I(r_i < d_{i,j} \le r_o).}
Where \eqn{d_{i,j} \ge 0} correspond to the pairwise distances between coordinates, \eqn{x(s_i)} are the \code{p} random field values at location \eqn{s_i} (which is the i-th row of the argument \code{x} and the location corresponds to the i-th row of the argument \code{coords}) and \eqn{\bar{x}} is the sample mean vector. The function argument \code{kernel_parameters} determines the parameters of the used ring kernel functions or alternatively a list of kernel matrices can be given with the argument \code{kernel_list}, see \code{\link{sbss}} for details. 

The null hypothesis specified with the argument \code{q} states that the last \eqn{p - q} components of the estimated latent field are white noise. The method orders the components of the latent field by the order of the decreasing sums of squares of the corresponding (pseudo-)eigenvalues of the local covariance matrices produced by the joint diagonalization algorithm (or the eigendecomposition if only one local covariance matrix is used). Under the null the lower right \eqn{(p - q) * (p - q)} block matrices of the jointly diagonalized local covariance matrices equal zero matrices. Therefore, the sum of their squared norms \eqn{m} is used as test statistic.  

This function conducts the hypothesis test using the asymptotic null distribution of \eqn{m}, a chi-squared distribution with \eqn{k(p - q)(p - q + 1)/2} degrees of freedom (\eqn{k} is the number jointly diagonalized local covariance matrices).

If more than one local covariance matrix is used \code{sbss_asymp} jointly diagonalizes these matrices with the function \code{\link[JADE]{frjd}}. \code{\dots} provides arguments for \code{frjd}, useful arguments might be:
\itemize{
  \item \code{eps}: tolerance for convergence.
  \item \code{maxiter}: maximum number of iterations.
}

}

\value{
\code{sbss_asymp} returns a list of class \code{'sbss_test'} inheriting from the classes \code{'htest'} and \code{'sbss'} with the following entries: 
  \item{alternative}{
  a string containing the alternative hypothesis.
} 

  \item{method}{
  a string which indicates which test methods was used.
} 

  \item{data.name}{
  a string specifying the name of the used data.
} 
  \item{statistic}{
  the value of the test statistic.
} 

  \item{parameters}{
  degrees of freedom for the asymptotic chi-squared distribution of the test statistic under the null hypothesis.
} 

  \item{p.value}{
  the p-value of the test.
} 

  \item{s}{
  object of \code{class(x)} containing the estimated source random field.
}  

  \item{coords}{
  coordinates of the observations. Is \code{NULL} if \code{x} was a matrix and the argument \code{kernel_list} was not \code{NULL} at the \code{sbss_asymp} call.
}

  \item{w}{
  estimated unmixing matrix.
}

  \item{w_inv}{
  inverse of the estimated unmixing matrix.
}

  \item{d}{
  matrix of stacked (jointly) diagonalized local covariance matrices with dimension \code{c((length(kernel_parameters) / 2)*p,p)}.
}

  \item{x_mu}{
  columnmeans of \code{x}.
}

  \item{cov_inv_sqrt}{
  square root of the inverse sample covariance matrix of \code{x}.
}

}

\references{
Muehlmann, C., Bachoc, F., Nordhausen, K. and Yi, M. (2021), \emph{Test of the Latent Dimension of a Spatial Blind Source Separation Model}, Submitted for publication. Preprint available at \url{https://arxiv.org/abs/2011.01711}.

}

\seealso{
\code{\link{sbss}}, \code{\link{spatial_kernel_matrix}}, \code{\link{local_covariance_matrix}}, \code{\link[sp]{sp}}, 
\code{\link[sf]{sf}}, \code{\link[JADE]{frjd}}
}

\examples{
# simulate coordinates
n <- 1000
coords <- runif(n * 2) * 20
dim(coords) <- c(n, 2)

# simulate random field
if (!requireNamespace('RandomFields', quietly = TRUE)) {
  message('Please install the package RandomFields to run the example code.')
} else {
  RandomFields::RFoptions(spConform = FALSE)
  field_1 <- RandomFields::RFsimulate(model = RandomFields::RMexp(), 
                                      x = coords)
  field_2 <- RandomFields::RFsimulate(model = RandomFields::RMspheric(), 
                                      x = coords)
  field_3 <- rnorm(n)
  field_4 <- rnorm(n)
  
  latent_field <- cbind(field_1, field_2, field_3, field_4)
  mixing_matrix <- matrix(rnorm(16), 4, 4)
  observed_field <- latent_field \%*\% t(mixing_matrix)
  
  # apply the asymptotic test for a hypothetical latent white noise dimension of q
  # q can lie between 0 and 3 in this case
  # using one ring kernel function and the null hypothesis q = 1
  asymp_res_1 <- 
    sbss_asymp(observed_field, coords, q = 1, kernel_parameters = c(0, 1))
  
  # using two ring kernel functions and the null hypothesis q = 3
  asymp_res_2 <- 
    sbss_asymp(observed_field, coords, q = 3, kernel_parameters = c(0, 1, 1, 2))
    
  # the result is of class sbss_test which is inherited from htest and sbss
  # print object (print method for an object of class htest)
  print(asymp_res_1)
  print(asymp_res_2)
  
  # plot latent field (plot method for an object of class sbss)
  plot(asymp_res_1, colorkey = TRUE, as.table = TRUE, cex = 1)
  
  # predict latent fields on grid (predict method for an object of class sbss)
  predict(asymp_res_1, colorkey = TRUE, as.table = TRUE, cex = 1)
  
  # unmixing matrix (coef method for an object of class sbss)
  w_unmix <- coef(asymp_res_1)
}

}

\keyword{ htest }
\keyword{ multivariate }
\keyword{ spatial }
