\name{sbss}
\alias{sbss}
\alias{sbss.default}
\alias{sbss.SpatialPointsDataFrame}
\alias{sbss.sf}

\title{
Spatial Blind Source Separation
}

\description{
\code{sbss} estimates the unmixing matrix assuming a spatial blind source separation model by simultaneous/jointly diagonalizing the covariance matrix and one/many local covariance matrices. These local covariance matrices are determined by spatial kernel functions. Three types of such kernel functions are supported.
}

\usage{
sbss(x, ...)

\method{sbss}{default}(x, coords, kernel_type = c('ring', 'ball', 'gauss'), 
     kernel_parameters, lcov = c('lcov', 'ldiff', 'lcov_norm'), ordered = TRUE, 
     kernel_list = NULL, rob_whitening = FALSE, \dots)
\method{sbss}{SpatialPointsDataFrame}(x, \dots)
\method{sbss}{sf}(x, \dots)
}

\arguments{
  \item{x}{
either a numeric matrix of dimension \code{c(n, p)} where the \code{p} columns correspond to the entries of the random field and the \code{n} rows are the observations, an object of class \code{\link[sp]{SpatialPointsDataFrame}} or an object of class \code{\link[sf]{sf}}. 
}
  \item{coords}{
a numeric matrix of dimension \code{c(n,2)} where each row represents the coordinates of a point in the spatial domain. Only needed if \code{x} is a matrix and the argument \code{kernel_list} is \code{NULL}.
}
  \item{kernel_type}{
a string indicating which kernel function to use. Either \code{'ring'} (default), \code{'ball'} or \code{'gauss'}.
}
  \item{kernel_parameters}{
a numeric vector that gives the parameters for the kernel function. At least length of one for \code{'ball'} and \code{'gauss'} or two for \code{'ring'} kernel, see details.
}
  \item{lcov}{
a string indicating which type of local covariance matrix to use. Either \code{'lcov'} (default), \code{'ldiff'} or \code{'lcov_norm'}. See \code{\link{sbss_asymp}} for details on the latter option.
}
  \item{ordered}{
logical. If \code{TRUE} the entries of the latent field are ordered by the sum of squared (pseudo-)eigenvalues of the diagonalized local covariance matrix/matrices. Default is \code{TRUE}.
}
  \item{kernel_list}{
a list of spatial kernel matrices with dimension \code{c(n,n)}, see details. Usually computed by the function \code{\link{spatial_kernel_matrix}}.
}
  \item{rob_whitening}{
logical. If \code{TRUE} whitening is carried out with respect to the first spatial scatter matrix and not the sample covariance matrix, see details. Default is \code{FALSE}.
}
  \item{\dots}{
further arguments for the fast real joint diagonalization algorithm that jointly diagonalizes the local covariance matrices. See details and \code{\link[JADE]{frjd}}.
}
}

\details{
Three versions of local covariance matrices are implemented, the argument \code{lcov} determines which version is used:
\itemize{
  \item \code{'lcov'}: \deqn{ LCov(f) = 1/n \sum_{i,j} f(d_{i,j}) (x(s_i)-\bar{x}) (x(s_j)-\bar{x})' ,}
  \item \code{'ldiff'}: \deqn{ LDiff(f) = 1/n \sum_{i,j} f(d_{i,j}) (x(s_i)-x(s_j)) (x(s_i)-x(s_j))',}
  \item \code{'lcov_norm'}: \deqn{ LCov^*(f) = 1/(n F^{1/2}_{f,n}) \sum_{i,j} f(d_{i,j}) (x(s_i)-\bar{x}) (x(s_j)-\bar{x})',} with \deqn{ F_{f,n} = 1 / n \sum_{i,j} f^2(d_{i,j}).}
}

Where \eqn{d_{i,j} \ge 0} correspond to the pairwise distances between coordinates, \eqn{x(s_i)} are the \code{p} random field values at location \eqn{s_i}, \eqn{\bar{x}} is the sample mean vector, and the kernel function \eqn{f(d)} determines the locality. The choice \code{'lcov_norm'} is useful when testing for the actual signal dimension of the latent field, see \code{\link{sbss_asymp}} and \code{\link{sbss_boot}}. LDiff matrices are supposed to be more robust when the random field shows a smooth trend. The following kernel functions are implemented and chosen with the argument \code{kernel_type}:
\itemize{
  \item \code{'ring'}: parameters are inner radius \eqn{r_{in}} and outer radius \eqn{r_{out}}, with \eqn{r_{in} < r_{out}}, and \eqn{r_{in}, r_{out} \ge 0}: \deqn{f(d;r_{in}, r_{out}) = I(r_{in} < d \le r_{out})}
  \item \code{'ball'}: parameter is the radius \eqn{r}, with \eqn{r \ge 0}: \deqn{f(d;r) = I(d \le r)}
  \item \code{'gauss'}: Gaussian function where 95\% of the mass is inside the parameter \eqn{r}, with \eqn{r \ge 0}: \deqn{f(d;r) = exp(-0.5 (\Phi^{-1}(0.95) d/r)^2)}
}

The argument \code{kernel_type} determines the used kernel function as presented above, the argument \code{kernel_parameters} gives the corresponding parameters for the kernel function. Specifically, if \code{kernel_type} equals \code{'ball'} or \code{'gauss'} then \code{kernel_parameters} is a numeric vector where each entry corresponds to one parameter. Hence, \code{length(kernel_parameters)} local covariance matrices are used. Whereas, if \code{kernel_type} equals \code{'ring'}, then \code{kernel_parameters} must be a numeric vector of even length where subsequently the inner and outer radii must be given (informally: \code{c(r_in1, r_out1, r_in2, r_out2, \dots)}). In that case \code{length(kernel_parameters) / 2} local covariance matrices are used.

Internally, \code{sbss} calls \code{\link{spatial_kernel_matrix}} to compute a list of \code{c(n,n)} kernel matrices based on the parameters given, where each entry of those matrices corresponds to \eqn{f(d_{i,j})}. Alternatively, such a list of kernel matrices can be given directly to the function \code{sbss} via the \code{kernel_list} argument. This is useful when \code{sbss} is called numerous times with the same coordinates/kernel functions as the computation of the kernel matrices is then done only once prior the actual \code{sbss} calls. For details see also \code{\link{spatial_kernel_matrix}}.

\code{rob_whitening} determines which scatter is used for the whitening step. If \code{TRUE}, whitening is carried out with respect to the scatter matrix defined by the \code{lcov} argument, where the kernel function is given by the argument \code{kernel_type} and the parameters correspond to the first occuring in the argument \code{kernel_parameters}. Therefore, at least two different kernel parameters need to be given. Note that only \eqn{LDiff(f)} matrices are positive definite, hence whitening with \code{'lcov'} is likely to produce an error. If the argument is \code{FALSE}, whitening is carried out with respect to the usual sample covariance matrix. \code{sbss} internally calls \code{\link{white_data}}.

If more than one local covariance matrix is used \code{sbss} jointly diagonalizes these matrices with the function \code{\link[JADE]{frjd}}. \code{\dots} provides arguments for \code{frjd}, useful arguments might be:
\itemize{
  \item \code{eps}: tolerance for convergence.
  \item \code{maxiter}: maximum number of iterations.
}

}

\value{
\code{sbss} returns a list of class \code{'sbss'} with the following entries: 
  \item{s}{
  object of \code{class(x)} containing the estimated source random field.
}  
  \item{coords}{
  coordinates of the observations. Is \code{NULL} if \code{x} was a matrix and the argument \code{kernel_list} was not \code{NULL} at the \code{sbss} call.
}

  \item{w}{
  estimated unmixing matrix.
}

  \item{w_inv}{
  inverse of the estimated unmixing matrix.
}
  \item{pevals}{
  (pseudo-)eigenvalues for each latent field entry.
}  
\item{d}{
  matrix of stacked (jointly) diagonalized local covariance matrices with dimension \code{c(length(kernel_parameters)*p,p)} for \code{'ball'} and \code{'gauss'} kernel or \code{c(}
  \code{(length(kernel_parameters) / 2)*p,p)} for \code{'ring'} kernel.
}
  \item{diags}{
  matrix of dimension \code{c(length(kernel_parameters),p)} where the rows contain the diagonal of the diagonalized local autocovariance matrices.
}
  \item{x_mu}{
  columnmeans of \code{x}.
}

  \item{cov_inv_sqrt}{
  square root of the inverse sample covariance matrix of \code{x}.
}

}

\references{
Muehlmann, C., Filzmoser, P. and Nordhausen, K. (2021), \emph{Spatial Blind Source Separation in the Presence of a Drift}, Submitted for publication. Preprint available at \url{https://arxiv.org/abs/2108.13813}.

Bachoc, F., Genton, M. G, Nordhausen, K., Ruiz-Gazen, A. and Virta, J. (2020), \emph{Spatial Blind Source Separation}, Biometrika, 107, 627-646, \doi{10.1093/biomet/asz079}.

Nordhausen, K., Oja, H., Filzmoser, P., Reimann, C. (2015), \emph{Blind Source Separation for Spatial Compositional Data}, Mathematical Geosciences 47, 753-770, \doi{10.1007/s11004-014-9559-5}.

Muehlmann, C., Bachoc, F., Nordhausen, K. and Yi, M. (2022), \emph{Test of the Latent Dimension of a Spatial Blind Source Separation Model}, to appear in Statistica Sinica, \doi{10.5705/ss.202021.0326}.v
}

\seealso{
\code{\link{spatial_kernel_matrix}}, \code{\link{local_covariance_matrix}}, \code{\link[sp]{sp}}, \code{\link[sf]{sf}}, \code{\link[JADE]{frjd}}
}

\examples{
# simulate coordinates
coords <- runif(1000 * 2) * 20
dim(coords) <- c(1000, 2)
coords_df <- as.data.frame(coords)
names(coords_df) <- c("x", "y")

# simulate random field
if (!requireNamespace('gstat', quietly = TRUE)) {
  message('Please install the package gstat to run the example code.')
} else {
  library(gstat)
  model_1 <- gstat(formula = z ~ 1, locations = ~ x + y, dummy = TRUE, beta = 0, 
                   model = vgm(psill = 0.025, range = 1, model = 'Exp'), nmax = 20)
  model_2 <- gstat(formula = z ~ 1, locations = ~ x + y, dummy = TRUE, beta = 0, 
                   model = vgm(psill = 0.025, range = 1, kappa = 2, model = 'Mat'), 
                   nmax = 20)
  model_3 <- gstat(formula = z ~ 1, locations = ~ x + y, dummy = TRUE, beta = 0, 
                   model = vgm(psill = 0.025, range = 1, model = 'Gau'), nmax = 20)
  field_1 <- predict(model_1, newdata = coords_df, nsim = 1)$sim1
  field_2 <- predict(model_2, newdata = coords_df, nsim = 1)$sim1
  field_3 <- predict(model_3, newdata = coords_df, nsim = 1)$sim1
  field <- as.matrix(cbind(field_1, field_2, field_3))

  # apply sbss with three ring kernels
  kernel_parameters <- c(0, 1, 1, 2, 2, 3)
  sbss_result <- 
    sbss(field, coords, kernel_type = 'ring', kernel_parameters = kernel_parameters)
  
  # print object
  print(sbss_result)
  
  # plot latent field
  plot(sbss_result, colorkey = TRUE, as.table = TRUE, cex = 1)
  
  # predict latent fields on grid
  predict(sbss_result, colorkey = TRUE, as.table = TRUE, cex = 1)
  
  # unmixing matrix
  w_unmix <- coef(sbss_result)
  
  # apply the same sbss with a kernel list
  kernel_list <- spatial_kernel_matrix(coords, kernel_type = 'ring', kernel_parameters)
  sbss_result_k <- sbss(field, kernel_list = kernel_list)
  
  # apply sbss with three ring kernels and local difference matrices
  sbss_result_ldiff <- 
    sbss(field, coords, kernel_type = 'ring', 
         kernel_parameters = kernel_parameters, lcov = 'ldiff')
}

}

\keyword{ multivariate }
\keyword{ spatial }
