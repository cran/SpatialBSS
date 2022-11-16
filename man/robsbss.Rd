\name{robsbss}
\alias{robsbss}
\alias{robsbss.default}
\alias{robsbss.SpatialPointsDataFrame}
\alias{robsbss.sf}

\title{
Robust Spatial Blind Source Separation
}

\description{
\code{robsbss} is a robust variant of \code{sbss}. It estimates the unmixing matrix assuming a spatial blind source separation model by jointly diagonalizing the Hettmansperger-Randles scatter matrix and one/many generalized local sign covariance matrices. These local generalized sign covariance matrices are determined by spatial kernel functions and radial functions. Three types of such kernel functions and three types of radial functions are supported.
}

\usage{
robsbss(x, ...)

\method{robsbss}{default}(x, coords, kernel_type = c('ring', 'ball', 'gauss'), 
     kernel_parameters, lcov = c('norm', 'winsor', 'qwinsor'), 
     ordered = TRUE, kernel_list = NULL, \dots)
\method{robsbss}{SpatialPointsDataFrame}(x, \dots)
\method{robsbss}{sf}(x, \dots)
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
a string indicating which radial function or type of robust local covariance matrix to use. Either \code{'norm'} (default), \code{'winsor'} or \code{'qwinsor'}. See also \cr 
\code{\link{local_gss_covariance_matrix}} for details.
}
  \item{ordered}{
logical. If \code{TRUE} the entries of the latent field are ordered by the sum of squared (pseudo-)eigenvalues of the diagonalized local covariance matrix/matrices. Default is \code{TRUE}.
}
  \item{kernel_list}{
a list of spatial kernel matrices with dimension \code{c(n,n)}, see details. Usually computed by the function \code{\link{spatial_kernel_matrix}}.
}
  \item{\dots}{
further arguments for the fast real joint diagonalization algorithm that jointly diagonalizes the local covariance matrices. See details and \code{\link[JADE]{frjd}}.
}
}

\details{
\code{robsbss} is a robust variant of \code{sbss} which uses Hettmansperger-Randles (HR) location and scatter estimates (Hettmansperger & Randles, 2002) for whitening (see \code{\link{white_data}} for details) and jointly diagonalizes HR scatter matrix and generalized local sign matrices to estimate the unmixing matrix. The generalized local sign matrices are determined by radial functions \eqn{w(l_i)}, where \eqn{l_i = ||x(s_i)-T(x)||} and \eqn{T(x)} is HR location estimator, and kernel functions \eqn{f(d_{i,j})}, where \eqn{d_{i,j}=||s_i - s_j||}. Generalized local sign covariance (gLSCM) matrix is then calculated as

\deqn{gLSCM(f,w) = 1/(n F^{1/2}_{f,n}) \sum_{i,j} f(d_{i,j})  w(l_i)w(l_j)(x(s_i)-T(x)) (x(s_j)-T(x))'} with \deqn{ F_{f,n} = 1 / n \sum_{i,j} f^2(d_{i,j}).}
Three radial functions (Raymaekers & Rousseeuw, 2019) \eqn{w(l_i)} are implemented, the parameter \code{lcov} defines which is used:
\itemize{
  \item \code{'norm'}: \deqn{ w(l_i) = 1/l_i }
  \item \code{'winsor'}: \deqn{ w(l_i) = Q/l_i }
  \item \code{'qwinsor'}: \deqn{ w(l_i) = Q^2/l_i^2.}
}
The cutoff \eqn{Q} is defined as \eqn{Q = l_{(h)}}, where \eqn{l_{(h)}} is \eqn{h}th order statistic of \eqn{\{l_1, ..., l_n\}} and \eqn{h = (n + p + 1)/2}.
In addition, three kernel functions \eqn{f(d)} are implemented, the parameter \code{kernel_type} defines which is used:
\itemize{
  \item \code{'ring'}: parameters are inner radius \eqn{r_{in}} and outer radius \eqn{r_{out}}, with \eqn{r_{in} < r_{out}}, and \eqn{r_{in}, r_{out} \ge 0}: \deqn{f(d;r_{in}, r_{out}) = I(r_{in} < d \le r_{out})}
  \item \code{'ball'}: parameter is the radius \eqn{r}, with \eqn{r \ge 0}: \deqn{f(d;r) = I(d \le r)}
  \item \code{'gauss'}: Gaussian function where 95\% of the mass is inside the parameter \eqn{r}, with \eqn{r \ge 0}: \deqn{f(d;r) = exp(-0.5 (\Phi^{-1}(0.95) d/r)^2).}
}
The argument \code{kernel_type} determines the used kernel function as presented above, the argument \code{kernel_parameters} gives the corresponding parameters for the kernel function. Specifically, if \code{kernel_type} equals \code{'ball'} or \code{'gauss'} then \code{kernel_parameters} is a numeric vector where each entry corresponds to one parameter. Hence, \code{length(kernel_parameters)} local covariance matrices are used. Whereas, if \code{kernel_type} equals \code{'ring'}, then \code{kernel_parameters} must be a numeric vector of even length where subsequently the inner and outer radii must be given (informally: \code{c(r_in1, r_out1, r_in2, r_out2, \dots)}). In that case \code{length(kernel_parameters) / 2} local covariance matrices are used.

\code{robsbss} calls \code{\link{spatial_kernel_matrix}} internally to compute a list of \code{c(n,n)} kernel matrices based on the parameters given, where each entry of those matrices corresponds to \eqn{f(d_{i,j})}. Alternatively, such a list of kernel matrices can be given directly to the function \code{robsbss} via the \code{kernel_list} argument. This is useful when \code{robsbss} is called numerous times with the same coordinates/kernel functions as the computation of the kernel matrices is then done only once prior the actual \code{robsbss} calls. For details see also \code{\link{spatial_kernel_matrix}}.


If more than one generalized local sign covariance matrix is used \code{robsbss} jointly diagonalizes these matrices with the function \code{\link[JADE]{frjd}}. \code{\dots} provides arguments for \code{frjd}, useful arguments might be:
\itemize{
  \item \code{eps}: tolerance for convergence.
  \item \code{maxiter}: maximum number of iterations.
}

}


\value{
\code{robsbss} returns a list of class \code{'sbss'} with the following entries: 
  \item{s}{
  object of \code{class(x)} containing the estimated source random field.
}  
  \item{coords}{
  coordinates of the observations. Is \code{NULL} if \code{x} was a matrix and the argument \code{kernel_list} was not \code{NULL} at the \code{robsbss} call.
}
  \item{w}{
  estimated unmixing matrix.
}
  \item{weights}{
  numeric vector of \code{length(n)} giving the weights for each observation for the robust local covariance estimation.
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
  robustly estimated columnmeans of \code{x}.
}

  \item{cov_inv_sqrt}{
  square root of the inverse sample covariance matrix of \code{x}.
}

}

\references{
Hettmansperger, T. P., & Randles, R. H. (2002). \emph{A practical affine equivariant multivariate median}. Biometrika, 89 , 851-860. \doi{10.1093/biomet/89.4.851}.

Raymaekers, J., & Rousseeuw, P. (2019). \emph{A generalized spatial sign covariance matrix. Journal of Multivariate Analysis}, 171 , 94-111. \doi{10.1016/j.jmva.2018.11.010}.

Sipila, M., Muehlmann, C. Nordhausen, K. & Taskinen, S. (2022). \emph{Robust second order stationary spatial blind source separation using generalized sign matrices}. Manuscript.
}

\seealso{
\code{\link{spatial_kernel_matrix}}, \code{\link{local_gss_covariance_matrix}}, \code{\link[sp]{sp}}, \code{\link[sf]{sf}}, \code{\link[JADE]{frjd}}
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
  field <- cbind(field_1, field_2, field_3)
  # Generate 5 \% local outliers to data
  field_cont <- gen_loc_outl(field, coords, radius = 2, 
                             swap_order = "regular")[,1:3]
  X <- as.matrix(field_cont)
  
  # apply sbss with three ring kernels
  kernel_parameters <- c(0, 1, 1, 2, 2, 3)
  robsbss_result <- 
    robsbss(X, coords, kernel_type = 'ring', kernel_parameters = kernel_parameters)
  
  # print object
  print(robsbss_result)
  
  # plot latent field
  plot(robsbss_result, colorkey = TRUE, as.table = TRUE, cex = 1)
  
  # predict latent fields on grid
  predict(robsbss_result, colorkey = TRUE, as.table = TRUE, cex = 1)
  
  # unmixing matrix
  w_unmix <- coef(robsbss_result)
}
}

\keyword{ multivariate }
\keyword{ spatial }
\keyword{ robust }
