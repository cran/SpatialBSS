\name{local_gss_covariance_matrix}
\alias{local_gss_covariance_matrix}

\title{
Computation of Robust Local Covariance Matrices
}

\description{
\code{local_gss_covariance_matrix} computes generalized local sign covariance matrices for a random field based on a given set of spatial kernel matrices. 
}

\usage{
local_gss_covariance_matrix(x, kernel_list, lcov = c('norm', 'winsor', 'qwinsor'), 
                        center = TRUE)
}
\arguments{
  \item{x}{
a numeric matrix of dimension \code{c(n, p)} where the \code{p} columns correspond to the entries of the random field and the \code{n} rows are the observations.
}
  \item{kernel_list}{
a list with spatial kernel matrices of dimension \code{c(n, n)}. This list is usually computed with the function \code{\link{spatial_kernel_matrix}}.
}
  \item{lcov}{
a string indicating which type of robust local covariance matrix to use. Either \code{'norm'} (default), \code{'winsor'} or \code{'qwinsor'}.
}
  \item{center}{
logical. If \code{TRUE} the data \code{x} is robustly centered prior computing the local covariance matrices. Default is \code{TRUE}. See also \code{\link{white_data}}.
}
}

\details{
Generalized local sign matrices are determined by radial functions \eqn{w(l_i)}, where \eqn{l_i = ||x(s_i)-T(x)||} and \eqn{T(x)} is Hettmansperger Randles location estimator (Hettmansperger & Randles, 2002), and kernel functions \eqn{f(d_{i,j})}, where \eqn{d_{i,j}=||s_i - s_j||}. Generalized local sign covariance (gLSCM) matrix is then calculated as

\deqn{gLSCM(f,w) = 1/(n F^{1/2}_{f,n}) \sum_{i,j} f(d_{i,j})  w(l_i)w(l_j)(x(s_i)-T(x)) (x(s_j)-T(x))'} with \deqn{ F_{f,n} = 1 / n \sum_{i,j} f^2(d_{i,j}).}

Three radial functions \eqn{w(l_i)} (Raymaekers & Rousseeuw, 2019) are implemented, the parameter \code{lcov} defines which is used:
\itemize{
  \item \code{'norm'}: \deqn{ w(l_i) = 1/l_i }
  \item \code{'winsor'}: \deqn{ w(l_i) = Q/l_i }
  \item \code{'qwinsor'}: \deqn{ w(l_i) = Q^2/l_i^2.}
}
The cutoff \eqn{Q} is defined as \eqn{Q = l_{(h)}}, where \eqn{l_{(h)}} is \eqn{h}th order statistic of \eqn{\{l_1, ..., l_n\}} and \eqn{h = (n + p + 1)/2}. If the argument \code{center} equals \code{FALSE} then the centering in the above formula for \eqn{gLSCM(f,w)} is not carried out. See also \code{\link{spatial_kernel_matrix}} for details. 
}

\value{
\code{local_gss_covariance_matrix} returns a list with two entries:
  \item{cov_sp_list}{
  List of equal length as the argument \code{kernel_list}. Each list entry is a numeric matrix of dimension \code{c(p, p)} corresponding to a robust local covariance matrix. The list has the attribute \code{'lcov'} which equals the function argument \code{lcov}.
}  
  \item{weights}{
  numeric vector of \code{length(n)} giving the weights for each observation for the robust local covariance estimation.
}
}

\references{
Hettmansperger, T. P., & Randles, R. H. (2002). \emph{A practical affine equivariant multivariate median}. Biometrika, 89 , 851-860. \doi{10.1093/biomet/89.4.851}.

Raymaekers, J., & Rousseeuw, P. (2019). \emph{A generalized spatial sign covariance matrix. Journal of Multivariate Analysis}, 171 , 94-111. \doi{10.1016/j.jmva.2018.11.010}.

Sipila, M., Muehlmann, C. Nordhausen, K. & Taskinen, S. (2024). \emph{Robust second order stationary spatial blind source separation using generalized sign matrices}, Spatial Statistics, 59, 100803, \doi{10.1016/j.spasta.2023.100803}.
}

\seealso{
\code{\link{spatial_kernel_matrix}}, \code{\link{robsbss}}
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
  
  # computing two ring kernel matrices and corresponding 
  # robust local covariance matrices using 'norm' radial function:
  kernel_params_ring <- c(0, 0.5, 0.5, 2)
  ring_kernel_list <- 
    spatial_kernel_matrix(coords, 'ring', kernel_params_ring)
  loc_cov_ring <- 
    local_gss_covariance_matrix(x = field, kernel_list = ring_kernel_list, 
                               lcov = 'norm')
    
  # computing three ball kernel matrices and corresponding 
  # robust local covariance matrices using 'winsor' radial function:
  kernel_params_ball <- c(0.5, 1, 2)
  ball_kernel_list <- 
    spatial_kernel_matrix(coords, 'ball', kernel_params_ball)
  loc_cov_ball <- 
    local_gss_covariance_matrix(x = field, kernel_list = ball_kernel_list, 
                               lcov = 'winsor')
    
  # computing three gauss kernel matrices and corresponding 
  # robust local covariance matrices using 'qwinsor' radial function:
  kernel_params_gauss <- c(0.5, 1, 2)
  gauss_kernel_list <- 
    spatial_kernel_matrix(coords, 'gauss', kernel_params_gauss)
  loc_cov_gauss <- 
    local_gss_covariance_matrix(x = field, kernel_list = gauss_kernel_list, 
                               lcov = 'qwinsor')
}
}

