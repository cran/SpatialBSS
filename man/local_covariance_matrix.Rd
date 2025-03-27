\name{local_covariance_matrix}
\alias{local_covariance_matrix}

\title{
Computation of Local Covariance Matrices
}

\description{
\code{local_covariance_matrix} computes local covariance matrices for a random field based on a given set of spatial kernel matrices.
}

\usage{
local_covariance_matrix(x, kernel_list, lcov = c('lcov', 'ldiff', 'lcov_norm'), 
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
a string indicating which type of local covariance matrix to use. Either \code{'lcov'} (default) or \code{'ldiff'}.
}
  \item{center}{
logical. If \code{TRUE} the data \code{x} is centered prior computing the local covariance matrices. Default is \code{TRUE}.
}
}

\details{
Two versions of local covariance matrices are implemented, the argument \code{lcov} determines which version is used:
\itemize{
  \item \code{'lcov'}: \deqn{ LCov(f) = 1/n \sum_{i,j} f(d_{i,j}) (x(s_i)-\bar{x}) (x(s_j)-\bar{x})' ,}
  \item \code{'ldiff'}: \deqn{ LDiff(f) = 1/n \sum_{i,j} f(d_{i,j}) (x(s_i)-x(s_j)) (x(s_i)-x(s_j))',}
    \item \code{'lcov_norm'}: \deqn{ LCov^*(f) = 1/(n F^{1/2}_{f,n}) \sum_{i,j} f(d_{i,j}) (x(s_i)-\bar{x}) (x(s_j)-\bar{x})',} with \deqn{ F_{f,n} = 1 / n \sum_{i,j} f^2(d_{i,j}).}
}

Where \eqn{d_{i,j} \ge 0} correspond to the pairwise distances between coordinates, \eqn{x(s_i)} are the \code{p} random field values at location \eqn{s_i}, \eqn{\bar{x}} is the sample mean vector, and the kernel function \eqn{f(d)} determines the locality. The choice \code{'lcov_norm'} is useful when testing for the actual signal dimension of the latent field, see \code{\link{sbss_asymp}} and \code{\link{sbss_boot}}. The function \code{local_covariance_matrix} computes local covariance matrices for a given random field and given spatial kernel matrices, the type of computed local covariance matrices is determined by the argument \code{'lcov'}. If the argument \code{center} equals \code{FALSE} then the centering in the above formula for \eqn{LCov(f)} is not carried out. See also \code{\link{spatial_kernel_matrix}} for details. 
}

\value{
\code{local_covariance_matrix} returns a list of equal length as the argument \code{kernel_list}. Each list entry is a numeric matrix of dimension \code{c(p, p)} corresponding to a local covariance matrix. The list has the attribute \code{'lcov'} which equals the function argument \code{lcov}.
}

\references{
Muehlmann, C., Filzmoser, P. and Nordhausen, K. (2024), \emph{Spatial Blind Source Separation in the Presence of a Drift}, Austrian Journal of Statistics, 53, 48-68, \doi{10.17713/ajs.v53i2.1668}.

Bachoc, F., Genton, M. G, Nordhausen, K., Ruiz-Gazen, A. and Virta, J. (2020), \emph{Spatial Blind Source Separation}, Biometrika, 107, 627-646, \doi{10.1093/biomet/asz079}.

}

\seealso{
\code{\link{spatial_kernel_matrix}}, \code{\link{sbss}}
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

  # computing two ring kernel matrices and corresponding local covariance matrices
  kernel_params_ring <- c(0, 0.5, 0.5, 2)
  ring_kernel_list <- 
    spatial_kernel_matrix(coords, 'ring', kernel_params_ring)
  loc_cov_ring <- 
    local_covariance_matrix(x = field, kernel_list = ring_kernel_list)
    
  # computing two ring kernel matrices and corresponding local difference matrices
  kernel_params_ring <- c(0, 0.5, 0.5, 2)
  ring_kernel_list <- 
    spatial_kernel_matrix(coords, 'ring', kernel_params_ring)
  loc_cov_ring <- 
    local_covariance_matrix(x = field, kernel_list = ring_kernel_list, lcov = 'ldiff')
  
  # computing three ball kernel matrices and corresponding local covariance matrices
  kernel_params_ball <- c(0.5, 1, 2)
  ball_kernel_list <- 
    spatial_kernel_matrix(coords, 'ball', kernel_params_ball)
  loc_cov_ball <- 
    local_covariance_matrix(x = field, kernel_list = ball_kernel_list)
  
  # computing three gauss kernel matrices and corresponding local covariance matrices
  kernel_params_gauss <- c(0.5, 1, 2)
  gauss_kernel_list <- 
    spatial_kernel_matrix(coords, 'gauss', kernel_params_gauss)
  loc_cov_gauss <- 
    local_covariance_matrix(x = field, kernel_list = gauss_kernel_list)
}
}

