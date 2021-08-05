\name{white_data}
\alias{white_data}

\title{
Different Approaches of Data Whitening
}

\description{
\code{white_data} whites the data with respect to the sample covariance matrix, or different spatial scatter matrices.
}

\usage{
white_data(x, rob_whitening = FALSE, lcov = c('lcov', 'ldiff', 'lcov_norm'), 
           kernel_mat = numeric(0))
}

\arguments{
  \item{x}{
a numeric matrix of dimension \code{c(n, p)} where the \code{p} columns correspond to the entries of the random field and the \code{n} rows are the observations.
}
  \item{rob_whitening}{
logical. If \code{TRUE} whitening is carried out with respect to the first spatial scatter matrix and not the sample covariance matrix, see details. Default is \code{FALSE}.
}
  \item{lcov}{
a string indicating which type of local covariance matrix is used for whitening. Either \code{'lcov'} (default) or \code{'ldiff'}.
}
  \item{kernel_mat}{
a spatial kernel matrix with dimension \code{c(n,n)}, see details. Usually computed by the function \code{\link{spatial_kernel_matrix}}.
}
}

\details{
The inverse square root of a positive definite matrix \eqn{M} with eigenvalue decomposition \eqn{UDU'} is defined as \eqn{M^{-1/2} = UD^{-1/2}U'}. \code{white_data} whitens the data by \eqn{M^{-1/2}(x - mu)} where \eqn{mu} are the column means of \code{x} and the matrix \eqn{M} is as follows. If the argument \code{rob_whitening} is \code{FALSE} then \eqn{M} is the sample covariance matrix. If the argument \code{rob_whitening} is \code{TRUE}, then the argument \code{lcov} determines the matrix \eqn{M} to be one of the following local scatter matrices:
\itemize{
  \item \code{'lcov'}: \deqn{ LCov(f) = 1/n \sum_{i,j} f(d_{i,j}) (x(s_i)-\bar{x}) (x(s_j)-\bar{x})' ,}
  \item \code{'ldiff'}: \deqn{ LDiff(f) = 1/n \sum_{i,j} f(d_{i,j}) (x(s_i)-x(s_j)) (x(s_i)-x(s_j))',}
  \item \code{'lcov_norm'}: \deqn{ LCov^*(f) = 1/(n F^{1/2}_{f,n}) \sum_{i,j} f(d_{i,j}) (x(s_i)-\bar{x}) (x(s_j)-\bar{x})',} with \deqn{ F_{f,n} = 1 / n \sum_{i,j} f^2(d_{i,j}).}
}

Where \eqn{d_{i,j} \ge 0} correspond to the pairwise distances between coordinates, \eqn{x(s_i)} are the \code{p} random field values at location \eqn{s_i}, \eqn{\bar{x}} is the sample mean vector, and the kernel function \eqn{f(d)} determines the locality. The choice \code{'lcov_norm'} is useful when testing for the actual signal dimension of the latent field, see \code{\link{sbss_asymp}} and \code{\link{sbss_boot}}. See also \code{\link{sbss}} for details.

Note that \eqn{LCov(f)} are usually not positive definite, therefore in that case the matrix cannot be inverted and an error is produced. Whitening with \eqn{LCov(f)} matrices might be favourable in the presence of spatially uncorrelated noise, and whitening with \eqn{LDiff(f)} might be favourable when a non-constant smooth drift is present in the data.

The argument \code{kernel_mat} is a matrix of dimension \code{c(n,n)} where each entry corresponds to the spatial kernel function evaluated at the distance between two sample locations, mathematically the entry ij of each kernel matrix is \eqn{f(d_{i,j})}. This matrix is usually computed with the function \code{\link{spatial_kernel_matrix}}.
}

\value{
\code{white_data} returns a list with the following entries: 
  \item{mu}{
  a numeric vector of length \code{ncol(x)} containing the column means of the data matrix \code{x}.
}  
  \item{x_0}{
  a numeric matrix of dimension \code{c(n, p)} containing the columns centerd data of \code{x}.
}

  \item{x_w}{
  a numeric matrix of dimension \code{c(n, p)} containing the whitened data of \code{x}.
}

  \item{s_inv_sqrt}{
  a numeric matrix of dimension \code{c(p, p)} which equals the inverse square root of the scatter matrix \eqn{M} used for whitening.
}

  \item{s_sqrt}{
  a numeric matrix of dimension \code{c(p, p)} which equals the square root of the scatter matrix \eqn{M}.
}

}

\references{
Muehlmann, C., Filzmoser, P. and Nordhausen, K. (2021), \emph{Spatial Blind Source Separation in the Presence of a Drift}, Submitted for publication.

Bachoc, F., Genton, M. G, Nordhausen, K., Ruiz-Gazen, A. and Virta, J. (2020), \emph{Spatial Blind Source Separation}, Biometrika, \doi{10.1093/biomet/asz079}.
}

\seealso{
\code{\link{sbss}}, \code{\link{spatial_kernel_matrix}}
}

\examples{
# simulate a set of coordinates
coords <- rnorm(100 * 2)
dim(coords) <- c(100, 2)

# simulate random field
if (!requireNamespace('RandomFields', quietly = TRUE)) {
  stop('Please install the package RandomFields to run the example code.')
} else {
  RandomFields::RFoptions(spConform = FALSE)
  field_1 <- RandomFields::RFsimulate(model = RandomFields::RMexp(), 
                                      x = coords)
  field_2 <- RandomFields::RFsimulate(model = RandomFields::RMspheric(), 
                                      x = coords)
  field_3 <- RandomFields::RFsimulate(model = RandomFields::RMwhittle(nu = 2),  
                                      x = coords)
  field <-  cbind(field_1, field_2, field_3) %*% matrix(rnorm(9), 3, 3) 
}

# white the data with the usual sample covariance 
x_w_1 <- white_data(field)

# white the data with a ldiff matrix and ring kernel
kernel_params_ring <- c(0, 1)
ring_kernel_list <- 
  spatial_kernel_matrix(coords, 'ring', kernel_params_ring)
x_w_2 <- white_data(field, rob_whitening = TRUE,
  lcov = 'ldiff', kernel_mat = ring_kernel_list[[1]])
}

\keyword{ array }

