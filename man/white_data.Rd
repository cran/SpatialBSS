\name{white_data}
\alias{white_data}

\title{
Different Approaches of Data Whitening
}

\description{
\code{white_data} whites the data with respect to the sample covariance matrix, or different spatial scatter matrices.
}

\usage{
white_data(x, whitening = c("standard", "rob", "hr"), 
           lcov = c('lcov', 'ldiff', 'lcov_norm'), 
           kernel_mat = numeric(0))
}

\arguments{
  \item{x}{
a numeric matrix of dimension \code{c(n, p)} where the \code{p} columns correspond to the entries of the random field and the \code{n} rows are the observations.
}
  \item{whitening}{
a string indicating the whitening method. If \code{'standard'} then the whitening is carried out with respect to sample covariance matrix, if \code{'rob'} then the first spatial scatter matrix is used instead of sample the covariance matrix and if \code{'hr'} then the Hettmansperger-Randles location and scatter estimates are used for whitening. See details for more. Default is \code{'standard'}.
}
  \item{lcov}{
a string indicating which type of local covariance matrix is used for whitening, when the whitening method \code{'rob'} is used. Either \code{'lcov'} (default) or \code{'ldiff'}.
}
  \item{kernel_mat}{
a spatial kernel matrix with dimension \code{c(n,n)}, see details. Usually computed by the function \code{\link{spatial_kernel_matrix}}.
}
}

\details{
The inverse square root of a positive definite matrix \eqn{M(x)} with eigenvalue decomposition \eqn{UDU'} is defined as \eqn{M(x)^{-1/2} = UD^{-1/2}U'}. \code{white_data} whitens the data by \eqn{M(x)^{-1/2}(x - T(x))} where \eqn{T(x)} is a location functional of \eqn{x} and the matrix \eqn{M(x)} is a scatter functional. If the argument \code{whitening} is \code{'standard'}, \eqn{M(x)} is the sample covariance matrix and \eqn{T(x)} is a vector of column means of \eqn{x}. If the argument \code{whitening} is \code{'hr'}, the Hettmansperger-Randles location and scatter estimates (Hettmansperger & Randles, 2002) are used as location functional \eqn{T(x)} and scatter functional \eqn{M(x)}. The Hettmansperger-Randles location and scatter estimates are robust variants of sample mean and covariance matrices, that are used for whitening in \code{\link{robsbss}}. If the argument \code{whitening} is \code{'rob'}, the argument \code{lcov} determines the scatter functional \eqn{M(x)} to be one of the following local scatter matrices:
\itemize{
  \item \code{'lcov'}: \deqn{ LCov(f) = 1/n \sum_{i,j} f(d_{i,j}) (x(s_i)-\bar{x}) (x(s_j)-\bar{x})' ,}
  \item \code{'ldiff'}: \deqn{ LDiff(f) = 1/n \sum_{i,j} f(d_{i,j}) (x(s_i)-x(s_j)) (x(s_i)-x(s_j))',}
  \item \code{'lcov_norm'}: \deqn{ LCov^*(f) = 1/(n F^{1/2}_{f,n}) \sum_{i,j} f(d_{i,j}) (x(s_i)-\bar{x}) (x(s_j)-\bar{x})',} with \deqn{ F_{f,n} = 1 / n \sum_{i,j} f^2(d_{i,j}),}
}
where \eqn{d_{i,j} \ge 0} correspond to the pairwise distances between coordinates, \eqn{x(s_i)} are the \code{p} random field values at location \eqn{s_i}, \eqn{\bar{x}} is the sample mean vector, and the kernel function \eqn{f(d)} determines the locality. The choice \code{'lcov_norm'} is useful when testing for the actual signal dimension of the latent field, see \code{\link{sbss_asymp}} and \code{\link{sbss_boot}}. See also \code{\link{sbss}} for details.

Note that \eqn{LCov(f)} are usually not positive definite, therefore in that case the matrix cannot be inverted and an error is produced. Whitening with \eqn{LCov(f)} matrices might be favorable in the presence of spatially uncorrelated noise, and whitening with \eqn{LDiff(f)} might be favorable when a non-constant smooth drift is present in the data.

The argument \code{kernel_mat} is a matrix of dimension \code{c(n,n)} where each entry corresponds to the spatial kernel function evaluated at the distance between two sample locations, mathematically the entry ij of each kernel matrix is \eqn{f(d_{i,j})}. This matrix is usually computed with the function \code{\link{spatial_kernel_matrix}}.
}

\value{
\code{white_data} returns a list with the following entries: 
  \item{mu}{
  a numeric vector of length \code{ncol(x)} containing the column means of the data matrix \code{x}.
}  
  \item{x_0}{
  a numeric matrix of dimension \code{c(n, p)} containing the columns centered data of \code{x}.
}

  \item{x_w}{
  a numeric matrix of dimension \code{c(n, p)} containing the whitened data of \code{x}.
}

  \item{s}{
  a numeric matrix of dimension \code{c(p, p)} which is the scatter matrix \eqn{M}.
}

  \item{s_inv_sqrt}{
  a numeric matrix of dimension \code{c(p, p)} which equals the inverse square root of the scatter matrix \eqn{M} used for whitening.
}

  \item{s_sqrt}{
  a numeric matrix of dimension \code{c(p, p)} which equals the square root of the scatter matrix \eqn{M}.
}

}

\references{
Muehlmann, C., Filzmoser, P. and Nordhausen, K. (2021), \emph{Spatial Blind Source Separation in the Presence of a Drift}, Submitted for publication. Preprint available at \url{https://arxiv.org/abs/2108.13813}.

Bachoc, F., Genton, M. G, Nordhausen, K., Ruiz-Gazen, A. and Virta, J. (2020), \emph{Spatial Blind Source Separation}, Biometrika, 107, 627-646, \doi{10.1093/biomet/asz079}.

Hettmansperger, T. P., & Randles, R. H. (2002). \emph{A practical affine equivariant multivariate median}. Biometrika, 89 , 851-860. \doi{10.1093/biomet/89.4.851}.
}

\seealso{
\code{\link{sbss}}, \code{\link{spatial_kernel_matrix}}
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
  X <- as.matrix(field)

  # white the data with the usual sample covariance 
  x_w_1 <- white_data(X)
  
  # white the data with a ldiff matrix and ring kernel
  kernel_params_ring <- c(0, 1)
  ring_kernel_list <- 
    spatial_kernel_matrix(coords, 'ring', kernel_params_ring)
  x_w_2 <- white_data(field, whitening = 'rob',
    lcov = 'ldiff', kernel_mat = ring_kernel_list[[1]])
  
  # Generate 5 \% of global outliers to data
  field_cont <- gen_glob_outl(field)[,1:3]
  X <- as.matrix(field_cont)
  # white the data using Hettmansperger-Randles location and scatter estimates
  x_w_3 <- white_data(X, whitening = 'hr')
}
}

\keyword{ array }

