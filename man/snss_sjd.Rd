\name{snss_sjd}
\alias{snss_sjd}
\alias{snss_sjd.default}
\alias{snss_sjd.list}
\alias{snss_sjd.SpatialPointsDataFrame}
\alias{snss_sjd.sf}

\title{
Spatial Non-Stationary Source Separation Spatial Joint Diagonalization
}

\description{
\code{snss_sjd} estimates the unmixing matrix assuming a spatial non-stationary source separation model implying non-constant (spatial) covariance by jointly diagonalizing several covariance and/or spatial covariance matrices computed for a subdivision of the spatial domain into at least two sub-domains.
}

\usage{
snss_sjd(x, \dots)

\method{snss_sjd}{default}(x, coords, n_block, kernel_type = c('ring', 'ball', 'gauss'), 
     kernel_parameters, with_cov = TRUE, lcov = c('lcov', 'ldiff', 'lcov_norm'), 
     ordered = TRUE, \dots)
\method{snss_sjd}{list}(x, coords, kernel_type = c('ring', 'ball', 'gauss'), 
     kernel_parameters, with_cov = TRUE, lcov = c('lcov', 'ldiff', 'lcov_norm'), 
     ordered = TRUE, \dots)
\method{snss_sjd}{SpatialPointsDataFrame}(x, \dots)
\method{snss_sjd}{sf}(x, \dots)
}

\arguments{
  \item{x}{
either a numeric matrix of dimension \code{c(n, p)} where the \code{p} columns correspond to the entries of the random field and the \code{n} rows are the observations, a list of length \code{K} defining the subdivision of the domain, an object of class \code{\link[sf]{sf}} or an object of class \code{\link[sp]{SpatialPointsDataFrame}}. 
}
  \item{coords}{
a numeric matrix of dimension \code{c(n,2)} when \code{x} is a matrix where each row represents the sample location of a point in the spatial domain or a list of length \code{K} if \code{x} is a list which defines the subdivision of the domain. Not needed otherwise.
}
  \item{n_block}{
either be an integer defining the subdivision of the domain, \code{'x'} or \code{'y'}. See details.
}
  \item{kernel_type}{
a string indicating which kernel function to use. Either \code{'ring'} (default), \code{'ball'} or \code{'gauss'}.
}
  \item{kernel_parameters}{
a numeric vector that gives the parameters for the kernel function. At least length of one for \code{'ball'} and \code{'gauss'} or two for \code{'ring'} kernel, see details.
}
  \item{with_cov}{
logical. If \code{TRUE} not only spatial covariance matrices but also the sample covariances matrices for each sub-domain are considered in the joint diagonalization procedure. Default is \code{TRUE}.
}
  \item{lcov}{
a string indicating which type of local covariance matrix to use. Either \code{'lcov'} (default), \code{'ldiff'} or \code{'lcov_norm'}. See \code{\link{sbss_asymp}} for details on the latter option.
}
  \item{ordered}{
logical. If \code{TRUE} the entries of the latent field are ordered by the sum of squared pseudo-eigenvalues of the diagonalized sub-domain (local) covariance matrices. Default is \code{TRUE}.
}
  \item{\dots}{
further arguments for the fast real joint diagonalization algorithm that jointly diagonalizes the sub-domain covariance matrices. See details and \code{\link[JADE]{frjd}}.
}
}

\details{
This function assumes that the random field \eqn{x} is formed by \deqn{ x(t) = A s(t) + b, } where \eqn{A} is the deterministic \eqn{p \times p} mixing matrix, \eqn{b} is the \eqn{p}-dimensional location vector, \eqn{x} is the observable \eqn{p}-variate random field given by the argument \code{x}, \eqn{t} are the spatial locations given by the argument \code{coords} and \eqn{s} is the latent \eqn{p}-variate random field assumed to consist of uncorrelated entries that have zero mean but non-constant (spatial) second order dependence. This function aims to recover \eqn{s} by \deqn{ W(x(t) - \bar{x}), } where \eqn{W} is the \eqn{p \times p} unmixing matrix and \eqn{\bar{x}} is the sample mean. The function does this by splitting the given spatial domain into \code{n_block^2} equally sized rectangular sub-domains and jointly diagonalizing the corresponding spatial covariance matrices for all sub-domains. If the argument \code{with_cov} equals \code{TRUE} (default) then additionally also the sample covariance matrices for each sub-domain are included in the joint diagonalization procedure.

The arguments \code{kernel_type}, \code{kernel_parameters} and \code{lcov} determine which spatial kernel functions and which type of local covariance matrices are used for each sub-domain. The usage is equal to the function \code{\link{sbss}}.

Alternatively the domain subdivision can be defined by providing lists of length \code{K} for the arguments \code{x} and \code{coords} where the first list entries correspond to the values and coordinates of the first sub-domain and the second entries to the values and coordinates of the second sub-domain, etc.. The argument \code{n_block} might be \code{'x'} or \code{'y'} indicating a split across the x or y coordinates similar as done by the function \code{\link{snss_sd}}.

\code{snss_sjd} jointly diagonalizes the covariance matrices for each sub-domain with the function \code{\link[JADE]{frjd}}. \code{\dots} provides arguments for \code{frjd}, useful arguments might be:
\itemize{
  \item \code{eps}: tolerance for convergence.
  \item \code{maxiter}: maximum number of iterations.
}

}

\value{
Similarly as \code{\link{sbss}} the function \code{snss_jd} returns a list of class \code{'snss'} and \code{'sbss'} with the following entries: 
  \item{s}{
  object of \code{class(x)} containing the estimated source random field.
}  
  \item{coords}{
  coordinates of the observations. Only given if \code{x} is a matrix or list.
}

  \item{w}{
  estimated unmixing matrix.
}

  \item{w_inv}{
  inverse of the estimated unmixing matrix.
}

  \item{d}{
  matrix of stacked (jointly) diagonalized sub-domain covariance and/or local covariance matrices.
}

  \item{x_mu}{
  columnmeans of \code{x}.
}

  \item{cov_inv_sqrt}{
  square root of the inverse sample covariance matrix of \code{x}.
}

}

\references{
Muehlmann, C., Bachoc, F. and Nordhausen, K. (2021), \emph{Blind Source Separation for Non-Stationary Random Fields}, Submitted for publication. Preprint available at \url{https://arxiv.org/abs/2107.01916}.
}

\seealso{
\code{\link{sbss}}, \code{\link[sp]{sp}}, \code{\link[sf]{sf}}
}

\examples{
# simulate coordinates
n <- 1000
coords <- runif(n * 2) * 20
dim(coords) <- c(n, 2)

# simulate random field
field_1 <- rnorm(n)
field_2 <- 2 * sin(pi / 20 * coords[, 1]) * rnorm(n)
field_3 <- rnorm(n) * (coords[, 1] < 10) + rnorm(n, 0, 3) * (coords[, 1] >= 10)

latent_field <- cbind(field_1, field_2, field_3)
mixing_matrix <- matrix(rnorm(9), 3, 3)
observed_field <- latent_field %*% t(mixing_matrix)

observed_field_sp <- sp::SpatialPointsDataFrame(coords = coords, 
                                                data = data.frame(observed_field))
sp::spplot(observed_field_sp, colorkey = TRUE, as.table = TRUE, cex = 1)

# apply snss_sjd with 4 sub-domains
# one ring kernel per sub-domain
# without covariances
res_4_ball <- snss_sjd(observed_field, coords, n_block = 2,
                  kernel_type = 'ball', kernel_parameters = c(0, 2), 
                  with_cov = TRUE)
JADE::MD(W.hat = coef(res_4_ball), A = mixing_matrix)

# apply snss_sjd with split across y
# one ring kernel per sub-domain
# without covariances
# should not work as field does not show spatial dependence
res_4_ring <- snss_sjd(observed_field, coords, n_block = 'y',
                       kernel_type = 'ring', kernel_parameters = c(0, 2), 
                       with_cov = FALSE)
JADE::MD(W.hat = coef(res_4_ring), A = mixing_matrix)

# print object
print(res_4_ball)

# plot latent field
plot(res_4_ball, colorkey = TRUE, as.table = TRUE, cex = 1)

# predict latent fields on grid
predict(res_4_ball, colorkey = TRUE, as.table = TRUE, cex = 1)

# unmixing matrix
w_unmix <- coef(res_4_ball)

# apply snss_jd with SpatialPointsDataFrame object 
res_4_ball_sp <- snss_sjd(observed_field_sp, n_block = 2,
                          kernel_type = 'ball', kernel_parameters = c(0, 2), 
                          with_cov = TRUE)

# apply with list arguments
# first axis split by 5
# second axis split by 10
# results in 4 sub-domains
flag_x <- coords[, 1] < 5
flag_y <- coords[, 2] < 10
coords_list <- list(coords[flag_x & flag_y, ],
                    coords[!flag_x & flag_y, ],
                    coords[flag_x & !flag_y, ],
                    coords[!flag_x & !flag_y, ])
field_list <- list(observed_field[flag_x & flag_y, ],
                   observed_field[!flag_x & flag_y, ],
                   observed_field[flag_x & !flag_y, ],
                   observed_field[!flag_x & !flag_y, ])
plot(coords, col = 1)
points(coords_list[[2]], col = 2)
points(coords_list[[3]], col = 3)
points(coords_list[[4]], col = 4)

res_list <- snss_sjd(x = field_list,
                    coords = coords_list,
                    kernel_type = 'ring', kernel_parameters = c(0, 2))
plot(res_list, colorkey = TRUE, as.table = TRUE, cex = 1)
JADE::MD(W.hat = coef(res_list), A = mixing_matrix)

}

\keyword{ multivariate }
\keyword{ spatial }
