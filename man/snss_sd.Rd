\name{snss_sd}
\alias{snss_sd}
\alias{snss_sd.default}
\alias{snss_sd.list}
\alias{snss_sd.SpatialPointsDataFrame}
\alias{snss_sd.sf}

\title{
Spatial Non-Stationary Source Separation Simultaneous Diagonalization
}

\description{
\code{snss_sd} estimates the unmixing matrix assuming a spatial non-stationary source separation model implying non-constant covariance by simultaneously diagonalizing two covariance matrices computed for two corresponding different sub-domains.
}

\usage{
snss_sd(x, \dots)

\method{snss_sd}{default}(x, coords, direction = c('x', 'y'), 
     ordered = TRUE, \dots)
\method{snss_sd}{list}(x, coords, ordered = TRUE, \dots)
\method{snss_sd}{SpatialPointsDataFrame}(x, \dots)
\method{snss_sd}{sf}(x, \dots)
}

\arguments{
  \item{x}{
either a numeric matrix of dimension \code{c(n, p)} where the \code{p} columns correspond to the entries of the random field and the \code{n} rows are the observations, a list of length two defining the subdivision of the domain, an object of class \code{\link[sf]{sf}} or an object of class \code{\link[sp]{SpatialPointsDataFrame}}. 
}
  \item{coords}{
a numeric matrix of dimension \code{c(n,2)} when \code{x} is a matrix where each row represents the sample location of a point in the spatial domain or a list of length two if \code{x} is a list which defines the subdivision of the domain. Not needed otherwise.
}
  \item{direction}{
a string indicating on which coordinate axis the domain is halved. Either \code{'x'} (default) or \code{'y'}.
}
  \item{ordered}{
logical. If \code{TRUE} the entries of the latent field are ordered according to the decreasingly ordered eigenvalues. Default is \code{TRUE}.
}
  \item{\dots}{
further arguments to be passed to or from methods.
}
}

\details{
This function assumes that the random field \eqn{x} is formed by \deqn{ x(t) = A s(t) + b, } where \eqn{A} is the deterministic \eqn{p \times p} mixing matrix, \eqn{b} is the \eqn{p}-dimensional location vector, \eqn{x} is the observable \eqn{p}-variate random field given by the argument \code{x}, \eqn{t} are the spatial locations given by the argument \code{coords} and \eqn{s} is the latent \eqn{p}-variate random field assumed to consist of uncorrelated entries that have zero mean but non-constant variances. This function aims to recover \eqn{s} by \deqn{ W(x(t) - \bar{x}), } where \eqn{W} is the \eqn{p \times p} unmixing matrix and \eqn{\bar{x}} is the sample mean. The function does this by splitting the given spatial domain in half according to the first coordinate (argument \code{direction} equals \code{'x'}) or the second coodinate (argument \code{direction} equals \code{'y'}) and simultaneously diagonalizing the sample covariance matrices for each of the two sub-domains.

Alternatively the domain subdivison can be defined by providing lists of length two for the arguments \code{x} and \code{coords} where the first list entries correspond to the values and coordinates of the first sub-domain and the second entries to the values and coordinates of the second sub-domain.

}

\value{
Similarly as \code{\link{sbss}} the function \code{snss_sd} returns a list of class \code{'snss'} and \code{'sbss'} with the following entries: 
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
  diagonal matrix containing the eigenvalues of the eigendecomposition.
}

  \item{x_mu}{
  columnmeans of \code{x}.
}

  \item{cov_inv_sqrt}{
  square root of the inverse sample covariance matrix for the first sub-domain.
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

# apply snss_sd with split in x 
res_x <- snss_sd(observed_field, coords, direction = 'x')
JADE::MD(W.hat = coef(res_x), A = mixing_matrix)

# apply snss_sd with split in y
# should be much worse as field shows only variation in x
res_y <- snss_sd(observed_field, coords, direction = 'y')
JADE::MD(W.hat = coef(res_y), A = mixing_matrix)

# print object
print(res_x)

# plot latent field
plot(res_x, colorkey = TRUE, as.table = TRUE, cex = 1)

# predict latent fields on grid
predict(res_x, colorkey = TRUE, as.table = TRUE, cex = 1)

# unmixing matrix
w_unmix <- coef(res_x)

# apply snss_sd with SpatialPointsDataFrame object 
res_x_sp <- snss_sd(observed_field_sp, direction = 'x')

# apply with list arguments
# first axis split by 5
flag_coords <- coords[, 1] < 5
coords_list <- list(coords[flag_coords, ],
                    coords[!flag_coords, ])
field_list <- list(observed_field[flag_coords, ],
                   observed_field[!flag_coords, ])
plot(coords, col = flag_coords + 1)

res_list <- snss_sd(x = field_list,
                    coords = coords_list)
plot(res_list, colorkey = TRUE, as.table = TRUE, cex = 1)
JADE::MD(W.hat = coef(res_list), A = mixing_matrix)

}

\keyword{ multivariate }
\keyword{ spatial }
