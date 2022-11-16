\name{gen_glob_outl}
\alias{gen_glob_outl}

\title{
Contamination with Global Outliers
}

\description{
Generates synthetic global outliers and contaminates a given p-variate random field
}

\usage{
gen_glob_outl(x, alpha = 0.05, h = 10, random_sign = FALSE)
}

\arguments{
  \item{x}{
a numeric matrix of dimension \code{c(n, p)} where the \code{p} columns correspond to the entries of the random field and the \code{n} rows are the observations.
}
  \item{alpha}{
a numerical value between 0 and 1 giving the proportion of observations to contaminate.
}
  \item{h}{
a numerical constant to determine how large the contaminated outliers are, see details.
}
  \item{random_sign}{
logical. If \code{TRUE}, the sign of each component of the outlier is randomly selected. Default is \code{FALSE}. See more in details.
}
}

\details{
\code{gen_glob_outl} generates outliers for a given field by selecting randomly \code{round(alpha * n)} observations \eqn{x_i} to be the outliers and contaminating them by setting \eqn{x^{out}_i = (c^i)'x_i}, where the elements \eqn{c^i_j} of vector \eqn{c^i} are determined by the parameter \code{random_sign}. If \code{random_sign = TRUE}, \eqn{c^i_j} is either \eqn{h} or \eqn{-h} with \eqn{P(c^i_j = h) = P(c^i_j = -h) = 0.5}. If \code{random_sign = FALSE}, \eqn{c^i_j=h} for all \eqn{j=1,...p}, \eqn{i=1,...,n}. The parameter \code{alpha} determines the contamination rate \eqn{\alpha} and the parameter \code{h} determines the size of the outliers.
}

\value{
\code{gen_glob_outl} returns a \code{data.frame} containing the contaminated fields as \eqn{p} first columns. The column \eqn{p + 1} contains a logical indicator whether the observation is outlier or not.
}

\seealso{
\code{\link{gen_loc_outl}}
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
  # Generate 10 \% global outliers to data, with size h=15.
  field_cont <- gen_glob_outl(field, alpha = 0.1, h = 15)
  
  # Generate 5 \% global outliers to data, with size h = 10 and random sign.
  field_cont2 <- gen_glob_outl(field, alpha = 0.05, h = 10, random_sign = TRUE)
}
}

\keyword{ array }

