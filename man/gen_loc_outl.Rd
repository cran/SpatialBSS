\name{gen_loc_outl}
\alias{gen_loc_outl}

\title{
Contamination with Local Outliers
}

\description{
Generates synthetic local outliers and contaminates a given p-variate random field by swapping observations based on the first principal component score.
}

\usage{
gen_loc_outl(x, coords, alpha = 0.05, 
             neighborhood_type = c("radius", "fixed_n"), 
             radius = NULL, 
             neighborhood_size = NULL, 
             swap_order = c("regular", "reverse", "random"))

}

\arguments{
  \item{x}{
a numeric matrix of dimension \code{c(n, p)} where the \code{p} columns correspond to the entries of the random field and the \code{n} rows are the observations.
}
  \item{coords}{
a numeric matrix or data frame with dimension \code{c(n,2)} containing the coordinates of the observations.
}
  \item{alpha}{
a numeric value between 0 and 1 determining the proportion of the contaminated observations.
}
  \item{neighborhood_type}{
a string determining the type of neighborhood. If \code{'radius'}, each neighborhood contains all points within the radius determined by the parameter \code{radius}. If \code{'fixed_n'}, each neighborhood contains a constant number of closest points, where the constant is determined by the parameter \code{neighborhood_size}. Default is \code{'radius'}.
}
\item{radius}{
a positive numeric value defining the size of the radius when the \cr
\code{neighborhood_type} is \code{'radius'}. If \code{NULL} the radius defaults as \code{0.01*n}.
}
  \item{neighborhood_size}{
a positive integer defining the number of points in each neighborhood when the \code{neighborhood_type} is \code{'fixed_n'}. If \code{NULL} the number of points defaults as \code{ceiling(0.01*n)}.
}
  \item{swap_order}{
a string to determine which swap order is used. Either \code{'regular'} (default), \code{'reverse'} or \code{'random'}. See details.
}
}

\details{
\code{gen_loc_outl} generates local outliers by swapping the most extreme and the least extreme observations based on the first principal component score under the condition that at most one outliers lies in each neighborhood. For each location \eqn{s_i}, the neighborhood \eqn{N_i} is defined based on the parameter 
\code{neighborhood_type}. When \code{neighborhood_type} is \code{'radius'}, the neighborhood \eqn{N_i} contains all locations \eqn{s_j} for which the Euclidean norm \eqn{||s_i - s_j|| < r}, where \eqn{r} is determined by the parameter \code{radius}. When \code{neighborhood_type} is \code{'fixed_n'}, the neighborhood \eqn{N_i} contains \eqn{m - 1} nearest locations of \eqn{s_i}, where \eqn{m} is determined by the parameter \code{neighborhood_size}. For more details see Ernst & Haesbroeck, (2017).

After calculating the neighborhoods, the local outliers are generated following Ernst & Haesbroeck, (2017) and Harris
et al. (2014) using the steps:
\enumerate{
  \item{
    Sort the observations from highest to lowest by their principle component analysis (PCA) scores of the first component (PC-1).
  }
  \item{
    Set \eqn{k} to be \eqn{\alpha N/2} rounded to nearest integer and select the set of local outlier points \eqn{S^{out}} by finding \eqn{k} observations with the highest PC-1 values and \eqn{k} observations with the lowest PC-1 values under the condition that for all \eqn{s_i, s_j \in S_{out}} it holds that \eqn{N_i \neq N_j}.
  }
  \item{
  Form sets \eqn{X_{large}}, which contains \eqn{k} observations with the largest PC-1 values of outlier points \eqn{S_{out}} and \eqn{X^{small}}, which contains \eqn{k} observations with the smallest PC-1 values of outlier points \eqn{S^{out}}. Generate the local outliers by swapping \eqn{X^{small,i}} with \eqn{X^{large, k + 1 - i}}, \eqn{i=1,...,k}. The parameter \code{swap_order} defines how the sets \eqn{X^{large}} and \eqn{X^{small}} are ordered.
  }
}
If the parameter \code{swap_order} is \code{'regular'}, \eqn{X^{small}} and \eqn{X^{large}} are sorted by PC-1 score from smallest to largest.
If the parameter \code{swap_order} is \code{'reverse'}, \eqn{X^{small}} is sorted from largest to smallest and \eqn{X^{large}} from smallest to largest.
If the parameter \code{swap_order} is \code{'random'}, \eqn{X^{small}} and \eqn{X^{large}} are in random order.
}

\value{
\code{gen_loc_outl} returns a \code{data.frame} containing the contaminated fields as \eqn{p} first columns. The column \eqn{p + 1} contains a logical indicator whether the observation is an outlier or not.
}

\note{
This function is a modified version of code originally
provided by M. Ernst and G. Haesbroeck.
}

\references{
Ernst, M., & Haesbroeck, G. (2017). \emph{Comparison of local outlier detection techniques in spatial multivariate data}. Data Mining and Knowledge Discovery, 31 , 371-399. \doi{10.1007/s10618-016-0471-0}

Harris, P., Brunsdon, C., Charlton, M., Juggins, S., & Clarke, A. (2014). \emph{Multivariate spatial outlier detection using robust geographically weighted methods}. Mathematical Geosciences, 46 , 1-31. \doi{10.1007/s11004-013-9491-0}
}

\seealso{
\code{\link{gen_glob_outl}}
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
  
  # Generate 5 \% local outliers to data using radius neighborhoods 
  # and regular swap_order.
  field_cont <- gen_loc_outl(field, coords, alpha = 0.05,
                             neighborhood_type = "radius", 
                             radius = 0.5, swap_order = "regular")

  # Generate 10 \% local outliers to data using fixed_n neighborhoods 
  # and reverse swap_order.
  field_cont2 <- gen_loc_outl(field, coords, alpha = 0.1, 
                              neighborhood_type = "fixed_n", 
                              neighborhood_size = 10, swap_order = "reverse")
}
}

\keyword{ array }

