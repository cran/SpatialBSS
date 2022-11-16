#-------------------------------------------#
# robsbss generic and methods
#-------------------------------------------#
robsbss <- function(x, ...) UseMethod("robsbss")

robsbss.default <- function(x, coords, kernel_type = c('ring', 'ball', 'gauss'), kernel_parameters, 
                  lcov = c('norm', 'winsor', 'qwinsor'), 
                  ordered = TRUE, kernel_list = NULL, ...) {
  # To avoid problems when x or coords are of class data.frame:
  x <- as.matrix(x)
  if (!missing(coords)) {
    coords <- as.matrix(coords)
  }
  # kernel matrix
  kernel_type <- match.arg(kernel_type)
  lcov <- match.arg(lcov)
  
  if (!missing(coords) && !missing(kernel_parameters) && is.vector(kernel_parameters)) {
    kernel_list <- spatial_kernel_matrix(coords, kernel_type = kernel_type, 
                                         kernel_parameters = kernel_parameters)
  } else if (!is.null(kernel_list) && is.list(kernel_list)) {
    if (missing(coords)) {
      coords <- NULL
    }
  } else {
    stop('Invalid input for kernels. Either coords (or a spatial object for the argument x), kernel_type and kernel_parameters (as vector) or kernel_list needs to be given.')
  }
  k <- length(kernel_list)
  
  # standardize data
  x_w <- white_data(x, whitening = "hr")
  
  # spatial covariance matrices
  cov_sp_list <- local_gss_covariance_matrix(x = x_w$x_w, kernel_list = kernel_list, 
                                            lcov = lcov, center = FALSE)
  
  # diagonalization
  cov_sp_d <- diag_scatters(cov_list = cov_sp_list$cov_sp_list, ordered = ordered, ...)
  
  # unmixing matrix
  w <- crossprod(cov_sp_d$u, x_w$s_inv_sqrt)
  w_inv <- crossprod(x_w$s_sqrt, cov_sp_d$u)
  s <- tcrossprod(x_w$x_0, w)
  colnames(s) <- paste0('IC.', 1:ncol(s))
  
  # results
  return(structure(list(s = s, coords = coords, w = w, weights = cov_sp_list$weights, w_inv = w_inv, 
                        pevals = cov_sp_d$pevals, d = cov_sp_d$d, diags = cov_sp_d$diags, 
                        x_mu = x_w$mu, cov_inv_sqrt = x_w$s_inv_sqrt), 
                   class = c("sbss")))
}


robsbss.SpatialPointsDataFrame <- function(x, ...) {
  result <- robsbss.default(x = as.matrix(x@data), coords = x@coords, ...)
  x@data <- data.frame(result$s)
  result$s <- x
  result$coords <- NULL
  return(result)
}

robsbss.sf <- function(x, ...) {
  if (!requireNamespace('sf', quietly = TRUE)) {
    stop('Please install the package sf to use this function.')
  } else {
    result <- robsbss.default(x = as.matrix(sf::st_drop_geometry(x)), coords = sf::st_coordinates(x), ...)
    result$s <- sf::st_set_geometry(x = data.frame(result$s), value = sf::st_geometry(x))
    result$coords <- NULL
    return(result)    
  }
}

