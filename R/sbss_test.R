#------------------------------------------------#
#             asymptotic test
#------------------------------------------------#
sbss_asymp <- function(x, ...) UseMethod("sbss_asymp")

sbss_asymp.default <- function(x, coords, q, kernel_parameters, 
                               kernel_list = NULL, ...) {
  # parameters
  data_name <- deparse(substitute(x))
  n_coords <- nrow(x)
  p <- ncol(x)  
  if(q >= p || q < 0) {
    stop('The number of hypothetical signals q must be between zero and ncol(x) - 1')
  }
  
  # apply sbss
  if (!missing(coords) && !missing(kernel_parameters) && is.vector(kernel_parameters)) {
    kernel_list <- spatial_kernel_matrix(coords, kernel_type = 'ring', kernel_parameters = kernel_parameters)
  } else if (!is.null(kernel_list) && is.list(kernel_list)) {
    if (missing(coords)) {
      coords <- NULL
    }
  } else {
    stop('Invalid input for kernels. Either coords (or a spatial object for the argument x) and kernel_parameters (as vector) or kernel_list needs to be given.')
  }
  k <- length(kernel_list)
  sbss_res <- sbss(x = x, coords = coords, kernel_list = kernel_list,
                   lcov = 'lcov_norm', ordered = TRUE, ...)

  # test
  t <- test_stat(sbss_res$d, q, n_coords)
  df <- k * (p - q) * (p - q + 1) / 2
  p_val <- stats::pchisq(t, df = df, lower.tail = FALSE, log.p = FALSE)
  
  # results
  names(t) <- 'T'
  names(df) <- 'df'
  res <- c(list(alternative = paste0("there are less than ", p - q, 
                                     " white noise components"), 
                method = c("SBSS asymptotic test for white noise processes"), 
                data.name = data_name,
                statistic = t,
                parameters = df,
                p.value = p_val), sbss_res)
  class(res) <- c('sbss_test', 'htest', 'sbss')
  
  return(res)
} 

sbss_asymp.SpatialPointsDataFrame <- function(x, ...) {
  result <- sbss_asymp.default(x = as.matrix(x@data), coords = x@coords, ...)
  x@data <- data.frame(result$s)
  result$s <- x
  result$coords <- NULL
  return(result)
}

sbss_asymp.sf <- function(x, ...) {
  if (!requireNamespace('sf', quietly = TRUE)) {
    stop('Please install the package sf to use this function.')
  } else {
    result <- sbss_asymp.default(x = as.matrix(sf::st_drop_geometry(x)), coords = sf::st_coordinates(x), ...)
    result$s <- sf::st_set_geometry(x = data.frame(result$s), value = sf::st_geometry(x))
    result$coords <- NULL
    return(result)    
  }
}

#------------------------------------------------#
#             bootstrap test
#------------------------------------------------#
sbss_boot <- function(x, ...) UseMethod("sbss_boot")

sbss_boot.default <- function(x, coords, q, kernel_parameters, boot_method = c('permute', 'parametric'), 
                              n_boot = 200, kernel_list = NULL, ...) {
  # parameters
  data_name <- deparse(substitute(x))
  n_coords <- nrow(x)
  boot_method <- match.arg(boot_method)
  p <- ncol(x)  
  if(q >= p || q < 0) {
    stop('The number of hypothetical signals q must be between zero and ncol(x) - 1')
  }
  
  # apply sbss
  if (!missing(coords) && !missing(kernel_parameters) && is.vector(kernel_parameters)) {
    kernel_list <- spatial_kernel_matrix(coords, kernel_type = 'ring', kernel_parameters = kernel_parameters)
  } else if (!is.null(kernel_list) && is.list(kernel_list)) {
    if (missing(coords)) {
      coords <- NULL
    }
  } else {
    stop('Invalid input for kernels. Either coords (or a spatial object for the argument x) and kernel_parameters (as vector) or kernel_list needs to be given.')
  }
  sbss_res <- sbss(x = x, coords = coords, kernel_list = kernel_list,
                   lcov = 'lcov_norm', ordered = TRUE, ...)
  t <- test_stat(sbss_res$d, q, n_coords) 

  # bootstrap loop
  t_boot <- vector(mode = 'numeric', length = n_boot)
  for (idx in 1:n_boot) {
    if (boot_method == 'permute') {
      noise <- sample(as.vector(sbss_res$s[, (q + 1):p]))
      dim(noise) <- c(n_coords, (p - q))
    } else if (boot_method == 'parametric') {
      noise <- stats::rnorm(n_coords * (p - q))
      dim(noise) <- c(n_coords, (p - q))
    } else {
      stop('Bootstrap method is not supported!')
    }
    x_bs <- tcrossprod(if (q == 0) noise else cbind(sbss_res$s[, 1:q], noise), 
                       sbss_res$w_inv)
    x_bs <- sweep(x_bs, MARGIN = 2, STATS = sbss_res$x_mu, FUN = '+')
    
    # apply sbss
    sbss_bs <- sbss(x = x_bs, kernel_list = kernel_list,
                    lcov = 'lcov_norm', ordered = TRUE, ...)
    
    # test statistic for bs sample
    t_boot[idx] <- test_stat(sbss_bs$d, q, n_coords)
  }
  p_val <- (sum(t_boot >= t) + 1) / (n_boot + 1) 
  
  # results
  names(t) <- 'T'
  names(n_boot) <- 'replications'
  res <- c(list(alternative = paste0("there are less than ", p - q, 
                                     " white noise components"), 
                method = paste0("SBSS ", boot_method, " bootstrap test for white noise processes"), 
                data.name = data_name,
                statistic = t,
                parameters = n_boot,
                p.value = p_val,
                t_boot = t_boot), sbss_res)
  class(res) <- c('sbss_test', 'htest', 'sbss')
  
  return(res)
}

sbss_boot.SpatialPointsDataFrame <- function(x, ...) {
  result <- sbss_boot.default(x = as.matrix(x@data), coords = x@coords, ...)
  x@data <- data.frame(result$s)
  result$s <- x
  result$coords <- NULL
  return(result)
}

sbss_boot.sf <- function(x, ...) {
  if (!requireNamespace('sf', quietly = TRUE)) {
    stop('Please install the package sf to use this function.')
  } else {
    result <- sbss_boot.default(x = as.matrix(sf::st_drop_geometry(x)), coords = sf::st_coordinates(x), ...)
    result$s <- sf::st_set_geometry(x = data.frame(result$s), value = sf::st_geometry(x))
    result$coords <- NULL
    return(result)    
  }
}
