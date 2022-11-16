#-------------------------------------------#
# snss_sd generic and methods
#-------------------------------------------#
snss_sd <- function(x, ...) UseMethod("snss_sd")

snss_sd.default <- function(x, coords, direction = c('x', 'y'), ordered = TRUE, ...) {
  direction = match.arg(direction)
  blocks <- make_blocks(x = x, coords = coords, n_block = direction)
  snss_sd.list(x = blocks$x_list, coords = blocks$coords_list, 
               ordered = ordered)
}

snss_sd.SpatialPointsDataFrame <- function(x, ...) {
  result <- snss_sd.default(x = as.matrix(x@data), coords = x@coords, ...)
  x@data <- data.frame(result$s)
  result$s <- x
  result$coords <- NULL
  return(result)
}

snss_sd.sf <- function(x, ...) {
  if (!requireNamespace('sf', quietly = TRUE)) {
    stop('Please install the package sf to use this function.')
  } else {
    result <- snss_sd.default(x = as.matrix(sf::st_drop_geometry(x)), 
                              coords = sf::st_coordinates(x), ...)
    result$s <- sf::st_set_geometry(x = data.frame(result$s), value = sf::st_geometry(x))
    result$coords <- NULL
    return(result)    
  }
}

snss_sd.list <- function(x, coords, ordered = TRUE, ...) {
  # white data
  x_w <- white_data(do.call(rbind, x))
  
  # scatters
  cov_1_inv_sqrt <- white_data(x[[1]])$s_inv_sqrt
  cov_2 <- white_data(x[[2]])$s
  
  # diagonalization
  cov_d <- diag_scatters(cov_list = list(cov_1_inv_sqrt %*% tcrossprod(cov_2, cov_1_inv_sqrt)), 
                         ordered = ordered)
  
  # unmixing matrix
  w <- crossprod(cov_d$u, cov_1_inv_sqrt)
  w_inv <- crossprod(cov_1_inv_sqrt, cov_d$u)
  s <- tcrossprod(x_w$x_0, w)
  colnames(s) <- paste0('IC.', 1:ncol(s))
  
  return(structure(list(s = s, coords = do.call(rbind, coords), w = w, w_inv = w_inv, d = cov_d$u, 
                        x_mu = x_w$mu, cov_inv_sqrt = cov_1_inv_sqrt), 
                   class = c("snss", "sbss")))
}

#-------------------------------------------#
# snss_jd generic and methods
#-------------------------------------------#
snss_jd <- function(x, ...) UseMethod("snss_jd")

snss_jd.default <- function(x, coords, n_block, ordered = TRUE, ...) {
  blocks <- make_blocks(x = x, coords = coords, n_block = n_block)
  snss_jd.list(x = blocks$x_list, coords = blocks$coords_list, 
               ordered = ordered, ...)
}

snss_jd.SpatialPointsDataFrame <- function(x, ...) {
  result <- snss_jd.default(x = as.matrix(x@data), coords = x@coords, ...)
  x@data <- data.frame(result$s)
  result$s <- x
  result$coords <- NULL
  return(result)
}

snss_jd.sf <- function(x, ...) {
  if (!requireNamespace('sf', quietly = TRUE)) {
    stop('Please install the package sf to use this function.')
  } else {
    result <- snss_jd.default(x = as.matrix(sf::st_drop_geometry(x)), 
                              coords = sf::st_coordinates(x), ...)
    result$s <- sf::st_set_geometry(x = data.frame(result$s), value = sf::st_geometry(x))
    result$coords <- NULL
    return(result)    
  }
}

snss_jd.list <- function(x, coords, ordered = TRUE, ...) {
  # white data
  x_w <- white_data(do.call(rbind, x))
  x_w_block <- lapply(x, function(x)
    tcrossprod(sweep(x, MARGIN = 2, STATS = x_w$mu, FUN = '-'),
               x_w$s_inv_sqrt))
  
  # scatters
  cov_list <- lapply(x_w_block, function(x) crossprod(x) / (nrow(x) - 1))
  
  # diagonalization
  cov_d <- diag_scatters(cov_list = cov_list, 
                         ordered = ordered, ...)
  
  # unmixing matrix
  w <- crossprod(cov_d$u, x_w$s_inv_sqrt)
  w_inv <- crossprod(x_w$s_inv_sqrt, cov_d$u)
  s <- tcrossprod(x_w$x_0, w)
  colnames(s) <- paste0('IC.', 1:ncol(s))
  
  return(structure(list(s = s, coords = do.call(rbind, coords), w = w, w_inv = w_inv, d = cov_d$u, 
                        x_mu = x_w$mu, cov_inv_sqrt = x_w$s_inv_sqrt), 
                   class = c("snss", "sbss")))  
}


#-------------------------------------------#
# snss_sjd generic and methods
#-------------------------------------------#
snss_sjd <- function(x, ...) UseMethod("snss_sjd")

snss_sjd.default <- function(x, coords, n_block, kernel_type = c('ring', 'ball', 'gauss'), 
                             kernel_parameters, with_cov = TRUE, 
                             lcov = c('lcov', 'ldiff', 'lcov_norm'), ordered = TRUE, ...) {
  kernel_type <- match.arg(kernel_type)
  lcov <- match.arg(lcov)
  blocks <- make_blocks(x = x, coords = coords, n_block = n_block)
  snss_sjd.list(x = blocks$x_list, coords = blocks$coords_list, 
                kernel_type = kernel_type, kernel_parameters = kernel_parameters, 
                lcov = lcov, with_cov = with_cov, ordered = ordered,  ...)
}

snss_sjd.SpatialPointsDataFrame <- function(x, ...) {
  result <- snss_sjd.default(x = as.matrix(x@data), coords = x@coords, ...)
  x@data <- data.frame(result$s)
  result$s <- x
  result$coords <- NULL
  return(result)
}

snss_sjd.sf <- function(x, ...) {
  if (!requireNamespace('sf', quietly = TRUE)) {
    stop('Please install the package sf to use this function.')
  } else {
    result <- snss_sjd.default(x = as.matrix(sf::st_drop_geometry(x)), 
                              coords = sf::st_coordinates(x), ...)
    result$s <- sf::st_set_geometry(x = data.frame(result$s), value = sf::st_geometry(x))
    result$coords <- NULL
    return(result)    
  }
}

snss_sjd.list <- function(x, coords, kernel_type = c('ring', 'ball', 'gauss'), 
                          kernel_parameters, with_cov = TRUE, 
                          lcov = c('lcov', 'ldiff', 'lcov_norm'), ordered = TRUE, ...) {
  # inputs
  kernel_type <- match.arg(kernel_type)
  lcov <- match.arg(lcov)  
  
  # white data
  x_w <- white_data(do.call(rbind, x))
  x_w_block <- lapply(x, function(x)
    tcrossprod(sweep(x, MARGIN = 2, STATS = x_w$mu, FUN = '-'),
               x_w$s_inv_sqrt))
  # scatters
  kernel_list <- lapply(coords, function(coords)
      spatial_kernel_matrix(coords = coords, kernel_type = kernel_type, kernel_parameters = kernel_parameters))

  cov_list <- lapply(1:length(kernel_list), function(idx)
    local_covariance_matrix(x = x_w_block[[idx]], kernel_list = kernel_list[[idx]],
                            lcov = lcov, center = TRUE))
  cov_list <- do.call(c, cov_list)
  if (with_cov) {
    cov_list <- c(cov_list, lapply(x_w_block, function(x) crossprod(x) / (nrow(x) - 1)))
  }    

  # diagonalization
  cov_d <- diag_scatters(cov_list = cov_list, 
                         ordered = ordered, ...)
  
  # unmixing matrix
  w <- crossprod(cov_d$u, x_w$s_inv_sqrt)
  w_inv <- crossprod(x_w$s_inv_sqrt, cov_d$u)
  s <- tcrossprod(x_w$x_0, w)
  colnames(s) <- paste0('IC.', 1:ncol(s))
  
  return(structure(list(s = s, coords = do.call(rbind, coords), w = w, w_inv = w_inv, d = cov_d$u, 
                        x_mu = x_w$mu, cov_inv_sqrt = x_w$s_inv_sqrt), 
                   class = c("snss", "sbss")))  
}
